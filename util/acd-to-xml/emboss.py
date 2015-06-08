# EMBOSS Parsing Library
from lxml import etree as ET
import shlex
import subprocess
import re
import sys
from bs4 import BeautifulSoup
import pprint
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()


class EPL(object):

    ACD_INPUTS = ('additional', 'advanced', 'input', 'required')
    ACD_OUTPUTS = ('output', )
    HTML_MAN_UNINTERESTING = (
        'Wiki',
        'Usage',
        'Command line arguments',
        'Input file format',
        'Output file format',
        'Data files',
        'References'
        'Diagnostic Error Messages',
        'Exit status',
        'See also',
        'Author(s)',
        'Target users'
    )

    TEST_CASE_DATA = (
        'Usage',
        'Input file format',
        'Output file format',
    )

    INPUT_TYPE_MAPPING = {
        'assembly': {},
        'codon': {},
        'cpdb': {},
        'datafile': {},
        'directory': {},
        'dirlist': {},
        'features': {},
        'filelist': {},
        'infile': {},
        'matrix': {},
        'matrixf': {},
        'obo': {'file_format': 'obo'},
        'pattern': {},
        'refseq': {},
        'regexp': {},
        'resource': {},
        'scop': {},
        'sequence': {'file_format': 'fasta'},
        'seqall': {'file_format': 'fasta'},
        'seqset': {'file_format': 'fasta'},
        'seqsetall': {'file_format': 'fasta'},
        'url': {'type': str},
        'variation': {'file_format': 'vcf'},
        'xml': {'file_format': 'xml'},
        'taxon': {'type': str},
        'text': {'type': str},
    }

    KNOWNTYPE_MAPPING = {
        'abi trace': 'abi'
    }

    SECTION_ORDERING = (
        'Standard (Mandatory) qualifiers',
        'Additional (Optional) qualifiers',
        'Advanced (Unprompted) qualifiers',
        # 'Associated qualifiers',  # Hidden, usually
        # 'General qualifiers',  # Also unused, -debug/-help/-version
    )

    def __init__(self, html_file, acd_file):
        self.html_file = html_file
        self.acd_file = acd_file
        with open(self.html_file, 'r') as handle:
            self.html = BeautifulSoup(handle.read())

        tree = ET.parse(self.acd_file)
        root = tree.getroot()

        tool = ET.Element("tool",
                          id='emboss_%s' % root.attrib['id'],
                          name=root.attrib['id'],
                          version='6.6.0')
        desc = ET.Element("description")
        desc.text = EPL.unbreak_strings(
            root.find('metadata').find('documentation').text)
        tool.append(desc)
        macros = ET.Element("macros")
        macros_import = ET.Element("import")
        macros_import.text = "macros.xml"
        macros.append(macros_import)
        tool.append(macros)
        tool.append(ET.Element("expand", macro="requirements"))
        tool.append(ET.Element("expand", macro="stdio"))
        version_command = ET.Element("version_command")
        version_command.text = "%s -version" % root.attrib['id']
        tool.append(version_command)
        command = ET.Element("command")
        command.text = "%s\n" % root.attrib['id']
        inputs = ET.Element("inputs")
        outputs = ET.Element("outputs")

        h2s = self.html.find_all('h2')
        param_table = None
        for i, h2 in enumerate(h2s):
            if h2.text.strip() == 'Command line arguments':
                param_table = self.parse_cli_table(
                    BeautifulSoup(
                        self.html_between_bounds(h2, h2s[i + 1])
                    )
                )

        if param_table is None:
            raise Exception("Couldn't parse out html command table")

        # Associates the "Associated (unprompted) qualifiers" with the
        # appropriate parameter.
        #
        # Sticking in self because we write HORRIBLY unclean code.
        self.reassociated = self._reassociate_associated_quals(param_table)

        for section in root.findall('section'):
            log.info('SECTION %s', pprint.pformat(section.attrib))

            # Get commands, inputs and outputs
            command_string, input_section, output_additions = self.build_section(section)
            # Command is just a string
            command.text += command_string
            # Inputs is wrapped up in a <section />
            if input_section is not None and not isinstance(input_section, list):
                inputs.append(input_section)
            # Outputs is a list however
            for output_addition in output_additions:
                outputs.append(output_addition)

        tool_help = ET.Element("help")
        testnode = ET.Element('tests')
        #test_cases = self.extract_test_cases(self.html)
        tool_help.text = ET.CDATA("\n" + self.clean_string(self.extract_useful_help(self.html)))

        tool.append(command)
        tool.append(inputs)
        tool.append(outputs)
        tool.append(testnode)
        citations = ET.Element("expand", macros="citation")

        tool.append(tool_help)
        tool.append(citations)

        print ET.tostring(tool)

    def parse_cli(self, command_line):
        # Ignore the first two parts, as they're useless (% and appname)
        parts = shlex.split(command_line)[2:]

        args = []
        kwargs = {}

        for i in range(len(parts)):
            if parts[i].startswith('-'):
                # Either a flag or a kwarg
                if i < len(parts) - 1:
                    # If we still have an item after this, treat naively
                    if parts[i + 1].startswith('-'):
                        # Flag case, dash and next is also a dash
                        kwargs[parts[i]] = "True"
                    else:
                        # Otherwise it's a flag with a string/argument after
                        kwargs[parts[i]] = parts[i + 1]
                        i += 1
                else:
                    kwargs[parts[i]] = "True"
        return args, kwargs

    @classmethod
    def html_between_bounds(cls, left_bound, right_bound):
        html = ''
        for tag in left_bound.next_siblings:
            if tag == right_bound:
                break
            else:
                html += str(tag)
        return html

    @classmethod
    def repeated_node_contents(cls, soup, node_tag, node_text):
        nodes = soup.find_all(node_tag)
        for i, node in enumerate(nodes):
            if node_text == node.text.strip():
                if i < len(nodes) - 1:
                    right_bound = nodes[i + 1]
                else:
                    right_bound = None

                return BeautifulSoup(EPL.html_between_bounds(node, right_bound))

    @classmethod
    def _reassociate_associated_quals(cls, clitable):
        keyed_params = {}
        for section in cls.SECTION_ORDERING:
            for param in clitable[section]:
                param_key = param['Qualifier']
                if '[-' in param_key:
                    param_key = param_key[1:param_key.index(']')]
                elif '-[' in param_key:
                    if ' ' in param_key:
                        param_key = param_key[0:param_key.index(' ')]
                else:
                    log.warn(param_key)

                if '[no]' in param_key:
                    param_key = param_key.replace('[no]', '')
                    # I have no idea what to do with this :(
                    param['boolean'] = '[no]'

                keyed_params[param_key] = param
                keyed_params[param_key]['section'] = section

        for key in clitable['Associated qualifiers']:
            real_param = key[1:key.rindex('"')]
            keyed_params[real_param]['extra'] = []

            for extra in clitable['Associated qualifiers'][key]:
                keyed_params[real_param]['extra'].append(extra)
        return keyed_params

    @classmethod
    def parse_cli_table(cls, soup):
        table = soup.find('table', attrs={'bgcolor': '#ccccff'})
        rows = table.find_all('tr')
        data = {}
        category = "None"
        subcategory = None
        col_types = []
        for i, row in enumerate(rows):

            if i == 0:
                cols = row.find_all('th')
                col_types = [x.text.strip() for x in cols]
            else:
                cols = row.find_all('td')
                if len(cols) == 1:
                    subcategory = cols[0].text.strip()
                    if isinstance(data[category], list):
                        data[category] = {}
                    data[category][subcategory] = []
                elif len(cols) == 0:
                    category = row.find_all('th')[0].text.strip()
                    data[category] = []
                    subcategory = None
                else:
                    col_vals = [e.text.strip() for e in cols]
                    col_dict = dict(zip(col_types, col_vals))
                    if subcategory is not None:
                        data[category][subcategory].append(col_dict)
                    else:
                        data[category].append(col_dict)
        return data

    def extract_useful_help(self, html):
        p = subprocess.Popen(['pandoc', '-f', 'html', '-t', 'rst'],
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        good_html = ""
        h2s = html.find_all('h2')
        for i, h2 in enumerate(h2s):
            # Get bounds for fetching HTML
            left_bound = h2
            if i < len(h2s) - 1:
                right_bound = h2s[i + 1]
            else:
                right_bound = None

            # Ignore uninteresting items
            if h2.text.strip() in self.HTML_MAN_UNINTERESTING:
                continue

            good_html += str(h2)
            good_html += str(self.html_between_bounds(left_bound, right_bound))

        # TODO References should be transferred to citations automatically
        out, err = p.communicate(input=good_html)
        return out

    @classmethod
    def unbreak_strings(cls, data):
        return re.sub('\s*[\r\n]\s*', ' ', data)

    def __parse_test_usage(self, html):
        code_statements = BeautifulSoup(html).find_all('pre')
        for code in code_statements:
            for line in code.text.strip().split('\n'):
                if line.startswith('%'):
                    args, kwargs = self.parse_cli(line)
        pass

    def __parse_test_inputs(self, html):
        pass

    def __parse_test_outputs(self, html):
        pass

    def extract_test_cases(self, html):
        h2s = html.find_all('h2')
        test_case = {}
        for i, h2 in enumerate(h2s):
            # Get bounds for fetching HTML
            left_bound = h2
            if i < len(h2s) - 1:
                right_bound = h2s[i + 1]
            else:
                right_bound = None

            # Ignore uninteresting items
            if h2.text.strip() in self.TEST_CASE_DATA:
                test_case[h2.text.strip()] = self.html_between_bounds(left_bound, right_bound)

        self.__parse_test_usage(test_case['Usage'])
        self.__parse_test_inputs(test_case['Input file format'])
        self.__parse_test_outputs(test_case['Output file format'])

    @classmethod
    def valid_xml_char_ordinal(cls, c):
        # http://stackoverflow.com/questions/8733233/filtering-out-certain-bytes-in-python
        codepoint = ord(c)
        # conditions ordered by presumed frequency
        return (
            0x20 <= codepoint <= 0xD7FF or
            codepoint in (0x9, 0xA, 0xD) or
            0xE000 <= codepoint <= 0xFFFD or
            0x10000 <= codepoint <= 0x10FFFF
        )

    @classmethod
    def clean_string(cls, data):
        return ''.join(c for c in data if EPL.valid_xml_char_ordinal(c)).decode('utf8')

    def __default_kwargs(self, child, assoc=None):
        kwargs = {
            'name': child.attrib['name'],
            'label': child.attrib['name'],
            'help': '',
        }

        if child.find('information') is not None:
            kwargs['label'] = self.unbreak_strings(child.find('information').text)

        if child.attrib['type'] in self.INPUT_TYPE_MAPPING.keys():
            kwargs['type'] = 'data'
        elif child.attrib['type'] in ('integer', 'float'):
            kwargs['type'] = child.attrib['type']
        elif child.attrib['type'] in ('boolean', 'toggle'):
            kwargs['type'] = 'boolean'
            if 'default' in child.attrib and child.attrib['default'] == 'Y':
                kwargs['selected'] = 'True'
        elif child.attrib['type'] in ('string',):
            kwargs['type'] = 'text'

        if 'Description' in assoc:
            if assoc['Description'] != kwargs['label']:
                kwargs['help'] = assoc['Description']

        if child.find('help') is not None:
            kwargs['help'] = self.unbreak_strings(child.find('help').text)

        # Add the (-argname)
        help_arg = ' (-%s)' % child.attrib['name']
        if len(kwargs['help']) == 0:
            help_arg = help_arg.strip()
        kwargs['help'] += help_arg
        #elif child.attrib['type'] == 'list':
        #elif child.attrib['type'] == 'selection':
        #elif child.attrib['type'] == 'range':
        return kwargs

    def __file_input(self, child, assoc=None):
        #<param label="Genome" name="positional_1" type="data" format="fasta"/>
        kwargs = self.__default_kwargs(child, assoc=assoc)

        if child.find('knowntype') is not None:
            kwargs['file_format'] = self.KNOWNTYPE_MAPPING.get(child.find('knowntype').text, 'UNKNOWN')

        # Unhandled: additional, default, nullok, relations, ... ?

        kwargs.update(self.INPUT_TYPE_MAPPING[child.attrib['type']])
        parameter = ET.Element('param', **kwargs)
        cmd = "-{0} ${0}".format(child.attrib['name'])
        return parameter, cmd

    def __numeric_input(self, child, assoc=None):
        num_min = None
        num_max = None
        if child.find('minimum') is not None:
            try:
                num_min = float(child.find('minimum').text)
            except Exception:
                pass
        if child.find('maximum') is not None:
            try:
                num_max = float(child.find('maximum').text)
            except Exception:
                pass

        try:
            num_default = float(child.find('default').text)
        except Exception:
            num_default = 0

        if child.attrib['type'] == 'integer':
            if num_min is not None:
                num_min = int(num_min)
            if num_max is not None:
                num_max = int(num_max)
            num_default = int(num_default)

        kwargs = self.__default_kwargs(child, assoc=assoc)
        kwargs.update({
            'default': str(num_default),
        })

        if num_min is not None:
            kwargs['minimum'] = str(num_min)

        if num_max is not None:
            kwargs['maximum'] = str(num_max)

        parameter = ET.Element('param', **kwargs)
        cmd = "-{0} ${0}".format(child.attrib['name'])
        return parameter, cmd

    def __text_input(self, child, assoc=None):
        kwargs = self.__default_kwargs(child, assoc=assoc)
        parameter = ET.Element('param', **kwargs)
        cmd = '-{0} "${0}"'.format(child.attrib['name'])
        return parameter, cmd

    def __boolean_input(self, child, assoc=None):
        kwargs = self.__default_kwargs(child, assoc=assoc)
        kwargs.update({
            'truevalue': '-%s' % child.attrib['name'],
            'falsevalue': '',
        })

        if 'boolean' in assoc and assoc['boolean'] == '[no]':
            kwargs['truevalue'] = '-%s' % child.attrib['name']
            kwargs['falsevalue'] = '-no%s' % child.attrib['name']

        parameter = ET.Element('param', **kwargs)
        cmd = "${0}".format(child.attrib['name'])
        return parameter, cmd

    def __list_input(self, child, assoc=None):
        kwargs = self.__default_kwargs(child, assoc=assoc)
        if child.find('delimiter') is not None:
            delimiter = child.find('delimiter').text
        else:
            delimiter = ';'

        if child.find('codedelimiter') is not None:
            codedelimiter = child.find('codedelimiter').text
        else:
            codedelimiter = ':'

        kv = {}
        for x in child.find('values').text.split(delimiter):
            if len(x.strip()) > 0:
                try:
                    key, value = x.split(codedelimiter)
                except:
                    tmp = x.split()
                    key = tmp[0]
                    value = ' '.join(tmp[1:])

                kv[key.strip()] = value.strip()

        kwargs['type'] = "select"
        #kwargs['choices'] = [kv[k] for k in kv.keys()]
        parameter = ET.Element('param', **kwargs)
        for key in kv:
            tmpkw = {'value': key}

            if child.find('default') is not None and child.find('default').text == key:
                tmpkw['checked'] = "True"

            option = ET.Element('option', **tmpkw)
            option.text = kv[key]
            parameter.append(option)
        cmd = '-{0} "${0}"'.format(child.attrib['name'])
        return parameter, cmd

    def __selection_input(self, child, assoc=None):
        kwargs = self.__default_kwargs(child, assoc=assoc)
        if child.find('delimiter') is not None:
            delimiter = child.find('delimiter').text
        else:
            delimiter = ';'

        kv = {}
        for i, x in enumerate(child.find('values').text.split(delimiter)):
            if len(x.strip()) > 0:
                kv[str(i)] = x.strip()

        kwargs['type'] = str
        #kwargs['choices'] = [kv[k] for k in kv.keys()]
        parameter = ET.Element('param', **kwargs)
        return parameter, None

    def handle_parameter(self, child, assoc=None):
        if child.attrib['type'] in self.INPUT_TYPE_MAPPING.keys():
            return self.__file_input(child, assoc=assoc)
        elif child.attrib['type'] in ('integer', 'float'):
            return self.__numeric_input(child, assoc=assoc)
        elif child.attrib['type'] in ('string'):
            return self.__text_input(child, assoc=assoc)
        elif child.attrib['type'] in ('boolean', 'toggle'):
            return self.__boolean_input(child, assoc=assoc)
        elif child.attrib['type'] == 'list':
            return self.__list_input(child, assoc=assoc)
        elif child.attrib['type'] == 'selection':
            return self.__selection_input(child, assoc=assoc)
        else:
            return ET.Comment("Unhanlded %s" % child.attrib['type']), None

    def _generate_extra_inputs(self, extra, parent_name=None):
        """
        extra is something like::

            {
                'Allowed values': 'Any string',
                'Default': '',
                'Description': 'Output seq format',
                'Qualifier': '-osformat2-osformat_outseq',
                'Type': 'string'
            }

        ``parent_name`` is used to pass the command line name of the parent
        output file
        """
        command = ''
        quals = extra['Qualifier'].split('-')[1:]
        if 'format' in quals[0]:
            # TODO: support more than just fasta
            command = '-%s "fasta"' % quals[0]
        elif 'extension' in quals[0]:
            command = '-%s "dat"' % quals[0]
        elif 'directory' in quals[0]:
            command = '-%s $(dirname $%s)' % (quals[0], parent_name)
        elif 'dbname' in quals[0]:
            pass
        elif 'name' in quals[0]:
            # This mega sucks. Sucks SO much.
            command = '-%s $(basename $%s dat)' % (quals[0], parent_name)
        elif 'single' in quals[0]:
            pass
        elif 'ufo' in quals[0]:
            pass
        else:
            log.info('Unhandled extra %s', quals[0])
        return command

    def _generate_outputs(self, section):
        assoc = self.reassociated.get('-' + section.attrib['name'], {})

        return ("\n## %s\n" % section.attrib['name'], [])

    def _generate_inputs(self, section):
        inputs = []

        assoc = self.reassociated.get('-' + section.attrib['name'], {})
        input_param, input_cmd = self.handle_parameter(section, assoc=assoc)
        inputs.append(input_param)
        if input_cmd is None:
            input_cmd = ''

        for extra in assoc.get('extra', []):
            id = section.attrib['name']
            extra_cmd = self._generate_extra_inputs(extra, parent_name=id)
            input_cmd += extra_cmd + "\n"

        return (input_cmd, inputs)

    def build_section(self, acd_xml_section):
        if acd_xml_section.tag == 'parameter':
            command_string = ''
            inputs = []
            outputs = []

            commands, input_additions = self._generate_inputs(acd_xml_section)
            command_string += ''.join(commands)
            for ia in input_additions:
                inputs.append(ia)

            commands, output_additions = self._generate_outputs(acd_xml_section)
            command_string += ''.join(commands)
            for oa in output_additions:
                outputs.append(oa)

            return command_string, inputs, outputs
        elif acd_xml_section.tag == 'section':
            # Get section information
            metadata = acd_xml_section.find('metadata')
            section_title = "None"
            if metadata is not None:
                section_title = metadata.find('information').text

            # Complete command for entire section
            section_cmd = ""
            # Some inputs
            inputs = ET.Element('section', name=acd_xml_section.attrib['id'],
                                title=section_title)
            # Some output datasets, maybe
            outputs = []

            # For all the kids
            for child in acd_xml_section.getchildren():
                # Skip metadata, we only want parameters and sections
                if child.tag == 'metadata':
                    continue

                command_string, input_extra, output = self.build_section(child)

                # Append command string
                section_cmd += command_string
                # Inputs are semi-special compared to outputs
                for i in input_extra:
                    inputs.append(i)
                # Outputs are just all lumped together
                outputs.extend(output)
            # Return all of them
            return section_cmd, inputs, outputs
        else:
            log.warn("Unhanlded section %s", acd_xml_section.tag)
            return None, None, None


if __name__ == '__main__':
    epl = EPL(sys.argv[1], sys.argv[2])
