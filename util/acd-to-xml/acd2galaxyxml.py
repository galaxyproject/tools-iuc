#!/usr/bin/env python
from lxml import etree as ET
import shlex
import subprocess
import argparse
import re
import sys
from bs4 import BeautifulSoup
import logging
import pprint
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

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


def unbreak_strings(data):
    return re.sub('\s*[\r\n]\s*', ' ', data)


def html_between_bounds(left_bound, right_bound):
    html = ''
    for tag in left_bound.next_siblings:
        if tag == right_bound:
            break
        else:
            html += str(tag)
    return html


def __parse_cli(command_line):
    # Ignore the first two parts, as they're useless (% and appname)
    parts = shlex.split(command_line)[2:]

    args = []
    kwargs = {}

    for i in range(len(parts)):
        if parts[i].startswith('-'):
            # Either a flag or a kwarg
            log.debug('%s -> %s', i, parts[i])

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


def __parse_test_usage(html):
    code_statements = BeautifulSoup(html).find_all('pre')
    for code in code_statements:
        for line in code.text.strip().split('\n'):
            if line.startswith('%'):
                args, kwargs = __parse_cli(line)
                log.debug(pprint.pformat(args))
                log.debug(pprint.pformat(kwargs))
    pass

def __parse_test_inputs(html):
    pass

def __parse_test_outputs(html):
    pass


def extract_test_cases(html):
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
        if h2.text.strip() in TEST_CASE_DATA:
            test_case[h2.text.strip()] = html_between_bounds(left_bound, right_bound)

    __parse_test_usage(test_case['Usage'])
    __parse_test_inputs(test_case['Input file format'])
    __parse_test_outputs(test_case['Output file format'])


def extract_useful_help(html):
    p = subprocess.Popen(['pandoc', '-f', 'html', '-t', 'rst'],
                         stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE
                         )
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
        if h2.text.strip() in HTML_MAN_UNINTERESTING:
            continue

        good_html += str(h2)
        good_html += str(html_between_bounds(left_bound, right_bound))

    # TODO References should be transferred to citations automatically
    out, err = p.communicate(input=good_html)
    return out


def valid_xml_char_ordinal(c):
    # http://stackoverflow.com/questions/8733233/filtering-out-certain-bytes-in-python
    codepoint = ord(c)
    # conditions ordered by presumed frequency
    return (
        0x20 <= codepoint <= 0xD7FF or
        codepoint in (0x9, 0xA, 0xD) or
        0xE000 <= codepoint <= 0xFFFD or
        0x10000 <= codepoint <= 0x10FFFF
    )


def clean_string(data):
    return ''.join(c for c in data if valid_xml_char_ordinal(c)).decode('utf8')


tree = ET.parse(sys.argv[1])
root = tree.getroot()

tool = ET.Element("tool",
                  id='emboss_%s' % root.attrib['id'],
                  name=root.attrib['id'],
                  version='6.6.0'
                  )
desc = ET.Element("description")
desc.text = unbreak_strings(root.find('metadata').find('documentation').text)
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


def __default_kwargs(child):
    kwargs = {
        'name': child.attrib['name'],
        'label': child.attrib['name'],
    }

    if child.find('information') is not None:
        kwargs['label'] = unbreak_strings(child.find('information').text)

    if child.attrib['type'] in INPUT_TYPE_MAPPING.keys():
        kwargs['type'] = 'data'
    elif child.attrib['type'] in ('integer', 'float'):
        kwargs['type'] = child.attrib['type']
    elif child.attrib['type'] in ('boolean', 'toggle'):
        kwargs['type'] = 'boolean'
        if 'default' in child.attrib and child.attrib['default'] == 'Y':
            kwargs['selected'] = 'True'
    elif child.attrib['type'] in ('string',):
        kwargs['type'] = 'text'
    #elif child.attrib['type'] == 'list':
    #elif child.attrib['type'] == 'selection':
    #elif child.attrib['type'] == 'range':
    return kwargs

def __file_input(child):
    #<param label="Genome" name="positional_1" type="data" format="fasta"/>
    kwargs = __default_kwargs(child)

    if child.find('knowntype') is not None:
        kwargs['file_format'] = child.find('knowntype').text

    if child.find('help') is not None:
        kwargs['help'] = unbreak_strings(child.find('help').text)
    # Unhandled: additional, default, nullok, relations, ... ?

    kwargs.update(INPUT_TYPE_MAPPING[child.attrib['type']])
    parameter = ET.Element('parameter', **kwargs)
    cmd = "-{0} ${0}".format(child.attrib['name'])
    return parameter, cmd


def __numeric_input(child):
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

    kwargs = __default_kwargs(child)
    kwargs.update({
        'default': str(num_default),
    })

    if num_min is not None:
        kwargs['minimum'] = str(num_min)

    if num_max is not None:
        kwargs['maximum'] = str(num_max)

    parameter = ET.Element('parameter', **kwargs)
    cmd = "-{0} ${0}".format(child.attrib['name'])
    return parameter, cmd

def __text_input(child):
    kwargs = __default_kwargs(child)
    parameter = ET.Element('parameter', **kwargs)
    cmd = '-{0} "${0}"'.format(child.attrib['name'])
    return parameter, cmd

def __boolean_input(child):
    kwargs = __default_kwargs(child)
    kwargs.update({
        'truevalue': '-%s' % child.attrib['name'],
        'falsevalue': '',
    })
    parameter = ET.Element('parameter', **kwargs)
    cmd = "${0}".format(child.attrib['name'])
    return parameter, cmd


def __list_input(child):
    kwargs = {}
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
    parameter = ET.Element('parameter', **kwargs)
    for key in kv:
        tmpkw = {'value': key}

        if child.find('default') is not None and child.find('default').text == key:
            tmpkw['checked'] = "True"

        option = ET.Element('option', **tmpkw)
        option.text = kv[key]
        parameter.append(option)
    cmd = '-{0} "${0}"'.format(child.attrib['name'])
    return parameter, cmd


def __selection_input(child):
    kwargs = {'name': child.attrib['name']}
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
    parameter = ET.Element('parameter', **kwargs)
    return parameter, None


def handle_parameter(child):
    if child.attrib['type'] in INPUT_TYPE_MAPPING.keys():
        return __file_input(child)
    elif child.attrib['type'] in ('integer', 'float'):
        return __numeric_input(child)
    elif child.attrib['type'] in ('string'):
        return __text_input(child)
    elif child.attrib['type'] in ('boolean', 'toggle'):
        return __boolean_input(child)
    elif child.attrib['type'] == 'list':
        return __list_input(child)
    elif child.attrib['type'] == 'selection':
        return __selection_input(child)
    else:
        return ET.Comment("Unhanlded %s" % child.attrib['type']), None


def build_section(acd_xml_section):
    if acd_xml_section.tag == 'parameter':
        parameter, command_string = handle_parameter(acd_xml_section)
        if command_string is None:
            command_string = "## $%s\n" % acd_xml_section.attrib['name']
        return command_string + "\n", parameter
    elif acd_xml_section.tag == 'section':
        # Get section information
        metadata = acd_xml_section.find('metadata')
        if metadata is not None:
            section_title = metadata.find('information').text
        else:
            section_title = "None"  # acd_xml_section.attrib['id']

        section = ET.Element('section',
                             name=acd_xml_section.attrib['id'],
                             title=section_title)
        section_cmd = ""
        for child in acd_xml_section.getchildren():
            if child.tag == 'metadata':
                continue
            else:
                command_string, parameter = build_section(child)
                if parameter is not None:
                    section_cmd += command_string
                    section.append(parameter)
        return section_cmd, section
    else:
        return None, None
        log.warn("Unhanlded section %s", acd_xml_section.tag)


command = ET.Element("command")
command.text = "%s\n" % root.attrib['id']
inputs = ET.Element("inputs")
outputs = ET.Element("outputs")
for section in root.findall('section'):
    log.info(section.attrib)
    # Top level attributes
    if section.attrib['id'] in ACD_INPUTS:
        command_string, input_additions = build_section(section)
        command.text += command_string + "\n"
        inputs.append(input_additions)
    else:
        for parameter in section.findall('parameter'):
            if parameter.attrib['type'] in ('toggle', ):
                log.warn("Unhandled output toggle")
            else:
                kwargs = {
                    'name': parameter.attrib['name']
                }
                if parameter.find('aformat') is not None:
                    kwargs['format'] = parameter.find('aformat').text
                elif parameter.find('extension') is not None:
                    kwargs['format'] = parameter.find('extension').text

                output_parameter = ET.Element('data', **kwargs)
                outputs.append(output_parameter)
                command.text += '-%s $%s\n' % (parameter.attrib['name'], parameter.attrib['name'])


tool_help = ET.Element("help")
with open(sys.argv[2], 'r') as handle:
    tool_help_unparsed_html = BeautifulSoup(handle.read())

testnode = ET.Element('tests')
test_cases = extract_test_cases(tool_help_unparsed_html)
tool_help.text = ET.CDATA("\n" + clean_string(extract_useful_help(tool_help_unparsed_html)))

tool.append(command)
tool.append(inputs)
tool.append(outputs)
tool.append(testnode)
citations = ET.Element("expand", macros="citation")

tool.append(tool_help)
tool.append(citations)


print ET.tostring(tool)
