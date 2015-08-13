import re
import logging
import fileinput
import xml.etree.ElementTree as etree
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

PARAMETER_REGEX = '(?P<option_param>[^ ]*[<\[][^ ]*[>\]][^ ]*)'
SPACES_THEN_HELP = '\s{2,}(?P<help_text>.*)$'

SINGLE_OPTION_PARAM = re.compile('^-(?P<option_short_a>.),\s+--(?P<option_long_a>[^ ]+) ' + PARAMETER_REGEX + SPACES_THEN_HELP)
SINGLE_OPTION_FLAG = re.compile('^-(?P<option_short_a>.),\s+--(?P<option_long_a>[^ ]+)' + SPACES_THEN_HELP)
ONLY_LONG_FLAG = re.compile('^--(?P<option_long_a>[^ ]+)' + SPACES_THEN_HELP)
ONLY_LONG_PARAM = re.compile('^--(?P<option_long_a>[^ ]+) ' + PARAMETER_REGEX + SPACES_THEN_HELP)

DUAL_OPTION_PARAM = re.compile('^-(?P<option_short_a>.)/(?P<option_short_b>.),\s+--(?P<option_long_a>[^ ]+)/--(?P<option_long_b>[^ ]+) ' + PARAMETER_REGEX + SPACES_THEN_HELP)
DUAL_OPTION_FLAG = re.compile('^-(?P<option_short_a>.)/(?P<option_short_b>.),\s+--(?P<option_long_a>[^ ]+)/--(?P<option_long_b>[^ ]+)' + SPACES_THEN_HELP)

DEFAULT_VALUE = re.compile('^(?P<help_text>.*) \[(?P<default_value>[^\]]*)\]$')

param_ds = {}
section = None
in_about = False
for line in fileinput.input():
    if line.startswith('About:   '):
        about_line = line.replace('About:   ', '').strip()
        in_about = True
    elif line.startswith('Usage:  '):
        usage_line = line.replace('Usage:   ', '').strip()
        in_about = False
    elif len(line.strip()) == 0:
        pass
    elif line.startswith(' '):
        if in_about:
            about_line += ' ' + line.strip()
        else:
            # Option or continuation of option text
            if section not in param_ds:
                param_ds[section] = []

            # Continuation of option text
            if line.startswith('                   '):
                param_ds[section][-1] += ' ' + line.strip()
            else:
                param_ds[section].append(line.strip())
    elif line.startswith('Examples:'):
        break
    else:
        section = line.replace('options:', '').replace('Options:', '').strip()
        if len(section) == 0:
            section = 'Default'
        log.debug("Saw section %s" % section)

import pprint
log.debug(pprint.pformat(param_ds))

tool = etree.Element("tool",
                     id="bcftools_@EXECUTABLE@",
                     name="bcftools @EXECUTABLE@",
                     version="@VERSION@.0")
desc = etree.SubElement(tool, 'description')
desc.text = about_line.split('. ')[0]
macros = etree.SubElement(tool, 'macros')
tool_id = usage_line.split()[1]
etree.SubElement(macros, 'token', name='@EXECUTABLE@').text = tool_id
etree.SubElement(macros, 'import').text = 'macros.xml'
etree.SubElement(tool, 'expand', macro="requirements")
etree.SubElement(tool, 'expand', macro="stdio")
etree.SubElement(tool, 'expand', macro="version_command")
command = etree.SubElement(tool, 'command')
command.text = 'bcftools @EXECUTABLE@\n'
inputs = etree.SubElement(tool, 'inputs')
outputs = etree.SubElement(tool, 'outputs')
tests = etree.SubElement(tool, 'tests')
help = etree.SubElement(tool, 'help')
help.text = about_line
etree.SubElement(tool, 'expand', macro='citations')

def parse_param(help_text):
    m = (re.match(SINGLE_OPTION_FLAG, help_text) or re.match(ONLY_LONG_FLAG, help_text))
    if m:
        return [{
            'type': 'flag',
            'long': m.group('option_long_a'),
            'help': m.group('help_text')
        }]

    m = (re.match(SINGLE_OPTION_PARAM, help_text) or re.match(ONLY_LONG_PARAM, help_text))
    if m:
        return [{
            'type': 'param',
            'long': m.group('option_long_a'),
            'help': m.group('help_text'),
            'param': m.group('option_param'),
        }]

    m = re.match(DUAL_OPTION_FLAG, help_text)
    if m:
        if 'exclude-' + m.group('option_long_a') == m.group('option_long_b'):
            return [{
                'type': 'flag/exclude',
                'truevalue': m.group('option_long_a'),
                'falsevalue': m.group('option_long_b'),
                'long': m.group('option_long_a') + '/' + m.group('option_long_b'),
                'help': m.group('help_text'),
            }]
        else:
            return [{
                'type': 'flag',
                'long': m.group('option_long_a'),
                'help': m.group('help_text')
            },{
                'type': 'flag',
                'long': m.group('option_long_b'),
                'help': m.group('help_text')
            }]

    m = re.match(DUAL_OPTION_PARAM, help_text)
    if m:
        if 'exclude-' + m.group('option_long_a') == m.group('option_long_b'):
            return [{
                'type': 'param/exclude',
                'long': m.group('option_long_a') + '_' + m.group('option_long_b'),
                'help': m.group('help_text'),
                'param': m.group('option_param'),
                'flag_choices': (m.group('option_long_a'), m.group('option_long_b')),
            }]
        else:
            return [{
                'type': 'param',
                'long': m.group('option_long_a'),
                'param': m.group('option_param'),
                'help': m.group('help_text')
            },{
                'type': 'param',
                'long': m.group('option_long_b'),
                'param': m.group('option_param'),
                'help': m.group('help_text')
            }]

    log.warn("Unhandled parameter: %s", help_text)



for section_name in sorted(param_ds.keys()):
    safe_sec_name = 'sec_' + section_name.replace(' ', '_').replace('/', '_').lower()
    section = etree.SubElement(inputs, 'section', name=safe_sec_name,
                               title=section_name + ' Options', expanded="true")

    command.text += "\n## %s section\n" % section_name
    for command_string in param_ds[section_name]:
        for parsed_command in parse_param(command_string):
            if parsed_command['type'] == 'flag':
                pkw = {
                    'name': parsed_command['long'].replace('-', '_'),
                    'label': parsed_command['long'].replace('-', ' ').title(),
                    'type': 'boolean',
                    'truevalue': '--' + parsed_command['long'],
                    'falsevalue': '',
                    'help': parsed_command['help'],
                }
                param = etree.SubElement(section, 'param', **pkw)
                # These we really don't have a choice about applying/not,
                # without turning them into ternaries.
                # However, it looks like bcftools is well behaved in this respect.
                #command.text += "#if $%s.%s:\n" % (safe_sec_name, pkw['name'])
                command.text += "${%s.%s}\n" % (safe_sec_name, pkw['name'])
                #command.text += "#end if\n"

            elif parsed_command['type'] == 'flag/exclude':
                pkw = {
                    'name': parsed_command['long'].replace('-', '_').replace('/', '_'),
                    'label': parsed_command['long'].replace('-', ' ').title(),
                    'type': 'boolean',
                    'truevalue': '--' + parsed_command['truevalue'],
                    'falsevalue': '--' + parsed_command['falsevalue'],
                    'help': parsed_command['help'],
                }
                param = etree.SubElement(section, 'param', **pkw)
                #command.text += "#if $%s.%s:\n" % (safe_sec_name, pkw['name'])
                command.text += "${%s.%s}\n" % (safe_sec_name, pkw['name'])
                #command.text += "#end if\n"

            elif parsed_command['type'].startswith('param'):
                pkw = {
                    'name': parsed_command['long'].replace('-', '_'),
                    'label': parsed_command['long'].replace('-', ' ').title(),
                    'help': parsed_command['help'],
                }

                # TODO: there are [default]s at the end of the help text.
                # Need to strip those...
                if '<file>' in parsed_command['param']:
                    pkw['type'] = 'data'
                    pkw['format'] = 'data'
                    pkw['optional'] = 'True'
                elif '<int>' in parsed_command['param']:
                    pkw['type'] = 'integer'
                    pkw['optional'] = 'True'
                    m = re.match(DEFAULT_VALUE, pkw['help'])
                    if m:
                        pkw['help'] = m.group('help_text')
                        pkw['default'] = m.group('default_value')
                elif '<str>' in parsed_command['param'] or '<name>' in parsed_command['param'] or '<string>' in parsed_command['param']:
                    pkw['type'] = 'string'
                    pkw['optional'] = 'True'
                elif '<float>' in parsed_command['param']:
                    pkw['type'] = 'float'
                    pkw['optional'] = 'True'
                    m = re.match(DEFAULT_VALUE, pkw['help'])
                    if m:
                        pkw['help'] = m.group('help_text')
                        pkw['default'] = m.group('default_value')
                elif '<' in parsed_command['param'] and '|' in parsed_command['param']:
                    parsed_command['select'] = True
                    help_text = pkw['help']
                    default = None
                    # Get the default if it exists
                    m = re.match(DEFAULT_VALUE, pkw['help'])
                    if m:
                        help_text = m.group('help_text')
                        default = m.group('default_value')


                    select = etree.SubElement(section, 'param',
                                              type='select',
                                              name='select_' + pkw['name'])

                    # Disabled until https://github.com/galaxyproject/galaxy/issues/599
                    #if default is None:
                        #pkw['optional'] = 'True'
                    if default is None:
                        etree.SubElement(select, 'option', value='__none__', selected="True").text = "No selection"


                    for kv in help_text.split(', '):
                        kvd = kv.split(': ')
                        k = kvd[0]
                        v = ': '.join(kvd[1:])

                        sokwd = {'value': k}
                        if k == default:
                            sokwd['selected'] = 'True'

                        etree.SubElement(select, 'option', **sokwd).text = v

                    command.text += "#if str($%s.%s) != \"__none__\":\n" % (safe_sec_name, 'select_' + pkw['name'])
                    command.text += "  --%s \"${%s.%s}\"\n" % (parsed_command['long'], safe_sec_name, 'select_' + pkw['name'])
                    command.text += "#end if\n"

                    if pkw['name'] == 'output_type':
                        out = etree.SubElement(outputs, 'data', format="vcf", name="output_file")
                        cf = etree.SubElement(out, 'change_format')
                        etree.SubElement(cf, 'when', input=safe_sec_name + "|select_output_type", value="b", format="bcf_bgz")
                        etree.SubElement(cf, 'when', input=safe_sec_name + "|select_output_type", value="u", format="bcf")
                        etree.SubElement(cf, 'when', input=safe_sec_name + "|select_output_type", value="z", format="vcf_bgz")
                        etree.SubElement(cf, 'when', input=safe_sec_name + "|select_output_type", value="v", format="vcf")

                elif '<region>' in parsed_command['param']:
                    repeat = etree.SubElement(section, 'repeat',
                                              name=pkw['name'] + '_repeat',
                                              title=pkw['label'])
                    etree.SubElement(repeat, 'param', **pkw)
                    command_id = 'values_%s_%s' % (safe_sec_name, pkw['name'])
                    command.text += """#set %s = '","'.join([str($value) for $value in $%s.%s])\n""" % (command_id, safe_sec_name, pkw['name'] + "_repeat")
                    command.text += "#if $%s:\n" % (command_id)
                    command.text += "  --%s \"${%s}\"\n" % (parsed_command['long'], command_id)
                    command.text += "#end if\n"

                elif '<list>' in parsed_command['param']:
                    repeat = etree.SubElement(section, 'repeat',
                                              name=pkw['name'] + '_repeat',
                                              title=pkw['label'])
                    etree.SubElement(repeat, 'param', **pkw)

                    command_id = 'values_%s_%s' % (safe_sec_name, pkw['name'])
                    command.text += """#set %s = '","'.join([str($value) for $value in $%s.%s])\n""" % (command_id, safe_sec_name, pkw['name'] + "_repeat")
                    command.text += "#if $%s:\n" % (command_id)
                    if '[^]' in parsed_command['param']:
                        mpkw = {
                            'name': 'invert_' + pkw['name'],
                            'label': 'Invert ' + pkw['label'],
                            'help': 'inverts the query/filtering applied by ' + pkw['label'],
                            'type': 'boolean',
                            'truevalue': '^',
                            'falsevalue': '',
                        }
                        etree.SubElement(section, 'param', **mpkw)
                        command.text += "  --%s \"${%s.%s}${%s}\"\n" % (parsed_command['long'], safe_sec_name, mpkw['name'], command_id)
                    else:
                        command.text += "  --%s \"${%s}\"\n" % (parsed_command['long'], command_id)
                    command.text += "#end if\n"
                else:
                    pkw['__TODO__'] = 'TODO'
                    log.warn("Unknown type: %s", parsed_command['param'])

                # Lists are repeats and handled elsewhere.
                if '<list>' not in parsed_command['param'] and 'select' not in parsed_command and '<region>' not in parsed_command['param']:
                    etree.SubElement(section, 'param', **pkw)

                    if '[^]' in parsed_command['param'] or parsed_command['type'] == 'param/exclude':
                        mpkw = {
                            'name': 'invert_' + pkw['name'],
                            'label': 'Invert ' + pkw['label'],
                            'help': 'inverts the query/filtering applied by ' + pkw['label'],
                            'type': 'boolean',
                            'truevalue': '^',
                            'falsevalue': '',
                        }

                        command.text += "#if $%s.%s:\n" % (safe_sec_name, pkw['name'])
                        if parsed_command['type'] == 'param/exclude':
                            # If it's an exclude instead of a carat, we've
                            # reduced the two flags into a single real option
                            # plus the inverter. This will be used in the template instead of the '^' or flag value.
                            mpkw.update({
                                'truevalue': parsed_command['flag_choices'][0],
                                'falsevalue': parsed_command['flag_choices'][1],
                            })
                            command.text += "  --{%s.%s} \"${%s.%s}\"\n" % (
                                safe_sec_name,
                                mpkw['name'],
                                safe_sec_name,
                                pkw['name']
                            )
                        else:
                            command.text += "  --%s \"${%s.%s}${%s.%s}\"\n" % (
                                parsed_command['long'],
                                safe_sec_name,
                                mpkw['name'],
                                safe_sec_name,
                                pkw['name']
                            )
                        command.text += "#end if\n"

                        etree.SubElement(section, 'param', **mpkw)

                    else:
                        command.text += "#if $%s.%s:\n" % (safe_sec_name, pkw['name'])
                        command.text += "  --%s \"${%s.%s}\"\n" % (
                            parsed_command['long'],
                            safe_sec_name,
                            pkw['name']
                        )
                        command.text += "#end if\n"
            else:
                print parsed_command

# Rest of command parsing
#         bcftools   annotate    [options]   <in.vcf.gz>
#         bcftools   call        [options]   <in.vcf.gz>
#         bcftools   filter      [options]   <in.vcf.gz>
#         bcftools   norm        [options]   <in.vcf.gz>
#         bcftools   reheader    [OPTIONS]   <in.vcf.gz>
#         bcftools   roh         [options]   <in.vcf.gz>
#         bcftools   view        [options]   <in.vcf.gz>   [region1              [...]]

#         bcftools   consensus   [OPTIONS]   <file.vcf>

#         bcftools   gtcheck     [options]   [-g           <genotypes.vcf.gz>]   <query.vcf.gz>

#         bcftools   isec        [options]   <A.vcf.gz>    <B.vcf.gz>            [...]
#         bcftools   merge       [options]   <A.vcf.gz>    <B.vcf.gz>            [...]

#         bcftools   concat      [options]   <A.vcf.gz>    [<B.vcf.gz>           [...]]
#         bcftools   query       [options]   <A.vcf.gz>    [<B.vcf.gz>           [...]]

#         bcftools   stats       [options]   <A.vcf.gz>    [<B.vcf.gz>]

command.text += "\n## Primary Input/Outputs\n\n"
elements = []
if '-g' in usage_line:
    # Special case, only happens once
    elements = [
        etree.Element('param', type="data", format="vcf,bcf", name="genotypes_file", label="Genotypes VCF/BCF Data", optional="True"),
        etree.Element('param', type="data", format="vcf,bcf", name="query_file", label="Query VCF/BCF Data")
    ]

    etree.SubElement(outputs, 'data', format="tabular", name="output_file")
    command.text += "#if $genotypes_file:\n  -g $genotypes_file\n#end if\n"
    command.text += "$query_file\n>\n$output_file"

else:
    if '<in.vcf.gz>' in usage_line:
        elements = [
            etree.Element('param', type="data", format="vcf,bcf", name="input_file", label="VCF/BCF Data")
        ]
        command.text += "$input_file\n>\n$output_file"
    elif '<file.vcf>' in usage_line:
        elements = [
            etree.Element('param', type="data", format="fasta", name="reference_fasta", label="Reference Fasta"),
            etree.Element('param', type="data", format="vcf_bgz,bcf_bgz", name="input_file", label="VCF/BCF Data")
        ]
        etree.SubElement(outputs, 'data', format="fasta", name="output_file")
        command.text = "cat $reference_fasta | " + command.text
        command.text += "$input_file\n>\n$output_file"
    elif '<A.vcf.gz> [<B.vcf.gz>]' in usage_line:
        elements = [
            etree.Element('param', type="data", format="vcf,bcf", name="input_file1", label="VCF/BCF Data"),
            etree.Element('param', type="data", format="vcf,bcf", name="input_file2", label="Second VCF/BCF Data", optional="True"),
        ]
    elif '<A.vcf.gz> [<B.vcf.gz [...]]' in usage_line:
        elements = [
            etree.Element('param', type="data", format="vcf,bcf", name="input_file1", label="VCF/BCF Data", multiple="True"),
        ]
    elif '<A.vcf.gz> <B.vcf.gz> [...]':
        elements = [
            etree.Element('param', type="data", format="vcf,bcf", name="input_file1", label="VCF/BCF Data"),
            etree.Element('param', type="data", format="vcf,bcf", name="input_file2", label="Other VCF/BCF Datasets", multiple="True"),
        ]
    else:
        log.warning("Unhandled input files from command line usage: %s", usage_line)


for elem in elements[::-1]:
    inputs.insert(0, elem)

starter = {
    'name': -100,
    'label': -99,
    'type': -80,
    'help': -50,
    'truevalue': -30,
    'falsevalue': -29,
}

def priority(key):
    if key in starter:
        return starter[key]
    else:
        return ord(key[0])

#http://stackoverflow.com/questions/2741480/can-elementtree-be-told-to-preserve-the-order-of-attributes
def _serialize_xml(write, elem, encoding, qnames, namespaces):
    tag = elem.tag
    text = elem.text
    if tag is etree.Comment:
        write("<!--%s-->" % etree._encode(text, encoding))
    elif tag is etree.ProcessingInstruction:
        write("<?%s?>" % etree._encode(text, encoding))
    else:
        tag = qnames[tag]
        if tag is None:
            if text:
                write(etree._escape_cdata(text, encoding))
            for e in elem:
                _serialize_xml(write, e, encoding, qnames, None)
        else:
            write("<" + tag)
            items = elem.items()
            if items or namespaces:
                if namespaces:
                    for v, k in sorted(namespaces.items(),
                                       key=lambda x: x[1]):  # sort on prefix
                        if k:
                            k = ":" + k
                        write(" xmlns%s=\"%s\"" % (
                            k.encode(encoding),
                            etree._escape_attrib(v, encoding)
                            ))
                #for k, v in sorted(items):  # lexical order
                #for k, v in items: # Monkey patch
                log.debug('LOOP')
                for k, v in sorted(items, key=lambda x: priority(x[0])):
                    log.debug(k)
                    if isinstance(k, etree.QName):
                        k = k.text
                    if isinstance(v, etree.QName):
                        v = qnames[v.text]
                    else:
                        v = etree._escape_attrib(v, encoding)
                    write(" %s=\"%s\"" % (qnames[k], v))
            if text or len(elem):
                write(">")
                if text:
                    write(etree._escape_cdata(text, encoding))
                for e in elem:
                    _serialize_xml(write, e, encoding, qnames, None)
                write("</" + tag + ">")
            else:
                write(" />")
    if elem.tail:
        write(etree._escape_cdata(elem.tail, encoding))

etree._serialize_xml = _serialize_xml

import sys
tree = etree.ElementTree(tool)
tree.write(sys.stdout)
