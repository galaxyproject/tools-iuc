import re
import logging
import fileinput
import xml.etree.cElementTree as etree
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
    else:
        section = line.replace(' options:', '').strip()
        section = 'Default'
        log.debug("Saw section %s" % section)


tool = etree.Element("tool",
                     id="bcftools_@EXECUTABLE@",
                     name="bcftools @EXECUTABLE@",
                     version="@VERSION@.0")
desc = etree.SubElement(tool, 'description')
desc.text = about_line.split('. ')[0]
etree.SubElement(tool, 'expand', macro="requirements")
etree.SubElement(tool, 'expand', macro="version_command")
etree.SubElement(tool, 'expand', macro="stdio")
macros = etree.SubElement(tool, 'macros')
etree.SubElement(macros, 'token', name='@EXECUTABLE@').text = 'view'
etree.SubElement(macros, 'import').text = 'bcftools_macros.xml'
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
    safe_sec_name = 'sec_' + section_name.replace(' ', '_').lower()
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
                command.text += "${%s.%s}\n" % (safe_sec_name, pkw['name'])
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
                command.text += "${%s.%s}\n" % (safe_sec_name, pkw['name'])
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
                elif '<' in parsed_command['param'] and '|' in parsed_command['param']:
                    parsed_command['select'] = True
                    help_text = pkw['help']
                    default = None
                    # Get the default if it exists
                    m = re.match(DEFAULT_VALUE, pkw['help'])
                    if m:
                        help_text = m.group('help_text')
                        default = m.group('default_value')

                    if default is None:
                        pkw['optional'] = 'True'

                    select = etree.SubElement(section, 'param',
                                              type='select',
                                              name='select_' + pkw['name'])

                    for kv in help_text.split(', '):
                        kvd = kv.split(': ')
                        k = kvd[0]
                        v = ': '.join(kvd[1:])

                        sokwd = {'value': k}
                        if k == default:
                            sokwd['selected'] = 'True'

                        etree.SubElement(select, 'option', **sokwd).text = v

                    command.text += "--%s \"${%s.%s}\"\n" % (parsed_command['long'], safe_sec_name, 'select_' + pkw['name'])

                elif '<region>' in parsed_command['param']:
                    repeat = etree.SubElement(section, 'repeat',
                                              name=pkw['name'] + '_repeat',
                                              title=pkw['label'])
                    etree.SubElement(repeat, 'param', **pkw)
                    command_id = 'values_%s_%s' % (safe_sec_name, pkw['name'])
                    command.text += """#set %s = '","'.join([str($value) for $value in $%s.%s])\n""" % (command_id, safe_sec_name, pkw['name'] + "_repeat")
                    command.text += "--%s \"${%s}\"\n" % (parsed_command['long'], command_id)

                elif '<list>' in parsed_command['param']:
                    repeat = etree.SubElement(section, 'repeat',
                                              name=pkw['name'] + '_repeat',
                                              title=pkw['label'])
                    etree.SubElement(repeat, 'param', **pkw)

                    command_id = 'values_%s_%s' % (safe_sec_name, pkw['name'])
                    command.text += """#set %s = '","'.join([str($value) for $value in $%s.%s])\n""" % (command_id, safe_sec_name, pkw['name'] + "_repeat")
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
                        command.text += "--%s \"${%s.%s}${%s}\"\n" % (parsed_command['long'], safe_sec_name, mpkw['name'], command_id)
                    else:
                        command.text += "--%s \"${%s}\"\n" % (parsed_command['long'], command_id)
                else:
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

                        if parsed_command['type'] == 'param/exclude':
                            # If it's an exclude instead of a carat, we've
                            # reduced the two flags into a single real option
                            # plus the inverter. This will be used in the template instead of the '^' or flag value.
                            mpkw.update({
                                'truevalue': parsed_command['flag_choices'][0],
                                'falsevalue': parsed_command['flag_choices'][1],
                            })
                            command.text += "--{%s.%s} \"${%s.%s}\"\n" % (
                                safe_sec_name,
                                mpkw['name'],
                                safe_sec_name,
                                pkw['name']
                            )
                        else:
                            command.text += "--%s \"${%s.%s}${%s.%s}\"\n" % (
                                parsed_command['long'],
                                safe_sec_name,
                                mpkw['name'],
                                safe_sec_name,
                                pkw['name']
                            )

                        etree.SubElement(section, 'param', **mpkw)

                    else:
                        command.text += "--%s \"${%s.%s}\"\n" % (
                            parsed_command['long'],
                            safe_sec_name,
                            pkw['name']
                        )
            else:
                print parsed_command




print etree.tostring(tool)
