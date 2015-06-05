#!/usr/bin/env python
from lxml import etree as ET
import re
import sys
import logging
import pprint
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

ACD_INPUTS = ('additional', 'advanced', 'input', 'required')
ACD_OUTPUTS = ('output', )

INPUT_TYPE_MAPPING = {
    'assembly': {},
    'codon': {},
    'cpdb': {},
    'datafile': {'file_format': 'data'},
    'directory': {},
    'dirlist': {},
    'features': {},
    'filelist': {},
    'infile': {'file_format': 'data'},
    'matrix': {},
    'matrixf': {},
    'obo': {'file_format': 'obo'},
    'pattern': {},
    'refseq': {},
    'regexp': {},
    'resource': {'file_format': 'data'},
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


def __input_file(child):
    #<param label="Genome" name="positional_1" type="data" format="fasta"/>
    kwargs = {
        'type': 'data'
    }
    if child.find('knowntype') is not None:
        kwargs['file_format'] = child.find('knowntype').text

    #if child.find('information') is not None:
        #kwargs['help'] = child.find('information').text

    if child.find('help') is not None:
        kwargs['help'] = unbreak_strings(child.find('help').text)
    # Unhandled: additional, default, nullok, relations, ... ?

    kwargs.update(INPUT_TYPE_MAPPING[child.attrib['type']])
    return kwargs


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

    kwargs = {
        'maximum': str(num_max),
        'default': str(num_default),
        'type': child.attrib['type']
    }

    if num_min is not None:
        kwargs['minimum'] = str(num_min)

    if num_max is not None:
        kwargs['maximum'] = str(num_max)

    return kwargs


def __boolean_input(child):
    kwargs = {
        'type': 'boolean',
        'name': child.attrib['name'],
        'truevalue': '-%s' % child.attrib['name'],
        'falsevalue': '',
    }
    return kwargs


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

    kwargs['type'] = "string"
    #kwargs['choices'] = [kv[k] for k in kv.keys()]
    return kwargs


def __selection_input(child):
    kwargs = {}
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
    return kwargs


def handle_parameter(child):
    kwargs = {'name': child.attrib['name']}
    if child.attrib['type'] in INPUT_TYPE_MAPPING.keys():
        kwargs.update(__input_file(child))
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] in ('integer', 'float'):
        kwargs.update(__numeric_input(child))
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] in ('boolean', 'toggle'):
        kwargs.update(__boolean_input(child))
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] == 'list':
        kwargs.update(__list_input(child))
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] == 'selection':
        kwargs.update(__selection_input(child))
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] == 'range':
        # TODO
        kwargs['type'] = 'string'
        parameter = ET.Element('parameter', **kwargs)
    elif child.attrib['type'] == 'array':
        # TODO
        kwargs['type'] = 'string'

    if child.find('information') is not None:
        kwargs['label'] = unbreak_strings(child.find('information').text)
        parameter = ET.Element('parameter', **kwargs)

    return parameter, None


def build_section(acd_xml_section):
    log.debug('build_section %s', acd_xml_section.tag)
    if acd_xml_section.tag == 'parameter':
        parameter, command_string = handle_parameter(acd_xml_section)
        command_string = "## $%s\n" % acd_xml_section.attrib['name']
        return command_string, parameter
    elif acd_xml_section.tag == 'section':
        # Get section information
        metadata = acd_xml_section.find('metadata')
        if metadata is not None:
            section_title = metadata.find('information').text
        else:
            section_title = "None"  # acd_xml_section.attrib['id']
        log.info('X %s', pprint.pformat(acd_xml_section.attrib))

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
        command.text += command_string
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
                command.text += '-%s $%s' % (parameter.attrib['name'], parameter.attrib['name'])

tool.append(command)
tool.append(inputs)
tool.append(outputs)
tool_help = ET.Element("help")
with open(sys.argv[2], 'r') as handle:
    tool_help.text = ET.CDATA("\n" + clean_string(handle.read()))
citations = ET.Element("expand", macros="citation")

tool.append(tool_help)
tool.append(citations)


print ET.tostring(tool)
