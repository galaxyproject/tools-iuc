#!/usr/bin/env python
from lxml import etree as ET
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
    'datafile': {'file_formats': ['data']},
    'directory': {},
    'dirlist': {},
    'features': {},
    'filelist': {},
    'infile': {'file_formats': ['data']},
    'matrix': {},
    'matrixf': {},
    'obo': {'file_formats': ['obo']},
    'pattern': {},
    'refseq': {},
    'regexp': {},
    'resource': {'file_formats': ['data']},
    'scop': {},
    'sequence': {'file_formats': ['fasta']},
    'seqall': {'file_formats': ['fasta']},
    'seqset': {'file_formats': ['fasta']},
    'seqsetall': {'file_formats': ['fasta']},
    'url': {'type': str},
    'variation': {'file_formats': ['vcf']},
    'xml': {'file_formats': ['xml']},
    'taxon': {'type': str},
    'text': {'type': str},
}


tree = ET.parse(sys.argv[1])
root = tree.getroot()

tool = ET.Element("tool",
                  id='emboss_%s' % root.attrib['id'],
                  name=root.attrib['id'],
                  version='6.6.0'
                  )
desc = ET.Element("description")
desc.text = root.find('metadata').find('documentation').text
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


def handle_parameter(child):
    kwargs = {
        'type': child.attrib['type']
    }
    if child.attrib['type'] in INPUT_TYPE_MAPPING.keys():
        kwargs['type'] = 'input-file'
        if child.find('knowntype') is not None:
            kwargs['file_formats'] = [child.find('knowntype').text]

        kwargs.update(INPUT_TYPE_MAPPING[child.attrib['type']])
    elif child.attrib['type'] in ('integer', 'float'):
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

        kwargs['num_range'] = (num_min, num_max)
        kwargs['default'] = num_default
        kwargs['type'] = int if child.attrib['type'] == 'integer' else float
    elif child.attrib['type'] in ('boolean', 'toggle'):
        kwargs['type'] = str
        kwargs['choices'] = ['true', 'false']
    elif child.attrib['type'] == 'list':
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

        kwargs['type'] = str
        kwargs['choices'] = [kv[k] for k in kv.keys()]
    elif child.attrib['type'] == 'selection':
        if child.find('delimiter') is not None:
            delimiter = child.find('delimiter').text
        else:
            delimiter = ';'

        kv = {}
        for i, x in enumerate(child.find('values').text.split(delimiter)):
            if len(x.strip()) > 0:
                kv[str(i)] = x.strip()

        kwargs['type'] = str
        kwargs['choices'] = [kv[k] for k in kv.keys()]
    elif child.attrib['type'] == 'range':
        # TODO
        kwargs['type'] = 'string'
    elif child.attrib['type'] == 'array':
        # TODO
        kwargs['type'] = 'string'

    if child.find('information') is not None:
        kwargs['description'] = child.find('information').text

    log.debug(pprint.pformat(kwargs))
    return kwargs


def build_section(acd_xml_section):
    log.debug('build_section %s', acd_xml_section.tag)
    if acd_xml_section.tag == 'parameter':
        kwargs, command_string = handle_parameter(acd_xml_section)
        log.debug('Param %s' % pprint.pformat(kwargs))
        command_string = "## $%s\n" % acd_xml_section.attrib['name']
        parameter = ET.Element('parameter')
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

tool.append(command)
tool.append(inputs)
tool.append(outputs)
tool_help = ET.Element("help")
tool_help.text = "TODO\n\n@ATTRIBUTION@\n  "
citations = ET.Element("expand", macros="citation")

tool.append(tool_help)
tool.append(citations)


print ET.tostring(tool)
