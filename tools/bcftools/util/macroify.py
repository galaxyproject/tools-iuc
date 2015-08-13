import xml.etree.ElementTree as ET
import copy
import os
import argparse

def hash_param(param):
    data = ''
    for key in sorted(param.attrib.keys()):
        data += '%s||%s++' % (key, param.attrib[key])
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Macro-ify some automatically generated wrappers')

    parser.add_argument('--macros', help='new/existing macros file', default='macros.xml')
    parser.add_argument('tools', type=file, nargs='+', help='Tool files to attempt to reduce')

    args = parser.parse_args()


    hashes = {}

    for xml_file in args.tools:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for param in root.findall('.//inputs/section/*'):
            h = hash_param(param)
            if h not in hashes:
                hashes[h] = {'param': param, 'files': []}

            hashes[h]['files'].append(xml_file)

    dupes = [key for key in hashes if len(hashes[key]['files']) > 1]

    try:
        macro_file = ET.parse(args.macros)
        macros = macro_file.getroot()
    except:
        macros = ET.Element("macros")

    macro_elements = []
    for dupe_key in dupes:
        # Append to macros
        macro_id = 'macro_%s' % hashes[dupe_key]['param'].attrib['name']
        macro_element = ET.Element('xml', name=macro_id)
        macro_element.append(copy.deepcopy(hashes[dupe_key]['param']))

        macro_elements.append(macro_element)
        # Store param that can be copied in as replacement.
        hashes[dupe_key]['mid'] = macro_id


    em = [hash_param(elem) for elem in macros.findall('.//macros/*')]

    for me in macro_elements:
        if hash_param(me) not in em:
            macros.append(copy.deepcopy(me))

    with open(args.macros, 'w') as handle:
        handle.write(ET.tostring(macros))

    for xml_file in args.tools:
        xml_file.seek(0)
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for section in root.findall('.//inputs/section'):
            children = section.getchildren()

            for idx in range(len(children)):
                child = children[idx]
                h = hash_param(child)
                if h in dupes:
                    section.insert(idx, ET.Element('expand', macro=hashes[h]['mid']))
                    section.remove(child)

        xml_file.seek(0)
        tree.write(xml_file.name)
