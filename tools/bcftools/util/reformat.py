import xml.etree.ElementTree as ET
import gtree
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reformat XML')
    parser.add_argument('tools', type=file, nargs='+', help='Tool files to attempt to reformat')
    args = parser.parse_args()

    for xml_file in args.tools:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        xml_file.seek(0)
        outtree = gtree.GElementTree(root)
        outtree.gwrite(xml_file.name)
