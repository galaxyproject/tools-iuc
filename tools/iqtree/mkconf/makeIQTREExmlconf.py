#!/usr/bin/env python3

from bs4 import BeautifulSoup as bs
from SectionHandler import *	   
from CommandParse import CommandParse


section_map = {} #  title -> section

class Document2Section:

    def __init__(self, html_file):
        
        with open(html_file, 'r') as f:
            html = f.read()
            self.inputs = Document2Section.parse(html)


    @staticmethod
    def parse(text):

        bsobj = bs(text, 'html.parser')
        main_div = bsobj.body.find('div', attrs={'class':'col-md-9'})
        tables = main_div.find_all('table')
        
        return_text=""

        inputs = doc.createElement('inputs')

        section_order = []
          
        for table in tables:
            # Previous H1
            h2 = table.find_previous_sibling('h2')    
            title = h2.text               

            #if tmp_section != None:
            #    sect_xml = tmp_section.getSection()
            #    if sect_xml != "":
            #        inputs.appendChild(sect_xml)
            #        h2_map[title] = True  # add to title map on non-empty sections

            tmp_section = None
            if title not in section_map:
                section_map[title] = Section(title)
                section_order.append(title)

            tmp_section = section_map[title]

            # Each table is a section
            thead = table.thead
            thead1 = thead.th.text
            thead2 = thead.th.nextSibling.text

            datatype_table = False

            if thead1 == "Option" and thead2 == "Usage and meaning":
                print("Section:", title, file=sys.stderr)
                #continue
            else:
                continue
                # Special tables for defining extra -m param types (RateType, DataType, FreqType)
                #datatype_table = True
                #print( "- subsection", thead1, thead2, file=sys.stderr)

            tbody = table.tbody

            try:
                helps = [tr.td.nextSibling.text for tr in tbody.findAll('tr')]
                codes = [tr.td.nextSibling.findAll('code') for tr in tbody.findAll('tr')]
                flags = []

                if not(datatype_table):       
                    flags = [tr.td.code.text for tr in tbody.findAll('tr')]
                else:
                    tmp_section.setSubSection(thead1)
                    flags = [tr.td.text for tr in tbody.findAll('tr')]

            except AttributeError:
                import pdb; pdb.set_trace()

            if len(flags) != len(helps):
                print("Table mismatch", file=sys.stderr)
                exit(-1)

            flag_helper = list(zip(flags, helps, codes))

            for flag, helper, codelist in flag_helper:    
                parts = flag.split(' ')
                flag = parts[0]
                flag_params = None

                if len(parts) > 1:
                    flag_params = parts[1:]

                #flag, flag_params, short, help
                tmp_section.insertFlag(flag, flag_params, helper, codelist)


        # Print final
        for section in section_order:
            sect_xml = section_map[section].getSection()
            if sect_xml != "":
                inputs.appendChild( sect_xml )
        
        return inputs



dd = Document2Section(sys.argv[1])
cc = CommandParse(flag_map, exclude_map)

with open('iqtree.macros.xml','w') as f:

    macros = doc.createElement('macros')
    
    xml_inp = doc.createElement('xml')
    xml_inp.setAttribute('name', 'inputs')

    xml_com = doc.createElement('xml')
    xml_com.setAttribute('name', 'command')

    command = doc.createElement('command')   
    command.appendChild(
        doc.createTextNode(cc.text)
    )


    xml_inp.appendChild(dd.inputs)
    xml_com.appendChild(command)

    macros.appendChild(xml_inp)
    macros.appendChild(xml_com)
    
    print(macros.toprettyxml(), file=f)
   
    f.close()

#import pdb;pdb.set_trace()
