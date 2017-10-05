#!/usr/bin/env python3

from bs4 import BeautifulSoup as bs
from SectionHandler import *	   
from CommandParse import CommandParse

class Document2Section:

    def __init__(self, html_file):
        
        with open(html_file, 'r') as f:
            html = f.read()

        self.parse(html)
        self.reorder()


    def reorder(self):
        # Replicate the section_nestmap whilst preserving orderings

        parent_ordering = [] # parent []
        child_ordering =  {} # parent : child []
        
        for title in self.section_order:

            if title not in section_nestmap:
                print("Unknown section:", title, file=sys.stderr)
                continue
                #exit(-1)

            expanded, parent_section_name = section_nestmap[title]

            # Unparented sections become parents
            if parent_section_name == "":
                parent_ordering.append( title )
                continue

            # Has parent, preserve ordering within group
            if parent_section_name not in child_ordering:
                child_ordering[parent_section_name] = []
                parent_ordering.append( parent_section_name )

            child_ordering[parent_section_name].append( title )

        # Now print out all

        # -- flat listing
        #for name in self.section_order:
        #    cxml = self.section_map[name].getSection()
        #    if cxml != "":
        #        self.inputs.appendChild(cxml)
        #
        #return 0

        self.inputs = doc.createElement('inputs')
        
        for parent in parent_ordering:
            # False parent, actual section
            if parent in section_nestmap:
                expanded = section_nestmap[parent][0]
                sxml = self.section_map[parent].getSection()
                
                if sxml == "":
                    continue
               
                sxml.setAttribute('expanded', str(expanded))
                self.inputs.appendChild( sxml )
                continue

            # True parent, create section and nest children within
            parent_section = doc.createElement('section')
            parent_section.setAttribute('name', '_'.join(parent.lower().split()))
            parent_section.setAttribute('title', parent)

            child_sections = child_ordering[parent]

            for child in child_sections:

                expanded = section_nestmap[child][0]
                cxml = self.section_map[child].getSection()

                if cxml == "":
                    continue

                cxml.setAttribute('expanded', str(expanded))
                parent_section.appendChild( cxml )

            self.inputs.appendChild(parent_section)

            

            
    def parse(self, text):

        bsobj = bs(text, 'html.parser')
        main_div = bsobj.body.find('div', attrs={'class':'col-md-9'})
        tables = main_div.find_all('table')
        
        self.section_map   = {}
        self.section_order = []
          
        for table in tables:
            # Previous H1
            h2 = table.find_previous_sibling('h2')    
            title = h2.text               

            if title not in self.section_map:               
                self.section_map[title] = Section(title)
                self.section_order.append(title)

            tmp_section = self.section_map[title]

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






dd = Document2Section(sys.argv[1])

with open('iqtree.inputs.xml','w') as f:

    macros = doc.createElement('macros')
    xml_inp = doc.createElement('xml')
    xml_inp.setAttribute('name', 'inputs')
    
    xml_inp.appendChild(dd.inputs)
    macros.appendChild(xml_inp)

    print(macros.toprettyxml(), file=f)  
    f.close()


CommandParse(
    "iqtree", dd.inputs, exclude_map,
    "iqtree.command.xml"
)

