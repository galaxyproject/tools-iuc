#!/usr/bin/env python3

from Names import Names
import sys

import defaults

class CommandParse:

    def __init__(self, tool, inputs, exclusion_map, filename):
        self.inputs = inputs
        self.exclude = exclusion_map

        with open(filename, 'w') as fout:

            print("<macros>",                                  file=fout)
            print("   <xml name='command' >",                  file=fout)
            print("     <command detect_errors='exit_code' >", file=fout, end="")
            print("<![CDATA[\n%s\n" % tool,                    file=fout)

            self.printDefaults(fout)
            self.printCheetah(fout)
            
            print("\n]]>\n",         file=fout)
            print("     </command>", file=fout)
            print("  </xml>",        file=fout)
            print("</macros>",       file=fout)



    def printDefaults(self, fout):
        
        for arg in defaults.exclude_map:
            val = defaults.exclude_map[arg]

            if val == True:
                print(arg, file=fout)
                continue

            if val == False:
                continue

            # use the default provided
            print("%s %s" % (arg, val), file=fout)

        print("", file=fout)
            
            
                
        

    def printCheetah(self, fout):
    
        def recurse(node, parent_list):
            if node.tagName == "param":
                arg = node.getAttribute("argument")

                if arg == "":
                    return

                argname = Names.sensibleCheetah(arg)
                fullname= '.'.join( parent_list + [argname] )

                print("Updated name from:", arg, "to", fullname, file=sys.stderr)

                if arg in self.exclude:
                    val = self.exclude[arg]
                    if val != False:
                        print("%s %s\n\n" % (flag, arg), file=fout)
                        
                    return


                ftype = node.getAttribute("type")            

                # text, integer, float, file, boolean
                #
                if ftype == "boolean":
                    # just print the var name, trueval and falseval fill the gap
                    print("$%s\n\n" % fullname, file=fout)

                elif ftype == "data":
                    print('''#if str( $%s ) != "None":\n%s $%s\n#end if''' % (
                        fullname, arg, fullname
                    ), file=fout)                    
                else:
                    print('''#if str( $%s ) != "":\n%s $%s\n#end if''' % (
                        fullname, arg, fullname
                    ), file=fout)
                   
                return


            if node.tagName == "section":
                name = node.getAttribute("name")
                
                for child in node.childNodes:
                    recurse(child, parent_list + [name])
                return

            
            if node.tagName == "conditional":
                name = node.getAttribute("name")               

                childs = node.childNodes
                
                bool_param = childs[0]
                true_when  = childs[1] if childs[1].getAttribute("value") == "true" else childs[2]
                false_when = childs[2] if childs[1].getAttribute("value") == "true" else childs[1]

                fullname = '.'.join(parent_list + [name, bool_param.getAttribute('name')])

                print("#if str( $%s ):\n" % fullname, file=fout)

                #import pdb;pdb.set_trace()

                for true_child in true_when.childNodes:
                    recurse(true_child, parent_list + [name]) # DONT skip conditional name

                print("#else:\n", file=fout)

                for false_child in false_when.childNodes:
                    recurse(false_child, parent_list + [name])

                print("#end if\n\n", file=fout)
                return


            for child in node.childNodes:
                recurse(child, parent_list) # don't add non-section names
            
            
            print("Unable to parse:", node, file=sys.stderr)

        recurse(self.inputs, [])
