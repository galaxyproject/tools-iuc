#!/usr/bin/env python

import sys

class Names:

    """--really-long-dashes  ->  --really_long_dashes 
    Planemo strips leading dashes, but throws errors for internal dashes
    """
    @staticmethod
    def sensibleFlag(arg):

        if arg.startswith("--"):
            return "--" + arg[2:].replace('-','_')
            
        if arg.startswith("-"):
            return "-" + arg[1:].replace('-','_')

        print("Unable to sanitize", arg, file=sys.stderr)
        exit(-1)


    @staticmethod
    def sensibleCheetah(arg):
        argname = Names.sensibleFlag(arg)

        if argname.startswith("--"):
            return argname[2:]

        if argname.startswith("-"):
            return argname[1:]

        print("Unable to parse sensible cheetah", file=sys.stderr)
        exit(-1)

    
