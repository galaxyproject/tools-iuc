# -*- coding: utf-8 -*-
import sys, getopt, os


class Args:

    def __init__(self):
        """
        Instanciate Files object
        """
        self.input = None
        self.output_dir = None
        self.pattern = None
        self.site_res_5 = None
        self.site_res_3 = None
        self.getargs()

    def usage(self, info):
        text = None
        text = "Fibronectin script.\n\n"
        if (info): text += info
        temp = "Option\t\t\t\tfile\t\t\tDescription\n"
        text += temp
        text += '-' * (len(temp) + 60)
        text += '\n'
        text += "-i, --input\t\t\tfile.fasta\t\tFasta file that contains the DNA sequences\n"
        text += "-o, --output_dir\t\t/path/for/output\tDirectory where output files will be written\n"
        text += "-p, --pattern\t\t\tstring\t\t\tPattern of the sequence bank\n"
        text += "-5, --restriction-site-5\tstring\t\t\tSequence of the restriction site in 5'\n"
        text += "-3, --restriction-site-3\tstring\t\t\tSequence of the restriction site in 3'\n"
        return text

    def case(self):
        # Test des fichiers et repertoires
        if not self.input:
            sys.exit(self.usage("input (-i,--input) : \"%s\" must be indicated\n" % (self.input)))
        if not self.output_dir:
            sys.exit(self.usage("output directory (-o,--output_dir) : \"%s\" must be indicated\n" % (self.output_dir)))
        if not self.pattern:
            sys.exit(
                self.usage("Pattern of the sequence bank (-p,--pattern) : \"%s\" must be indicated\n" % (self.pattern)))
        if not self.site_res_5:
            sys.exit(self.usage(
                "Sequence of the restriction site in 5' (-5,--restriction-site-5) : \"%s\" must be indicated\n" % (
                    self.site_res_5)))
        if not self.site_res_3:
            sys.exit(self.usage(
                "Sequence of the restriction site in 3' (-3,--restriction-site-3) : \"%s\" must be indicated\n" % (
                    self.site_res_3)))

    def data_format(self):
        """
        Check if information are correct
        """
        # Run without arguments
        if len(sys.argv) == 1:
            sys.exit(self.usage(None))
            # Test input file argument
        if self.input:
            if not os.path.isfile(self.input):
                print(self.input)
                print(os.path.isfile(self.input))
                #sys.exit(self.usage("Error with \"%s\" : -i required an input file\n" % self.multilist))

                # Determine les fichiers fournis en arguments

    def getargs(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """
        # Sans argument
        if len(sys.argv) <= 1: sys.exit("Do './fibronectin.py -h' for a usage summary")
        # test des option
        try:
            (opts, args) = getopt.getopt(sys.argv[1:], "i:o:p:5:3:h",
                                         ["input=", "output_dir=", "pattern=", "site_res_5=", "site_res_3="])
        except getopt.GetoptError as err:
            # print help information and exit:
            print(str(err))  # will print something like "option -a not recognized"
            sys.exit(self.usage(None))
        # Identification of options
        for (o, a) in opts:
            if o in ("-i", "--input"):
                self.input = a
            elif o in ("-o", "--output_dir"):
                self.output_dir = a
            elif o in ("-p", "--pattern"):
                self.pattern = a
            elif o in ("-5", "--restriction-site-5"):
                self.site_res_5 = a
            elif o in ("-3", "--restriction-site-3"):
                self.site_res_3 = a
            elif o in ("-h", "--help"):
                sys.exit(self.usage(None))
            else:
                assert False, "unhandled option"
            # Verification of cases
        self.case()
        self.data_format()
