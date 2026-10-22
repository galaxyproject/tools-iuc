# -*- coding: utf-8 -*-
import getopt
import os
import sys


def usage(info) -> str:
    text = "SynDivA script.\n\n"
    if info:
        text += info
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


def get_os_path_join(directory, filename):
    return os.path.join(directory, filename)


def get_os_path_name(input):
    return os.path.basename(input)


def check_pattern(pattern):
    authorized_pattern_letter = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                                 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', ':', '0',
                                 '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '*']
    return len([letter in authorized_pattern_letter for letter in pattern]) == len(pattern)


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

    def case(self):
        # Test des fichiers et repertoires
        if not self.input:
            sys.exit(usage("input (-i,--input) : \"%s\" must be indicated\n" % (self.input)))
        if not self.output_dir:
            sys.exit(usage("output directory (-o,--output_dir) : \"%s\" must be indicated\n" % (self.output_dir)))
        if not self.pattern:
            sys.exit(
                usage("Pattern of the sequence bank (-p,--pattern) : \"%s\" must be indicated\n" % (self.pattern)))
        if not self.site_res_5:
            sys.exit(usage(
                "Sequence of the restriction site in 5' (-5,--restriction-site-5) : \"%s\" must be indicated\n" % (
                    self.site_res_5)))
        if not self.site_res_3:
            sys.exit(usage(
                "Sequence of the restriction site in 3' (-3,--restriction-site-3) : \"%s\" must be indicated\n" % (
                    self.site_res_3)))

    def data_format(self):
        """
        Check if information are correct
        """
        # Run without arguments
        if len(sys.argv) == 1:
            sys.exit(usage(None))
            # Test input file argument
        if self.input and not os.path.isfile(self.input):
            print(self.input)
            print(os.path.isfile(self.input))

    def getargs(self):
        """
        Determine the files provided as arguments
        @return: Choosen options
        """
        # Sans argument
        if len(sys.argv) <= 1:
            sys.exit("Do './fibronectin.py -h' for a usage summary")
        # options test
        try:
            (opts, args) = getopt.getopt(sys.argv[1:], "i:o:p:5:3:h",
                                         ["input=", "output_dir=", "pattern=", "site_res_5=", "site_res_3="])
        except getopt.GetoptError as err:
            # print help information and exit:
            print(str(err))  # will print something like "option -a not recognized"
            sys.exit(usage(None))
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
                sys.exit(usage(None))
            else:
                assert False, "unhandled option"
            # Verification of cases
        self.case()
        self.data_format()
