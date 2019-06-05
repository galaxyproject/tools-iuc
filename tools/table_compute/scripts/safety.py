import math
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

# Removing linting error on Travis
# with junk function  declaration
__ = math.log, np.log, pd.DataFrame.sum


class Safety():
    """
    Class to safely evaluate mathematical expression on single
    or table data
    """

    __conditionals = ('(', ')', '()', '.', 'if', 'else')
    __basicops = ('+', '-', '*', '/', '%', '!=', '==')
    __allowed = ('np', 'math')
    __disallow = ('memory_usage', 'to_clipboard',
                  'to_csv', 'to_dense', 'to_dict', 'to_excel', 'to_feather',
                  'to_gbq', 'to_hdf', 'to_html', 'to_json', 'to_latex',
                  'to_msgpack', 'to_panel', 'to_parquet', 'to_period',
                  'to_pickle', 'to_records', 'to_sparse', 'to_sql',
                  'to_stata', 'to_string', 'to_timestamp', 'to_xarray',
                  # builtins
                  'open', 'write', 'os', 'exec', 'eval', 'return', 'import')

    def __init__(self, obj_name, custom_string, strict_functiondefs):
        if isinstance(obj_name, list):
            self.name = obj_name[0]
            self.allowed_names = obj_name + list(self.__allowed)
        else:
            self.name = obj_name
            self.allowed_names = list(self.__allowed)

        self.expr = custom_string
        self.parseXML(strict_functiondefs)
        self.__assertSafe()

    def parseXML(self, xml_file):
        root = ET.parse(xml_file).getroot()
        options = [option.attrib["value"]
                   for macro in root for option in macro
                   if macro.tag == "macro" and option.tag == "option"]

        self.allowed_names = self.allowed_names + options

    def generateMultiLine(self):
        """
        Generates a multi-line function to be evaluated outside the class
        """
        cust_fun = "def fun(%s):%s" % (
            self.name, "\n".join("\t" + x for x in self.expr.splitlines()))
        return cust_fun

    def generateFunction(self):
        """
        Generates a function to be evaluated outside the class
        """
        cust_fun = "def fun(%s):\n\treturn(%s)" % (self.name, self.expr)
        return cust_fun

    def __assertSafe(self):
        indeed = self.__isSafeStatement()
        if not indeed:
            print("Custom Expression is not safe.")
            exit(-1)

    @staticmethod
    def arrange(str_array):
        """
        Method to always arrange list set objects in
        the same way across multiple tests
        """
        return sorted(str_array, key=len, reverse=True)

    @staticmethod
    def detailedExcuse(word):
        """
        Gives a verbose statement for why users should not use some specific operators.
        """
        mess = None
        if word == "for":
            mess = "For loops are not permitted. Use numpy or pandas table operations instead."
        elif word == ":":
            mess = "Please do not use colons for expressions. Use inline Python if/else statements"
        elif word == "=":
            mess = "Please do not use variable assignment in expressions. Use the in-built functions for substituting values"
        elif word in (",", "[", "]"):
            mess = "Please do not subscript arrays in custom functions. Use the in-built functions for substituting values"
        else:
            mess = "Not a safe operation"

        print("( '%s' ) %s" % (word, mess))

    def __isSafeStatement(self):
        """
        Determines how safe a user-expression is.

        We wish only for operators from the math, pd, or np library
        to be used as well as the builtins of +/-* and inline if/else.
        """
        #  '-log(1 - elem/4096) * 4096 if elem != bn else elem - 0.5'
        # 'vec.median() +  vec.sum()'
        safe = True

        # 1. Split statement into keywords and operators
        keywords = Safety.arrange(set(re.findall(r'[a-zA-Z0-9]+', self.expr)))
        # ['log', '1', 'elem', '4096', 'if', 'else', '0', '5']
        # ['sum', 'median', 'vec']
        operators = Safety.arrange(set(re.findall(r'[^a-zA-Z0-9_.() ]+', self.expr)))
        # ['!=', '*', '-', '/']
        # ['+']

        # 2. Check operators for basic ops
        for opw in Safety.arrange(filter(
                lambda x: x not in self.__conditionals, Safety.arrange(operators))):
            if opw not in self.__basicops:
                Safety.detailedExcuse(opw)
                safe = False

        # 3. Check keywords for numeric, if, else, mathfns, and elem
        for keyw in Safety.arrange(keywords):
            try:
                # '.' isn't matched, so only ints expected
                int(keyw)
            except ValueError:
                is_good, is_bad = self.__checkKeyword(keyw)
                if is_bad or not is_good:
                    Safety.detailedExcuse(keyw)
                    safe = False

        # 4. Test the remaining string
        valid_ops = operators + keywords \
            + Safety.arrange(list(self.__conditionals)) \
            + Safety.arrange(list(self.allowed_names))

        remstring = self.expr

        for rem in Safety.arrange(valid_ops):
            remstring = remstring.replace(rem, "")

        remstring = remstring.strip()

        if remstring:
            print("Unable to parse leftover string:", remstring)
            safe = False

        return safe

    def __checkKeyword(self, keyw):
        """Tests a token keyword for malicious intent"""
        is_name = keyw == self.name
        is_ifel = keyw in self.__conditionals
        # is_math = hasattr(math, keyw)
        # is_pd = hasattr(pd.DataFrame, keyw)
        # is_np = hasattr(np, keyw)
        is_allow = keyw in self.allowed_names
        not_dis = keyw in self.__disallow

        # is_good = is_name or is_ifel or is_math or is_pd or is_np or is_allow
        is_good = is_name or is_ifel or is_allow
        is_bad = not_dis
        return is_good, is_bad
