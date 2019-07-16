import math
import operator

import pandas


class Utils():
    @staticmethod
    def getOneValueMathOp(op_name):
        """
        Returns a simple one value math operator such as log, sqrt, etc
        """
        return getattr(math, op_name)

    @staticmethod
    def getTwoValueBaseOp(op_name):
        """
        Returns a basic two value operator such as +, -, etc
        """
        return getattr(operator, "__" + op_name + "__")

    @staticmethod
    def getVectorPandaOp(op_name):
        """Returns a valid DataFrame vector operator"""
        return getattr(pandas.DataFrame, op_name)

    @staticmethod
    def getTwoValuePandaOp(op_name):
        """Returns a valid two value DataFrame operator"""
        if len(op_name) > 2:
            print(op_name, "Not a valid DataFrame op")
            exit(-1)
        return getattr(pandas.DataFrame, "__" + op_name + "__")

    @staticmethod
    def parseType(param):
        """Parse config parameter strings"""

        mod = param.strip("'")
        try:
            return int(mod)
        except ValueError:
            try:
                return float(mod)
            except ValueError:
                pass

        if mod in ('True', 'False'):
            return mod == 'True'

        if mod == "None":
            return None

        if mod[0] == '[' and mod[-1] == ']':
            if mod[1] in '0123456789':
                # integer array:
                #    Here we need to convert 1-base
                #     indexing to 0-base indexing
                return [x - 1 for x in map(int, mod[1:len(mod) - 1].split(','))]
            else:
                # string array
                return list(map(lambda x: x.strip(), mod[1:-1].split(',')))

        # String
        return mod
