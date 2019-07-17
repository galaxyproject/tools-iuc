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
    __allowedlibs = ('np', 'math')
    __allowed = (
        # Numpy
        'abs', 'absolute', 'add', 'alen', 'all', 'alltrue', 'amax', 'amin',
        'angle', 'any', 'append', 'arccos', 'arccosh', 'arcsin',
        'arcsinh', 'arctan', 'arctan2', 'arctanh', 'argmax', 'argmin',
        'argsort', 'around', 'array',
        'array2string',
        'array_equal', 'array_equiv', 'array_split', 'array_str', 'asmatrix',
        'asscalar', 'average', 'bincount', 'bitwise_and', 'bitwise_not', 'bitwise_or', 'bitwise_xor',
        'bmat', 'bool', 'byte', 'cdouble', 'ceil', 'cfloat', 'char', 'character',
        'choose', 'clip', 'complex', 'concatenate', 'conj', 'conjugate', 'convolve',
        'correlate', 'cos', 'cosh', 'count_nonzero', 'cov', 'cross', 'cumprod', 'cumproduct',
        'cumsum', 'datetime64', 'datetime_as_string', 'deg2rad', 'degrees', 'diag',
        'diagonal', 'digitize', 'divide', 'division', 'divmod', 'dot', 'double',
        'e', 'empty', 'equal', 'euler_gamma', 'exp', 'exp2', 'expm1',
        'fabs', 'fft', 'fill_diagonal', 'flip', 'fliplr', 'flipud', 'float',
        'float_power', 'floor', 'floor_divide', 'fmax', 'fmin', 'fmod', 'half', 'hamming',
        'hanning', 'heaviside', 'histogram', 'histogram2d', 'hsplit', 'hypot',
        'indices', 'inf', 'inner', 'insert', 'int', 'integer', 'invert',
        'iscomplex', 'isfinite', 'isin', 'isinf', 'isnan', 'isnat', 'isneginf', 'isposinf',
        'isreal', 'isscalar', 'kaiser', 'kron', 'lcm', 'ldexp', 'left_shift', 'less', 'less_equal',
        'lexsort', 'linalg', 'log', 'log10', 'log1p', 'log2', 'logaddexp', 'logaddexp2',
        'logical_and', 'logical_not', 'logical_or', 'logical_xor', 'logspace', 'long',
        'longcomplex', 'longdouble', 'longfloat', 'longlong',
        'mask_indices', 'mat', 'matmul', 'matrix', 'max', 'maximum', 'mean', 'median',
        'min', 'minimum', 'mod', 'modf', 'moveaxis', 'msort', 'multiply', 'nan', 'nan_to_num',
        'nanargmax', 'nanargmin', 'nancumprod', 'nancumsum', 'nanmax', 'nanmean', 'nanmedian',
        'nanmin', 'nanpercentile', 'nanprod', 'nanquantile', 'nanstd', 'nansum', 'nanvar',
        'nbytes', 'ndarray', 'ndim', 'negative', 'nonzero', 'not_equal', 'nper',
        'numarray', 'number', 'ones', 'outer', 'pad', 'partition', 'percentile', 'pi',
        'piecewise', 'poly', 'poly1d', 'polyadd', 'polyder', 'polydiv', 'polyfit', 'polyint',
        'polymul', 'polynomial', 'polysub', 'polyval', 'positive', 'power',
        'prod', 'product', 'pv', 'quantile', 'rad2deg', 'radians', 'random', 'rank', 'rate',
        'real', 'reciprocal', 'remainder', 'repeat', 'reshape', 'resize',
        'right_shift', 'rint', 'roll', 'rollaxis', 'roots', 'rot90', 'round',
        'short', 'signedinteger', 'sin', 'sinc', 'single', 'singlecomplex', 'sinh', 'size',
        'sometrue', 'sort', 'sort_complex', 'spacing', 'split', 'sqrt', 'square', 'squeeze',
        'stack', 'std', 'str', 'subtract', 'sum', 'swapaxes', 'take', 'take_along_axis',
        'tan', 'tanh', 'tensordot', 'tile', 'transpose', 'trapz', 'trim_zeros',
        'true_divide', 'ubyte', 'uint', 'ulonglong', 'unicode',
        'unique', 'unwrap', 'ushort', 'vander', 'var', 'vdot', 'vectorize',
        'vsplit', 'vstack', 'where', 'zeros',
        # Math
        'abs', 'all', 'any', 'ascii', 'bin', 'chr',
        'divmod', 'hash', 'hex', 'id', 'len', 'max',
        'min', 'oct', 'ord', 'pow', 'print', 'round',
        'sorted', 'sum', 'complex', 'dict', 'enumerate',
        'filter', 'float', 'int', 'list', 'map', 'range',
        'reversed', 'set', 'slice', 'str', 'type', 'zip',
        'acos', 'acosh', 'asin', 'asinh', 'atan', 'atan2',
        'atanh', 'ceil', 'cos', 'cosh', 'degrees', 'e',
        'exp', 'expm1', 'fabs', 'factorial', 'floor', 'fmod',
        'frexp', 'fsum', 'gamma', 'hypot', 'inf', 'isclose',
        'isfinite', 'isinf', 'isnan', 'ldexp', 'lgamma',
        'log', 'log10', 'log1p', 'log2', 'modf', 'nan',
        'pi', 'pow', 'radians', 'remainder', 'sin', 'sinh',
        'sqrt', 'tan', 'tanh', 'tau', 'trunc',
        # Pandas DataFrame
        'abs', 'add', 'agg', 'aggregate', 'align', 'all', 'any', 'append',
        'apply', 'applymap', 'as_matrix', 'asfreq', 'at', 'axes', 'bool',
        'boxplot', 'clip', 'clip_lower', 'clip_upper', 'columns', 'combine',
        'compound', 'corr', 'count', 'cov', 'cummax', 'cummin', 'cumprod',
        'cumsum', 'div', 'divide', 'dot', 'drop', 'drop_duplicates',
        'droplevel', 'dropna', 'duplicated', 'empty', 'eq', 'equals',
        'expanding', 'ffill', 'fillna', 'filter', 'first', 'first_valid_index',
        'floordiv', 'ge', 'groupby', 'gt', 'head', 'hist', 'iloc', 'index',
        'insert', 'interpolate', 'isin', 'isna', 'isnull', 'items', 'iteritems',
        'iterrows', 'itertuples', 'ix', 'join', 'keys', 'kurt', 'kurtosis',
        'last', 'last_valid_index', 'le', 'loc', 'lookup', 'lt', 'mad',
        'mask', 'max', 'mean', 'median', 'melt', 'merge', 'min', 'mod',
        'mode', 'mul', 'multiply', 'ndim', 'ne', 'nlargest', 'notna', 'notnull',
        'nsmallest', 'nunique', 'pct_change', 'pivot', 'pivot_table', 'plot',
        'pop', 'pow', 'prod', 'product', 'quantile', 'radd', 'rank', 'rdiv',
        'replace', 'resample', 'rfloordiv', 'rmod', 'rmul', 'rolling', 'round',
        'rpow', 'rsub', 'rtruediv', 'sample', 'select',
        'sem', 'shape', 'shift', 'size', 'skew', 'slice_shift',
        'squeeze', 'stack', 'std', 'sub', 'subtract', 'sum', 'swapaxes',
        'swaplevel', 'tail', 'transform', 'transpose', 'truediv', 'truncate',
        'tshift', 'style', 'take', 'unstack', 'var', 'where',
        # Pandas
        'DataFrame', 'array', 'concat', 'cut', 'date_range', 'datetime',
        'factorize', 'interval_range', 'isna', 'isnull', 'melt', 'merge',
        'notna', 'notnull', 'np', 'period_range', 'pivot', 'pivot_table',
        'plotting', 'unique', 'value_counts', 'wide_to_long'
    )

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
                is_good = self.__checkKeyword(keyw)
                if not is_good:
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
        is_lib = keyw in self.__allowedlibs
        # is_math = hasattr(math, keyw)
        # is_pd = hasattr(pd.DataFrame, keyw)
        # is_np = hasattr(np, keyw)
        is_allowedname = keyw in self.allowed_names
        is_allow = keyw in self.__allowed

        # is_good = is_name or is_ifel or is_math or is_pd or is_np or is_allow
        is_good = is_name or is_ifel or is_lib or is_allowedname or is_allow
        return is_good
