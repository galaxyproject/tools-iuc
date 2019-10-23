import re


class Safety():
    """
    Class to safely evaluate mathematical expression on single
    or table data
    """

    __allowed_tokens = (
        '(', ')', 'if', 'else', 'or', 'and', 'not', 'in',
        '+', '-', '*', '/', '%', ',', '!=', '==', '>', '>=', '<', '<=',
        'min', 'max', 'sum',
    )
    __allowed_ref_types = {
        'pd.DataFrame': {
            'abs', 'add', 'agg', 'aggregate', 'align', 'all', 'any', 'append',
            'apply', 'applymap', 'as_matrix', 'asfreq', 'at', 'axes', 'bool',
            'clip', 'clip_lower', 'clip_upper', 'columns', 'combine',
            'compound', 'corr', 'count', 'cov', 'cummax', 'cummin', 'cumprod',
            'cumsum', 'describe', 'div', 'divide', 'dot', 'drop',
            'drop_duplicates', 'droplevel', 'dropna', 'duplicated', 'empty',
            'eq', 'equals', 'expanding', 'ffill', 'fillna', 'filter', 'first',
            'first_valid_index', 'floordiv', 'ge', 'groupby', 'gt', 'head',
            'iat', 'iloc', 'index', 'insert', 'interpolate', 'isin', 'isna',
            'isnull', 'items', 'iteritems', 'iterrows', 'itertuples', 'ix',
            'join', 'keys', 'kurt', 'kurtosis', 'last', 'last_valid_index',
            'le', 'loc', 'lookup', 'lt', 'mad', 'mask', 'max', 'mean',
            'median', 'melt', 'merge', 'min', 'mod', 'mode', 'mul', 'multiply',
            'ndim', 'ne', 'nlargest', 'notna', 'notnull', 'nsmallest',
            'nunique', 'pct_change', 'pivot', 'pivot_table', 'pop', 'pow',
            'prod', 'product', 'quantile', 'radd', 'rank', 'rdiv', 'replace',
            'resample', 'rfloordiv', 'rmod', 'rmul', 'rolling', 'round',
            'rpow', 'rsub', 'rtruediv', 'sample', 'select',
            'sem', 'shape', 'shift', 'size', 'skew', 'slice_shift',
            'squeeze', 'stack', 'std', 'sub', 'subtract', 'sum', 'swapaxes',
            'swaplevel', 'T', 'tail', 'take', 'transform', 'transpose',
            'truediv', 'truncate', 'tshift', 'unstack', 'var', 'where',
        },
        'pd.Series': {
            'abs', 'add', 'agg', 'aggregate', 'align', 'all', 'any', 'append',
            'apply', 'argsort', 'as_matrix', 'asfreq', 'asof', 'astype', 'at',
            'at_time', 'autocorr', 'axes', 'between', 'between_time', 'bfill',
            'bool', 'cat', 'clip', 'clip_lower', 'clip_upper', 'combine',
            'combine_first', 'compound', 'corr', 'count', 'cov', 'cummax',
            'cummin', 'cumprod', 'cumsum', 'describe', 'diff', 'div', 'divide',
            'divmod', 'dot', 'drop', 'drop_duplicates', 'droplevel', 'dropna',
            'dt', 'dtype', 'dtypes', 'duplicated', 'empty', 'eq', 'equals',
            'ewm', 'expanding', 'factorize', 'ffill', 'fillna', 'filter',
            'first', 'first_valid_index', 'flags', 'floordiv', 'ge', 'groupby',
            'gt', 'hasnans', 'head', 'iat', 'idxmax', 'idxmin', 'iloc', 'imag',
            'index', 'interpolate', 'is_monotonic', 'is_monotonic_decreasing',
            'is_monotonic_increasing', 'is_unique', 'isin', 'isna', 'isnull',
            'item', 'items', 'iteritems', 'ix', 'keys', 'kurt', 'kurtosis',
            'last', 'last_valid_index', 'le', 'loc', 'lt', 'mad', 'map',
            'mask', 'max', 'mean', 'median', 'min', 'mod', 'mode', 'mul',
            'multiply', 'name', 'ndim', 'ne', 'nlargest', 'nonzero', 'notna',
            'notnull', 'nsmallest', 'nunique', 'pct_change', 'pop', 'pow',
            'prod', 'product', 'ptp', 'quantile', 'radd', 'rank', 'rdiv',
            'rdivmod', 'real', 'repeat', 'replace', 'resample', 'rfloordiv',
            'rmod', 'rmul', 'rolling', 'round', 'rpow', 'rsub', 'rtruediv',
            'sample', 'searchsorted', 'select', 'sem', 'shape', 'shift',
            'size', 'skew', 'slice_shift', 'sort_index', 'sort_values',
            'squeeze', 'std', 'sub', 'subtract', 'sum', 'swapaxes',
            'swaplevel', 'T', 'tail', 'take', 'transform', 'transpose',
            'truediv', 'truncate', 'tshift', 'unique', 'unstack',
            'value_counts', 'var', 'where', 'xs',
        },
    }

    __allowed_qualified = {
        # allowed numpy functionality
        'np': {
            'abs', 'add', 'all', 'any', 'append', 'array', 'bool', 'ceil',
            'complex', 'cos', 'cosh', 'cov', 'cumprod', 'cumsum', 'degrees',
            'divide', 'divmod', 'dot', 'e', 'empty', 'exp', 'float', 'floor',
            'hypot', 'inf', 'int', 'isfinite', 'isin', 'isinf', 'isnan', 'log',
            'log10', 'log2', 'max', 'mean', 'median', 'min', 'mod', 'multiply',
            'nan', 'ndim', 'pi', 'product', 'quantile', 'radians', 'rank',
            'remainder', 'round', 'sin', 'sinh', 'size', 'sqrt', 'squeeze',
            'stack', 'std', 'str', 'subtract', 'sum', 'swapaxes', 'take',
            'tan', 'tanh', 'transpose', 'unique', 'var', 'where',
        },
        # allowed math functionality
        'math': {
            'acos', 'acosh', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 'ceil',
            'copysign', 'cos', 'cosh', 'degrees', 'e', 'erf', 'erfc', 'exp',
            'expm1', 'fabs', 'factorial', 'floor', 'fmod', 'frexp', 'fsum',
            'gamma', 'gcd', 'hypot', 'inf', 'isclose', 'isfinite', 'isinf',
            'isnan', 'ldexp', 'lgamma', 'log', 'log10', 'log1p', 'log2',
            'modf', 'nan', 'pi', 'pow', 'radians', 'remainder', 'sin', 'sinh',
            'sqrt', 'tan', 'tanh', 'tau', 'trunc',
        },
        # allowed pd functionality
        'pd': {
            'DataFrame', 'array', 'concat', 'cut', 'date_range', 'factorize',
            'interval_range', 'isna', 'isnull', 'melt', 'merge', 'notna',
            'notnull', 'period_range', 'pivot', 'pivot_table', 'unique',
            'value_counts', 'wide_to_long',
        },
    }

    def __init__(self, expression,
                 ref_whitelist=None, ref_type=None,
                 custom_qualified=None):
        self.allowed_qualified = self.__allowed_qualified.copy()
        if ref_whitelist is None:
            self.these = []
        else:
            self.these = ref_whitelist
            if ref_type is None or ref_type not in self.__allowed_ref_types:
                self.allowed_qualified['_this'] = set()
            else:
                self.allowed_qualified[
                    '_this'
                ] = self.__allowed_ref_types[ref_type]
        if custom_qualified is not None:
            self.allowed_qualified.update(custom_qualified)
        self.expr = expression
        self.__assertSafe()

    def generateFunction(self):
        "Generates a function to be evaluated outside the class"
        cust_fun = "def fun(%s):\n\treturn(%s)" % (self.these[0], self.expr)
        return cust_fun

    def __assertSafe(self):
        indeed, problematic_token = self.__isSafeStatement()
        if not indeed:
            self.detailedExcuse(problematic_token)
            raise ValueError("Custom Expression is not safe.")

    @staticmethod
    def detailedExcuse(word):
        "Gives a verbose statement for why users should not use some specific operators."
        mess = None
        if word == "for":
            mess = "for loops and comprehensions are not allowed. Use numpy or pandas table operations instead."
        elif word == ":":
            mess = "Colons are not allowed. Use inline Python if/else statements."
        elif word == "=":
            mess = "Variable assignment is not allowed. Use object methods to substitute values."
        elif word in ("[", "]"):
            mess = "Direct indexing of arrays is not allowed. Use numpy or pandas functions/methods to address specific parts of tables."
        else:
            mess = "Not an allowed token in this operation"
        print("( '%s' ) %s" % (word, mess))

    def __isSafeStatement(self):
        """
        Determines if a user-expression is safe to evaluate.

        To be considered safe an expression may contain only:
        - standard Python operators and numbers
        - inline conditional expressions
        - select functions and objects
          by default, these come from the math, numpy and pandas
          libraries, and must be qualified with the modules' conventional
          names math, np, pd; can be overridden at the instance level
        - references to a whitelist of objects (pd.DataFrames by default)
          and their methods
        """

        safe = True
        # examples of user-expressions
        # '-math.log(1 - elem/4096) * 4096 if elem != 1 else elem - 0.5'
        # 'vec.median() +  vec.sum()'

        # 1. Break expressions into tokens
        # e.g.,
        # [
        #     '-', 'math.log', '(', '1', '-', 'elem', '/', '4096', ')', '*',
        #     '4096', 'if', 'elem', '!=', '1', 'else', 'elem', '-', '0.5'
        # ]
        # or
        # ['vec.median', '(', ')', '+', 'vec.sum', '(', ')']
        tokens = [
            e for e in re.split(
                r'([a-zA-Z0-9_.]+|[^a-zA-Z0-9_.() ]+|[()])', self.expr
            ) if e.strip()
        ]

        # 2. Subtract allowed standard tokens
        rem = [e for e in tokens if e not in self.__allowed_tokens]

        # 3. Subtract allowed qualified objects from allowed modules
        #    and whitelisted references and their attributes
        rem2 = []
        for e in rem:
            parts = e.split('.')
            if len(parts) == 1:
                if parts[0] in self.these:
                    continue
            if len(parts) == 2:
                if parts[0] in self.these:
                    parts[0] = '_this'
                if parts[0] in self.allowed_qualified:
                    if parts[1] in self.allowed_qualified[parts[0]]:
                        continue
            rem2.append(e)

        # 4. Assert that rest are real numbers or strings
        e = ''
        for e in rem2:
            try:
                _ = float(e)
            except ValueError:
                safe = False
                break

        return safe, e
