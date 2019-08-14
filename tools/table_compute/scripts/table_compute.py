#!/usr/bin/env python3
"""
Table Compute

This script is an all-in-one reproduction of the table_compute.R counterpart
which was re-written in Python due to sandboxing concerns.
"""

__version__ = "0.7"


import csv
import math
import operator
from sys import argv

import numpy as np
import pandas as pd
import userconfig as uc
from safety import Safety

# This should be generated in the same directory

if len(argv) == 2 and argv[1] == "--version":
    print(__version__)
    exit(-1)


class Utils:
    @staticmethod
    def getOneValueMathOp(op_name):
        "Returns a simple one value math operator such as log, sqrt, etc"
        return getattr(math, op_name)

    @staticmethod
    def getTwoValueBaseOp(op_name):
        "Returns a basic two value operator such as +, -, etc"
        return getattr(operator, "__" + op_name + "__")

    @staticmethod
    def getVectorPandaOp(op_name):
        "Returns a valid DataFrame vector operator"
        return getattr(pd.DataFrame, op_name)

    @staticmethod
    def getTwoValuePandaOp(op_name):
        "Returns a valid two value DataFrame operator"
        return getattr(pd.DataFrame, "__" + op_name + "__")


# Math is imported but not directly used because users
# may specify a "math.<function>" when inserting a custom
# function. To remove linting errors, which break Travis
# we will just use an arbitrary math statement here.
__ = math.log

strict_xml = uc.Default["function_defs"]

# Set decimal precision
pd.options.display.precision = uc.Default["precision"]

user_mode = uc.Default["user_mode"]
user_mode_single = None
out_table = None

if user_mode == "single":
    # Read in TSV file
    data = pd.read_csv(
        uc.SingleTableOps["reader_file"],
        header=uc.SingleTableOps["reader_header"],
        index_col=uc.SingleTableOps["reader_row_col"],
        keep_default_na=uc.Default["narm"],
        sep='\t'
    )
    user_mode_single = uc.SingleTableOps["user_mode_single"]

    if user_mode_single == "precision":
        # Useful for changing decimal precision on write out
        out_table = data

    elif user_mode_single == "select":
        sto_set = uc.SingleTableOps["SELECT"]

        cols_specified = sto_set["select_cols_wanted"]
        rows_specified = sto_set["select_rows_wanted"]

        # Select all indexes if empty array of values
        if not cols_specified:
            cols_specified = list(range(len(data.columns)))
        if not rows_specified:
            rows_specified = list(range(len(data)))

        # do not use duplicate indexes
        # e.g. [2,3,2,5,5,4,2] to [2,3,5,4]
        nodupes_col = not sto_set["select_cols_unique"]
        nodupes_row = not sto_set["select_rows_unique"]

        if nodupes_col:
            tab = cols_specified
            cols_specified = [x for i, x in enumerate(tab)
                              if x not in tab[:i]]
        if nodupes_row:
            tab = rows_specified
            rows_specified = [x for i, x in enumerate(tab)
                              if x not in tab[:i]]

        out_table = data.iloc[rows_specified, cols_specified]

    elif user_mode_single == "filtersumval":
        sto_set = uc.SingleTableOps["FILTERSUMVAL"]

        mode = sto_set["filtersumval_mode"]
        axis = sto_set["filtersumval_axis"]
        operation = sto_set["filtersumval_op"]
        compare_operation = sto_set["filtersumval_compare"]
        value = sto_set["filtersumval_against"]
        minmatch = sto_set["filtersumval_minmatch"]

        if mode == "operation":
            # Perform axis operation
            summary_op = Utils.getVectorPandaOp(operation)
            axis_summary = summary_op(data, axis=axis)
            # Perform vector comparison
            compare_op = Utils.getTwoValueBaseOp(compare_operation)
            axis_bool = compare_op(axis_summary, value)

        elif mode == "element":
            if operation.startswith("str_"):
                data = data.astype("str")
                value = str(value)
                # Convert str_eq to eq
                operation = operation[4:]
            else:
                value = float(value)

            op = Utils.getTwoValueBaseOp(operation)
            bool_mat = op(data, value)
            axis_bool = np.sum(bool_mat, axis=axis) >= minmatch

        out_table = data.loc[:, axis_bool] if axis == 0 else data.loc[axis_bool, :]

    elif user_mode_single == "matrixapply":
        sto_set = uc.SingleTableOps["MATRIXAPPLY"]
        # 0 - column, 1 - row
        axis = sto_set["matrixapply_dimension"]
        # sd, mean, max, min, sum, median, summary
        operation = sto_set["matrixapply_op"]

        if operation is None:
            use_custom = sto_set["matrixapply_custom"]
            if use_custom is True:
                custom_func = sto_set["matrixapply_custom_func"]

                def fun(vec):
                    """Dummy Function"""
                    return vec

                ss = Safety(custom_func, ['vec'], 'pd.Series')
                fun_string = ss.generateFunction()
                exec(fun_string)  # SUPER DUPER SAFE...

                out_table = data.apply(fun, axis)

            else:
                print("No operation given")
                exit(-1)
        else:
            op = getattr(pd.DataFrame, operation)
            out_table = op(data, axis)

    elif user_mode_single == "element":
        sto_set = uc.SingleTableOps["ELEMENT"]

        # Add None values for missing keys
        # - we check if operation exists later
        for co in ('value', 'mode', 'replace', 'modify_op',
                   'scale_op', 'scale_value', 'customop'):
            if 'element_' + co not in sto_set.keys():
                sto_set["element_" + co] = 'None'

        # lt, gt, ge, etc.
        operation = sto_set["element_op"]
        if operation == "None":
            operation = None
        # Here we first create a boolean matrix of all the
        # elements we wish to replace
        #
        # First we define a filter matrix of identical size to
        # the original, full of True values
        out_table = data.copy()
        bool_mat = data.copy()
        bool_mat.loc[:, :] = True

        if operation:
            op = Utils.getTwoValuePandaOp(operation)
            value = sto_set["element_value"]
            try:
                # Could be numeric
                value = float(value)
            except ValueError:
                pass
            # Otherwise we subset the data to be
            # replaced using T/F values
            bool_mat = op(data, value)

        # Get the main processing mode
        mode = sto_set["element_mode"]

        if mode == "replace":
            replacement_val = sto_set["element_replace"]
            out_table[bool_mat] = replacement_val

        elif mode == "modify":
            mod_op = Utils.getOneValueMathOp(sto_set["element_modify_op"])
            out_table[bool_mat] = out_table[bool_mat].applymap(mod_op)

        elif mode == "scale":
            scale_op = Utils.getTwoValueBaseOp(sto_set["element_scale_op"])
            scale_value = sto_set["element_scale_value"]
            out_table[bool_mat] = scale_op(out_table[bool_mat], scale_value)

        elif mode == "custom":
            element_customop = sto_set["element_customop"]

            def fun(elem):
                """Dummy Function"""
                return elem

            ss = Safety(element_customop, ['elem'])
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table[bool_mat] = out_table[bool_mat].applymap(fun)

        else:
            print("No such element mode!", mode)
            exit(-1)

    elif user_mode_single == "fulltable":

        general_mode = uc.SingleTableOps["FULLTABLE"]["mode"]

        if general_mode == "melt":
            sto_set = uc.SingleTableOps["FULLTABLE"]["MELT"]
            melt_ids = sto_set["melt_ids"]
            melt_values = sto_set["melt_values"]

            out_table = pd.melt(data, id_vars=melt_ids, value_vars=melt_values)

        elif general_mode == "pivot":
            sto_set = uc.SingleTableOps["FULLTABLE"]["PIVOT"]
            pivot_index = sto_set["pivot_index"]
            pivot_column = sto_set["pivot_column"]
            pivot_values = sto_set["pivot_values"]

            out_table = data.pivot(index=pivot_index, columns=pivot_column, values=pivot_values)

        elif general_mode == "custom":
            sto_set = uc.SingleTableOps["FULLTABLE"]["FULLTABLE_CUSTOM"]
            custom_func = sto_set["fulltable_customop"]

            def fun(tableau):
                """Dummy Function"""
                return tableau

            ss = Safety(custom_func, ['table'], 'pd.DataFrame')
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table = fun(data)

    # elif user_mode_single == "sort":
    #     sto_set = uc.SingleTableOps["SELECT"]
    else:
        print("No such mode!", user_mode_single)
        exit(-1)


elif user_mode == "multiple":

    table_sections = uc.MultipleTableOps["TABLES"]

    if not table_sections:
        print("Multiple table sets not given!")
        exit(-1)

    reader_skip = uc.Default["reader_skip"]

    # Data
    table = [None]
    # Handlers for users "table1", "table2", etc.
    table_names = []
    # Actual references "table[1]", "table[2]", etc.
    # - as a result we need a DummyTable to fill the table[0] slot
    table_names_real = []

    # Read and populate tables
    for x, t_sect in enumerate(table_sections):
        if x == 0:
            # Skip padding table
            continue

        tmp = pd.read_csv(
            t_sect["file"],
            header=t_sect["header"],
            index_col=t_sect["row_names"],
            keep_default_na=uc.Default["narm"],
            sep="\t"
        )
        table.append(tmp)
        table_names.append("table" + str(x))
        table_names_real.append("table[" + str(x) + "]")

    custom_op = uc.MultipleTableOps["fulltable_customop"]
    ss = Safety(custom_op, table_names, 'pd.DataFrame')
    fun_string = ss.generateFunction()
    # Change the argument to table
    fun_string = fun_string.replace("fun(table1):", "fun():")
    # table1 to table[1]
    for name, name_real in zip(table_names, table_names_real):
        fun_string = fun_string.replace(name, name_real)

    fun_string = fun_string.replace("fun():", "fun(table):")
    exec(fun_string)  # SUPER DUPER SAFE...
    out_table = fun(table)

else:
    print("No such mode!", user_mode)
    exit(-1)

if not isinstance(out_table, pd.DataFrame):
    print('The specified operation did not result in a table to return.')
    raise RuntimeError(
        'The operation did not result in a pd.DataFrame to return.'
    )
out_parameters = {
    "sep": "\t",
    "float_format": "%%.%df" % pd.options.display.precision,
    "header": uc.Default["out_headers_col"],
    "index": uc.Default["out_headers_row"]
}
if user_mode_single not in ('matrixapply', None):
    out_parameters["quoting"] = csv.QUOTE_NONE

out_table.to_csv(uc.Default["outtable"], **out_parameters)
