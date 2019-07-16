#!/usr/bin/env python3
"""
Table Compute

This script is an all-in-one reproduction of the table_compute.R counterpart
which was re-written in Python due to sandboxing concerns.
"""

__version__ = "0.2"


import csv
import math
from configparser import ConfigParser
from sys import argv

import numpy as np
import pandas as pd
from safety import Safety
from utils import Utils

if len(argv) == 1 or argv[1] == "--version":
    print(__version__)
    exit(-1)


# Math is imported but not directly used because users
# may specify a "math.<function>" when inserting a custom
# function. To remove linting errors, which break Travis
# we will just use an arbitrary math statement here.
__ = math.log

conf = argv[1]

config = ConfigParser()
config.read(conf)

wean = Utils.parseType
strict_xml = wean(config["Internal"]["function_defs"])

# Set decimal precision
pd.options.display.precision = wean(config["Default"]["precision"])

user_mode = wean(config["Default"]["user_mode"])
user_mode_single = None
out_table = None

if user_mode == "single":
    # Read in TSV file
    data = pd.read_csv(
        wean(config["SingleTableOps"]["reader_file"]),
        header=wean(config["SingleTableOps"]["reader_header"]),
        index_col=wean(config["SingleTableOps"]["reader_row_col"]),
        keep_default_na=wean(config["Default"]["narm"]),
        sep='\t'
    )
    user_mode_single = wean(config["SingleTableOps"]["user_mode_single"])

    if user_mode_single == "precision":
        # Useful for changing decimal precision on write out
        out_table = data

    elif user_mode_single == "select":
        sto_set = config["SingleTableOps.SELECT"]

        # do not use duplicates?
        unique_col = not wean(sto_set["select_cols_unique"])
        unique_row = not wean(sto_set["select_rows_unique"])

        unique_col_wanted = wean(sto_set["select_cols_wanted"])
        unique_row_wanted = wean(sto_set["select_rows_wanted"])

        # Select all indexes if True
        if unique_col_wanted is True:
            unique_col_wanted = list(range(len(data.columns)))
        if unique_row_wanted is True:
            unique_row_wanted = list(range(len(data)))

        if unique_col:
            tab = unique_col_wanted
            unique_col_wanted = [x for i, x in enumerate(tab)
                                 if x not in tab[:i]]
        if unique_row:
            tab = unique_row_wanted
            unique_row_wanted = [x for i, x in enumerate(tab)
                                 if x not in tab[:i]]

        out_table = data.iloc[unique_row_wanted, unique_col_wanted]

    elif user_mode_single == "filtersumval":
        sto_set = config["SingleTableOps.FILTERSUMVAL"]

        mode = wean(sto_set["filtersumval_mode"])
        axis = wean(sto_set["filtersumval_axis"])
        operation = wean(sto_set["filtersumval_op"])
        compare_operation = wean(sto_set["filtersumval_compare"])
        value = wean(sto_set["filtersumval_against"])
        minmatch = wean(sto_set["filtersumval_minmatch"])

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

            op = Utils.getTwoValueBaseOp(operation)
            bool_mat = op(data, value)
            axis_bool = np.sum(bool_mat, axis=axis) >= minmatch

        out_table = data.loc[:, axis_bool] if axis == 0 else data.loc[axis_bool, :]

    elif user_mode_single == "matrixapply":
        sto_set = config["SingleTableOps.MATRIXAPPLY"]
        # 0 - column, 1 - row
        axis = wean(sto_set["matrixapply_dimension"])
        # sd, mean, max, min, sum, median, summary
        operation = wean(sto_set["matrixapply_op"])

        if operation is None:
            use_custom = wean(sto_set["matrixapply_custom"])
            if use_custom is True:
                custom_func = wean(sto_set["matrixapply_custom_func"])

                def fun(vec):
                    """Dummy Function"""
                    return vec

                ss = Safety("vec", custom_func, strict_xml)
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
        sto_set = config["SingleTableOps.ELEMENT"]

        # Add None values for missing keys
        for co in ('value', 'mode', 'replace', 'modify_op',
                   'scale_op', 'scale_value', 'customop'):
            if 'element_' + co not in sto_set.keys():
                sto_set["element_" + co] = 'None'

        # lt, gt, ge, etc.
        operation = wean(sto_set["element_op"])
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
            value = wean(sto_set["element_value"])
            # Otherwise we subset the data to be
            # replaced using T/F values
            bool_mat = op(data, value)

        # Get the main processing mode
        mode = wean(sto_set["element_mode"])

        if mode == "replace":
            replacement_val = wean(sto_set["element_replace"])
            out_table[bool_mat] = replacement_val

        elif mode == "modify":
            mod_op = Utils.getOneValueMathOp(wean(sto_set["element_modify_op"]))
            out_table[bool_mat] = out_table[bool_mat].applymap(mod_op)

        elif mode == "scale":
            scale_op = Utils.getTwoValueBaseOp(wean(sto_set["element_scale_op"]))
            scale_value = wean(sto_set["element_scale_value"])
            out_table[bool_mat] = scale_op(out_table[bool_mat], scale_value)

        elif mode == "custom":
            element_customop = wean(sto_set["element_customop"])

            def fun(elem):
                """Dummy Function"""
                return elem

            ss = Safety("elem", element_customop, strict_xml)
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table[bool_mat] = out_table[bool_mat].applymap(fun)

        else:
            print("No such element mode!", mode)
            exit(-1)

    elif user_mode_single == "fulltable":
        general_mode = wean(config["SingleTableOps.FULLTABLE"]["mode"])

        if general_mode == "melt":
            sto_set = config["SingleTableOps.MELT"]
            melt_ids = wean(sto_set["melt_ids"])
            melt_values = wean(sto_set["melt_values"])



        elif general_mode == "pivot":
            sto_set = config["SingleTableOps.PIVOT"]
            pivot_index = wean(sto_set["pivot_index"])
            pivot_column = wean(sto_set["pivot_column"])
            pivot_values = wean(sto_set["pivot_values"])
        elif general_mode == "custom":
            sto_set = config["SingleTableOps.FULLTABLE_CUSTOM"]
            custom_func = wean(sto_set["fulltable_customop"])

            def fun(tableau):
                """Dummy Function"""
                return tableau

            ss = Safety("table", custom_func, strict_xml)
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table = fun(data)

    # elif user_mode_single == "sort":
    #     sto_set = config["SingleTableOps.SELECT"]
    else:
        print("No such mode!", user_mode_single)
        exit(-1)


elif user_mode == "multiple":

    table_sections = list(filter(
        lambda x: x.startswith("MultipleTableOps.TABLE"), config.sections()
    ))

    if not table_sections:
        print("Multiple table sets not given!")
        exit(-1)

    reader_skip = wean(config["Default"]["reader_skip"])

    table = ["Dummy Table"]

    # Read and populate tables
    for sect in table_sections:
        tmp = pd.read_csv(
            wean(config[sect]["file"]),
            header=wean(config[sect]["header"]),
            index_col=wean(config[sect]["row_names"]),
            keep_default_na=wean(config["Default"]["narm"]),
            sep="\t"
        )
        table.append(tmp)

    custom_op = wean(config["MultipleTableOps"]["fulltable_customop"])
    table_names = list(map(lambda x: "table" + str(x), range(1, len(table))))
    table_real_names = list(map(lambda x: "table[" + str(x) + "]", range(1, len(table))))

    ss = Safety(table_names, custom_op, strict_xml)
    fun_string = ss.generateFunction()
    # Change the argument to table
    fun_string = fun_string.replace("fun(table1):", "fun():")
    # table1 to table[1]
    for i, tab in enumerate(table_names):
        oldtab = tab
        newtab = table_real_names[i]
        fun_string = fun_string.replace(oldtab, newtab, 1)

    fun_string = fun_string.replace("fun():", "fun(table):")
    exec(fun_string)  # SUPER DUPER SAFE...
    out_table = fun(table)

else:
    print("No such mode!", user_mode)
    exit(-1)

out_parameters = {
    "sep": "\t",
    "float_format": "%%.%df" % pd.options.display.precision,
    "header": wean(config["Default"]["out_headers_col"]),
    "index": wean(config["Default"]["out_headers_row"])
}
if user_mode_single not in ('matrixapply', None):
    out_parameters["quoting"] = csv.QUOTE_NONE

out_table.to_csv(wean(config["Default"]["outtable"]), **out_parameters)
