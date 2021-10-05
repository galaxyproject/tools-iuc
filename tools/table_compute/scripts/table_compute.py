#!/usr/bin/env python3
"""
Table Compute tool - a wrapper around pandas with parameter input validation.
"""


__version__ = "0.9.2"

import csv
import math
from sys import argv

import numpy as np
import pandas as pd
from safety import Safety

if len(argv) == 2 and argv[1] == "--version":
    print(__version__)
    exit(-1)

# The import below should be generated in the same directory as
# the table_compute.py script.
# It is placed here so that the --version switch does not fail
import userconfig as uc  # noqa: I100,I202


class Utils:
    @staticmethod
    def getOneValueMathOp(op_name):
        "Returns a simple one value math operator such as log, sqrt, etc"
        return getattr(math, op_name)

    @staticmethod
    def getVectorPandaOp(op_name):
        "Returns a valid DataFrame vector operator"
        return getattr(pd.DataFrame, op_name)

    @staticmethod
    def getTwoValuePandaOp(op_name, pd_obj):
        "Returns a valid two value DataFrame or Series operator"
        return getattr(type(pd_obj), "__" + op_name + "__")

    @staticmethod
    def readcsv(filedict, narm):
        data = pd.read_csv(
            filedict["file"],
            header=filedict["header"],
            index_col=filedict["row_names"],
            keep_default_na=narm,
            nrows=filedict["nrows"],
            skipfooter=filedict["skipfooter"],
            skip_blank_lines=filedict["skip_blank_lines"],
            sep='\t'
        )
        # Fix whitespace issues in index or column names
        data.columns = [col.strip() if type(col) is str else col
                        for col in data.columns]
        data.index = [row.strip() if type(row) is str else row
                      for row in data.index]
        return(data)

    @staticmethod
    def rangemaker(tab):
        # e.g. "1:3,2:-2" specifies "1,2,3,2,1,0,-1,-2" to give [0,1,2,1,0,-1,-2]
        # Positive indices are decremented by 1 to reference 0-base numbering
        # Negative indices are unaltered, so that -1 refers to the last column
        out = []
        err_mess = None
        for ranges in tab.split(","):
            nums = ranges.split(":")
            if len(nums) == 1:
                numb = int(nums[0])
                # Positive numbers get decremented.
                # i.e. column "3" refers to index 2
                #      column "-1" still refers to index -1
                if numb != 0:
                    out.append(numb if (numb < 0) else (numb - 1))
                else:
                    err_mess = "Please do not use 0 as an index"
            elif len(nums) == 2:
                left, right = map(int, nums)
                if 0 in (left, right):
                    err_mess = "Please do not use 0 as an index"
                elif left < right:
                    if left > 0:  # and right > 0 too
                        # 1:3 to 0,1,2
                        out.extend(range(left - 1, right))
                    elif right < 0:  # and left < 0 too
                        # -3:-1 to -3,-2,-1
                        out.extend(range(left, right + 1))
                    elif left < 0 and right > 0:
                        # -2:2 to -2,-1,0,1
                        out.extend(range(left, 0))
                        out.extend(range(0, right))
                elif right < left:
                    if right > 0:  # and left > 0
                        # 3:1 to 2,1,0
                        out.extend(range(left - 1, right - 2, -1))
                    elif left < 0:  # and right < 0
                        # -1:-3 to -1,-2,-3
                        out.extend(range(left, right - 1, -1))
                    elif right < 0 and left > 0:
                        # 2:-2 to 1,0,-1,-2
                        out.extend(range(left - 1, right - 1, -1))
                else:
                    err_mess = "%s should not be equal or contain a zero" % nums
            if err_mess:
                print(err_mess)
                return(None)
        return(out)


# Set decimal precision
pd.options.display.precision = uc.Default["precision"]

user_mode = uc.Default["user_mode"]
user_mode_single = None
out_table = None
params = uc.Data["params"]

if user_mode == "single":
    # Read in TSV file
    data = Utils.readcsv(uc.Data["tables"][0], uc.Default["narm"])
    user_mode_single = params["user_mode_single"]

    if user_mode_single == "precision":
        # Useful for changing decimal precision on write out
        out_table = data

    elif user_mode_single == "select":
        cols_specified = params["select_cols_wanted"]
        rows_specified = params["select_rows_wanted"]

        # Select all indexes if empty array of values
        if cols_specified:
            cols_specified = Utils.rangemaker(cols_specified)
        else:
            cols_specified = range(len(data.columns))
        if rows_specified:
            rows_specified = Utils.rangemaker(rows_specified)
        else:
            rows_specified = range(len(data))

        # do not use duplicate indexes
        # e.g. [2,3,2,5,5,4,2] to [2,3,5,4]
        nodupes_col = not params["select_cols_unique"]
        nodupes_row = not params["select_rows_unique"]

        if nodupes_col:
            cols_specified = [x for i, x in enumerate(cols_specified)
                              if x not in cols_specified[:i]]
        if nodupes_row:
            rows_specified = [x for i, x in enumerate(rows_specified)
                              if x not in rows_specified[:i]]

        out_table = data.iloc[rows_specified, cols_specified]

    elif user_mode_single == "filtersumval":
        mode = params["filtersumval_mode"]
        axis = params["filtersumval_axis"]
        operation = params["filtersumval_op"]
        compare_operation = params["filtersumval_compare"]
        value = params["filtersumval_against"]
        minmatch = params["filtersumval_minmatch"]

        if mode == "operation":
            # Perform axis operation
            summary_op = Utils.getVectorPandaOp(operation)
            axis_summary = summary_op(data, axis=axis)
            # Perform vector comparison
            compare_op = Utils.getTwoValuePandaOp(
                compare_operation, axis_summary
            )
            axis_bool = compare_op(axis_summary, value)

        elif mode == "element":
            if operation.startswith("str_"):
                data = data.astype("str")
                value = str(value)
                # Convert str_eq to eq
                operation = operation[4:]
            else:
                value = float(value)

            op = Utils.getTwoValuePandaOp(operation, data)
            bool_mat = op(data, value)
            axis_bool = np.sum(bool_mat, axis=axis) >= minmatch

        out_table = data.loc[:, axis_bool] if axis == 0 else data.loc[axis_bool, :]

    elif user_mode_single == "matrixapply":
        # 0 - column, 1 - row
        axis = params["matrixapply_dimension"]
        # sd, mean, max, min, sum, median, summary
        operation = params["matrixapply_op"]

        if operation is None:
            use_custom = params["matrixapply_custom"]
            if use_custom:
                custom_func = params["matrixapply_custom_func"]

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
        # lt, gt, ge, etc.
        operation = params["element_op"]
        bool_mat = None
        if operation is not None:
            if operation == "rowcol":
                # Select all indexes if empty array of values
                if "element_cols" in params:
                    cols_specified = Utils.rangemaker(params["element_cols"])
                else:
                    cols_specified = range(len(data.columns))
                if "element_rows" in params:
                    rows_specified = Utils.rangemaker(params["element_rows"])
                else:
                    rows_specified = range(len(data))

                # Inclusive selection:
                # - True: Giving a row or column will match all elements in that row or column
                # - False: Give a row or column will match only elements in both those rows or columns
                inclusive = params["element_inclusive"]

                # Create a bool matrix (intialised to False) with selected
                # rows and columns set to True
                bool_mat = data.copy()
                bool_mat[:] = False
                if inclusive:
                    bool_mat.iloc[rows_specified, :] = True
                    bool_mat.iloc[:, cols_specified] = True
                else:
                    bool_mat.iloc[rows_specified, cols_specified] = True

            else:
                op = Utils.getTwoValuePandaOp(operation, data)
                value = params["element_value"]
                try:
                    # Could be numeric
                    value = float(value)
                except ValueError:
                    pass
                # generate filter matrix of True/False values
                bool_mat = op(data, value)
        else:
            # implement no filtering through a filter matrix filled with
            # True values.
            bool_mat = np.full(data.shape, True)

        # Get the main processing mode
        mode = params["element_mode"]
        if mode == "replace":
            replacement_val = params["element_replace"]
            out_table = data.mask(
                bool_mat,
                data.where(bool_mat).applymap(
                    lambda x: replacement_val.format(elem=x)
                )
            )
        elif mode == "modify":
            mod_op = Utils.getOneValueMathOp(params["element_modify_op"])
            out_table = data.mask(
                bool_mat, data.where(bool_mat).applymap(mod_op)
            )
        elif mode == "scale":
            scale_op = Utils.getTwoValuePandaOp(
                params["element_scale_op"], data
            )
            scale_value = params["element_scale_value"]
            out_table = data.mask(
                bool_mat, scale_op(data.where(bool_mat), scale_value)
            )
        elif mode == "custom":
            element_customop = params["element_customop"]

            def fun(elem):
                """Dummy Function"""
                return elem

            ss = Safety(element_customop, ['elem'])
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table = data.mask(
                bool_mat, data.where(bool_mat).applymap(fun)
            )
        else:
            print("No such element mode!", mode)
            exit(-1)

    elif user_mode_single == "fulltable":
        general_mode = params["mode"]

        if general_mode == "transpose":
            out_table = data.T
        elif general_mode == "melt":
            melt_ids = params["MELT"]["melt_ids"]
            melt_values = params["MELT"]["melt_values"]

            out_table = pd.melt(data, id_vars=melt_ids, value_vars=melt_values)
        elif general_mode == "pivot":
            pivot_index = params["PIVOT"]["pivot_index"]
            pivot_column = params["PIVOT"]["pivot_columns"]
            pivot_values = params["PIVOT"]["pivot_values"]
            pivot_aggfunc = params["PIVOT"]["pivot_aggfunc"]

            if not(pivot_aggfunc):
                out_table = data.pivot(
                    index=pivot_index, columns=pivot_column, values=pivot_values
                )
            else:
                out_table = data.pivot_table(
                    index=pivot_index, columns=pivot_column, values=pivot_values,
                    aggfunc=pivot_aggfunc
                )

        elif general_mode == "custom":
            custom_func = params["fulltable_customop"]

            def fun(tableau):
                """Dummy Function"""
                return tableau

            ss = Safety(custom_func, ['table'], 'pd.DataFrame')
            fun_string = ss.generateFunction()
            exec(fun_string)  # SUPER DUPER SAFE...

            out_table = fun(data)

    else:
        print("No such mode!", user_mode_single)
        exit(-1)


elif user_mode == "multiple":

    table_sections = uc.Data["tables"]

    if not table_sections:
        print("Multiple table sets not given!")
        exit(-1)

    reader_skip = uc.Default["reader_skip"]

    # Data
    table = []
    # 1-based handlers for users "table1", "table2", etc.
    table_names = []
    # Actual 0-based references "table[0]", "table[1]", etc.
    table_names_real = []

    # Read and populate tables
    for x, t_sect in enumerate(table_sections):
        tmp = Utils.readcsv(t_sect, uc.Default["narm"])
        table.append(tmp)
        table_names.append("table" + str(x + 1))
        table_names_real.append("table[" + str(x) + "]")

    custom_op = params["fulltable_customop"]
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

if not isinstance(out_table, (pd.DataFrame, pd.Series)):
    print('The specified operation did not result in a table to return.')
    raise RuntimeError(
        'The operation did not generate a pd.DataFrame or pd.Series to return.'
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
