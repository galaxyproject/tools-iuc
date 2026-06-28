#!/usr/bin/env python3

# Written by Thiseas C. Lamnidis and released under the MIT license.
# See git repository (https://github.com/nf-core/eager) for full license text.

import json
import re
import sys
from collections import OrderedDict

jsonOut = OrderedDict()
data = OrderedDict()


def make_float(x):
    """Convert elements to float, else leave as is."""
    output = [None for i in range(len(x))]
    # If value for an estimate/error is -nan, replace with "NA".
    # JSON does not accept NaN as a valid field.
    for i in range(len(x)):
        if x[i] == "-nan" or x[i] == "nan":
            output[i] = "N/A"
            continue
        try:
            output[i] = float(x[i])
        except (ValueError, TypeError):
            output[i] = x[i]

    return tuple(output)


input_files = sys.argv[1:]

with open("nuclear_contamination.txt", 'w') as output:
    print(
        "Individual", "Num_SNPs", "Method1_MOM_estimate", "Method1_MOM_SE",
        "Method1_ML_estimate", "Method1_ML_SE", "Method2_MOM_estimate",
        "Method2_MOM_SE", "Method2_ML_estimate", "Method2_ML_SE",
        sep="\t", file=output
    )
    for fn in input_files:
        # Reset values for each file to "N/A"
        mom1, err_mom1 = "N/A", "N/A"
        ml1, err_ml1 = "N/A", "N/A"
        mom2, err_mom2 = "N/A", "N/A"
        ml2, err_ml2 = "N/A", "N/A"
        nSNPs = "0"
        with open(fn, 'r') as f:
            ind = re.sub(r'\.X.contamination.out$', '', fn).split("/")[-1]
            for line in f:
                fields = line.strip().split()
                if not fields:
                    continue
                if line.strip().startswith("We have nSNP sites:"):
                    nSNPs = fields[4].rstrip(",")
                elif line.startswith("Method1") and 'new_llh' in fields[1]:
                    mom1 = fields[3].split(":")[1]
                    err_mom1 = fields[4].split(":")[1]
                    ml1 = fields[5].split(":")[1]
                    err_ml1 = fields[6].split(":")[1]
                    if err_ml1.endswith("contamination"):
                        err_ml1 = err_ml1[:-13]
                elif line.startswith("Method2") and 'new_llh' in fields[1]:
                    mom2 = fields[3].split(":")[1]
                    err_mom2 = fields[4].split(":")[1]
                    ml2 = fields[5].split(":")[1]
                    err_ml2 = fields[6].split(":")[1]

            # Convert estimates and errors to floating point numbers
            estimates = make_float((
                ml1, err_ml1, mom1, err_mom1, ml2, err_ml2, mom2, err_mom2
            ))
            (ml1, err_ml1, mom1, err_mom1, ml2, err_ml2, mom2, err_mom2) = estimates
            data[ind] = {
                "Num_SNPs": int(nSNPs),
                "Method1_MOM_estimate": mom1,
                "Method1_MOM_SE": err_mom1,
                "Method1_ML_estimate": ml1,
                "Method1_ML_SE": err_ml1,
                "Method2_MOM_estimate": mom2,
                "Method2_MOM_SE": err_mom2,
                "Method2_ML_estimate": ml2,
                "Method2_ML_SE": err_ml2
            }
            print(
                ind, nSNPs, mom1, err_mom1, ml1, err_ml1,
                mom2, err_mom2, ml2, err_ml2,
                sep="\t", file=output
            )


jsonOut = {
    "plot_type": "generalstats",
    "id": "nuclear_contamination",
    "pconfig": {
        "Num_SNPs": {"title": "Number of SNPs"},
        "Method1_MOM_estimate": {"title": "Contamination Estimate (Method1_MOM)"},
        "Method1_MOM_SE": {"title": "Estimate Error (Method1_MOM)"},
        "Method1_ML_estimate": {"title": "Contamination Estimate (Method1_ML)"},
        "Method1_ML_SE": {"title": "Estimate Error (Method1_ML)"},
        "Method2_MOM_estimate": {"title": "Contamination Estimate (Method2_MOM)"},
        "Method2_MOM_SE": {"title": "Estimate Error (Method2_MOM)"},
        "Method2_ML_estimate": {"title": "Contamination Estimate (Method2_ML)"},
        "Method2_ML_SE": {"title": "Estimate Error (Method2_ML)"}
    },
    "data": data
}

with open('nuclear_contamination_mqc.json', 'w') as outfile:
    json.dump(jsonOut, outfile)
