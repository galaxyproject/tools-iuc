import argparse
import cobra
import pandas as pd

def __main__():
    parser = argparse.ArgumentParser(
        prog = "FluxVariabilityAnalysis",
        description = "This program performs flux variability analysis on a GEM",
        epilog = "Adding an epilog, but doubt it's needed.",
    )
    parser.add_argument(
        "-m", "--cb_model_location", dest="cb_model_location", action="store", type=str, default=None, required=True, help="The model to use."
    )
    parser.add_argument(
        "-output", "--output", dest="out_file", action="store", type=str, default=None, required=True, help="The output file."
    )
    parser.add_argument(
        "-f", "--fraction", dest="fraction_of_optimum", action="store", type=float, default=None, required=True, help="The fraction of optimum the FVA solutions should come within."
    )
    parser.add_argument(
        "-u", "--uptake_constraints_file", dest="uptake_constraints_file", action="store", type=str, default=None, required=False, help="File containing new uptake constraits."
    )

    args = parser.parse_args()

    cb_model = cobra.io.read_sbml_model(args.cb_model_location)

    if args.uptake_constraints_file is not None and args.uptake_constraints_file != "None":
        constraints_df = pd.read_csv(args.uptake_constraints_file, sep=";", header=0, index_col=False)
        for index, row in constraints_df.iterrows():
            cb_model.reactions.get_by_id(row["reaction_id"]).lower_bound = row["lower_bound"]
            cb_model.reactions.get_by_id(row["reaction_id"]).upper_bound = row["upper_bound"]

    fraction_of_optimum = args.fraction_of_optimum

    # perform fva
    fva_result = cobra.flux_analysis.flux_variability_analysis(cb_model, fraction_of_optimum=fraction_of_optimum)

    # add reaction names and ids to the dataframe
    fva_result["reaction_id"] = fva_result.index
    fva_result["reaction_name"] = fva_result["reaction_id"].apply(lambda x: cb_model.reactions.get_by_id(x).name)

    # reorder the columns
    fva_result = fva_result[["reaction_id", "reaction_name", "minimum", "maximum"]]
    
    fva_result.to_csv(args.out_file, sep=";", index=False, header=True)
    
if __name__ == "__main__":
    __main__()