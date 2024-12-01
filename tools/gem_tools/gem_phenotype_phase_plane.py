import argparse
import cobra
import pandas as pd
import numpy as np

def __main__():
    parser = argparse.ArgumentParser(
        prog = "ExpectedGrowthRate",
        description = "This program calculates the expected growth rate of a GEM",
        epilog = "Adding an epilog, but doubt it's needed.",
    )
    parser.add_argument(
        "-m", "--cb_model_location", dest="cb_model_location", action="store", type=str, default=None, required=True, help="The model to use."
    )
    parser.add_argument(
        "-output_csv", "--output_csv", dest="out_file_csv", action="store", type=str, default=None, required=True, help="The output csv file name."
    )
    parser.add_argument(
        "-r1", "--reaction1", dest="reaction1", action="store", type=str, default=None, required=True, help="The first reaction to scan."
    )
    parser.add_argument(
        "-r2", "--reaction2", dest="reaction2", action="store", type=str, default=None, required=True, help="The second reaction to scan."
    )
    parser.add_argument(
        "-p", "--points", dest="points", action="store", type=int, default=10, required=False, help="The number of points to scan."
    )
    parser.add_argument(
        "-c", "--objective", dest="objective", action="store", type=str, default=None, required=False, help="The reaction to use as objective."
    )
    parser.add_argument(
        "-u", "--uptake_constraints_file", dest="uptake_constraints_file", action="store", type=str, default=None, required=False, help="The file containing the uptake constraints."
    )

    args = parser.parse_args()
    try:
        assert(len(vars(args)) == 7)
    except:
        raise Exception("{} arguments were received. 7 were expected.".format(len(vars(args))))

    try:
        cb_model = cobra.io.read_sbml_model(args.cb_model_location)
    except:
        raise Exception("The model could not be read. Ensure it is in correct SBML format.")

    # set the uptake constraints if provided
    if args.uptake_constraints_file is not None and args.uptake_constraints_file != "None":
        constraints_df = pd.read_csv(args.uptake_constraints_file, sep=";", header=0, index_col=False)
        for index, row in constraints_df.iterrows():
            cb_model.reactions.get_by_id(row["reaction_id"]).lower_bound = row["lower_bound"]
            cb_model.reactions.get_by_id(row["reaction_id"]).upper_bound = row["upper_bound"]
    

    # get the reactions
    reactions = [args.reaction1, args.reaction2]

    # checking if reactions are in model
    for reaction in reactions:
        if reaction not in cb_model.reactions:
            raise Exception(f"Reaction {reaction} not found in model {args.cb_model_location.split('/')[-1]}")

    # get the points
    points = args.points
    if not isinstance(points, int):
        raise Exception("Points must be an integer")
    if points < 1:
        raise Exception("Must have at least one point in the phase plane")

    # perform phenotype phase plane analysis
    objective = None
    if args.objective is not None and args.objective != "None":
        objective = args.objective

    results = cobra.flux_analysis.phenotype_phase_plane.production_envelope(
        model = cb_model,
        reactions = reactions,
        points = points,
        objective=None,
    )

    # save the results
    results.to_csv(args.out_file_csv, sep=";", header=True, index=False)
    
if __name__ == "__main__":
    __main__()