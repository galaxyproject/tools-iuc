import argparse

import cobra
import pandas as pd


def __main__():
    parser = argparse.ArgumentParser(
        prog="FluxDistribution",
        description="This program calculates the flux distribution of a GEM",
        epilog="Adding an epilog, but doubt it's needed.",
    )
    parser.add_argument(
        "-m",
        "--cb_model_location",
        dest="cb_model_location",
        action="store",
        type=str,
        default=None,
        required=True,
        help="The model to use."
    )
    parser.add_argument(
        "-output",
        "--output",
        dest="out_file",
        action="store",
        type=str,
        default=None,
        required=True,
        help="The output file."
    )
    parser.add_argument(
        "-u",
        "--uptake_constraints_file",
        dest="uptake_constraints_file",
        action="store",
        type=str,
        default=None,
        required=False,
        help="File containing new uptake constraits."
    )

    args = parser.parse_args()

    try:
        cb_model = cobra.io.read_sbml_model(args.cb_model_location)
    except Exception as e:
        raise Exception(
            "The model could not be read. Ensure "
            "it is in correct SBML format."
        ) from e

    if args.uptake_constraints_file is not None\
            and args.uptake_constraints_file != "None":
        constraints_df = pd.read_csv(
            args.uptake_constraints_file,
            sep=";",
            header=0,
            index_col=False
        )
        for _, row in constraints_df.iterrows():
            cb_model.reactions.get_by_id(
                row["reaction_id"]
            ).lower_bound = row["lower_bound"]
            cb_model.reactions.get_by_id(
                row["reaction_id"]
            ).upper_bound = row["upper_bound"]

    # do pFBA
    solution = cobra.flux_analysis.pfba(cb_model)

    # make a dataframe with the reaction names,
    # reaction ids, and flux distribution
    flux_distribution = pd.DataFrame(
        columns=["reaction_name", "reaction_id", "flux"]
    )

    flux_distribution["reaction_name"] = \
        [reaction.name for reaction in cb_model.reactions]
    flux_distribution["reaction_id"] = \
        [reaction.id for reaction in cb_model.reactions]
    flux_distribution["flux"] = \
        [solution.fluxes[reaction.id] for reaction in cb_model.reactions]

    flux_distribution.to_csv(args.out_file, sep=";", index=False)


if __name__ == "__main__":
    __main__()
