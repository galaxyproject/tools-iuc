import argparse

import cobra
import pandas as pd
from escher import Builder


def __main__():
    parser = argparse.ArgumentParser(
        prog="EscherVisualization",
        description="This program visualizes an Escher map",
        epilog="Adding an epilog, but doubt it's needed.",
    )
    parser.add_argument(
        "-m",
        "--cb_model_location",
        dest="cb_model_location",
        action="store",
        type=str,
        default=None,
        required=False,
        help="The model to use."
    )
    parser.add_argument(
        "-f",
        "--flux_distribution_location",
        dest="flux_distribution_location",
        action="store",
        type=str,
        default=None,
        required=False,
        help="The flux distribution to visualize."
    )
    parser.add_argument(
        "-e",
        "--expect_map",
        dest="expect_map",
        action="store",
        type=str,
        default=None,
        required=True,
        help="Is a map expected to be uploaded?"
    )
    parser.add_argument(
        "-l",
        "--model_to_download",
        dest="model_to_download",
        action="store",
        type=str,
        default=None,
        required=False,
        help="The model to download."
    )
    parser.add_argument(
        "--map_load_name",
        dest="map_load_name",
        action="store",
        type=str,
        default=None,
        required=False,
        help="The name of the map to use."
    )
    parser.add_argument(
        "--map_upload_name",
        dest="map_upload_name",
        action="store",
        type=str,
        default=None,
        required=False,
        help="The name of the map to use."
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

    args = parser.parse_args()

    if args.expect_map not in ["True", "False"]:
        raise Exception("The expect_map argument must be either True or False.")
    if args.expect_map == "True" and args.map_load_name is None and \
            args.map_upload_name is None:
        raise Exception(
            "You must specify a map name if a map is expected to be uploaded."
        )

    cb_model = None
    model_name = None
    if args.model_to_download is not None and args.model_to_download != "None":
        if args.cb_model_location is not None \
                and args.cb_model_location != "None":
            raise Exception(
                "You cannot specify both a model to "
                "download and a model to use."
            )
        model_name = args.model_to_download
    elif args.cb_model_location is not None\
            and args.cb_model_location != "None":
        try:
            cb_model = cobra.io.read_sbml_model(args.cb_model_location)
        except Exception as e:
            raise Exception(
                "The model could not be read. "
                "Ensure it is in correct SBML format."
            ) from e

    map_name = None
    map_location = None
    if args.map_upload_name is not None and args.map_upload_name != "None":
        if args.map_load_name is not None and args.map_load_name != "None":
            raise Exception(
                "You cannot specify both a map to upload and a map to load."
            )
        map_location = args.map_upload_name
    elif args.map_load_name is not None and args.map_load_name != "None":
        map_name = args.map_load_name

    if args.uptake_constraints_file is not None and \
            args.uptake_constraints_file != "None":
        if cb_model is None:
            raise Exception(
                "You cannot specify uptake constraints "
                "without uploading a model."
            )
        else:
            constraints_df = pd.read_csv(
                args.uptake_constraints_file,
                sep=";",
                header=0,
                index_col=False
            )
            for index, row in constraints_df.iterrows():
                rxn_id = row["reaction_id"]
                cb_model.reactions.get_by_id(rxn_id).lower_bound = \
                    row["lower_bound"]
                cb_model.reactions.get_by_id(rxn_id).upper_bound = \
                    row["upper_bound"]

    flux_dict = None
    if args.flux_distribution_location is not None and \
            args.flux_distribution_location != "None":
        if cb_model is None:
            raise Exception(
                "You cannot specify a flux distribution "
                "without uploading a model."
            )
        if args.uptake_constraints_file is not None and \
                args.uptake_constraints_file != "None":
            raise Exception(
                "You cannot specify both uptake constraints and a flux "
                "distribution."
            )
        try:
            flux_df = pd.read_csv(
                args.flux_distribution_location,
                sep=";",
                header=0,
                index_col=False
            )
            flux_dict = {
                key: value for key, value in zip(
                    flux_df['reaction_name'],
                    flux_df['flux']
                )
            }
        except Exception as e:
            raise Exception(
                "The flux distribution file could not be read. "
                "Ensure the file has semicolon-separated "
                "columns and a header row."
            ) from e

    if cb_model is not None and flux_dict is None:
        solution = cobra.flux_analysis.pfba(cb_model)

        # make a dataframe with the reaction names, reaction ids, and flux
        flux_distribution = pd.DataFrame(
            columns=["reaction_name", "reaction_id", "flux"]
        )
        flux_distribution["reaction_name"] = [
            reaction.name for reaction in cb_model.reactions
        ]
        flux_distribution["reaction_id"] = [
            reaction.id for reaction in cb_model.reactions
        ]
        flux_distribution["flux"] = [
            solution.fluxes[reaction.id] for reaction in cb_model.reactions
        ]
        flux_dict = {
            key: value for key, value in zip(
                flux_distribution['reaction_name'],
                flux_distribution['flux']
            )
        }

    builder = Builder()
    if map_name is not None:
        builder.map_name = map_name
        print("Downloading map...")
    if map_location is not None:
        builder.map_json = map_location
        print("Uploading map...")
    if model_name is not None:
        builder.model_name = model_name
        print("Downloading model...")
    if cb_model is not None:
        builder.model = cb_model
        print("Uploading model...")

    if flux_dict is not None:
        builder.reaction_data = flux_dict

    builder.save_html(args.out_file)


if __name__ == "__main__":
    __main__()
