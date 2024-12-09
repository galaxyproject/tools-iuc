import argparse

import cobra


def read_model(model_location):
    model = cobra.io.read_sbml_model(model_location)
    return model


def get_exchange_reactions_info(model):
    exchange_reactions = model.exchanges
    exchange_reactions_info = []
    for reaction in exchange_reactions:
        exchange_reactions_info.append(
            [
                reaction.id,
                reaction.name,
                reaction.reaction,
                reaction.lower_bound,
                reaction.upper_bound
            ])
    txt_object = (
        "reaction_id;reaction_name;reaction_stoichiometry;"
        "lower_bound;upper_bound\n"
    )
    for reaction in exchange_reactions_info:
        txt_object += ";".join([str(x) for x in reaction]) + "\n"
    return txt_object


def __main__():

    # Parsing arguments
    parser = argparse.ArgumentParser(
        prog="GEM ",
        description="This program retrieves the exchange fluxes "
        "of a GEM model to be used in Galaxy.",
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
    args = parser.parse_args()

    # Reading model from file
    try:
        cb_model = read_model(args.cb_model_location)
    except Exception as e:
        raise Exception(
            "The model could not be read. Ensure it is in correct SBML format."
        ) from e

    # Getting exchange reactions info
    answer = get_exchange_reactions_info(cb_model)

    # Writing exchange reactions info to file
    with open(args.out_file, "w") as outfile:
        outfile.write(str(answer))


if __name__ == "__main__":
    __main__()
