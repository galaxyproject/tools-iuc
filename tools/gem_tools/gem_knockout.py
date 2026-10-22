import argparse

import cobra
import pandas as pd


def __main__():
    parser = argparse.ArgumentParser(
        prog="FluxDistribution",
        description="Performs FBA knockout analysis on a GEM.",
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
        "-k",
        "--knockout_type",
        dest="knockout_type",
        action="store",
        type=str,
        default="single",
        required=False,
        help="Type of knockout to perform, single or double"
    )
    parser.add_argument(
        "-g",
        "--gene_knockouts",
        dest="gene_knockouts",
        action="store",
        type=str, default=None,
        required=False,
        help="List of genes to knock out. Defaults to all."
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

    # Reading the model
    try:
        cb_model = cobra.io.read_sbml_model(args.cb_model_location)
    except Exception as e:
        raise Exception(
            "The model could not be read. "
            "Ensure it is in correct SBML format."
        ) from e

    # Verifying the genes are present in the model
    gene_ids = [gene.id for gene in cb_model.genes]

    genes_to_knockout_1 = args.gene_knockouts.split(',')\
        if args.gene_knockouts is not None else []
    gene_bool = [
        True if gene in gene_ids else False for gene in genes_to_knockout_1
    ]
    if not all(gene_bool):
        print(
            f'Found {sum(gene_bool)} of {len(genes_to_knockout_1)} genes '
            'in the model.'
        )
        raise Exception(
            "One or more of the genes to knockout are not present "
            "in the model."
        )

    # Adding all genes to knockout if none are specified
    if genes_to_knockout_1 is None or len(genes_to_knockout_1) == 0:
        genes_to_knockout_1 = [gene.id for gene in cb_model.genes]
    # Applying uptake constraints
    if (args.uptake_constraints_file is not None
            and args.uptake_constraints_file != "None"):
        constraints_df = pd.read_csv(
            args.uptake_constraints_file,
            sep=";",
            header=0,
            index_col=False
        )
        for index, row in constraints_df.iterrows():
            reaction = cb_model.reactions.get_by_id(row["reaction_id"])
            reaction.lower_bound = row["lower_bound"]
            reaction.upper_bound = row["upper_bound"]

    result = pd.DataFrame(columns=[
        "reaction_id", "ko_gene_id_1", "ko_gene_id_2",
        "reaction", "wildtype_flux", "knockout_flux"
    ])

    if args.knockout_type == "single":
        genes_to_knockout_2 = [0]
    elif args.knockout_type == "double":
        genes_to_knockout_2 = genes_to_knockout_1.copy()
    else:
        raise Exception(
            f"Invalid knockout type {args.knockout_type}. "
            "Only single and double are allowed."
        )

    # Wildtype pFBA
    with cb_model as model:
        wildtype_solution = model.optimize()

    # Performing gene knockouts
    for gene1 in genes_to_knockout_1:
        for gene2 in genes_to_knockout_2:
            with cb_model as model:
                model.genes.get_by_id(gene1).knock_out()
                if args.knockout_type == "double":
                    model.genes.get_by_id(gene2).knock_out()
                solution = model.optimize()
                for reaction in model.reactions:
                    result = pd.concat([result, pd.DataFrame([{
                        "reaction_id": reaction.id,
                        "ko_gene_id_1": gene1,
                        "ko_gene_id_2": gene2
                        if args.knockout_type == "double" else None,
                        "reaction": reaction.reaction,
                        "wildtype_flux": wildtype_solution.fluxes[reaction.id],
                        "knockout_flux": solution.fluxes[reaction.id],
                    }])], ignore_index=True)

    # Writing the results to file
    result.to_csv(args.out_file, sep=";", index=False)


if __name__ == "__main__":
    __main__()
