import click

import run


@click.command(context_settings={"max_content_width": 120}, name="fit")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help="""Path to TSV format file with copy number and allele count information for all samples. """
    """See the examples directory in the GitHub repository for format.""",
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(resolve_path=True),
    help="""Path to where results will be written in HDF5 format.""",
)
@click.option(
    "-a",
    "--num-annealing-steps",
    default=1,
    type=int,
    help="""Number of simulated annealing steps to use. """
    """Default is one step i.e. not to use simulated annealing.""",
)
@click.option(
    "-c",
    "--num-clusters",
    default=10,
    type=int,
    help="""Number of clusters to use in variational approximation distribution. """
    """Note that not all clusters may not be assigned data points, so the final number of clusters could be lower. """
    """Default is 10.""",
)
@click.option(
    "-d",
    "--density",
    default="binomial",
    type=click.Choice(["beta-binomial", "binomial"]),
    help="""Allele count density in the PyClone model. Use beta-binomial for high coverage sequencing. """
    """Default binomial.""",
)
@click.option(
    "-g",
    "--num-grid-points",
    default=100,
    type=int,
    help="""Number of points used to approximate CCF values. Default is 100.""",
)
@click.option(
    "-r",
    "--num-restarts",
    default=1,
    type=int,
    help="""Number of random restarts of variational inference. Default is 1.""",
)
@click.option(
    "-t",
    "--num-threads",
    default=1,
    type=int,
    help="""Number of threads to use. Default is 1.""",
)
@click.option(
    "--annealing-power",
    default=1.0,
    type=float,
    help="""Exponent of entries in the annealing ladder.""" """Default is 1.0.""",
)
@click.option(
    "--convergence-threshold",
    default=1e-6,
    type=float,
    help="""Maximum relative ELBO difference between iterations to decide on convergence. """
    """Default is 10^-6.""",
)
@click.option(
    "--max-iters",
    default=int(1e4),
    type=int,
    help="""Maximum number of ELBO optimization iterations."""
    """Default is 10,0000.""",
)
@click.option(
    "--mix-weight-prior",
    default=1.0,
    type=float,
    help="""Parameter value of symmetric Dirichlet prior distribution on mixture weights. Higher values will produce more clusters. """
    """Default is 1.0 which is the uniform prior.""",
)
@click.option(
    "--precision",
    default=200,
    type=float,
    help="""Precision for Beta-Binomial density. Has no effect when using Binomial. """
    """Default is 200.""",
)
@click.option(
    "--print-freq",
    default=100,
    type=int,
    help="""How often to print information about optimization. """
    """Default is every 100 iteration.""",
)
@click.option(
    "--seed",
    default=None,
    type=int,
    help="""Set random seed so results can be reproduced. """
    """By default a random seed is chosen.""",
)
def fit(**kwargs):
    """Fit PyClone-VI model to data."""
    run.fit(**kwargs)


@click.command(context_settings={"max_content_width": 120}, name="write-results-file")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help="""Path to HDF5 format file produced by the `fit` command.""",
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(resolve_path=True),
    help="""Path to where results will be written in tsv format.""",
)
@click.option(
    "-c",
    "--compress",
    is_flag=True,
    help="""If set the output file will be compressed using gzip.""",
)
def write_results_file(**kwargs):
    """Write the results of a fitted model to file."""
    run.write_results_file(**kwargs)


@click.group(name="pyclone_vi-vi")
def main():
    pass


main.add_command(fit)
main.add_command(write_results_file)
