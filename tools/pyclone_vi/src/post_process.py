import h5py
import numpy as np
import pandas as pd


def fix_cluster_ids(df):
    """Renumber clusters so they are sequential from 0."""
    cluster_map = dict(
        zip(df["cluster_id"].unique(), np.arange(df["cluster_id"].nunique()))
    )

    df["cluster_id"] = df["cluster_id"].map(cluster_map)

    return df


def load_results_df(file_name):
    cluster_df = load_cluster_df(file_name)

    loci_df = load_loci_df(file_name)

    df = pd.merge(cluster_df, loci_df, on=["cluster_id"])

    df = df[
        [
            "mutation_id",
            "sample_id",
            "cluster_id",
            "cellular_prevalence",
            "cellular_prevalence_std",
            "cluster_assignment_prob",
        ]
    ]

    df = df.sort_values(by=["cluster_id", "mutation_id", "sample_id"])

    return df


def load_cluster_df(file_name):
    with h5py.File(file_name, "r") as fh:
        samples = fh["/data/samples"].asstr()[()]

        theta = fh["/var_params/theta"][()]

        z = fh["/var_params/z"][()]

    labels = np.argmax(z, axis=1)

    clusters = np.unique(labels)

    df = []

    x = np.linspace(0, 1, theta.shape[2])

    for i, sample_id in enumerate(samples):
        for cluster_id in clusters:
            q = theta[cluster_id, i]

            mean = np.sum(x * q)

            var = np.sum(x**2 * q) - mean**2

            std = np.sqrt(var)

            df.append(
                {
                    "sample_id": sample_id,
                    "cluster_id": cluster_id,
                    "cellular_prevalence": mean,
                    "cellular_prevalence_std": std,
                }
            )

    df = pd.DataFrame(df)

    df = df[
        ["sample_id", "cluster_id", "cellular_prevalence", "cellular_prevalence_std"]
    ]

    return df


def load_loci_df(file_name):
    with h5py.File(file_name, "r") as fh:
        mutations = fh["/data/mutations"].asstr()[()]

        z = fh["/var_params/z"][()]

    labels = np.argmax(z, axis=1)

    probs = np.max(z, axis=1)

    df = pd.DataFrame(
        {
            "mutation_id": mutations,
            "cluster_id": labels,
            "cluster_assignment_prob": probs,
        }
    )

    df = df[["mutation_id", "cluster_id", "cluster_assignment_prob"]]

    return df
