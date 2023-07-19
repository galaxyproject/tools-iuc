import logging
import os
from pathlib import Path
from typing import Dict, Literal, Optional, Set

from constants import molecular_types_for_metabologram

import helpers

from hydra.core.hydra_config import HydraConfig

from omegaconf import DictConfig, ListConfig

import pandas as pd

from pydantic import BaseModel as PydanticBaseModel


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


logger = logging.getLogger(__name__)


class DatasetConfig(BaseModel):
    label: str
    name: str
    subfolder: str
    metadata: str
    # We allow for some default values for the following files
    # will be ignored if they do not exist

    # First condition is the reference condition (control)
    # Conditions should be a subset of the medata corresponding column
    conditions: ListConfig
    abundances: str = "AbundanceCorrected"
    mean_enrichment: str = "MeanEnrichment13C"
    isotopologue_proportions: str = "IsotopologuesProportions"
    isotopologues: str = "Isotopologues"  # isotopologue absolute values

    def build(self) -> "Dataset":
        return Dataset(config=self)


class Dataset(BaseModel):
    config: DatasetConfig
    raw_data_folder: str = None
    processed_data_folder: str = None
    sub_folder_absolute: str = None
    metadata_df: Optional[pd.DataFrame] = None
    abundances_df: Optional[pd.DataFrame] = None
    mean_enrichment_df: Optional[pd.DataFrame] = None
    isotopologue_proportions_df: Optional[pd.DataFrame] = None
    isotopologues_df: Optional[pd.DataFrame] = None
    available_datasets: Set[
        Literal[
            "metadata", "abundances", "mean_enrichment",
            "isotopologue_proportions", "isotopologues"]
    ] = set()
    compartmentalized_dfs: Dict[str, Dict[str, pd.DataFrame]] = {}

    def preload(self):
        # check if we have a relative or absolute path, compute the absolute
        # path then load the data using pandas
        # if the path is relative, we assume it is relative to the original
        # CWD (befre hydra changed it)
        # store the data in self.metadata
        original_cwd = HydraConfig.get().runtime.cwd
        logger.info("Current config directory is %s", original_cwd)
        if not self.config.subfolder.startswith("/"):
            self.sub_folder_absolute = os.path.join(Path(original_cwd),
                                                    "data",
                                                    self.config.subfolder)
            logger.info("looking for data in %s", self.sub_folder_absolute)
        else:
            self.sub_folder_absolute = self.config.subfolder
        self.raw_data_folder = os.path.join(self.sub_folder_absolute, "raw")
        self.processed_data_folder = os.path.join(self.sub_folder_absolute,
                                                  "processed")
        # start loading the dataframes
        file_paths = [
            ("metadata", os.path.join(self.raw_data_folder,
                                      self.config.metadata + ".csv")),
            ("abundances", os.path.join(self.raw_data_folder,
                                        self.config.abundances + ".csv")),
            ("mean_enrichment", os.path.join(
                self.raw_data_folder,
                self.config.mean_enrichment + ".csv")),
            ("isotopologue_proportions", os.path.join(
                self.raw_data_folder,
                self.config.isotopologue_proportions + ".csv")),
            ("isotopologues", os.path.join(
                self.raw_data_folder, self.config.isotopologues + ".csv")),
        ]
        dfs = []
        for label, file_path in file_paths:
            try:
                if label != "metadata":
                    # the quantifications dfs take first column as index
                    # (metabolites), regardless the name of that column
                    dfs.append(pd.read_csv(file_path, sep="\t",
                                           header=0, index_col=0))
                else:
                    dfs.append(pd.read_csv(file_path, sep="\t", header=0))
                self.available_datasets.add(label)
            except FileNotFoundError:
                logger.critical(
                    "File %s not found, continuing, "
                    "but this might fail miserably",
                    file_path)
                dfs.append(None)
            except Exception as e:
                logger.error(
                    "Failed to load file %s during preload, aborting",
                    file_path)
                raise e

        (
            self.metadata_df,
            self.abundances_df,
            self.mean_enrichment_df,
            self.isotopologue_proportions_df,
            self.isotopologues_df,
        ) = dfs

        # log the first 5 rows of the metadata
        logger.info("Loaded metadata: \n%s", self.metadata_df.head())
        logger.info(
            "Finished loading raw dataset %s, available dataframes are : %s",
            self.config.label, self.available_datasets
        )
        self.check_expectations()

    def check_expectations(self):
        # conditions should be a subset of the metadata corresponding column
        if not set(self.config.conditions).issubset(
                set(self.metadata_df["condition"].unique())):
            logger.error("Conditions %s are not a subset of "
                         "the metadata declared conditions",
                         self.config.conditions)
            raise ValueError(
                f"Conditions {self.config.conditions} are not a subset "
                f"of the metadata declared conditions"
            )

        helpers.verify_metadata_sample_not_duplicated(self.metadata_df)

    def split_datafiles_by_compartment(self) -> None:
        frames_dict = {}
        for data_file_label in self.available_datasets:
            if 'metadata' in data_file_label:
                continue
            dataframe_label = data_file_label + "_df"  # TODO: this is fragile!
            tmp_co_dict = helpers.df_to_dict_by_compartment(
                getattr(self, dataframe_label),
                self.metadata_df)  # split by compartment
            frames_dict[data_file_label] = tmp_co_dict

        frames_dict = helpers.drop_all_nan_metabolites_on_comp_frames(
            frames_dict, self.metadata_df)
        frames_dict = helpers.set_samples_names(frames_dict, self.metadata_df)
        self.compartmentalized_dfs = frames_dict

    def save_datafiles_split_by_compartment(self) -> None:
        os.makedirs(self.processed_data_folder, exist_ok=True)
        out_data_path = self.processed_data_folder
        for file_name in self.compartmentalized_dfs.keys():
            for compartment in self.compartmentalized_dfs[file_name].keys():
                df = self.compartmentalized_dfs[file_name][compartment]
                tmp_file_name = self.get_file_for_label(file_name)
                output_file_name = f"{tmp_file_name}-{compartment}.csv"
                df.to_csv(os.path.join(out_data_path, output_file_name),
                          sep="\t", header=True, index=False)
                logger.info(
                    f"Saved the {compartment} compartment version "
                    f"of {file_name} in {out_data_path}")

    def get_file_for_label(self, label):
        if label == "abundances":
            return self.config.abundances
        elif label == "mean_enrichment":
            return self.config.mean_enrichment
        elif label == "isotopologue_proportions":
            return self.config.isotopologue_proportions
        elif label == "isotopologues":
            return self.config.isotopologues
        else:
            raise ValueError(f"Unknown label {label}")


class DataIntegrationConfig(DatasetConfig):
    transcripts: ListConfig
    pathways: DictConfig


class DataIntegration(Dataset):
    config: DataIntegrationConfig
    integration_files_folder_absolute: str = None  # absolute path directory
    deg_dfs: Dict[int, pd.DataFrame] = {}
    pathways_dfs: Dict[str, pd.DataFrame] = {}

    def set_dataset_integration_config(self):
        self.preload()
        self.split_datafiles_by_compartment()
        self.save_datafiles_split_by_compartment()

        self.integration_files_folder_absolute = os.path.join(
            self.sub_folder_absolute, "integration_files")

        self.check_expectations()  # of the Dataset class
        self.check_expectations_integration_data()  # of this child class

    def check_expectations_integration_data(self):
        if not len(set(self.config.transcripts)) == \
               len(self.config.transcripts):
            logger.error(
                f"Duplicated names: {self.config.transcripts}"
                f" in transcripts")
            raise ValueError(
                f"Duplicated names: {self.config.transcripts}"
                f" in transcripts")

        if not len(set(molecular_types_for_metabologram).difference(
                set(self.config.pathways.keys())
        )) == 0:
            logger.error("Unrecognized pathways configuration"
                         " in dataset yaml file")
            raise ValueError(
                "Unrecognized pathways_file_names configuration"
                " in dataset yaml file"
            )

    def load_deg_dfs(self):
        # generate dictionary of transcripts dataframes (DEGs) :
        # the keys are integers, with the order of files in the dataset yml
        for i, file_name in enumerate(self.config.transcripts):
            try:
                path_deg_file = os.path.join(
                    self.integration_files_folder_absolute,
                    f"{file_name}.csv")
                deg_df = pd.read_csv(path_deg_file, sep='\t', header=0)
                self.deg_dfs[i] = deg_df
            except FileNotFoundError:
                logger.info(f"{file_name}.csv: file not found")
            except Exception as e:
                logger.info(f'Error while opening file {file_name}.csv {e}')

        logger.info("Finished loading transcripts dataframes: "
                    "%s", self.config.transcripts)

    def load_pathways_dfs(self):
        for k in self.config.pathways.keys():
            try:
                path_file = os.path.join(
                    self.integration_files_folder_absolute,
                    f"{self.config.pathways[k]}.csv")
                pathway_df = pd.read_csv(path_file, sep='\t', header=0)
                self.pathways_dfs[k] = pathway_df
            except FileNotFoundError:
                logger.info(f"{self.config.pathways[k]}.csv: file not found")
            except Exception as e:
                logger.info(f'Error while opening file '
                            f'{self.config.pathways[k]}.csv {e}')

        logger.info("Finished loading pathways dataframes: "
                    "%s", self.config.pathways)

    def get_names_transcripts_files(self) -> Dict[int, str]:
        out_dict = {}
        for i, name in enumerate(self.config.transcripts):
            out_dict[i] = name
        return out_dict
