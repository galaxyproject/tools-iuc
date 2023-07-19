import logging
import os
import sys
from typing import Union

from constants import (
    assert_literal,
    availtest_methods,
    columns_transcripts_config_keys,
    data_files_keys_type,
    data_types_suitable_for_metabologram,
    metabolites_values_for_metabologram,
)

from data import DataIntegration, Dataset

from helpers import flatten

import hydra

from omegaconf import DictConfig, ListConfig, OmegaConf, open_dict
from omegaconf.errors import ConfigAttributeError

from processing.differential_analysis import (
    differential_comparison,
    multi_group_compairson,
    time_course_analysis
)
from processing.pca_analysis import run_pca_analysis

from pydantic import BaseModel as PydanticBaseModel

from visualization.abundance_bars import run_plot_abundance_bars
from visualization.distr_fit_plot import run_distr_fit_plot
from visualization.isotopologue_proportions import (
    run_isotopologue_proportions_plot)
from visualization.mean_enrichment_line_plot import (
    run_mean_enrichment_line_plot)
from visualization.metabologram import run_metabologram
from visualization.pca_plot import run_pca_plot


logger = logging.getLogger(__name__)


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class MethodConfig(BaseModel):
    label: str
    name: str

    def build(self) -> "Method":
        raise NotImplementedError


class AbundancePlotConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """

    barcolor: str = "timepoint"
    axisx: str = "condition"
    axisx_labeltilt: int = 20  # 0 is no tilt
    height_each_subfig: float = 5.4
    as_grid: Union[bool, None] = False

    def build(self) -> "AbundancePlot":
        return AbundancePlot(config=self)


class DifferentialAnalysisConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """

    grouping: ListConfig = ["condition", "timepoint"]
    qualityDistanceOverSpan: float
    correction_method: str = "bonferroni"
    impute_values: DictConfig

    def build(self) -> "DifferentialAnalysis":
        return DifferentialAnalysis(config=self)


class MultiGroupComparisonConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """

    grouping: ListConfig = ["condition", "timepoint"]
    correction_method: str = "bonferroni"
    impute_values: DictConfig

    def build(self) -> "MultiGroupComparison":
        return MultiGroupComparison(config=self)


class IsotopologueProportionsPlotConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """
    max_nb_carbons_possible: int = 24
    appearance_separated_time: bool = True
    separated_plots_by_condition: bool = False
    height_each_stack: float = 4.9
    sharey: bool = False  # share y axis across subplots
    x_ticks_text_size: int = 18
    y_ticks_text_size: int = 19
    as_grid: Union[bool, None] = False

    def build(self) -> "IsotopologueProportionsPlot":
        return IsotopologueProportionsPlot(config=self)


class MeanEnrichmentLinePlotConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """
    alpha: float = 1
    xaxis_title: str = "Time"
    color_lines_by: str = "condition"  # or  "metabolite"
    palette_condition: str = "paired"  # seaborn/matplotlib pals
    palette_metabolite: str = "auto_multi_color"  # or .csv path
    height_subplot: float = 4.3
    as_grid: Union[bool, None] = False

    def build(self) -> "MeanEnrichmentLinePlot":
        return MeanEnrichmentLinePlot(config=self)


class PcaAnalysisConfig(MethodConfig):
    pca_split_further: Union[ListConfig, None] = ["timepoint"]

    def build(self) -> "PcaAnalysis":
        return PcaAnalysis(config=self)


class PcaPlotConfig(MethodConfig):
    pca_split_further: Union[ListConfig, None] = ["timepoint"]
    draw_ellipses: Union[str, None] = "condition"
    run_iris_demo: bool = False

    def build(self) -> "PcaPlot":
        return PcaPlot(config=self)


class TimeCourseAnalysisConfig(MethodConfig):
    grouping: ListConfig = ["condition", "timepoint"]
    qualityDistanceOverSpan: float
    correction_method: str = "bonferroni"
    impute_values: DictConfig

    def build(self) -> "TimeCourseAnalysis":
        return TimeCourseAnalysis(config=self)


class DistrFitPlotConfig(MethodConfig):
    grouping: ListConfig = ["condition", "timepoint"]
    qualityDistanceOverSpan: float
    impute_values: DictConfig

    def build(self) -> "DistrFitPlot":
        return DistrFitPlot(config=self)


class MetabologramIntegrationConfig(MethodConfig):
    """
    Sets default values or fills them from the method yaml file
    """
    grouping: ListConfig = ["condition", "timepoint"]
    qualityDistanceOverSpan: float
    correction_method: str = "bonferroni"
    impute_values: DictConfig
    abs_values_scale_color_bar: DictConfig
    colors_divergent_palette: ListConfig = ["darkcyan", "white", "orangered"]
    fig_width: float = 7
    fig_height: float = 8
    edge_color: ListConfig = ['#000000', '#000000']
    line_width: ListConfig = [1, 1]
    font_size: float = 8
    display_label_and_value: bool = False

    def build(self) -> "MetabologramIntegration":
        return MetabologramIntegration(config=self)


class Method(BaseModel):
    config: MethodConfig

    def plot(self):
        logger.info("Will plot the method, with the following config: %s",
                    self.config)

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info("Not instantialted in the parent class.")
        raise NotImplementedError

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info("Not instantialted in the parent class.")
        raise NotImplementedError


class AbundancePlot(Method):
    config: AbundancePlotConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(
            "Will plot the abundance plot, with the following config: %s",
            self.config)
        if not (
                "metabolites" in cfg.analysis.keys()
        ):  # plotting for _all_ metabolites
            logger.warning(
                "No selected metabolites provided, plotting for all")
            with open_dict(cfg):
                cfg.analysis["metabolites"] = {}
                for c in set(dataset.metadata_df["short_comp"]):
                    metabolites_compartment = \
                        dataset.compartmentalized_dfs[
                            'abundances'][c].index.to_list()
                    cfg.analysis["metabolites"][c] = metabolites_compartment

        self.check_expectations(cfg, dataset)
        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        run_plot_abundance_bars(dataset, out_plot_dir, cfg)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        # check that necessary information is provided in the analysis config
        try:
            if not set(cfg.analysis.metabolites.keys()).issubset(
                    dataset.metadata_df["short_comp"]):
                raise ValueError(
                    "[Analysis > Metabolites > compartments] are missing "
                    "from [Metadata > Compartments]"
                )
            if not set(cfg.analysis.timepoints).issubset(
                    set(dataset.metadata_df["timepoint"])):
                raise ValueError(
                    "[Analysis > Timepoints] time points provided in the "
                    "config file are not present in [Metadata > timepoint]"
                )
            if cfg.analysis.width_each_subfig is not None:
                try:
                    float(cfg.analysis.width_each_subfig)
                except TypeError:
                    logger.error("Unrecognized value for width_each_subfig"
                                 " in the config file")
                assert cfg.analysis.width_each_subfig >= 0, logger.error(
                    "width_each_subfig must be a positive number"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency:{e}")
            sys.exit(1)


class DifferentialAnalysis(Method):
    config: DifferentialAnalysisConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(f"The current working directory is {os.getcwd()}")
        logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))
        logger.info(
            "Will perform differential analysis, with the following "
            "config: %s", self.config)
        out_table_dir = os.path.join(os.getcwd(), cfg.table_path)
        os.makedirs(out_table_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)
        for file_name, test in cfg.analysis.statistical_test.items():
            if test is None:
                continue
            tmp = dataset.get_file_for_label(file_name)  # current file name
            logger.info(
                f"Running differential analysis of {tmp} using {test} test")
            differential_comparison(file_name, dataset, cfg, test,
                                    out_table_dir=out_table_dir)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        # check that necessary information is provided in the analysis config
        try:
            if not (all(len(c) == 2 for c in cfg.analysis.comparisons)):
                raise ValueError(
                    "Number of conditions has to be 2 for a pairwise "
                    "comparison, see config file"
                )

            if not all(any(
                    item[0] in set(dataset.metadata_df["condition"]) for item
                    in sublist)
                       for sublist in cfg.analysis.comparisons):
                raise ValueError(
                    "Conditions provided for comparisons in the config file "
                    "are not present in the metadata file, aborting"
                )
            if not all(any(
                    item[1] in set(dataset.metadata_df["timepoint"]) for item
                    in sublist)
                       for sublist in cfg.analysis.comparisons):
                raise ValueError(
                    "Timepoints provided for comparisons in the config file "
                    "are not present in the metadata file, aborting"
                )
            # all values in the comparisons arrays of arrays are in the
            # metadata, either as conditions or timepoints
            conditions_and_tp = set(dataset.metadata_df["condition"]).union(
                set(dataset.metadata_df["timepoint"]))
            comparisons = set(flatten(cfg.analysis.comparisons))
            diff = comparisons.difference(conditions_and_tp)
            if not comparisons.issubset(conditions_and_tp):
                raise ValueError(
                    f"Comparisons > Conditions or timepoints provided in the "
                    f"config file {diff} are not present in the metadata "
                    f"file, aborting"
                )
            # comparison_mode is one of the constant values
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency:{e}")
            sys.exit(1)


class MultiGroupComparison(Method):
    config: MultiGroupComparisonConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(f"The current working directory is {os.getcwd()}")
        logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))
        logger.info(
            "Will perform multi group analysis, with the following "
            "config: %s", self.config)
        out_table_dir = os.path.join(os.getcwd(), cfg.table_path)
        os.makedirs(out_table_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)

        for file_name in cfg.analysis.datatypes:
            logger.info("Running multi group analysis of %s",
                        dataset.get_file_for_label(file_name))
            multi_group_compairson(file_name, dataset, cfg,
                                   out_table_dir=out_table_dir)

    # TODO: add expectations on the compartementalised dfs?
    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        try:
            [assert_literal(dt, data_files_keys_type, "datatype") for dt in
             cfg.analysis.datatypes]

            if not set(
                    [sublist[0] for sublist in
                        cfg.analysis.conditions]
            ).issubset(set(dataset.metadata_df["condition"])):
                raise ValueError(
                    "Conditions provided for comparisons in the config file "
                    "are not present in the metadata file, aborting"
                )
            if not set(
                    [sublist[1] for sublist in cfg.analysis.conditions]
            ).issubset(set(dataset.metadata_df["timepoint"])):
                raise ValueError(
                    "Timepoints provided for comparisons in the config file"
                    " are not present in the metadata file, aborting"
                )

        except ValueError as e:
            logger.error(f"Data inconsistency:{e}")


class IsotopologueProportionsPlot(Method):
    config: IsotopologueProportionsPlotConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(
            "Will perform isotopologue proportions stacked-bar plots, "
            "with the following config: %s", self.config)

        if not (
                "metabolites" in cfg.analysis.keys()):
            # plotting for _all_ metabolites
            logger.warning(
                "No selected metabolites provided, plotting for all")
            with open_dict(cfg):
                compartments = list(set(dataset.metadata_df["short_comp"]))
                cfg.analysis["metabolites"] = dict()
                for c in compartments:
                    isotopologues_names = \
                        dataset.compartmentalized_dfs[
                            "isotopologue_proportions"][c].index.to_list()
                    metabolites_c = set(
                        [i.split("_m+")[0] for i in isotopologues_names]
                    )
                    cfg.analysis["metabolites"][c] = list(metabolites_c)

        self.check_expectations(cfg, dataset)
        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        run_isotopologue_proportions_plot(dataset, out_plot_dir, cfg)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        try:
            if not set(cfg.analysis.metabolites.keys()).issubset(
                    dataset.metadata_df['short_comp']):
                raise ValueError(
                    "[Analysis > Metabolites > compartments] "
                    "are missing from [Metadata > Compartments]"
                )
            if not set(cfg.analysis.timepoints).issubset(
                    set(dataset.metadata_df["timepoint"])):
                raise ValueError(
                    "[Analysis > Timepoints] time points provided in the "
                    "config file are not present in [Metadata > timepoint]"
                )
            if not cfg.analysis.width_each_stack >= float(0.0):
                raise ValueError(
                    "width_each_stack in the config file "
                    "must be a positive number"
                )
            if not cfg.analysis.inner_numbers_size >= float(0.0):
                raise ValueError(
                    "inner_numbers_size in the config file"
                    " must be a positive number"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:"
                f"{e}, aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)


class MeanEnrichmentLinePlot(Method):
    config: MeanEnrichmentLinePlotConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info("Will perform Mean Enrichment (syn. Fractional "
                    "Contributions) line-plot "
                    "with the following config: %s", self.config)
        if not ("metabolites" in cfg.analysis.keys()):
            logger.warning(
                "No selected metabolites provided, plotting for all")
            with open_dict(cfg):
                cfg.analysis["metabolites"] = {}
                for c in set(dataset.metadata_df["short_comp"]):
                    cfg.analysis["metabolites"][c] = \
                        dataset.compartmentalized_dfs[
                            'mean_enrichment'][c].index.to_list()

        self.check_expectations(cfg, dataset)
        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        run_mean_enrichment_line_plot(dataset, out_plot_dir, cfg)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        try:
            if not set(cfg.analysis.metabolites.keys()).issubset(
                    dataset.metadata_df['short_comp']):
                raise ValueError(
                    "[Analysis > Metabolites > compartments] are missing "
                    "from [Metadata > Compartments]"
                )
            if cfg.analysis.method.color_lines_by not \
                    in ["metabolite", "condition"]:
                raise ValueError(
                    "[config > analysis > method > color_lines_by] "
                    "must be metabolite or condition"
                )
            if cfg.analysis.width_subplot is not None:
                try:
                    float(cfg.analysis.width_subplot)
                except TypeError:
                    logger.error("Unrecognized value for width_subplot"
                                 " in the config file")
                assert cfg.analysis.width_subplot >= 0, logger.error(
                    "width_subplot must be a positive number"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided "
                f"in the config file: {e}, aborting"
            )
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)


class PcaAnalysis(Method):
    config: PcaAnalysisConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info("Will perform PCA analysis and save tables, "
                    "with the following config: %s", self.config)
        out_table_dir = os.path.join(os.getcwd(), cfg.table_path)
        os.makedirs(out_table_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)
        available_pca_suitable_datatypes = set(
            ['abundances', 'mean_enrichment']
        ).intersection(dataset.available_datasets)
        for file_name in available_pca_suitable_datatypes:
            logger.info("Running pca analysis of %s,",
                        dataset.get_file_for_label(file_name))
            run_pca_analysis(file_name, dataset, cfg,
                             out_table_dir, mode="save_tables")

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        try:
            if (dataset.abundances_df is None) and (
                    dataset.mean_enrichment_df is None
            ):
                raise ValueError(
                    "abundances and mean_enrichment are missing in [Dataset]"
                )
            if (cfg.analysis.method.pca_split_further is not None) and not (
                    set(cfg.analysis.method.pca_split_further).issubset(
                        set(["condition", "timepoint"]))):
                raise ValueError(
                    "Unknown parameters in [config > analysis > method]"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided "
                f"in the config file: {e}, aborting"
            )
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)


class PcaPlot(Method):
    config: PcaPlotConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info("Will perform PCA plots and save figures, "
                    "with the following config: %s", self.config)
        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)
        available_pca_suitable_datatypes = set(
            ['abundances', 'mean_enrichment']
        ).intersection(dataset.available_datasets)
        for file_name in available_pca_suitable_datatypes:
            # call analysis:
            pca_results_dict = run_pca_analysis(
                file_name, dataset, cfg,
                out_table_dir="",  # no writing tables in PcaPlot Method
                mode="return_results_dict")
            # plot:
            logger.info("Running pca plot(s) of %s ",
                        dataset.get_file_for_label(file_name))
            run_pca_plot(pca_results_dict, cfg, out_plot_dir)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        try:
            if (dataset.abundances_df is None) and \
                    (dataset.mean_enrichment_df is None):
                raise ValueError(
                    "abundances and mean_enrichment are missing in [Dataset]"
                )
            if (cfg.analysis.method.pca_split_further is not None) and not (
                    set(cfg.analysis.method.pca_split_further).issubset(
                        set(["condition", "timepoint"]))):
                raise ValueError(
                    "Unknown parameters in [config > analysis > method > "
                    "pca_plot > pca_split_further]"
                )
            if (cfg.analysis.method.draw_ellipses is not None) and not (
                    cfg.analysis.method.draw_ellipses in
                    ["condition", "timepoint"]):
                raise ValueError(
                    "Unknown parameters in [config > analysis > method > "
                    "pca_plot > draw_ellipses]"
                )
        except ConfigAttributeError:
            logger.error(
                "Mandatory parameter not provided "
                "in the config file: {e}, aborting"
            )
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)


class TimeCourseAnalysis(Method):
    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(f"The current working directory is {os.getcwd()}")
        logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))
        logger.info(
            "Will perform time-course analysis, with the following config:"
            " %s", self.config)
        out_table_dir = os.path.join(os.getcwd(), cfg.table_path)
        os.makedirs(out_table_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)
        for file_name, test in cfg.analysis.statistical_test.items():
            if test is None:
                continue
            tmp = dataset.get_file_for_label(file_name)  # current file
            logger.info(
                f"Running time-course analysis of {tmp} using {test} test")
            time_course_analysis(file_name, dataset, cfg, test,
                                 out_table_dir=out_table_dir)

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        user_tests = [t[1] for t in cfg.analysis.statistical_test.items()]

        try:
            if not set(user_tests).issubset(
                    set(availtest_methods)):
                raise ValueError(
                    "Statistical tests provided in the config file not "
                    "recognized, aborting"
                )
            if not (
                    len(dataset.metadata_df['timepoint'].unique()) ==
                    len(dataset.metadata_df['timenum'].unique())
            ) and (
                    len(dataset.metadata_df['timenum'].unique()) ==
                    len(list(set(list(zip(dataset.metadata_df['timenum'],
                                          dataset.metadata_df['timenum'])))))
            ):
                raise ValueError(
                    "Inconsistent timepoint and timenum columns in metadata"
                )

        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency:{e}")
            sys.exit(1)


class DistrFitPlot(Method):
    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(f"The current working directory is {os.getcwd()}")
        logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))
        logger.info(
            "Will perform distribution fitting plots, with the following "
            "config: %s", self.config)
        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        self.check_expectations(cfg, dataset)
        for file_name, test in cfg.analysis.statistical_test.items():
            if test is None:
                continue
            elif test == "disfit":
                logger.info("Running plot of %s ",
                            dataset.get_file_for_label(file_name))
                if cfg.analysis.comparisons is not None:  # pairwise
                    run_distr_fit_plot(file_name, dataset, cfg, test,
                                       out_plot_dir=out_plot_dir,
                                       mode="pairwise")
                else:
                    continue
                    # not set for time_course, reason: too many plots

    def check_expectations(self, cfg: DictConfig, dataset: Dataset) -> None:
        user_tests = [t[1] for t in cfg.analysis.statistical_test.items()]
        try:
            if "disfit" not in user_tests:
                raise ValueError(
                    "No distribution fitting in config file, aborting"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency:{e}")
            sys.exit(1)


class MetabologramIntegration(Method):
    config: MetabologramIntegrationConfig

    def run(self, cfg: DictConfig, dataset: Dataset) -> None:
        logger.info(f"The current working directory is {os.getcwd()}")
        logger.info("Current configugration is %s", OmegaConf.to_yaml(cfg))
        logger.info("Will perform metabologram, "
                    "with the following config: %s", self.config)
        # 'data_integration' (data for integration) inherits from dataset
        data_integration: DataIntegration = DataIntegration(
            config=hydra.utils.instantiate(cfg.analysis.dataset))
        data_integration.set_dataset_integration_config()
        data_integration.load_deg_dfs()
        data_integration.load_pathways_dfs()

        self.check_expectations(cfg, data_integration)
        self.check_expectations_config_transcripto(cfg, data_integration)
        self.check_expectations_config_metabo(cfg, data_integration)

        out_plot_dir = os.path.join(os.getcwd(), cfg.figure_path)
        os.makedirs(out_plot_dir, exist_ok=True)
        file_name = list(cfg.analysis.statistical_test.keys())[0]
        test = cfg.analysis.statistical_test[file_name]

        try:
            #  'ID' is facultative in config 'columns_metabolites'
            str(cfg.analysis.columns_metabolites['ID'])
        except KeyError:
            with open_dict(cfg):  # if not set, assign "metabolite"
                cfg.analysis.columns_metabolites['ID'] = "metabolite"

        metabolome_file_name = data_integration.get_file_for_label(file_name)
        logger.info(
            f"Running metabologram of {metabolome_file_name}"
            f" (differences computed using {test} test)")
        run_metabologram(file_name, data_integration,
                         cfg, test, out_plot_dir=out_plot_dir)

    def check_expectations(self, cfg: DictConfig,
                           data_integration: DataIntegration) -> None:
        # expectations proper to the metabologram
        try:
            if 'values' not in set(cfg.analysis.columns_metabolites.keys()):
                raise ValueError("Missing 'values' in columns_metabolites")
            if cfg.analysis.columns_metabolites['values'] not in \
                    metabolites_values_for_metabologram:
                raise ValueError(
                    f"column '{cfg.analysis.columns_metabolites['values']}'"
                    f" in metabolites not recognized , currently supported:"
                    f" {metabolites_values_for_metabologram}"
                )
            if not isinstance(cfg.analysis.compartment, str):
                raise ValueError("compartment must be string in config file")
            if cfg.analysis.compartment not in \
                    set(data_integration.metadata_df['short_comp']):
                raise ValueError(
                    f"the compartment '{cfg.analysis.compartment}' "
                    f"in the config file does not exist. Must be one of: "
                    f"{set(data_integration.metadata_df['short_comp'])}"
                )
            if not len(cfg.analysis.statistical_test.keys()) == 1:
                raise ValueError(
                    f"More than one type of data in config file: "
                    f"{list(cfg.analysis.statistical_test.keys())}, "
                    f"must be only one: abundances OR mean_enrichment"
                )
            type_of_data = list(cfg.analysis.statistical_test.keys())
            if type_of_data[0] not in \
                    data_types_suitable_for_metabologram:
                raise ValueError(
                    f"Unauthorized data types '{type_of_data[0]}' "
                    f" in config file for the Metabologram. Must be one of: "
                    f"{data_types_suitable_for_metabologram}"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)

    def check_expectations_config_transcripto(
            self, cfg: DictConfig, data_integration: DataIntegration) -> None:
        try:
            if not set(
                    cfg.analysis.columns_transcripts.keys()
            ) == set(columns_transcripts_config_keys):
                raise ValueError(
                    f"Unrecognized 'columns_transcripts' parameters "
                    f"{cfg.analysis.columns_transcripts.keys()}. Must be: "
                    f"{columns_transcripts_config_keys}"
                )

            if not all(
                    set(cfg.analysis.columns_transcripts.values()).issubset(
                        set(df.columns)) for df in
                    data_integration.deg_dfs.values()):
                raise ValueError(
                    f"{set(cfg.analysis.columns_transcripts.values())} absent"
                    f" in one or more of the transcripts dataframes"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)

    def check_expectations_config_metabo(
            self, cfg: DictConfig, data_integration: DataIntegration) -> None:
        try:

            if not len(cfg.analysis.comparisons) == len(
                    data_integration.deg_dfs):
                raise ValueError(
                    "the number of comparisons, and the number of transcripts"
                    " dataframes, are unequal"
                )
            if not (all(len(c) == 2 for c in cfg.analysis.comparisons)):
                raise ValueError(
                    "Number of conditions has to be 2 "
                    "for a pairwise comparison, see config file"
                )
            if not all(any(
                    item[0] in set(data_integration.metadata_df["condition"])
                    for item
                    in sublist)
                       for sublist in cfg.analysis.comparisons):
                raise ValueError(
                    "Conditions provided for comparisons in the config file "
                    "are not present in the metadata file, aborting"
                )
            if not all(any(
                    item[1] in set(data_integration.metadata_df["timepoint"])
                    for item
                    in sublist)
                       for sublist in cfg.analysis.comparisons):
                raise ValueError(
                    "Timepoints provided for comparisons in the config file "
                    "are not present in the metadata file, aborting"
                )
            # all values in the comparisons arrays of arrays are in the
            # metadata, either as conditions or timepoints
            conditions_and_tp = set(
                data_integration.metadata_df["condition"]).union(
                set(data_integration.metadata_df["timepoint"]))

            comparisons = set(flatten(cfg.analysis.comparisons))
            diff = comparisons.difference(conditions_and_tp)
            if not comparisons.issubset(conditions_and_tp):
                raise ValueError(
                    f"Comparisons > Conditions or timepoints provided in the "
                    f"config file {diff} are not present in the metadata "
                    f"file, aborting"
                )
        except ConfigAttributeError as e:
            logger.error(
                f"Mandatory parameter not provided in the config file:{e}, "
                f"aborting")
            sys.exit(1)
        except ValueError as e:
            logger.error(f"Data inconsistency: {e}")
            sys.exit(1)
