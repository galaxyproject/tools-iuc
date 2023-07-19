#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

import logging
import os

from data import Dataset

import hydra

from method import Method


from omegaconf import DictConfig, OmegaConf


logger = logging.getLogger(__name__)


@hydra.main(config_path="../config", config_name="config")
def main_run_analysis(cfg: DictConfig) -> None:
    logger.info(f"The current working directory is {os.getcwd()}")
    logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))

    dataset: Dataset = Dataset(
        config=hydra.utils.instantiate(cfg.analysis.dataset))
    dataset.preload()
    dataset.split_datafiles_by_compartment()
    dataset.save_datafiles_split_by_compartment()
    method: Method = hydra.utils.instantiate(
        cfg.analysis.method).build()  # method factory

    method.run(cfg, dataset)


if __name__ == "__main__":
    main_run_analysis()
