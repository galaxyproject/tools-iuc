#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
from typing import Literal, Optional, get_args


def assert_literal(value: str, lit_type, check: Optional[str] = None):
    assert value in get_args(lit_type), \
        (check if check else "") + f"Value {value} is not in {lit_type}"


data_files_keys = [
    "abundances",
    "mean_enrichment",
    "isotopologue_proportions",
    "isotopologues",
]

data_files_keys_type = Literal[
    "abundances",
    "mean_enrichment",
    "isotopologue_proportions",
    "isotopologues",
]

availtest_methods = ["MW", "KW", "ranksum", "Wcox", "Tt", "BrMu", "prm-scipy",
                     "disfit", "none"]

availtest_methods_type = Literal[
    "MW", "KW", "ranksum", "Wcox", "Tt", "BrMu", "prm-scipy", "disfit", "none"]

correction_methods = [
    "bonferroni",
    "sidak",
    "holm-sidak",
    "holm",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
]

correction_methods_type = Literal[
    "bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg", "hommel",
    "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"
]

comparison_modes = ["pairwise", "multigroup"]  # TODO: verify if ever used

comparison_modes_types = Literal[
    "pairwise", "multigroup"]  # TODO: verify if ever used

overlap_methods = ["symmetric", "asymmetric"]

overlap_methods_types = Literal["symmetric", "asymmetric"]

data_types_suitable_for_metabologram = ['abundances', 'mean_enrichment']

molecular_types_for_metabologram = ["transcripts", "metabolites"]

columns_transcripts_config_keys = ['ID', 'values']

metabolites_values_for_metabologram = ['log2FC', 'FC']
