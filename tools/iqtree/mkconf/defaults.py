#!/usr/bin/env python3

# Section nesting

section_nestmap = {
    # section : (expanded, parent)
    'General options' :                       (True, ""),
    'Likelihood mapping analysis' :           (False,""),
    'Checkpointing to resume stopped run' :   (False,""),
    'Automatic model selection' :             (False, "Modelling Parameters"),
    'Specifying substitution models' :        (False, "Modelling Parameters"),
    'Rate heterogeneity' :                    (False, "Modelling Parameters"),
    'Partition model options' :               (False, "Modelling Parameters"),
    'Site-specific frequency model options' : (False, "Modelling Parameters"),
    'Tree search parameters' :                (False, "Tree Parameters"),
    'Ultrafast bootstrap parameters' :        (False, "Bootstrap Parameters"),
    'Nonparametric bootstrap' :               (False, "Bootstrap Parameters",),
    'Single branch tests' :                   (False, "Tree Parameters"),
    'Tree topology tests' :                   (False, "Tree Parameters"),
    'Constructing consensus tree' :           (False, "Tree Parameters"),
    'Computing Robinson-Foulds distance' :    (False, "Tree Parameters"),
    'Generating random trees' :               (False, "Tree Parameters"),
    'Miscellaneous options' :                 (False, "")
}

# Args to exclude from user customatization
# Give a value to override default
exclude_map = {
    # General
    '-pre' : 'PREF',   # prefix
    '-o'   : False,   # outfile
    '-h'   : False,   # help
    '-v'   : False,   # verbose
    '-mem' : False,   # max ram spec
    '-quiet' : False, # suppress printing to screen
    '-nt' : "AUTO",    # number of cores

    # checkpointing
    '-redo' : True,
    '-cptime' : False, # default 20s   

    # Likelihood Mapping Analysis
    # Automatic Model Selection
    '-mredo' : False,

    # Rate heterogeneity
    # Partition Models
    # Site specific
    # Tree search params
    # Ultrafast bootstrap
    # nonparametric bootstrap
    # single branch tests
    # tree topology tests
    # constructing consensus tree
    # computing robinson-foulds distance
    # Generating random trees
    # miscellaneous options
    '-wt' : False,
    '-fixbr' : False,
    '-wsl' : False,    # throws errors if "true"
    '-wslg' : False
}

# Args to override, and show to user.
override_defaults = {
#    '-seed': '1547',
}
