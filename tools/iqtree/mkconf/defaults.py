#!/usr/bin/env python3

# Args to exclude from user customatization
# Give a value to override default
exclude_map = {
    # General
    '-pre' : 'PRE',   # prefix
    '-o'   : 'OUT',   # outfile
    '-h'   : False,   # help
    '-v'   : False,   # verbose
    '-mem' : False,   # max ram spec
    '-quiet' : False, # suppress printing to screen
    '-nt' : "AUTO",    # number of cores

    # checkpointing
    '-redo' : False,
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
    '-wt' : "true",
    '-fixbr' : "true",
    '-wsl' : "true",
    '-wslg' : "true"   
}

# Args to override, and show to user.
override_defaults = {
    '-seed': '1547',
}
