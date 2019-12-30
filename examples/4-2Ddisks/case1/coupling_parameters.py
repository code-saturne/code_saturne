# -*- coding: utf-8 -*-

#===============================================================================
# User variable settings to specify a coupling computation environnement.

# A coupling case is defined by a dictionnary, containing the following:

# Solver type ('Code_Saturne', 'SYRTHES', 'NEPTUNE_CFD' or 'code_aster')
# Domain directory name
# Run parameter setting file
# Number of processors (or None for automatic setting)
# Optional command line parameters. If not useful = None
#===============================================================================

# Ensure the correct SYRTHES install is used.
sys.path.insert(1, '<code_saturne.cfg/syrthes>')

# Define coupled domains

domains = [

    {'solver': 'Code_Saturne',
     'domain': 'fluid',
     'script': 'runcase',
     'n_procs_weight': None,
     'n_procs_min': 1,
     'n_procs_max': None}

    ,
    {'solver': 'SYRTHES',
     'domain': 'solid',
     'script': 'solid-coupling.syd',
     'n_procs_weight': None,
     'n_procs_min': 1,
     'n_procs_max': None,
     'opt' : '-v ens'}         # Additional SYRTHES options
                               # (ex.: postprocessing with '-v ens' or '-v med')

    ]

