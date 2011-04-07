#!/usr/bin/env python
#
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009-2011  EDF
#
#   Code_Saturne is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License,
#   or (at your option) any later version.
#
#   Code_Saturne is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import os
import os.path
import sys

import cs_config

from cs_exec_environment import *
from cs_case_domain import *
from cs_case import *

#-------------------------------------------------------------------------------
# Extract a parameters file name from a shell run sript
#-------------------------------------------------------------------------------

def get_param(path):
    """
    Extract a parameters file name from a shell run sript
    """
    f = open(path, 'r')
    lines = f.readlines()
    f.close()

    param = None

    for line in lines:
        if line[0] == '#':
            continue
        i = line.find('#')
        if i > -1:
            line = line[0:i]
        line = line.strip()
        i = string.find(line, '--param')
        if i >= 0:
            i += len('--param')
            line = line[i:].strip()
            if not line:
                continue
            # Find file name, possibly protected by quotes
            # (protection by escape character not handled)
            sep = line[0]
            if sep == '"' or sep == "'":
                param = line.split(sep)[1]
            else:
                param = line.split()[0]

    return param

#===============================================================================
# Main function for code coupling execution
#===============================================================================

def coupling(package,
             domains,
             casedir,
             exec_preprocess = True,
             exec_partition = True,
             exec_solver = True):

    use_saturne = False
    use_syrthes = False
    use_neptune = False

    # Initialize code domains
    sat_domains = []
    syr_domains = []
    nep_domains = []

    if domains == None:
        raise RunCaseError('No domains defined.')

    for d in domains:

        if (d.get('script') == None or d.get('domain') == None):
            msg = 'Check your coupling definition.\n'
            msg += 'script or domain key is missing.'
            raise RunCaseError(msg)

        if (d.get('solver') == 'Code_Saturne' or d.get('solver') == 'Saturne'):

            try:
                runcase = os.path.join(os.getcwd(),
                                       d.get('domain'),
                                       'SCRIPTS',
                                       d.get('script'))
                param = get_param(runcase)

            except Exception:
                err_str = 'Cannot read Code_Saturne script: ' + runcase
                raise RunCaseError(err_str)

            dom = domain(package,
                         name = d.get('domain'),
                         param = param,
                         n_procs_weight = d.get('n_procs_weight'),
                         n_procs_min = d.get('n_procs_min'),
                         n_procs_max = d.get('n_procs_max'),
                         n_procs_partition = None)

            use_saturne = True
            sat_domains.append(dom)

        elif (d.get('solver') == 'SYRTHES 3' or d.get('solver') == 'SYRTHES3'):

            try:
                dom = syrthes3_domain(package,
                                      name = d.get('domain'),
                                      echo_comm = d.get('echo_comm'),
                                      coupling_mode = 'MPI',  # 'MPI' or 'Sockets'
                                      coupled_apps = d.get('coupled_apps'))

            except Exception:
                err_str = 'Cannot create SYRTHES 3 domain.\n'
                err_str += ' case = ' + d.get('domain')
                raise RunCaseError(err_str)

            use_syrthes = True
            syr_domains.append(dom)

        elif (d.get('solver') == 'SYRTHES'):

            try:
                dom = syrthes_domain(package,
                                     cmd_line = d.get('opt'),
                                     name = d.get('domain'),
                                     param = d.get('script'),
                                     n_procs_weight = d.get('n_procs_weight'),
                                     n_procs_min = d.get('n_procs_min'),
                                     n_procs_max = d.get('n_procs_max'))

            except Exception:
                err_str = 'Cannot create SYRTHES domain. Opt = ' + d.get('opt') + '\n'
                err_str += ' domain = ' + d.get('domain')
                err_str += ' script = ' + d.get('script') + '\n'
                err_str += ' n_procs_weight = ' + str(d.get('n_procs_weight')) + '\n'
                raise RunCaseError(err_str)

            use_syrthes = True
            syr_domains.append(dom)

        elif (d.get('solver') == 'NEPTUNE_CFD'):

            try:
                runcase = os.path.join(os.getcwd(),
                                       d.get('domain'),
                                       'SCRIPTS',
                                       d.get('script'))
                param = get_param(runcase)

            except Exception:
                err_str = 'Cannot read NEPTUNE_CFD script: ' + runcase
                raise RunCaseError(err_str)

            dom = domain(package,
                         name = d.get('domain'),
                         param = param,
                         n_procs_weight = d.get('n_procs_weight'),
                         n_procs_min = d.get('n_procs_min'),
                         n_procs_max = d.get('n_procs_max'),
                         n_procs_partition = None)

            use_neptune = True
            nep_domains.append(dom)

        elif (d.get('solver') == 'Code_Aster' or d.get('solver') == 'Aster'):
            err_str = 'Code_Aster code coupling not handled yet.\n'
            raise RunCaseError(err_str)

        else:
            err_str = 'Unknown code type : ' + d.get('solver') + '.\n'
            raise RunCaseError(err_str)

    # Now handle case for the corresponding calculation domain(s).

    c = case(package,
             casedir,
             sat_domains + nep_domains,
             syr_domains,
             exec_preprocess = exec_preprocess,
             exec_partition = exec_partition,
             exec_solver = exec_solver)

    msg = ' Coupling execution between: \n'
    if use_saturne == True:
        msg += '   o Code_Saturne [' + str(len(sat_domains)) + ' domain(s)];\n'
    if use_syrthes == True:
        msg += '   o SYRTHES      [' + str(len(syr_domains)) + ' domain(s)];\n'
    if use_neptune == True:
        msg += '   o NEPTUNE_CFD  [' + str(len(nep_domains)) + ' domain(s)];\n'
    sys.stdout.write(msg+'\n')

    return c

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
