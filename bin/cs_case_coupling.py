#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------
try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3
import os
import os.path
import sys

from cs_case_domain import *
from cs_case import *
import cs_exec_environment
import cs_runcase

#===============================================================================
# Main function for code coupling execution
#===============================================================================

def coupling(package,
             domains,
             casedir,
             staging_dir = None,
             verbose = True,
             package_compute = None):

    use_saturne = False
    use_syrthes = False
    use_neptune = False
    use_cathare = False
    use_py_code = False

    # Use alternate compute (back-end) package if defined

    config = configparser.ConfigParser()
    config.read(package.get_global_configfile())

    if package_compute == None:
        package_compute = package

    # Initialize code domains
    sat_domains = []
    syr_domains = []
    nep_domains = []
    cat_domains = []
    ast_domain = []
    fsi_coupler = None
    py_domains = []

    if domains == None:
        raise RunCaseError('No domains defined.')

    for d in domains:

        if ((d.get('script') == None or d.get('domain') == None) \
            and d.get('coupler') == None):
            msg = 'Check your coupling definition.\n'
            msg += 'script or domain key is missing.'
            raise RunCaseError(msg)

        if (d.get('solver') == 'Code_Saturne' or d.get('solver') == 'NEPTUNE_CFD'):

            script = d.get('script')

            if script[-4:] == '.xml':
                param = script

            else:
                runcase_path = os.path.join(os.getcwd(),
                                            d.get('domain'),
                                            'SCRIPTS',
                                            script)
                try:
                    runcase = cs_runcase.runcase(runcase_path)
                    param = runcase.get_parameters()

                except Exception:
                    err_str = 'Cannot read ' + d.get('solver') \
                              + ' script: ' + runcase_path
                    raise RunCaseError(err_str)
                d['script'] = param

            dom = domain(package,
                         package_compute = package_compute,
                         name = d.get('domain'),
                         param = param,
                         n_procs_weight = d.get('n_procs_weight'),
                         n_procs_min = d.get('n_procs_min'),
                         n_procs_max = d.get('n_procs_max'))

            if (d.get('solver')) == 'Code_Saturne':
                use_saturne = True
                sat_domains.append(dom)
            elif (d.get('solver')) == 'NEPTUNE_CFD':
                use_neptune = True
                nep_domains.append(dom)

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

        elif (d.get('solver') == 'Code_Aster' or d.get('solver') == 'Aster'):

            if len(ast_domain) > 0:
                err_str = 'Only 1 Code_Aster domain is currently handled\n'
                raise RunCaseError(err_str)

            try:
                dom = aster_domain(package,
                                   name = d.get('domain'),
                                   param = d.get('script'))

            except Exception:
                err_str = 'Cannot create Code_Aster domain.\n'
                err_str += ' domain = ' + d.get('domain') + '\n'
                err_str += ' script = ' + d.get('script') + '\n'
                raise RunCaseError(err_str)

            ast_domain.append(dom)

        elif (d.get('coupler') == 'FSI_coupler'):

            if fsi_coupler:
                err_str = 'Only 1 FSI coupler is currently handled\n'
                raise RunCaseError(err_str)

            try:
                fsi_coupler = {'max_time_steps' : d.get('max_time_steps'),
                               'n_sub_iterations' : d.get('n_sub_iterations'),
                               'time_step' : d.get('time_step'),
                               'start_time' : d.get('start_time'),
                               'epsilon' : d.get('epsilon')}

            except Exception:
                err_str = 'Cannot create FSI coupler\n'
                err_str += '  max_time_steps = ' + d.get('max_time_steps') + '\n'
                err_str += '  n_sub_iterations = ' + d.get('n_sub_iterations' + '\n')
                err_str += '  time_step = ' + d.get('time_step') + '\n'
                err_str += '  start_time = ' + d.get('start_time') + '\n'
                err_str += '  epsilon = ' + d.get('epsilon') + '\n'
                raise RunCaseError(err_str)

        elif (d.get('solver') == 'CATHARE'):
            # Current version using Cathare2: the cathare case is converted to a
            # .so library which is opened and launched by a NEPTUNE_CFD executable

            script = d.get('script')

            if script[-4:] == '.xml':
                param = script

            else:
                runcase_path = os.path.join(os.getcwd(),
                                            d.get('domain'),
                                            'SCRIPTS',
                                            script)
                try:
                    runcase = cs_runcase.runcase(runcase_path)
                    param = runcase.get_parameters()

                except Exception:
                    err_str = 'Cannot read ' + d.get('solver') \
                              + ' script: ' + runcase_path
                    raise RunCaseError(err_str)
                d['script'] = param

            dom = cathare_domain(package,
                                 package_compute = package_compute,
                                 name = d.get('domain'),
                                 param = d.get('paramfile'),
                                 n_procs_weight = d.get('n_procs_weight'),
                                 n_procs_min = d.get('n_procs_min'),
                                 n_procs_max = d.get('n_procs_max'),
                                 cathare_case_file=d.get('cathare_case_file'),
                                 neptune_cfd_dom=d.get('neptune_cfd_domain'))

            use_cathare = True
            cat_domains.append(dom)

        elif (d.get('solver') == 'PYTHON_CODE'):
            # Generic Code_Saturne/Python Script coupling
            # The python script can contain any MPI compatible code or supevisor

            try:
                dom = python_domain(package,
                                    name = d.get('domain'),
                                    cmd_line = d.get('command_line'),
                                    script_name = d.get('py_code'))

            except Exception:
                err_str = 'Cannot create Python code domain.\n'
                err_str += ' domain = ' + d.get('domain') + '\n'
                err_str += ' script = ' + d.get('script') + '\n'
                raise RunCaseError(err_str)

            use_py_code = True
            py_domains.append(dom)

        else:
            err_str = 'Unknown code type : ' + d.get('solver') + '.\n'
            raise RunCaseError(err_str)

    # Now handle case for the corresponding calculation domain(s).

    coupling_parameters = \
"""
# -*- coding: utf-8 -*-"

domains = [

"""

    d_first = True
    for d in domains:
        if not d_first:
            coupling_parameters += ",\n\n"
        d_first = False
        k_first = True
        for k in ('solver', 'domain', 'script', 'n_procs_weight', 'n_procs_min',
                  'coupler', 'max_time_steps', 'n_sub_iterations', 'time_step',
                  'start_time', 'epsilon'):
            v = d.get(k)
            if v:
                if type(v) == str:
                    s_v = "'" + v + "'"
                else:
                    s_v = str(v)
                if k_first:
                    coupling_parameters += "    {"
                    k_first = False
                else:
                    coupling_parameters += ",\n     "
                coupling_parameters += "'" + k + "': " + s_v
        coupling_parameters += "}"
    coupling_parameters += "\n\n    ]\n\n"

    c = case(package,
             package_compute = package_compute,
             case_dir = casedir,
             staging_dir = staging_dir,
             coupling_parameters = coupling_parameters,
             domains = sat_domains + nep_domains + cat_domains,
             syr_domains = syr_domains,
             ast_domain = ast_domain,
             fsi_coupler = fsi_coupler,
             py_domains = py_domains)

    if verbose:
        msg = ' Coupling execution between: \n'
        if use_saturne == True:
            msg += '   o Code_Saturne [' + str(len(sat_domains)) + ' domain(s)];\n'
        if use_syrthes == True:
            msg += '   o SYRTHES      [' + str(len(syr_domains)) + ' domain(s)];\n'
        if use_neptune == True:
            msg += '   o NEPTUNE_CFD  [' + str(len(nep_domains)) + ' domain(s)];\n'
        if ast_domain or fsi_coupler:
            msg += '   o Code_Aster   [1 domain(s)];\n'
            msg += '                  [1 coupler(s)];\n'
        if use_cathare == True:
            msg += '   o CATHARE2     [' + str(len(cat_domains)) + ' domain(s)];\n'
        if use_py_code == True:
            msg += '   o Python Script  [' + str(len(py_domains)) + ' domain(s)];\n'
        sys.stdout.write(msg+'\n')

    return c

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
