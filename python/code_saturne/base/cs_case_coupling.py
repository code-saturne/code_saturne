#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

import configparser
import os
import os.path
import sys

from code_saturne.base.cs_case_domain import *
from code_saturne.base.cs_case import *
from code_saturne.base import cs_exec_environment
from code_saturne.base import cs_run_conf
from code_saturne.base import cs_runcase

#===============================================================================
# Main function for code coupling execution
#===============================================================================

def coupling(package,
             domains,
             casedir,
             dest_dir = None,
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

    if package_compute is None:
        package_compute = package

    # Initialize code domains
    sat_domains = []
    syr_domains = []
    nep_domains = []
    cat_domains = []
    py_domains = []

    if domains is None:
        raise RunCaseError('No domains defined.')

    for d in domains:

        domain_s = d.get('domain')
        solver_s = d.get('solver').lower()
        script_s = None
        param_s = None

        if (domain_s is None):
            msg = 'Check your coupling definition.\n'
            msg += 'domain key is missing.'
            raise RunCaseError(msg)

        # First, determine parameter file to use for code_saturne
        # or associated modules (ensuring backwards compatibiliy)

        if solver_s in package.config.solver_modules.keys() \
           or solver_s == 'cathare':

            param = None

            script = d.get('script')          # v6.1 and older structure
            if script != None:
                if script[-4:] == '.xml':
                    param = script
            else:
                param = d.get('param')
                if not param:
                    param = d.get('paramfile') # older_name

            s_dir = staging_dir
            if not s_dir:
                s_dir = casedir
            if not s_dir:
                s_dir = os.getcwd()

            if param is None:
                run_conf_path = os.path.join(s_dir,
                                             domain_s,
                                             'DATA',
                                             'run.cfg')
                if os.path.isfile(run_conf_path):
                    run_conf = cs_run_conf.run_conf(run_conf_path)
                    if 'setup' in run_conf.sections:
                        if 'param' in run_conf.sections['setup']:
                            param = run_conf.sections['setup']['param']

            if script and not param: #  v6.1 and older case structure
                runcase_path = os.path.join(s_dir,
                                            domain_s,
                                            'SCRIPTS',
                                            script)
                if os.path.isfile(runcase_path):
                    try:
                        runcase = cs_runcase.runcase(runcase_path)
                        param = runcase.get_parameters()
                    except Exception:
                        err_str = 'Cannot read ' + d.get('solver') \
                                  + ' script: ' + runcase_path
                        raise RunCaseError(err_str)

            # Remark: if param is undefined, the code_saturne domain will
            # default to 'setup.xml' if present.

            d['param'] = param

        # Now build case domain for the different solvers:

        if solver_s in package.config.solver_modules.keys():

            dom = domain(package,
                         package_compute = package_compute,
                         name = domain_s,
                         param = d.get('param'),
                         n_procs_weight = d.get('n_procs_weight'),
                         n_procs_min = d.get('n_procs_min'),
                         n_procs_max = d.get('n_procs_max'))

            if solver_s == 'code_saturne':
                use_saturne = True
                sat_domains.append(dom)
            elif solver_s == 'neptune_cfd':
                use_neptune = True
                nep_domains.append(dom)

        elif solver_s == 'syrthes':

            param_s = d.get('param')
            if param_s is None:
                param_s = d.get('script') # older name

            if (param_s is None):
                msg = 'Check your coupling definition.\n'
                msg += 'parameters file selection is missing for domain: '
                msg += domain_s + '.\n'
                raise RunCaseError(msg)

            try:
                dom = syrthes_domain(package,
                                     cmd_line = d.get('opt'),
                                     name = domain_s,
                                     param = d.get('param'),
                                     n_procs_weight = d.get('n_procs_weight'),
                                     n_procs_min = d.get('n_procs_min'),
                                     n_procs_max = d.get('n_procs_max'),
                                     verbose = verbose)

            except Exception:
                err_str = 'Cannot create SYRTHES domain. Opt = ' + d.get('opt') + '\n'
                err_str += ' domain = ' + domain_s + '\n'
                err_str += ' n_procs_weight = ' + str(d.get('n_procs_weight')) + '\n'
                raise RunCaseError(err_str)

            use_syrthes = True
            syr_domains.append(dom)

        elif solver_s == 'cathare':

            # Current version using Cathare2: the cathare case is converted to a
            # .so library which is opened and launched by a neptune_cfd executable

            dom = cathare_domain(package,
                                 package_compute = package_compute,
                                 name = domain_s,
                                 param = d.get('param'),
                                 n_procs_weight = None,
                                 n_procs_min = 1,
                                 n_procs_max = 1,
                                 cathare_case_file = d.get('cathare_case_file'),
                                 neptune_cfd_dom = d.get('neptune_cfd_domain'))

            use_cathare = True
            cat_domains.append(dom)

        elif solver_s == 'python_code':

            script_s = d.get('script')
            if (script_s is None):
                msg = 'Check your coupling definition.\n'
                msg += 'Python script file selection is missing for domain: '
                msg += domain_s + '.\n'
                raise RunCaseError(msg)

            # Generic code_saturne/Python Script coupling
            # The python script can contain any MPI compatible code or supervisor

            try:
                dom = python_domain(package,
                                    name = domain_s,
                                    cmd_line = d.get('command_line'),
                                    script_name = script_s)

            except Exception:
                err_str = 'Cannot create Python code domain.\n'
                err_str += ' domain = ' + domain_s + '\n'
                err_str += ' script = ' + str(d.get('script')) + '\n'
                raise RunCaseError(err_str)

            use_py_code = True
            py_domains.append(dom)

        else:
            err_str = 'Unknown code type : ' + d.get('solver') + '.\n'
            raise RunCaseError(err_str)

    # Now handle case for the corresponding calculation domain(s).

    c = case(package,
             package_compute = package_compute,
             case_dir = casedir,
             dest_dir = dest_dir,
             staging_dir = staging_dir,
             domains = sat_domains + nep_domains + cat_domains,
             syr_domains = syr_domains,
             py_domains = py_domains)

    if verbose:
        msg = ' Coupling execution between: \n'
        if use_saturne == True:
            msg += '   o code_saturne [' + str(len(sat_domains)) + ' domain(s)];\n'
        if use_syrthes == True:
            msg += '   o SYRTHES      [' + str(len(syr_domains)) + ' domain(s)];\n'
        if use_neptune == True:
            msg += '   o neptune_cfd  [' + str(len(nep_domains)) + ' domain(s)];\n'
        if use_cathare == True:
            msg += '   o CATHARE2     [' + str(len(cat_domains)) + ' domain(s)];\n'
        if use_py_code == True:
            msg += '   o Python Script  [' + str(len(py_domains)) + ' domain(s)];\n'
        sys.stdout.write(msg+'\n')

    return c

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
