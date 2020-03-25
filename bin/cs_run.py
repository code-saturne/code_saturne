#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

"""
This module describes the script used to run a study/case for Code_Saturne.

This module defines the following functions:
- process_cmd_line
- main
"""

#===============================================================================
# Import required Python modules
#===============================================================================

from __future__ import print_function

import os, sys
import types, string, re, fnmatch
from optparse import OptionParser
try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

from code_saturne import cs_exec_environment
from code_saturne import cs_case_domain
from code_saturne import cs_case
from code_saturne import cs_run_conf

#-------------------------------------------------------------------------------
# Update run steps based on run_conf object
#-------------------------------------------------------------------------------

def update_run_steps(s_c, run_conf, final=False):
    """
    Process the passed command line arguments.
    """

    filter_stages = False
    for k in s_c:
        if s_c[k] != None:
            filter_stages = True

    if run_conf and not filter_stages:
        if 'run' in run_conf.sections:
            for kw in s_c:
                s_c[kw] = run_conf.get_bool('run', kw)

    filter_stages = False
    for k in s_c:
        if s_c[k] != None:
            filter_stages = True

    # Default if nothing provided, ensure range is filled otherwise

    if filter_stages:
        i_s = -1
        i_f = -1
        for i, k in enumerate(s_c):
            if s_c[k] == True:
                if i_s < 0:
                    i_s = i
                i_f = i + 1
        for i, k in enumerate(s_c):
            if i < i_s:
                s_c[k] = False
            elif i < i_f:
                s_c[k] = True
            else:
                s_c[k] = False

    elif final:
        for i, k in enumerate(s_c):
            s_c[k] = True

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog run [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("--compute-build", dest="compute_build", type="string",
                      metavar="<build>",
                      help="base name or full path to the compute build")

    parser.add_option("-n", "--nprocs", dest="nprocs", type="int",
                      metavar="<nprocs>",
                      help="number of MPI processes for the computation")

    parser.add_option("--nt", "--threads-per-task", dest="nthreads", type="int",
                      help="number of OpenMP threads per task")

    parser.add_option("-p", "--param", dest="param", type="string",
                      metavar="<param>",
                      help="path or name of the parameters file")

    parser.add_option("--case", dest="case", type="string",
                      metavar="<case>",
                      help="path to the case's directory")

    parser.add_option("--coupling", dest="coupling", type="string",
                      metavar="<coupling>",
                      help="path or name of the coupling descriptor file")

    parser.add_option("--id", dest="id", type="string",
                      metavar="<id>",
                      help="use the given run id")

    parser.add_option("--id-prefix", dest="id_prefix", type="string",
                      metavar="<prefix>",
                      help="prefix the run id with the given string")

    parser.add_option("--id-suffix", dest="id_suffix", type="string",
                      metavar="<suffix>",
                      help="suffix the run id with the given string")

    parser.add_option("--suggest-id", dest="suggest_id",
                      action="store_true",
                      help="suggest a run id for the next run")

    parser.add_option("--force", dest="force",
                      action="store_true",
                      help="run the data preparation stage even if " \
                           + "the matching execution directory exists")

    parser.add_option("--stage", dest="stage",
                      action="store_true",
                      help="stage data prior to preparation and execution")

    parser.add_option("--initialize", "--preprocess", dest="initialize",
                      action="store_true",
                      help="run the data preparation stage")

    parser.add_option("--compute", "--execute", dest="compute",
                      action="store_true",
                      help="run the compute stage")

    parser.add_option("--finalize", dest="finalize",
                      action="store_true",
                      help="run the results copy/cleanup stage")

    parser.set_defaults(compute_build=False)
    parser.set_defaults(suggest_id=False)
    parser.set_defaults(stage=None)
    parser.set_defaults(initialize=None)
    parser.set_defaults(compute=None)
    parser.set_defaults(finalize=None)
    parser.set_defaults(param=None)
    parser.set_defaults(coupling=None)
    parser.set_defaults(domain=None)
    parser.set_defaults(id=None)
    parser.set_defaults(nprocs=None)
    parser.set_defaults(nthreads=None)

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    (options, args) = parser.parse_args(argv)

    # Stages to run (if no filter given, all are done).

    s_c = {'stage': options.stage,
           'initialize': options.initialize,
           'compute': options.compute,
           'finalize': options.finalize}

    filter_stages = False
    for k in s_c:
        if s_c[k]:
            filter_stages = True

    # Try to determine case directory

    casedir = None
    staging_dir = None
    param = None
    compute_build = None

    if options.coupling and options.param:
        # Multiple domain case
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        err_str = 'Error:' + os.linesep
        err_str += cmd_line + os.linesep
        err_str += '--coupling and -p/--param options are incompatible.'
        raise RunCaseError(err_str)

    # Also check for possible settings file

    coupling = options.coupling
    run_id = options.id
    run_conf = None

    run_config_path = os.path.join(os.getcwd(), 'run.cfg')
    if os.path.isfile(run_config_path):
        run_conf = cs_run_conf.run_conf(run_config_path, package=pkg)
        if not coupling:
            if 'setup' in run_conf.sections:
                if 'coupling' in run_conf.sections['setup']:
                    coupling = run_conf.sections['setup']['coupling']
        if not run_id and not filter_stages:
            if 'run' in run_conf.sections and not filter_stages:
                update_run_steps(s_c, run_conf)
                if s_c['stage'] == False:
                    if 'id' in run_conf.sections['run']:
                        run_id = run_conf.sections['run']['id']

    casedir, staging_dir = cs_case.get_case_dir(case=options.case,
                                                param=options.param,
                                                coupling=coupling,
                                                id=run_id)

    if casedir == None:
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        print('Error:', file = sys.stderr)
        print(cmd_line, file = sys.stderr)
        print('run from directory \"' + str(os.getcwd()) + '\",', file = sys.stderr)
        print('which does not seem to be inside a case directory.', file = sys.stderr)

    param = options.param

    # If no parameter file passed, and a setup.xml is present in DATA, run it
    if param is None:
        has_setup = os.path.isfile(os.path.join(casedir, 'DATA', 'setup.xml'))
        if has_setup:
            param = "setup.xml"

    compute_build = options.compute_build

    if not options.force:
        force_id = False
    else:
        force_id = True

    # Return associated dictionary (also force number of ranks and threads)

    r_c = {'casedir': casedir,
           'staging_dir': staging_dir,
           'run_id': run_id,
           'param': param,
           'coupling': coupling,
           'id_prefix': options.id_prefix,
           'id_suffix': options.id_suffix,
           'suggest_id': options.suggest_id,
           'force_id': force_id,
           'n_procs': options.nprocs,
           'n_threads': options.nthreads,
           'compute_build': compute_build}

    return r_c, s_c, run_conf

#-------------------------------------------------------------------------------
# Read the run configuration file
#-------------------------------------------------------------------------------

def read_run_config_file(i_c, r_c, s_c, pkg, run_conf=None):
    """
    Process the passed command line arguments.
    """

    casedir = r_c['casedir']
    run_config_path = ""
    setup_default_path = ""

    if r_c['coupling']:
        run_config_path = os.path.join(casedir, 'run.cfg')
    else:
        run_config_path = os.path.join(casedir, 'DATA', 'run.cfg')
        setup_default_path = os.path.join(casedir, 'DATA', 'setup.xml')

    if run_conf == None:
        if not os.path.isfile(run_config_path):
            print('Warning:', file = sys.stderr)
            print('  \'run.cfg\' not found in case directory; case update recommended.',
                  file = sys.stderr)
            print('', file = sys.stderr)
            return

    # Only load run.cfg if not already done

    if run_conf and s_c['stage'] != False:
        if not run_conf.path == run_config_path:
            run_conf = None

    if not run_conf:
        run_conf = cs_run_conf.run_conf(run_config_path, package=pkg)

    # Parameters file

    for kw in ('param',):
        if not r_c[kw]:
            r_c[kw] = run_conf.get('setup', kw)

    if not r_c['param'] and setup_default_path:
        if os.path.isfile(setup_default_path):
            r_c['param'] = setup_default_path

    # Run id

    if not r_c['run_id']:
        r_c['run_id'] = run_conf.get('run', 'id')

    if not r_c['force_id']:
        r_c['force_id'] = run_conf.get_bool('run', 'force_id')

    # Compute stages

    update_run_steps(s_c, run_conf)

    # Resources: try to find a matching section, using
    # resource_name, batch, and job_defaults in decreasing priority.

    resource_name = i_c['resource_name']
    if not resource_name or not resource_name in run_conf.sections:
        resource_name = i_c['batch']
        if resource_name:
            resource_name = os.path.basename(resource_name).lower()
    if not resource_name or not resource_name in run_conf.sections:
        resource_name = 'job_defaults'

    run_conf_r = None
    if resource_name in run_conf.sections:
        run_conf_r = run_conf.sections[resource_name]

    if run_conf_r:
        for kw in ('n_procs', 'n_threads', 'time_limit'):
            if kw in r_c:
                if r_c[kw] != None:
                    continue
            r_c[kw] = None
            v = run_conf.get_int(resource_name, kw)
            if v:
                r_c[kw] = v

    run_conf_kw= ('job_parameters', 'job_header',
                  'run_prologue', 'run_epilogue',
                  'compute_prologue', 'compute_epilogue')

    for kw in run_conf_kw:
        if not kw in r_c:
            r_c[kw] = None

    if run_conf_r:
        for kw in run_conf_kw:
            if r_c[kw] != None:
                continue
            if kw in run_conf_r:
                r_c[kw] = run_conf_r[kw]

    # Handle case where files are used

    if not (r_c['job_parameters'] or r_c['job_header']):
        kw = 'job_header_file'
        f_path = None
        if kw in run_conf_kw:
            f_path = run_conf_kw[kw]
            if f_path:
                if not os.path.isabs(f_path):
                    f_prefix = os.path.basename(run_config_path)
                    f_path= os.path.join(f_prefix, f_path)
                if os.path.isfile(f_path):
                    f = file.open(f_path)
                    r_c['job_header'] = f.read()
                    f.close
                else:
                    err_str = """warning in run.cfg: [{0}] {1} = {2}
  "{3}" not present (use defaults)"""
                    print(err_str.format(resource_name, kw, r_c[kw], f_path),
                          file = sys.stderr)
                    r_c['job_header'] = None
        elif 'jobmanager' in r_c:
            err_str = 'warning in run.cfg: [{0}] {1} = {2}; not currently handled (ignored)'
            print(err_str.format(resource_name, kw, r_c[kw]),
                  file = sys.stderr)
            r_c[kw] = None

#-------------------------------------------------------------------------------
# Generate the run configuration file for further steps
#-------------------------------------------------------------------------------

def generate_run_config_file(path, resource_name, r_c, s_c, pkg):
    """
    Generate a minimalist run configuration file in the execution directory
    for successive run steps.

    Returns job submission info dictionnary
    """

    sections = {}

    if 'coupling' in r_c:
        if r_c['coupling']:
            sections['setup'] = {'coupling': r_c['coupling']}

    sections['run'] = {'id': r_c['run_id'],
                       'stage': False,
                       'initialize': s_c['initialize'],
                       'compute': s_c['run_solver'],
                       'finalize': s_c['save_results']}

    r_d = {}
    for kw in ('n_procs', 'n_threads', 'time_limit'):
        if r_c[kw]:
            r_d[kw] = r_c[kw]

    for kw in ('job_parameters', 'job_header',
               'compute_prologue', 'compute_epilogue'):
        if r_c[kw]:
            r_d[kw] = r_c[kw]

    sections[resource_name] = r_d

    run_conf = cs_run_conf.run_conf(None)
    run_conf.sections = sections
    run_conf.save(path, new=True)

#===============================================================================
# Run the calculation
#===============================================================================

def run(argv=[], pkg=None, submit_args=None):
    """
    Run calculation;
    returns return code, run id, and results directory path when created.
    """

    if not pkg:
        from code_saturne.cs_package import package
        pkg = package()

    i_c = cs_run_conf.get_install_config_info(pkg)

    r_c, s_c, run_conf = process_cmd_line(argv, pkg)

    casedir = r_c['casedir']

    if not casedir:
        return 1, None, None

    # Read run configuration

    read_run_config_file(i_c, r_c, s_c, pkg, run_conf=run_conf)

    # Determine compute stages

    update_run_steps(s_c, None, final=True)

    stages = {'prepare_data': s_c['stage'],
              'initialize': s_c['initialize'],
              'run_solver': s_c['compute'],
              'save_results': s_c['finalize']}

    # Use alternate compute (back-end) package if defined

    pkg_compute = None
    if not r_c['compute_build']:
        if i_c['compute_builds']:
            r_c['compute_build'] = i_c['compute_builds'][0]
    if r_c['compute_build']:
        pkg_compute = pkg.get_alternate_version(r_c['compute_build'])

    # Specific case for coupling

    if r_c['coupling']:

        from code_saturne import cs_case_coupling

        coupling =r_c['coupling']
        if os.path.isfile(coupling):
            try:
                c_locals = {}
                exec(compile(open(coupling).read(), coupling, 'exec'),
                     globals(),
                     c_locals)
                domains = c_locals['domains']
            except Exception:
                execfile(coupling)

        verbose = True
        if r_c['suggest_id']:
            verbose = False

        c = cs_case_coupling.coupling(pkg,
                                      domains,
                                      casedir,
                                      staging_dir=r_c['staging_dir'],
                                      verbose=verbose,
                                      package_compute=pkg_compute)

    else:
        # Values in case and associated domain set from parameters
        d = cs_case_domain.domain(pkg,
                                  package_compute=pkg_compute,
                                  param=r_c['param'])

        # Now handle case for the corresponding calculation domain(s).
        c = cs_case.case(pkg,
                         package_compute=pkg_compute,
                         case_dir=casedir,
                         staging_dir=r_c['staging_dir'],
                         domains=d)

    # Determine run id if not forced

    if not r_c['run_id']:
        r_c['run_id'] = c.suggest_id(r_c['id_prefix'], r_c['id_suffix'])

    if r_c['suggest_id']:
        print(r_c['run_id'])
        return 0, r_c['run_id'], None

    if submit_args != None:
        submit_stages = {'prepare_data': False}
        for s in ('initialize', 'run_solver', 'save_results'):
            submit_stages[s] = stages[s]
            stages[s] = False

    # Now run case

    retval = c.run(n_procs=r_c['n_procs'],
                   n_threads=r_c['n_threads'],
                   run_id=r_c['run_id'],
                   force_id=r_c['force_id'],
                   stages=stages)

    if submit_args != None:
        resource_name = cs_run_conf.get_resource_name(i_c)
        run_cfg_path = os.path.join(c.result_dir, 'run.cfg')
        if len(submit_args) > 1:
            job_parameters = cs_exec_environment.assemble_args(submit_args)
            r_c['job_parameters'] = job_parameters
        generate_run_config_file(run_cfg_path, resource_name,
                                 r_c, submit_stages, pkg)

    return retval, c.result_dir, r_c

#===============================================================================
# Main function
#===============================================================================

def main(argv, pkg):
    """
    Main function.
    """
    return run(argv, pkg)[0]

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    retval = main(argv=sys.argv[1:])

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
