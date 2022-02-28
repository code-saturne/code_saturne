#!/usr/bin/env python3

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

"""
This module describes the script used to run a study/case for code_saturne.
"""

#===============================================================================
# Import required Python modules
#===============================================================================

import os, sys
import types, string, re, fnmatch
import argparse

from code_saturne.base import cs_exec_environment
from code_saturne.base import cs_case_domain
from code_saturne.base import cs_case
from code_saturne.base import cs_run_conf

#-------------------------------------------------------------------------------
# Update run steps based on run_conf object
#-------------------------------------------------------------------------------

def update_run_steps(s_c, run_conf, final=False):
    """
    Process the passed command line arguments.
    """

    filter_stages = False
    for k in s_c:
        if s_c[k] != None and s_c[k] != '':
            filter_stages = True

    if run_conf and not filter_stages:
        if 'run' in run_conf.sections:
            for kw in s_c:
                s_c[kw] = run_conf.get_bool('run', kw)

    filter_stages = False
    for k in s_c:
        if s_c[k] != None and s_c[k] != '':
            filter_stages = True

    # Default if nothing provided, ensure range is filled otherwise

    if filter_stages:

        stages = ('stage', 'initialize', 'compute', 'finalize')

        stage_ini = s_c['stage']
        i_s = -1
        i_f = -1
        for i, k in enumerate(stages):
            if s_c[k] == True:
                if i_s < 0:
                    i_s = i
                i_f = i + 1
            elif s_c[k] == False:
               if i_s < 0 or i_s == i:
                    i_s = i + 1
        if i_f < 0:
            i_f = len(stages)
        for i, k in enumerate(stages):
            if i < i_s:
                s_c[k] = False
            elif i < i_f:
                s_c[k] = True
            else:
                s_c[k] = False
        # Special handling of stage defaults, as it is a substage
        # Of initialization
        if stage_ini is None:
            if s_c['initialize'] == True:
                s_c['stage'] = True
            else:
                s_c['stage'] = False

    elif final:
        for i, k in enumerate(s_c):
            if s_c[k] != False:
                s_c[k] = True

#-------------------------------------------------------------------------------
# Clean append for command-line arguments parser
#-------------------------------------------------------------------------------

class multi_append(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, list())
        for value in values:
            v_args = cs_exec_environment.separate_args(value)
            for sv in v_args:
                getattr(namespace, self.dest).append(sv)

#-------------------------------------------------------------------------------

class multi_append_kv(argparse.Action):
    """
    Parse and append key-value pairs.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest) is None:
            setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

#-------------------------------------------------------------------------------
# Build command-line arguments parser
#-------------------------------------------------------------------------------

def arg_parser(argv):
    """
    Process the passed command line arguments.
    """

    prog = os.path.basename(sys.argv[0]) + " " + sys.argv[1]
    parser = argparse.ArgumentParser(description="Run a case or specified run stages.",
                                     prog=prog)

    parser.add_argument("--compute-build", dest="compute_build", type=str,
                        metavar="<build>",
                        help="base name or full path to the compute build")

    parser.add_argument("-n", "--nprocs" , "--n-procs", dest="nprocs", type=int,
                        metavar="<nprocs>",
                        help="number of MPI processes for the computation")

    parser.add_argument("--nt", "--threads-per-task", dest="nthreads", type=int,
                        help="number of OpenMP threads per task")

    parser.add_argument("--notebook-args", nargs='*', action = multi_append_kv,
                        help="key=value pairs to pass to user scripts and notebook")

    parser.add_argument("--parametric-args", nargs='*', action = multi_append,
                        help="key=value pairs to pass to cs_parametric filter")

    parser.add_argument("--kw-args", nargs='*', action = multi_append,
                        help="additional keywords to pass to user scripts")

    parser.add_argument("--debug-args", nargs='*', action = multi_append,
                        help="additional commands and keywords for " \
                        + "debug wrapper (see developer guide for options)")

    parser.add_argument("--tool-args", nargs='*', action = multi_append,
                        help="additional commands and keywords to run tool " \
                        + "with solver (such as profiler)")

    parser.add_argument("--mpi-tool-args", nargs='*', action = multi_append,
                        help="additional commands and keywords to run tool " \
                        + "(such as debugger or profiler) prefixing MPI launch")

    parser.add_argument("-p", "--param", dest="param", type=str,
                        metavar="<param>",
                        help="path or name of the parameters file")

    parser.add_argument("--case", dest="case", type=str,
                        metavar="<case>",
                        help="path to the case's directory")

    parser.add_argument("--dest", dest="dest", type=str,
                        metavar="<dest>",
                        help="path to the destination top directory")

    parser.add_argument("--id", dest="id", type=str,
                        metavar="<id>",
                        help="use the given run id")

    parser.add_argument("--id-prefix", dest="id_prefix", type=str,
                        metavar="<prefix>",
                        help="prefix the run id with the given string")

    parser.add_argument("--id-suffix", dest="id_suffix", type=str,
                        metavar="<suffix>",
                        help="suffix the run id with the given string")

    parser.add_argument("--suggest-id", dest="suggest_id",
                        action="store_true",
                        help="suggest a run id for the next run")

    parser.add_argument("--force", dest="force",
                        action="store_true",

                        help="run the data preparation stage even if " \
                              + "the matching execution directory exists")

    parser.add_argument("--stage", dest="stage",
                        action="store_true",
                        help="stage data prior to preparation and execution")

    parser.add_argument("--no-stage", dest="stage",
                        action="store_false",
                        help="do not stage data prior to preparation and execution")

    parser.add_argument("--initialize", dest="initialize",
                        action="store_true",
                        help="run the data preparation step")

    parser.add_argument("--compute", "--execute", dest="compute",
                        action="store_true",
                        help="run the compute step")

    parser.add_argument("--finalize", dest="finalize",
                        action="store_true",
                        help="run the results copy/cleanup step")

    parser.add_argument("--with-resource", dest="resource_name", type=str,
                        metavar="<resource>",
                        help="use resource settings based on given name")

    parser.set_defaults(compute_build=False)
    parser.set_defaults(suggest_id=False)
    parser.set_defaults(stage=None)
    parser.set_defaults(initialize=None)
    parser.set_defaults(compute=None)
    parser.set_defaults(finalize=None)
    parser.set_defaults(param=None)
    parser.set_defaults(domain=None)
    parser.set_defaults(id=None)
    parser.set_defaults(nprocs=None)
    parser.set_defaults(nthreads=None)
    parser.set_defaults(resource=None)

    return parser

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def parse_cmd_line(argv):
    """
    Process the passed command line arguments.
    """

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    parser = arg_parser(argv)
    options = parser.parse_args(argv)

    return options

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_options(options, pkg):
    """
    Process the passed command line arguments.
    """

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

    # Also check for possible settings file

    run_id = options.id
    run_conf = None

    run_config_path = os.path.join(os.getcwd(), 'run.cfg')
    if os.path.isfile(run_config_path):
        run_conf = cs_run_conf.run_conf(run_config_path, package=pkg)
        if not run_id and not filter_stages:
            if 'run' in run_conf.sections:
                update_run_steps(s_c, run_conf)
                if s_c['stage'] == False:
                    if 'id' in run_conf.sections['run']:
                        run_id = run_conf.sections['run']['id']

    if options.stage == False:    # Handle --no-stage option
        s_c['stage'] = False

    if s_c['stage'] == False and not run_id:
        err_str = os.linesep + os.linesep + 'Error:' + os.linesep
        err_str += 'Incompatible options in the run.cfg file or command arguments'
        err_str += os.linesep
        err_str += 'When the "stage" step is set to False, a run id is required.'
        raise cs_case_domain.RunCaseError(err_str)

    # Check for multiple domain case
    # Kept the 'coupling' file def for the definition of the case_dir function
    coupling = None
    coupled_domains =[]

    if run_conf:
        coupled_domains = run_conf.get_coupling_parameters()

    if coupled_domains != []:
        coupling = run_config_path

    if coupling and options.param:
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        err_str = os.linesep + os.linesep + 'Error:' + os.linesep
        err_str += cmd_line + os.linesep
        err_str += '-p/--param option is incompatible with '
        err_str += '"coupled_domains" option defined within the run.cfg file.'
        raise cs_case_domain.RunCaseError(err_str)

    casedir, staging_dir = cs_case.get_case_dir(case=options.case,
                                                param=options.param,
                                                coupling=coupling,
                                                id=run_id)

    if casedir is None and staging_dir is None:
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        print('Error:', file = sys.stderr)
        print(cmd_line, file = sys.stderr)
        print('run from directory \"' + str(os.getcwd()) + '\",', file = sys.stderr)
        print('which does not seem to be inside a case directory.', file = sys.stderr)

    param = options.param

    compute_build = options.compute_build

    if not options.force:
        force_id = False
    else:
        force_id = True

    # Return associated dictionary (also force number of ranks and threads)

    r_c = {'casedir': casedir,
           'dest_dir': options.dest,
           'staging_dir': staging_dir,
           'run_id': run_id,
           'param': param,
           'coupled_domains': coupled_domains,
           'id_prefix': options.id_prefix,
           'id_suffix': options.id_suffix,
           'suggest_id': options.suggest_id,
           'force_id': force_id,
           'n_procs': options.nprocs,
           'n_threads': options.nthreads,
           'time_limit': None,
           'compute_build': compute_build,
           'debug_args': None,
           'tool_args': None,
           'mpi_tool_args': None}

    if options.debug_args:
        r_c['debug_args'] = cs_exec_environment.assemble_args(options.debug_args)
    if options.tool_args:
        r_c['tool_args'] = cs_exec_environment.assemble_args(options.tool_args)
    if options.mpi_tool_args:
        r_c['mpi_tool_args'] = cs_exec_environment.assemble_args(options.mpi_tool_args)

    return r_c, s_c, run_conf

#-------------------------------------------------------------------------------
# Read the run configuration file
#-------------------------------------------------------------------------------

def read_run_config_file(i_c, r_c, s_c, pkg, run_conf=None):
    """
    Process the passed command line arguments.
    """

    run_config_path = ""

    if r_c['staging_dir']:
        run_config_path = os.path.join(r_c['staging_dir'], 'run.cfg')

    elif r_c['casedir']:
        casedir = r_c['casedir']
        if r_c['coupled_domains'] != []:
            run_config_path = os.path.join(casedir, 'run.cfg')
        else:
            run_config_path = os.path.join(casedir, 'DATA', 'run.cfg')

    # Ensure some keys are set in all cases to simplify future tests

    run_conf_kw = ('job_parameters', 'job_header',
                   'run_prologue', 'run_epilogue',
                   'compute_prologue', 'compute_epilogue',
                   'debug_args', 'tool_args', 'mpi_tool_args')

    for kw in run_conf_kw:
        if not kw in r_c:
            r_c[kw] = None

    if run_conf is None:
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

    # Case path if not determined yet
    # (when destination or staging directory is outside case directory)

    if 'paths' in run_conf.sections:
        if not r_c['casedir']:
            if 'case' in run_conf.sections['paths']:
                r_c['casedir'] = run_conf.sections['paths']['case']
        if not r_c['dest_dir']:
            if 'top_results_directory' in run_conf.sections['paths']:
                r_c['dest_dir'] = run_conf.sections['paths']['top_results_directory']

    # Parameters file

    for kw in ('param',):
        if run_conf.get('setup', kw):
            r_c[kw] = run_conf.get('setup', kw)

    # Run id

    if not r_c['run_id']:
        r_c['run_id'] = run_conf.get('run', 'id')

    if not r_c['force_id']:
        r_c['force_id'] = run_conf.get_bool('run', 'force_id')

    # Compute stages

    update_run_steps(s_c, run_conf)

    # Resources: try to find a matching section, using
    # resource_name, batch, and job_defaults in decreasing priority.

    if not r_c['compute_build']:
        r_c['compute_build'] = run_conf.get_bool('run', 'compute_build')

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

    if path and r_c['casedir']:
        in_case = cs_case.is_exec_dir_in_case(r_c['casedir'],
                                              os.path.dirname(path))
        if not in_case:
            sections['paths'] = {'case': r_c['casedir']}
            if r_c['dest_dir']:
                sections['paths']['top_results_directory'] = r_c['dest_dir']

    if 'coupled_domains' in r_c:
        if r_c['coupled_domains'] != []:
            dom_str=''
            for i, d in enumerate(r_c['coupled_domains']):
                dom_name = d['domain']
                if i != 0:
                    dom_str += ':'
                dom_str += dom_name

                # Add the domain section
                sections[dom_name] = {key:str(d[key]) for key in d.keys()}

            sections['setup'] = {'coupled_domains': dom_str}

    sections['run'] = {'compute_build': r_c['compute_build'],
                       'id': r_c['run_id'],
                       'stage': False,
                       'initialize': s_c['initialize'],
                       'compute': s_c['run_solver'],
                       'finalize': s_c['save_results']}

    r_d = {}
    for kw in ('n_procs', 'n_threads', 'time_limit'):
        if r_c[kw]:
            r_d[kw] = r_c[kw]

    for kw in ('job_parameters', 'job_header',
               'compute_prologue', 'compute_epilogue',
               'debug_args', 'tool_args', 'mpi_tool_args'):
        if r_c[kw]:
            r_d[kw] = r_c[kw]

    sections[resource_name] = r_d

    run_conf = cs_run_conf.run_conf(None)
    run_conf.sections = sections
    run_conf.save(path, new=True)

#===============================================================================
# Run the calculation
#===============================================================================

def run(argv=[], pkg=None, run_args=None, submit_args=None):
    """
    Run calculation;
    returns return code, run id, and results directory path when created.
    """

    if not pkg:
        from code_saturne.base.cs_package import package
        pkg = package()

    if run_args is None:
        options = parse_cmd_line(argv)
    else:
        options = run_args

    i_c = cs_run_conf.get_install_config_info(pkg)

    if options.resource_name != None:
        i_c['resource_name'] = options.resource_name.lower()

    r_c, s_c, run_conf = process_options(options, pkg)

    if not r_c['casedir'] and not r_c['staging_dir']:
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

    if r_c['coupled_domains'] != []:

        from code_saturne.base import cs_case_coupling

        domains = r_c['coupled_domains']

        verbose = True
        if r_c['suggest_id']:
            verbose = False

        c = cs_case_coupling.coupling(pkg,
                                      domains,
                                      r_c['casedir'],
                                      r_c['dest_dir'],
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
                         case_dir=r_c['casedir'],
                         dest_dir=r_c['dest_dir'],
                         staging_dir=r_c['staging_dir'],
                         domains=d)

    # Determine run id if not forced

    if not r_c['run_id']:
        r_c['run_id'] = c.suggest_id(r_c['id_prefix'], r_c['id_suffix'])

    if r_c['suggest_id']:
        print(r_c['run_id'])
        return 0, r_c['run_id'], None

    # Check stages

    if submit_args != None:
        submit_stages = {'prepare_data': False}
        for s in ('initialize', 'run_solver', 'save_results'):
            submit_stages[s] = stages[s]
            stages[s] = False

    # Set script hooks

    c.run_prologue = r_c['run_prologue']
    c.run_epilogue = r_c['run_epilogue']
    c.compute_prologue = r_c['compute_prologue']
    c.compute_epilogue = r_c['compute_epilogue']

    c.debug_args = r_c['debug_args']
    c.tool_args = r_c['tool_args']
    c.mpi_tool_args = r_c['mpi_tool_args']

    # Now run case

    retval = c.run(n_procs=r_c['n_procs'],
                   n_threads=r_c['n_threads'],
                   time_limit=r_c['time_limit'],
                   run_id=r_c['run_id'],
                   force_id=r_c['force_id'],
                   stages=stages,
                   notebook_args=options.notebook_args,
                   parametric_args=options.parametric_args,
                   kw_args=options.kw_args)

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
