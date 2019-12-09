#!/usr/bin/env python

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

"""
This module describes the script used to run a study/case for Code_Saturne.

This module defines the following functions:
- process_cmd_line
- main
"""

#===============================================================================
# Import required Python modules
#===============================================================================

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

    parser.add_option("--initialize", dest="initialize",
                      action="store_true",
                      help="run the data preparation stage")

    parser.add_option("--execute", dest="execute",
                      action="store_true",
                      help="run the execution stage")

    parser.add_option("--finalize", dest="finalize",
                      action="store_true",
                      help="run the results copy/cleanup stage")

    parser.set_defaults(compute_build=False)
    parser.set_defaults(suggest_id=False)
    parser.set_defaults(stage=False)
    parser.set_defaults(initialize=False)
    parser.set_defaults(execute=False)
    parser.set_defaults(finalize=False)
    parser.set_defaults(param=None)
    parser.set_defaults(coupling=None)
    parser.set_defaults(domain=None)
    parser.set_defaults(id=None)
    parser.set_defaults(nprocs=None)
    parser.set_defaults(nthreads=None)

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    (options, args) = parser.parse_args(argv)

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
        err_str = 'Error:\n' + cmd_line + '\n' \
                  '--coupling and -p/--param options are incompatible.\n'
        sys.stderr.write(err_str)
        sys.exit(1)

    casedir, staging_dir = cs_case.get_case_dir(case=options.case,
                                                param=options.param,
                                                coupling=options.coupling,
                                                id=options.id)

    if casedir == None:
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        err_str = 'Error:\n' + cmd_line + '\n' \
                  'run from directory \"' + str(os.getcwd()) + '\",\n' \
                  'which does not seem to be inside a case directory.\n'
        sys.stderr.write(err_str)

    param = options.param

    # If no parameter file passed, and a setup.xml is present in DATA, run it
    if param is None:
        has_setup = os.path.isfile(os.path.join(casedir, 'DATA', 'setup.xml'))
        if has_setup:
            param = "setup.xml"

    # Stages to run (if no filter given, all are done).

    compute_build = options.compute_build

    stages = {'prepare_data':options.stage,
              'initialize':options.initialize,
              'run_solver':options.execute,
              'save_results':options.finalize}

    if not options.force:
        force_id = False
    else:
        force_id = True

    # Stages to run (if no filter given, all are run; specific stages
    # are given, all stages in thar range are run).

    ordered_stages = ['prepare_data',
                      'initialize',
                      'run_solver',
                      'save_results']

    stages = {'prepare_data':options.stage,
              'initialize':options.initialize,
              'run_solver':options.execute,
              'save_results':options.finalize}

    stages_start = len(ordered_stages)
    stages_end = -1

    i = 0
    for k in ordered_stages:
        if stages[k] == True and stages_start > i:
            stages_start = i
        if stages[k] and stages_end < i+1:
            stages_end = i + 1
        i += 1

    # Default if nothing provided, ensure range is filled otherwise

    if stages_end < 0:
        for k in ordered_stages:
            stages[k] = True

    else:
        for i in range(stages_start + 1, stages_end -1):
            stages[ordered_stages[i]] = True

    # Forced number of ranks and threads

    n_procs = options.nprocs
    n_threads = options.nthreads

    return  (casedir, staging_dir, options.id, param, options.coupling,
             options.id_prefix, options.id_suffix, options.suggest_id, force_id,
             n_procs, n_threads, stages, compute_build)

#===============================================================================
# Run the calculation
#===============================================================================

def run(argv, pkg):
    """
    Run calculation;
    returns return code, run id, and results directory path when created.
    """

    (casedir, staging_dir, run_id, param, coupling,
     id_prefix, id_suffix, suggest_id, force, n_procs, n_threads,
     stages, compute_build) = process_cmd_line(argv, pkg)

    if not casedir:
        return 1, None, None

    # Use alternate compute (back-end) package if defined

    config = configparser.ConfigParser()
    config.read(pkg.get_global_configfile())

    pkg_compute = None
    if not compute_build:
        if config.has_option('install', 'compute_versions'):
            compute_versions = config.get('install', 'compute_versions').split(':')
            if compute_versions[0]:
                pkg_compute = pkg.get_alternate_version(compute_versions[0])
    else:
        pkg_compute = pkg.get_alternate_version(compute_build)

    if coupling:

        # Specific case for coupling
        from code_saturne import cs_case_coupling

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
        if suggest_id:
            verbose = False

        c = cs_case_coupling.coupling(pkg,
                                      domains,
                                      casedir,
                                      staging_dir=staging_dir,
                                      verbose=verbose,
                                      package_compute=pkg_compute)

    else:
        # Values in case and associated domain set from parameters
        d = cs_case_domain.domain(pkg, package_compute=pkg_compute, param=param)

        # Now handle case for the corresponding calculation domain(s).
        c = cs_case.case(pkg,
                         package_compute=pkg_compute,
                         case_dir=casedir,
                         staging_dir=staging_dir,
                         domains=d)

    # Determine run id if not forced

    if not run_id or suggest_id:
        run_id = c.suggest_id(id_prefix, id_suffix)

        if suggest_id:
            print(run_id)
            return 0, run_id, None

    # Now run case

    retval = c.run(n_procs=n_procs,
                   n_threads=n_threads,
                   run_id=run_id,
                   force_id=force,
                   stages=stages)

    return retval, c.run_id, c.result_dir

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

    # Run package
    from code_saturne.cs_package import package
    pkg = package()

    retval = main(sys.argv[1:], pkg)

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
