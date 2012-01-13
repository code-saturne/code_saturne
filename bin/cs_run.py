#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

import datetime
import os, sys, pwd
import types, string, re, fnmatch
from optparse import OptionParser
import ConfigParser

import cs_exec_environment
import cs_case_domain
import cs_case

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-p", "--param", dest="param", type="string",
                      metavar="<param>",
                      help="path or name of the parameters file")

    parser.add_option("--case", dest="case", type="string",
                      metavar="<case>",
                      help="path to the case's directory")

    parser.add_option("--id", dest="id", type="string",
                      metavar="<id>",
                      help="use the given run id")

    parser.add_option("--suggest-id", dest="suggest_id",
                      action="store_true",
                      help="suggest a run id for the next run")

    parser.add_option("--initialize", dest="initialize",
                      action="store_true",
                      help="run the data preparation stage")

    parser.add_option("--execute", dest="execute",
                      action="store_true",
                      help="run the execution stage")

    parser.add_option("--finalize", dest="finalize",
                      action="store_true",
                      help="run the results copy/cleanup stage")

    parser.set_defaults(suggest_id=False)
    parser.set_defaults(initialize=False)
    parser.set_defaults(execute=False)
    parser.set_defaults(finalize=False)
    parser.set_defaults(param=None)
    parser.set_defaults(domain=None)
    parser.set_defaults(id=None)

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    (options, args) = parser.parse_args(argv)

    # Try to determine case directory

    casedir = None
    param = None

    if options.param:
        param = os.path.basename(options.param)
        if param != options.param:
            datadir = os.path.split(os.path.realpath(options.param))[0]
            (casedir, data) = os.path.split(datadir)
            if data != 'DATA': # inconsistent paramaters location.
                casedir = None

    if options.case:
        casedir = os.path.realpath(options.case)

    if not casedir:
        casedir = os.getcwd()
        while os.path.basename(casedir):
            data = os.path.join(casedir, 'DATA')
            src = os.path.join(casedir, 'SRC')
            if os.path.isdir(data) and os.path.isdir(src):
                break
            casedir = os.path.split(casedir)[0]

    if not (os.path.isdir(data) and os.path.isdir(src)):
        casedir = None
        cmd_line = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd_line += ' ' + arg
        err_str = 'Error:\n' + cmd_line + '\n' \
            'run from directory \"' + str(os.getcwd()) + '\",\n' \
            'which does not seem to be inside a case directory.\n'
        sys.stderr.write(err_str)

    # Stages to run (if no filter given, all are done).

    run_id = options.id

    suggest_id = options.suggest_id

    prepare_data = options.initialize
    run_solver = options.execute
    save_results = options.finalize

    if not (prepare_data or run_solver or save_results):
        prepare_data = True
        run_solver = True
        save_results = True

    return  (casedir, run_id, param, suggest_id,
             prepare_data, run_solver, save_results)

#===============================================================================
# Run the calculation
#===============================================================================

def main(argv, pkg):
    """
    Main function.
    """

    (casedir, run_id, param, suggest_id,
     prepare_data, run_solver, save_results) = process_cmd_line(argv, pkg)

    if not casedir:
        return 1

    if suggest_id:
        now = datetime.datetime.now()
        run_id = now.strftime('%Y%m%d-%H%M')
        print(run_id)
        return 0

    # Use alternate compute (back-end) package if defined

    config = ConfigParser.ConfigParser()
    config.read([pkg.get_configfile()])

    pkg_compute = None
    if config.has_option('install', 'compute_versions'):
        compute_versions = config.get('install', 'compute_versions').split(':')
        if compute_versions[0]:
            pkg_compute = pkg.get_alternate_version(compute_versions[0])

    # Values in case and associated domain set from parameters

    d = cs_case_domain.domain(pkg, package_compute=pkg_compute, param=param)

    # Now handle case for the corresponding calculation domain(s).

    c = cs_case.case(pkg,
                     package_compute=pkg_compute,
                     case_dir=casedir,
                     domains=d)

    # Now run case

    retval = c.run(mpi_environment=None,
                   run_id=run_id,
                   prepare_data=prepare_data,
                   run_solver=run_solver,
                   save_results=save_results)

    return retval

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Retrieve package information (name, version, installation dirs, ...)
    from cs_package import package
    pkg = package()

    retval = main(sys.argv[1:], pkg)

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

