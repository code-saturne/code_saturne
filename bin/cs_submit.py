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

import os, sys, subprocess
import types, string

from code_saturne import cs_batch
from code_saturne.cs_exec_environment import enquote_arg, get_shell_type
from code_saturne import cs_run
from code_saturne import cs_runcase

#-------------------------------------------------------------------------------
# Print a help page.
#-------------------------------------------------------------------------------

def print_help(submit_cmd, pkg):
    """
    Print a help page.
    """

    help_string = \
"""
Usage: %s [batch options] <runcase_script> [script args]

Options:

  Any option provided by the batch submission command (%s)

Runcase script:

  Shell script containing "%s run" command.
"""
    print(help_string % (sys.argv[0], submit_cmd, pkg.name))

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, submit_cmd, pkg):
    """
    Process the passed command line arguments.
    """

    runcase_path = None
    need_help = False

    if "-h" in argv or "--help" in argv:
        need_help = True
    elif len(argv) < 1:
        need_help = True
    else:
        for a in argv:
            if os.path.isfile(a):
                runcase_path = a
                break

    if not runcase_path:
        need_help = True

    if need_help:
        print_help(submit_cmd, pkg)

    return runcase_path

#===============================================================================
# Run the calculation
#===============================================================================

def main(argv, pkg):
    """
    Main function.
    """

    # Use alternate compute (back-end) package if defined

    batch = cs_batch.batch(pkg)

    submit_cmd = batch.submit_command_prefix()

    if not submit_cmd:
        submit_cmd = get_shell_type()

    runcase_path = process_cmd_line(argv, submit_cmd, pkg)

    if not runcase_path:
        return 1

    scripts_dir = os.path.abspath(os.path.dirname(runcase_path))

    runcase = cs_runcase.runcase(runcase_path)

    # Adjust run command for staging only

    run_args = runcase.get_run_args()

    for a in ['--stage', '--initialize', '--execute', '--finalize']:
        while a in run_args:
            run_args.remove(a)

    run_args.append('--stage')

    retcode, run_id, result_path = cs_run.run(run_args, pkg)

    # Now prepare runcase for submission

    if retcode != 0 or not result_path:
        err_str = ' run not staged due to prior error.'
        raise Exception(err_str)

    run_args = runcase.get_run_args()

    while '--stage' in run_args:
        run_args.remove(a)

    stages = []
    for a in ['--initialize', '--execute', '--finalize']:
        if a in run_args:
            stages.append(a)

    if not stages:
        run_args.append('--initialize')
        run_args.append('--finalize')

    runcase.set_run_args(run_args)
    runcase.set_run_id(run_id=run_id)

    runcase_path = os.path.join(result_path,
                                os.path.basename(runcase_path))

    runcase.save(runcase_path)

    # Now submit case

    save_wd = os.getcwd()
    os.chdir(result_path)

    for a in argv:
        if os.path.isfile(a) and runcase_path:
            submit_cmd += ' ' + enquote_arg(runcase_path)
            runcase_path = None
        else:
            submit_cmd += ' ' + enquote_arg(a)

    return subprocess.call(submit_cmd, shell=True)

    os.chdir(save_wd)

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

