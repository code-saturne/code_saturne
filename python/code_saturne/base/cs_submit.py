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

This module defines the following functions:
- process_cmd_line
- main
"""

#===============================================================================
# Import required Python modules
#===============================================================================

import os, sys, subprocess
import types, string
from argparse import ArgumentParser

from code_saturne.base import cs_batch
from code_saturne.base.cs_exec_environment import enquote_arg, get_shell_type
from code_saturne.base import cs_run
from code_saturne.base import cs_runcase

#-------------------------------------------------------------------------------
# Print a help page.
#-------------------------------------------------------------------------------

def print_help(submit_cmd, pkg):
    """
    Print a help page.
    """

    help_string = \
"""
Usage: %s [batch options]

Options:
   Any option provided by the batch submission command (%s)
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

    if need_help:
        print_help(submit_cmd, pkg)

    return

#===============================================================================
# Sumbit the calculation
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

    epilog = ("Options not listed above are passed to the batch "
              "submission commands; see your batch documentation for this.")

    run_parser = cs_run.arg_parser(argv)
    prog = os.path.basename(sys.argv[0]) + " " + sys.argv[1]

    parser = ArgumentParser(parents=[run_parser],
                            prog=prog,
                            description="Submit a case or specified run stages.",
                            usage='%(prog)s [run options] [batch options]',
                            epilog=epilog,
                            conflict_handler='resolve')

    run_args, submit_args = parser.parse_known_args(argv)

    retcode, result_path, r_c = cs_run.run(pkg=pkg,
                                           run_args=run_args,
                                           submit_args=submit_args)

    # Now prepare runcase for submission

    if retcode != 0 or not result_path:
        err_str = ' run not staged due to prior error.'
        raise Exception(err_str)

    runcase_path = os.path.join(result_path, 'runcase')
    runcase = cs_runcase.runcase(runcase_path,
                                 package=pkg,
                                 submit=True,
                                 job_header = r_c['job_header'],
                                 prologue = r_c['run_prologue'],
                                 epilogue = r_c['run_epilogue'])
    runcase.save()

    # Now submit case

    save_wd = os.getcwd()
    os.chdir(result_path)

    if r_c['job_parameters']:
        submit_cmd += ' ' + r_c['job_parameters']

    submit_cmd += ' ./runcase'

    return subprocess.call(submit_cmd, shell=True)

    os.chdir(save_wd)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Run package
    from code_saturne.base.cs_package import package
    pkg = package()

    retval = main(sys.argv[1:], pkg)

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
