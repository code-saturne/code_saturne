#!/usr/bin/env python
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009  EDF
#
#   The Code_Saturne Preprocessor is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Preprocessor is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import os
import sys
import tempfile
from cs_exec_environment import run_command
from optparse import OptionParser

#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string the
             string is split to create a list of arguments.
    """

    parser = OptionParser(usage="Usage: %prog <monitoring file>")

    (options, args) = parser.parse_args(argv)

    if len(args) == 1:
        file = args[0]
    else:
        parser.print_help()
        sys.exit(1)

    return file

#-------------------------------------------------------------------------------

def plot(f):
    """
    Plot a monitoring file with XmGrace.
    """
    retval = 0

    # Create a temporary file and correctly format the monitoring file

    temp_file = tempfile.mkstemp(suffix=".hst")

    cmd = "grep -v '#'" + " " + str(f) + "|cut -c9- >" + " " + temp_file[1]
    run_command(cmd)

    # Run XmGrace -nxy command

    cmd = "xmgrace -geometry 1100x850 -nxy" + " " + temp_file[1]
    run_command(cmd)

    # Cleanup and exit

    os.remove(temp_file[1])

    return retval

#-------------------------------------------------------------------------------

def main(argv, main):
    """
    Main function.
    """

    file = process_cmd_line(argv, pkg)

    retcode = plot(file)

    sys.exit(retcode)

if __name__ == '__main__':
    main(sys.argv[1:], None)
