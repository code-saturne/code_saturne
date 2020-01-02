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

import datetime
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
        usage = "usage: %prog [options] <file_name>"
    else:
        usage = "usage: %prog bdump [options] <file_name>"

    usage += """

Dump headers and optionnaly content of a """ + pkg.code_name + """ Preprocessor,
Partitioner, or restart file.
"""

    parser = OptionParser(usage=usage)

    parser.add_option("-e", "--extract", dest="extract",
                      action="store_true",
                      help="extract mode (extract full section data, with " \
                          + "no metadata).")

    parser.add_option("--f-format", dest="f_format", type="string",
                      metavar="<fmt>",
                      help="define format for floating-point numbers (default: " \
                          + "\"15.9e\" for floats, \"22.15e\" for doubles).")

    parser.add_option("--location", dest="location", type="string",
                      metavar="<id>",
                      help="only output section(s) with given location id.")

    parser.add_option("-n", dest="level", type="int",
                      metavar="<level>",
                      help="number of first and last elements of each section " \
                          + "to output (default: print headers only).")

    parser.add_option("--section", dest="section", type="string",
                      metavar="<name>",
                      help="only consider section matching given criteria.")

    parser.set_defaults(extract=False)
    parser.set_defaults(f_format=None)
    parser.set_defaults(location=None)
    parser.set_defaults(level=None)
    parser.set_defaults(section=None)

    (options, args) = parser.parse_args(argv)

    if len(args) == 1:
        if args[0] != argv[len(argv) - 1]:
            args = None
    else:
            args = None

    if not args:
        parser.print_help()

    return  options, args

#===============================================================================
# Run the utility
#===============================================================================

def main(argv, pkg):
    """
    Main function.
    """

    options, args = process_cmd_line(argv, pkg)

    if not args:
        return 1

    cmd = [pkg.get_io_dump()]
    cmd += argv

    # Set environment modules if present

    cs_exec_environment.set_modules(pkg)
    cs_exec_environment.source_rcfile(pkg)

    # Run executable

    return cs_exec_environment.run_command(cmd)

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

