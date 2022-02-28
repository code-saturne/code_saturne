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
This modules describes the script used to translate runtime symbols
obtained from a stack into the corresponding file and line.
"""

#===============================================================================
# Import required Python modules
#===============================================================================

import sys
from optparse import OptionParser

from code_saturne.base.cs_exec_environment import get_command_output

#-------------------------------------------------------------------------------
# Print a help page.
#-------------------------------------------------------------------------------

def print_help(pkg):
    """
    Print a help page.
    """

    help_string = \
"""
This is a debug symbol to file/line utility. Usage:

%s symbol2line -s symbol1 [-s symbol2, .. -s symbolN] [optional arguments]
Example:
If the stack in the listing/error.log file contains the following line :
1: 0x55834c6102d7 <cs_function+0x17>             (cs_solver)

then the command to launch is : '%s symbol2line -s cs_function+0x17 [optional arguments]'
"""

    print(help_string % (pkg.name, pkg.name))

#-------------------------------------------------------------------------------
# Process command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.
    """

    for arg in argv:
        if arg == "-h" or arg == "--help":
            print_help(pkg)

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog compile [options]"

    parser = OptionParser(usage=usage)

    parser.add_option('-s', '--symbol', dest='symbols', type='string',
                      action='append',
                      help="Symbols to translate into file:line")

    parser.add_option('-p', '--path', dest='path', type='string',
                      help="Optional path to executable")

    parser.add_option('-e', '--executable', dest='solver', type='string',
                      help="Optional name of executable.")

    parser.set_defaults(path=None)
    parser.set_defaults(solver=None)
    parser.set_defaults(symbols=[])

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)

    return options


#===============================================================================
# Main class
#===============================================================================

class cs_debug_symbol_translator:
    """
    Class with utility functions to translate debug symbols into
    file:line form using addr2line.
    """

    #---------------------------------------------------------------------------

    def __init__(self, pkg, path=None, solver=None):
        """
        Constructor
        """

        self.pkg    = pkg
        self.path   = path
        self.solver = solver
        if self.solver is None:
            self.solver = pkg.solver

    #---------------------------------------------------------------------------

    def split_symbol(self, symbol):
        """
        Split stack symbol into name and offset
        """

        name, offset = symbol.split("+")

        return name, offset

    #---------------------------------------------------------------------------

    def get_runtime_adresse(self, symbol_name):
        """
        Get runtime address of a symbol
        """

        solver_path = "."
        if self.path:
            solver_path = self.path

        solver_path = "/".join([solver_path, self.solver])


        cmd = "nm -D %s | grep %s" % (solver_path, symbol_name)

        out = get_command_output(cmd)

        if out == "":
            sys.exit(1)

        run_time_addr = "0x"+out.split(" ")[0]


        return run_time_addr

    #---------------------------------------------------------------------------

    def get_line_from_address(self, addr, offset):
        """
        Get file and line using the runtime addresse and offset of symbol
        """

        solver_path = "."
        if self.path:
            solver_path = self.path

        solver_path = "/".join([solver_path, self.solver])

        addr_int = int(addr, 0)
        off_int  = int(offset, 0)
        hex_addr = hex(addr_int + off_int)

        cmd = "addr2line -e %s %s" % (solver_path, hex_addr)

        out = get_command_output(cmd)

        if out is None or out == "":
            sys.exit(1)

        return out

    def symbol2line(self, sym_offset):

        sym, offset  = self.split_symbol(sym_offset)
        runtime_addr = self.get_runtime_adresse(sym)
        sym_line     = self.get_line_from_address(runtime_addr, offset)

        return sym_line


#===============================================================================
# Get symbols' lines
#===============================================================================

def symbols_to_lines(argv, pkg):
    """
    For each debug symbol, get the corresponding file and line.
    """

    if not pkg:
        from code_saturne.base.cs_package import package
        pkg = package()

    options = process_cmd_line(argv, pkg)

    translator = cs_debug_symbol_translator(pkg,
                                            path=options.path,
                                            solver=options.solver)

    for sym in options.symbols:
        file_line = translator.symbol2line(sym)

        if "??:?" in file_line:
            sys.stdout.write("Symbol '%s' yielded no file/line.\n" % sym)
            sys.stdout.write("Please rerun case with an instance compiled with --debug option.\n")

        else:
            res = "%s -> %s" % (sym, file_line)
            sys.stdout.write(res)

    return 0

#===============================================================================
# Main function
#===============================================================================

def main(argv, pkg):
    """
    Main function
    """
    return symbols_to_lines(argv, pkg)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    retval = main(argv=sys.argv[1:])

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
