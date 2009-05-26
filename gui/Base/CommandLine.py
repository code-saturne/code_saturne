# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module describes what to do when the line command is parsed.

This module defines the following functions:
- usage
- process_cmd_line
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types, string, getopt

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Usage/help message
#-------------------------------------------------------------------------------

def usage():
    """
    Usage of the GUI.
    """
    txt="""
Usage:\n\ncs_gui [options]... [xml file name]

Valid options are one or more of the following:

-f [xml file name]
--file [xml file name]

     Upload a previous case at the interface start.

-n
--new

     Open a new case. None effect if '-f' option is used.

-v
--version

     Prints the GUI version.

-h
--help

     Prints this help message.

-z
--no-splash

     Without splash screen.

-m
--matisse

    Load matisse version.

-b [batchfile]
--batch (batchfile]

    Set batchrunning window with batch file
    (-f or --file option is mandatory).

-t
--no-tree

    No tree window loaded.


Examples:

     cs_gui -f name_of_file.xml

     cs_gui --file name_of_file.xml

     cs_gui -n

----------------------------------------------------------------------
"""
    return txt


#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------


def process_cmd_line (arg):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string the
             string is split to create a list of arguments.
    Returned Values:
      case -- Name of the file. If the variable case is set 
              to "new case" GUI is open with a new case.
    """

    # The "-n" option has none effect when "-f" is used
    #
    if ( ('-n' in arg) or ('--new' in arg) ) and \
       ( ('-f' in arg) or ('--file' in arg) ):
        try:
            del arg[arg.index('-n')]
        except:
            del arg[arg.index('--new')]

    # When xml file name is the only one argument
    #
    if type(arg) is types.StringType:
        arg = "-f " + arg
        arg = string.split(arg)

    options = "f:b:mntrz"
    long_opts = ['file=', 'batch=', 'matisse', 'new', 'no-tree', 'read-only', 'no-splash']

    try:
        opts, args = getopt.getopt(arg, options, long_opts)
    except getopt.error, msg:
        print "\nWarning:\n", msg
        print usage()
        print "Warning:\n", msg, "\n"
        sys.exit(1)

    if (opts and args) or len(args) > 1: 
        print usage()
        sys.exit(1)

    case = ""
    matisse = False
    batch_window = False
    batch_file = 'lance'
    tree_window = True
    read_only = False
    splash = True

    if not opts and args :
        case = args[0]

    for o, a in opts:

        if o in ('-f', '--file'):
            case = a

        if o in ('-n', '--new'):
            case = "new case"

        if o in ('-m', '--matisse'):
            matisse = True

        if o in ('-b', '--batch'):
            batch_window = True
            batch_file = os.path.basename(a)
            if not case:
                print usage()
                sys.exit(1)

        if o in ('-z', '--no-splash'):
            splash = False

        if o in ('-t', '--no-tree'):
            tree_window = False

        if o in ('-r', '--read-only'):
            read_only = True


    return case, splash, matisse, batch_window, batch_file, tree_window, read_only


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
