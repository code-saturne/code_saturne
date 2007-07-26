#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2007 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne Kernel is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne Kernel is distributed in the hope that it will be
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
#============================================================================
#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, getopt

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------
import Autovalidation.Common as Common


def usage():
    """
    Usage of autovalid
    """
    txt="""
Usage:\n\nautovalid -f [xml file name] [-d [local directory]]

Valid options are one or more of the following:

-f [xml file name]
--file [xml file name]

     Upload a studies list.

-d [directory]
--tmpdirectory [directory]

     Define the local directory.

-h
--help

     Prints this help message.
----------------------------------------------------------------------
"""
    return txt

def process_cmd_line (arg):
    
    options = "f:hd:"
    long_opts = ['file=', 'help', 'tmpdirectory=']

    try:
        opts, args = getopt.getopt(arg, options, long_opts)
    except getopt.error, msg:
        print "\nWarning:\n", msg
        print usage()
        print "Warning:\n", msg, "\n"
        sys.exit(1)

    if (opts and args) or len(args)>1 : 
        print usage()
        sys.exit(1)

    if not opts and args :
        print usage()
        sys.exit(1)

    for o, a in opts:
        if o in ('-f', '--file'):
            Common.XMLFileName = a
        if o in ('-d', '--tmpdirectory'):
            Common.tmpDirectory = a
        if o in ('-h', '--help'):
            print usage()
            sys.exit(1)

    if Common.XMLFileName == "":
        print usage()
        sys.exit(1)
