#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
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
This module describes the script used to display information on Code_Saturne.

This module defines the following functions:
- process_cmd_line
- print_version
- launch_manual
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import os, sys, pwd, shutil, stat
from optparse import OptionParser

import cs_config


#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------


def process_cmd_line(argv):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string the
             string is split to create a list of arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-r", "--reader", dest="pdf_reader", type="string",
                      metavar="<pdfreader>",
                      help="define a pdf reader")

    parser.add_option("-g", "--guide", dest="guides", type="string",
                      metavar="<guide>", action="append",
                      help="open a manual [refcard, user, theory, tutorial]")

    parser.add_option("--version", dest="version",
                      action="store_true",
                      help="print Code_Saturne version number")

    parser.set_defaults(pdf_reader=None)
    parser.set_defaults(guides=[])
    parser.set_defaults(version=False)

    (options, args) = parser.parse_args(argv)

    if options.version:
        print_version()
        sys.exit(0)

    if len(args) > 0 or len(options.guides) == 0:
        parser.print_help()
        sys.exit(1)

    return options.guides, options.pdf_reader


#-------------------------------------------------------------------------------
# Print Code_Saturne version
#-------------------------------------------------------------------------------


def print_version():
    """
    Print Code_Saturne version.
    """

    print("Code_Saturne version: %s" % cs_config.package.version)


#-------------------------------------------------------------------------------
# Launch the PDF manual
#-------------------------------------------------------------------------------


def launch_manual(reader, m):
    """
    Launch the corresponding PDF manual.
    """

    sys_tools = ["xdg-open",       # Generic
                 "gnome-open",     # Gnome
                 "kfmclient exec", # KDE
                 "kde-open",       # KDE 4
                 "exo-open"]       # Xfce

    readers = ["evince", "gpdf", "kpdf", "xpdf", "acroread"]

    manual = os.path.join(cs_config.dirs.pdfdir, m) + '.pdf'

    if not os.path.isfile(manual):
        print("File %s not found." % manual)
        return

    # First:  use the reader specified by the user, if any
    # Second: try the different tool that open with the default system reader
    # Last:   try some classical pdf viewers

    if reader is not None:
        cmd = reader + ' ' + manual + ' 2>/dev/null &'
        os.system(cmd)

    else:
        if os.name == "nt":
            os.filestart(manual)

        elif os.name == "posix":

            for r in (sys_tools + readers):
                cmd = r + ' ' + manual + ' 2>/dev/null &'
                if os.system(cmd) == 0:
                    return


#-------------------------------------------------------------------------------
# Creation of the study directory
#-------------------------------------------------------------------------------


def main(argv):
    """
    Main function.
    """

    manuals, reader = process_cmd_line(argv)

    for m in manuals:
        launch_manual(reader, m)


if __name__ == "__main__":
    main(sys.argv[1:])

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
