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
This module describes the script used to display information on Code_Saturne.

This module defines the following functions:
- process_cmd_line
- print_version
- launch_manual
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, shutil, stat, fnmatch
from optparse import OptionParser

#-------------------------------------------------------------------------------
# Licence text
#-------------------------------------------------------------------------------

licence_text = \
"""Copyright (C) 1998-2020 EDF S.A.

This program is free software; you can redistribute it and/or modify it under \
the terms of the GNU General Public License as published by the Free Software \
Foundation; either version 2 of the License, or (at your option) any later \
version.

This program is distributed in the hope that it will be useful, but WITHOUT \
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS \
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more \
details.

You should have received a copy of the GNU General Public License along with \
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin \
Street, Fifth Floor, Boston, MA 02110-1301, USA."""

#-------------------------------------------------------------------------------
# Processes the passed command line arguments
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

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog info [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--reader", dest="pdf_reader", type="string",
                      metavar="<pdfreader>",
                      help="define a pdf reader")

    parser.add_option("-g", "--guide", dest="guides", type="string",
                      metavar="<guide>", action="append",
                      help="open a manual " + str(get_docs(pkg)))

    parser.add_option("--modules", dest="modules",
                      action="store_true",
                      help="print evironment modules")

    parser.add_option("--version", dest="version",
                      action="store_true",
                      help="print version number")

    parser.set_defaults(pdf_reader=None)
    parser.set_defaults(guides=[])
    parser.set_defaults(version=False)

    (options, args) = parser.parse_args(argv)

    if options.version:
        print_version(pkg)
        sys.exit(0)

    if options.modules:
        print_modules(pkg)
        sys.exit(0)

    if len(args) > 0 or len(options.guides) == 0:
        parser.print_help()
        sys.exit(1)

    return options.guides, options.pdf_reader


#-------------------------------------------------------------------------------
# Print Code_Saturne version
#-------------------------------------------------------------------------------

def print_version(pkg):
    """
    Print Code_Saturne version.
    """

    print(pkg.code_name + " version: " + pkg.version_full)

#-------------------------------------------------------------------------------
# Print Environment modules info
#-------------------------------------------------------------------------------

def print_modules(pkg):
    """
    Print Code_Saturne environment modules info.
    """

    from code_saturne import cs_config
    c = cs_config.config()

    if c.env_modulecmd:
        print("Module command:      " + str(c.env_modulecmd))
        print("Environment modules: " + str(c.env_modules))

#-------------------------------------------------------------------------------
# Launch the PDF manual
#-------------------------------------------------------------------------------

def get_docs(pkg):
    """
    Return the list of available PDF and Doxygen manuals for the command line.
    """
    l = []
    if os.path.isdir(pkg.get_dir('docdir')):
        for docs in fnmatch.filter(os.listdir(pkg.get_dir('docdir')), '*.pdf'):
            l.append(docs[:-4])
        doxy_dir = os.path.join(pkg.get_dir('docdir'),'doxygen', 'src')
        if os.path.isdir(doxy_dir):
            for docs in fnmatch.filter(os.listdir(doxy_dir), 'index.html'):
                l.append('Doxygen')
    return l

def launch_manual(reader, m, pkg):
    """
    Launch the corresponding PDF manual.
    """

    sys_tools = ["xdg-open",       # Generic
                 "gnome-open",     # Gnome
                 "kfmclient exec", # KDE
                 "kde-open",       # KDE 4
                 "exo-open"]       # Xfce

    if not m == "Doxygen":
        readers = ["okular", "evince", "kpdf", "gpdf", "xpdf", "acroread"]
        manual = os.path.join(pkg.get_dir('docdir'), m) + '.pdf'
    else:
        readers = ["firefox"]
        manual = os.path.join(pkg.get_dir('docdir'),'doxygen', 'src', 'index.html')

    if not os.path.isfile(manual):
        print("File %s not found." % manual)
        return

    # On Windows platform, use the system standard launcher
    if sys.platform.startswith('win'):
        os.startfile(manual)

    # On Linux and Unix-like platforms,
    #  - first:  use the reader specified by the user, if any
    #  - second: try the different tool that open with the default system reader
    #  - last:   try some classical pdf viewers
    else:
        if reader is not None:
            cmd = reader + ' ' + manual + ' 2>/dev/null &'
            os.system(cmd)
        else:
            for r in (sys_tools + readers):
                cmd = r + ' ' + manual + ' 2>/dev/null &'
                if os.system(cmd) == 0:
                    return

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    manuals, reader = process_cmd_line(argv, pkg)

    for m in manuals:
        launch_manual(reader, m, pkg)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
