# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
Parse command line arguments and launch tool for tracking convergence
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string

from optparse import OptionParser

if not hasattr(sys, 'version_info') or sys.version_info <= (2, 4, 0, 'final'):
    raise SystemExit("Graphical user interface requires python 2.4 or later.")

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

# to use PyQt API 2
#if sys.version_info[0] == 2:
#    import sip
#    sip.setapi('QString', 2)
#    sip.setapi('QVariant', 2)
#    sip.setapi('QTime', 2)
#    sip.setapi('QUrl', 2)

try:
    from code_saturne.Base.QtCore    import *
    from code_saturne.Base.QtGui     import *
    from code_saturne.Base.QtWidgets import *
except ImportError:
    print("\n  Error: Unable to import QtCore or QtGui modules.")
    print("  Please check your PyQt4 or PyQt5 installation.\n")
    sys.exit(0)


if list(map(int, QT_VERSION_STR.split( "."))) < [4, 3, 0]:
    raise SystemExit("Graphical user interface requires Qt 4.3 or later "\
                     "(found %s)." % QT_VERSION_STR)


if list(map(int, PYQT_VERSION_STR.split("."))) < [4, 3, 0]:
    raise SystemExit("Graphical user interface requires PyQt 4.3 or later "\
                     "(found %s)." % PYQT_VERSION_STR)

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import cs_config

#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog trackcvg [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--resu", dest="file_name", type="string",
                      metavar="<directory>",
                      help="open a result directory at the interface start")

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        if options.file_name or len(args) > 1:
            parser.error("Multiple filenames are given")
        else:
            options.file_name = args[0]

    return options.file_name

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Start Qt and a session of the application.
    """

    from cs_exec_environment import set_modules, source_rcfile
    set_modules(pkg)
    source_rcfile(pkg)

    # Test the package name to know which modules have to be imported
    if pkg.name == 'code_saturne':
        images_path = os.path.join(pkg.get_dir('pkgdatadir'), 'images')
        sys.path.insert(1, os.path.join(pkg.get_dir('pkgpythondir'), 'Base'))
    else:
        images_path = os.path.join(pkg.get_dir('pkgpythondir'), 'core', 'icons')
        sys.path.insert(1, os.path.join(pkg.get_dir('pkgpythondir'), 'core'))
        sys.path.insert(1, pkg.get_dir('pythondir'))

    case = process_cmd_line(argv)

    app = QApplication(argv)
    app.setOrganizationName(pkg.code_name) # Defines the name of subdirectory under .config
    app.setOrganizationDomain(pkg.url)
    app.setApplicationName("tracking convergence") # Defines the name of the configuration file
    app.lastWindowClosed.connect(app.quit)

    # Locale detection
    locale = QLocale.system().name()
    translator = QTranslator(app)
    if translator.load("qt_" + locale,
                       QLibraryInfo.location(QLibraryInfo.TranslationsPath)):
        app.installTranslator(translator)

    from code_saturne.trackcvg.MainView import MainView
    mv = MainView(cmd_package = pkg, cmd_case = case)

    try:
        mv.show()
    except:
        print("\n  Unable to display a Qt window.")
        print("  Please check your display environment.\n")
        sys.exit(0)

    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
