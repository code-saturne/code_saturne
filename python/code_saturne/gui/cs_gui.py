# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
Parse command line arguments and launch the GUI.
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string

from optparse import OptionParser

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

# Force PyQt API 2 for PyQt4

try:
    import sip
    for api in ('QDate', 'QDateTime', 'QString', 'QTextStream',
                'QTime', 'QUrl', 'QVariant'):
        sip.setapi(api, 2)
except Exception:
    pass

try:
    from code_saturne.gui.base.QtCore    import *
    from code_saturne.gui.base.QtGui     import *
    from code_saturne.gui.base.QtWidgets import *
except ImportError:
    print("\n  Error: Unable to import QtCore or QtGui modules.")
    print("  Please check your PyQt4 or PyQt5 installation.\n")
    sys.exit(0)


if list(map(int, QT_VERSION_STR.split( "."))) < [4, 3, 0]:
    raise SystemExit("Graphical user interface requires Qt 4.3 or later "\
                     "(found %s)." % QT_VERSION_STR)


if list(map(int, PYQT_VERSION_STR.split("."))) < [4, 5, 0]:
    raise SystemExit("Graphical user interface requires PyQt 4.5 or later "\
                     "(found %s)." % PYQT_VERSION_STR)

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
        usage = "usage: %prog gui [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--param", dest="file_name", type="string",
                      metavar="<file>",
                      help="upload a previous case at the interface start")

    parser.add_option("-n", "--new", dest="new",
                      action="store_true",
                      help="open a new case")

    parser.add_option("-z", "--no-splash", dest="splash_screen",
                      action="store_false",
                      help="deactivate splash screen")

    parser.set_defaults(splash_screen=True, enable_neptune_cfd=True)

    (options, args) = parser.parse_args(argv)

    if options.new and options.file_name:
        parser.error("Options --new and --param are mutually exclusive")

    if options.new:
        options.file_name = "new case"

    if len(args) > 0:
        if options.file_name or len(args) > 1:
            parser.error("Multiple filenames are given")
        else:
            options.file_name = args[0]

    # If no parameter file passed and a setup.xml is present, open it,
    # otherwise create a new case.
    has_setup = os.path.isfile(os.path.join(os.getcwd(), 'setup.xml'))
    if not options.file_name:
        if has_setup:
            options.file_name = "setup.xml"
        else:
            options.file_name = "new case"

    return options.file_name, options.splash_screen

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Start Qt and a session of the application.
    """

    from code_saturne.base.cs_exec_environment import set_modules, source_rcfile
    set_modules(pkg)
    source_rcfile(pkg)

    images_path = os.path.join(pkg.get_dir('pkgdatadir'), 'images')

    # Test if EOS modules could be imported
    cfg = pkg.config
    if cfg.libs['eos'].have:
        eosprefix = cfg.libs['eos'].prefix
        try:
            from distutils import sysconfig
            eospath = os.path.join(sysconfig.get_python_lib(0, 0, prefix=eosprefix), 'eos')
        except Exception:
            eospath = ''

        if sys.platform.startswith('win'):
            python_version = "%d.%d" % (sys.version_info.major, sys.version_info.minor)
            eospath = os.path.join(eosprefix,
                          'lib', 'python' + python_version, 'site-packages',
                          'eos')

        if eospath:
            if os.path.isdir(eospath) and not eospath in sys.path:
                sys.path.insert(0, eospath)

    case, spl = process_cmd_line(argv)

    app = QApplication(sys.argv)
    app.setOrganizationName(pkg.code_name) # Defines the name of subdirectory under .config
    app.setOrganizationDomain(pkg.url)
    app.setApplicationName(pkg.name)
    #app.setApplicationName("code_saturne") # Defines the name of the configuration file
    #app.setWindowIcon(QIcon(":/icon.png"))
    app.lastWindowClosed.connect(app.quit)

    # Locale detection
    locale = QLocale.system().name()
    translator = QTranslator(app)

    tr_file = "code_saturne"
    localedir = pkg.get_dir('localedir')
    tr_loaded = translator.load(QLocale(), tr_file, "_", localedir, ".qm")
    if tr_loaded:
        app.installTranslator(translator)

    if spl:
        app.setOverrideCursor(QCursor(Qt.WaitCursor))
        # Choose correct splahs screen based on solver at runtime
        if pkg.name == 'code_saturne':
            pixmap = QPixmap('%s/splashscreen.png' % images_path)
        else:
            pixmap = QPixmap('%s/logo_salome_cfd.png' % images_path)

        splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
        splash.setMask(pixmap.mask()) # this is useful if the splashscreen is not a regular rectangle...
        splash.show()
        app.processEvents()
        QTimer.singleShot(1500, splash.hide)

    if os.path.split(case)[1] == "run.cfg":
        from code_saturne.gui.base.QCouplingEditorView import QCouplingEditor
        mv = QCouplingEditor(cfgfile=case, standalone_mode=True)
    else:
        from code_saturne.gui.base.MainView import MainView
        mv = MainView(cmd_package = pkg, cmd_case = case)

    try:
        mv.show()
        if spl:
            app.processEvents()
            app.restoreOverrideCursor()
    except:
        print("\n  Unable to display a Qt window.")
        print("  Please check your display environment.\n")
        sys.exit(0)

    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
