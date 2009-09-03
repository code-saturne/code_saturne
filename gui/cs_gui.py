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
Parse command line arguments and launch the GUI.
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string

if not hasattr(sys, 'version_info') or sys.version_info <= (2, 4, 0, 'final'):
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires python 2.4 or later."

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui  import *
except ImportError:
    print "\n  Error: Unable to import PyQt4.QtCore or PyQt4.QtGui modules."
    print "  Please check your PyQt4 installation.\n"
    sys.exit(0)


if map(int, string.split(QT_VERSION_STR, ".")) < [4, 3, 0]:
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires Qt 4.3 or later (found %s)." % QT_VERSION_STR


if map(int, string.split(PYQT_VERSION_STR, ".")) < [4, 3, 0]:
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires PyQt 4.3 or later (found %s)." % PYQT_VERSION_STR

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

try:
    import ncs
except:
    pass

from Base.Common import icon_base_path
from Base.CommandLine import usage, process_cmd_line
from Base.MainView import MainView

import cs_config

#-------------------------------------------------------------------------------
# Help messages
#-------------------------------------------------------------------------------

if ('-h' in sys.argv[1:]) or ('--help' in sys.argv[1:]):
    print usage()
    sys.exit(0)


if ('-v' in sys.argv[1:]) or ('--version' in sys.argv[1:]):
    print "Graphical users interface of Code_Saturne %s" % cs_config.package.version
    sys.exit(0)

#-------------------------------------------------------------------------------
# Start point of the Graphical User Interface
#-------------------------------------------------------------------------------

def main(argv):
    """
    Start Qt and a session of the application.
    """
    case, spl, matisse, batch_window, batch_file, tree_window, read_only \
       = process_cmd_line(argv)

    app = QApplication(argv)
    app.setOrganizationName("EDF R&D")
    app.setOrganizationDomain("www.code_saturne.org")
    app.setApplicationName("Code_Saturne GUI")
    #app.setWindowIcon(QIcon(":/icon.png"))
    app.connect(app, SIGNAL("lastWindowClosed()"), app, SLOT("quit()"))

    if spl:
        app.setOverrideCursor(QCursor(Qt.WaitCursor))
        pixmap = QPixmap('%s/SplashScreen/logocs.png' % icon_base_path)
        splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
        splash.setMask(pixmap.mask()) # this is usefull if the splashscreen is not a regular ractangle...
        splash.show()
        splash.showMessage('GUI %s starting...' % cs_config.package.version,
                           Qt.AlignHCenter | Qt.AlignVCenter, Qt.black)
        app.processEvents()
        QTimer.singleShot(1500, splash.hide)

    main = MainView(cmd_case = case,
                    cmd_matisse = matisse,
                    cmd_batch_window = batch_window,
                    cmd_batch_file = batch_file,
                    cmd_tree_window = tree_window,
                    cmd_read_only = read_only)

    try:
        main.show()
        if spl:
            app.processEvents()
            app.restoreOverrideCursor()
    except:
        print "\n  Unable to display a Qt window."
        print "  Please check your display environment.\n"
        sys.exit(0)

    sys.exit(app.exec_())

if __name__ == '__main__':
    main(sys.argv[1:])

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
