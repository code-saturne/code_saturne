# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the main application classes for the Qt GUI.
This GUI provides a simple way to display independante pages, in order to put
informations in the XML document, which reflets the treated case.

This module defines the following classes:
- MainView

    @copyright: 1998-2017 EDF S.A., France
    @author: U{EDF<mailto:saturne-support@edf.fr>}
    @license: GNU GPL v2, see COPYING for details.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, shutil, signal, logging
import subprocess, platform

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import cs_info
from cs_exec_environment import \
    separate_args, update_command_single_value, assemble_args, enquote_arg
import cs_runcase

try:
    from code_saturne.studymanager_gui.MainForm import Ui_MainForm
except:
    sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Base"))
    from code_saturne.studymanager_gui.MainForm import Ui_MainForm

from code_saturne.Base import XMLengine
from code_saturne.Base.XMLmodel import *
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.Common import XML_DOC_VERSION
from code_saturne.studymanager_gui.Toolbox import displaySelectedPage
from code_saturne.studymanager_gui.BrowserView import BrowserView
from code_saturne.studymanager_gui.XMLinitialize import *

try:
    import code_saturne.Pages
except:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from code_saturne.Pages.WelcomeView import WelcomeView
from code_saturne.Pages.IdentityAndPathesModel import IdentityAndPathesModel
from code_saturne.Pages.XMLEditorView import XMLEditorView
from code_saturne.Base.QtPage import getexistingdirectory
from code_saturne.Base.QtPage import from_qvariant, to_text_string, getopenfilename, getsavefilename


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MainView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Base Main Window
#-------------------------------------------------------------------------------

class MainView(object):
    """
    Abstract class
    """
    NextId = 1
    Instances = set()

    def __new__(cls, cmd_package = None, cmd_case = ""):
        """
        Factory
        """
        return MainViewSaturne.__new__(MainViewSaturne, cmd_package, cmd_case)


    @staticmethod
    def updateInstances(qobj):
        """
        Overwrites the Instances set with a set that contains only those
        window instances that are still alive.
        """
        MainView.Instances = set([window for window \
                in MainView.Instances if isAlive(window)])


    def ui_initialize(self):
        self.setAttribute(Qt.WA_DeleteOnClose)
        MainView.Instances.add(self)

        iconpath = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
        iconpath = os.path.join(iconpath, "Base", "MONO-bulle-HD.png")
        icon = QIcon(QPixmap(iconpath))
        self.setWindowIcon(icon)

        self.setWindowTitle(self.package.code_name + " STUDYMANAGER GUI" + " - " + self.package.version)

        self.dockWidgetBrowser.setWidget(self.Browser)

        self.scrollArea = QScrollArea(self.frame)
        self.gridlayout1.addWidget(self.scrollArea,0,0,1,1)
        self.gridlayout1.setSpacing(0)
        self.gridlayout.addWidget(self.frame,0,0,1,1)

        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setFrameShape(QFrame.StyledPanel)
        self.scrollArea.setFrameShadow(QFrame.Raised)
        self.scrollArea.setFrameStyle(QFrame.NoFrame)

        # connections

        self.fileOpenAction.triggered.connect(self.fileOpen)
        self.fileNewAction.triggered.connect(self.fileNew)
        self.menuRecent.aboutToShow.connect(self.updateRecentFileMenu)
        self.fileSaveAction.triggered.connect(self.fileSave)
        self.fileSaveAsAction.triggered.connect(self.fileSaveAs)
        self.fileCloseAction.triggered.connect(self.close)
        self.fileQuitAction.triggered.connect(self.fileQuit)

        self.BrowserAction.toggled.connect(self.dockWidgetBrowserDisplay)

        self.displayAboutAction.triggered.connect(self.displayAbout)
        self.backgroundColorAction.triggered.connect(self.setColor)
        self.actionFont.triggered.connect(self.setFontSize)
        self.RestoreStyleDefaults.triggered.connect(self.restoreStyleDefaults)

        self.displayLicenceAction.triggered.connect(self.displayLicence)

        # connection for page layout

        self.Browser.treeView.pressed.connect(self.displayNewPage)
        self.destroyed.connect(MainView.updateInstances)

        # Ctrl+C signal handler (allow to shutdown the GUI with Ctrl+C)

        signal.signal(signal.SIGINT, signal.SIG_DFL)

        self.resize(800, 700)

        # restore system settings

        settings = QSettings()

        try:
            # API 2
            if settings.value("RecentFiles", []) is not None:
                try:
                    recentFiles = settings.value("RecentFiles").toStringList()
                    self.recentFiles = []
                    for f in recentFiles:
                        self.recentFiles.append(str(f))
                except:
                    self.recentFiles = list(settings.value("RecentFiles", []))
            else:
                self.recentFiles = []
            self.restoreGeometry(settings.value("MainWindow/Geometry", QByteArray()))
            self.restoreState(settings.value("MainWindow/State", QByteArray()))
        except:
            # API 1
            self.recentFiles = settings.value("RecentFiles").toStringList()
            self.restoreGeometry(settings.value("MainWindow/Geometry").toByteArray())
            self.restoreState(settings.value("MainWindow/State").toByteArray())

        app = QCoreApplication.instance()

        self.palette_default = None
        self.font_default = None

        if settings.contains("MainWindow/Color"):
            color = settings.value("MainWindow/Color",
                                   self.palette().color(QPalette.Window).name())
            color = QColor(color)
            if color.isValid():
                if not self.palette_default:
                    self.palette_default = QPalette().resolve(self.palette())
                self.setPalette(QPalette(color))
                app.setPalette(QPalette(color))

        if settings.contains("MainWindow/Font"):
            f = settings.value("MainWindow/Font", str(self.font()))
            if f:
                if not self.font_default:
                    self.font_default = self.font()
                font = QFont()
                if (font.fromString(from_qvariant(f, to_text_string))):
                    self.setFont(font)
                    app.setFont(font)

        self.updateRecentFileMenu()
        QTimer.singleShot(0, self.loadInitialFile)

        self.statusbar.setSizeGripEnabled(False)
        self.statusbar.showMessage(self.tr("Ready"), 5000)


    def loadInitialFile(self):
        """
        Private method.

        Checks the opening mode (from command line).
        """
        log.debug("loadInitialFile -> %s" % self.cmd_case)

        # 1) new case
        if self.cmd_case == "new case":
            MainView.NextId += 1
            self.fileNew()
            self.dockWidgetBrowserDisplay(True)

        # 2) existing case

        elif self.cmd_case:
            try:
                self.loadFile(self.cmd_case)
                self.dockWidgetBrowserDisplay(True)
            except:
                raise

        # 3) neutral point (default page layout)

        else:
            self.dockWidgetBrowserDisplay(False)


    def dockWidgetBrowserDisplay(self, bool=True):
        """
        Private slot.

        Show or hide the browser dock window.

        @type bool: C{True} or C{False}
        @param bool: if C{True}, shows the browser dock window
        """
        if bool:
            self.dockWidgetBrowser.show()
        else:
            self.dockWidgetBrowser.hide()


    def updateRecentFileMenu(self):
        """
        private method

        update the File menu with the recent files list
        """
        self.menuRecent.clear()

        if hasattr(self, 'case'):
            current = str(self.case['xmlfile'])
        else:
            current = ""

        recentFiles = []
        for f in self.recentFiles:
            if f != current and QFile.exists(f):
                recentFiles.append(f)

        if recentFiles:
            for i, f in enumerate(recentFiles):
                action = QAction(QIcon(":/icons/22x22/document-open.png"), "&%d %s" % (
                           i + 1, QFileInfo(f).fileName()), self)
                action.setData(f)
                action.triggered.connect(self.loadRecentFile)
                self.menuRecent.addAction(action)


    def addRecentFile(self, fname):
        """
        private method

        creates and update the list of the recent files

        @type fname: C{str}
        @param fname: filename to add in the recent files list
        """
        if fname is None:
            return
        if not fname in self.recentFiles:
            self.recentFiles.insert(0, str(fname))
            # l'idée est de garder les 8 premiers élements de la
            # liste. On pourrait donc remplacer le code ci-dessus par
            # :
            # self.recentFiles = self.recentFiles[:8]
            while len(self.recentFiles) > 9:
                try:
                    self.recentFiles.removeLast()
                except:
                    try:
                        self.recentFiles.pop()
                    except:
                        self.recentFiles.takeLast()


    def closeEvent(self, event):
        """
        public slot

        try to quit all the current MainWindow
        """
        if self.okToContinue():
            settings = QSettings()
            if self.recentFiles:
                recentFiles = self.recentFiles
            else:
                recentFiles = []

            settings.setValue("RecentFiles", recentFiles)
            settings.setValue("MainWindow/Geometry",
                              self.saveGeometry())
            settings.setValue("MainWindow/State",
                              self.saveState())

            event.accept()
            log.debug("closeEvent -> accept")

        else:
            event.ignore()
            log.debug("closeEvent -> ignore")


    def okToContinue(self):
        """
        private method

        ask for unsaved changes before quit

        @return: C{True} or C{False}
        """
        title = self.tr("Quit")
        msg   = self.tr("Save unsaved changes?")

        if hasattr(self, 'case'):
            log.debug("okToContinue -> %s" % self.case.isModified())

        if hasattr(self, 'case') and self.case.isModified():
            reply = QMessageBox.question(self,
                                         title,
                                         msg,
                                         QMessageBox.Yes|
                                         QMessageBox.No|
                                         QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()

        return True


    def fileQuit(self):
        """
        Public slot.

        try to quit all window
        """
        QApplication.closeAllWindows()


    def fileNew(self):
        """
        Public slot.

        create new Code_Saturne studymanager parameter file
        """
        if not hasattr(self, 'case'):
            self.case = XMLengine.Case(package=self.package, studymanager=True)
            self.case.root()['version'] = self.XML_DOC_VERSION
            self.initCase()
            title = self.tr("New parameters set") + \
                     " - " + self.tr(self.package.code_name) + self.tr(" STUDYMANAGER GUI") \
                     + " - " + self.package.version
            self.setWindowTitle(title)

            self.Browser.configureTree(self.case)
            self.dockWidgetBrowserDisplay(True)

            self.case['saved'] = "yes"
        else:
            MainView(cmd_package=self.package, cmd_case="new case").show()
        # TODO
        # faire le detect du nom study et charger tous les cases du repertoire par defaut
        # peut etre chrge les nouveau cas lorsqu'on re-ouvre avec status a off

    def fileAlreadyLoaded(self, f):
        """
        private method

        check if the file to load is not already loaded

        @type fname: C{str}
        @param fname: file name to load
        @return: C{True} or C{False}
        """
        for win in MainView.Instances:
            if isAlive(win) and hasattr(win, 'case') \
               and win.case['xmlfile'] == f:
                win.activateWindow()
                win.raise_()
                return True
        return False


    def loadRecentFile(self, file_name=None):
        """
        private slot

        reload an existing recent file

        @type fname: C{str}
        @param fname: file name to load
        """
        # reload  from File menu
        if file_name is None:
            action = self.sender()
            if isinstance(action, QAction):
                file_name = unicode(action.data().toString())
                if not self.okToContinue():
                    return
            else:
                return

        # check if the file to load is not already loaded
        if hasattr(self, 'case'):
            if not self.fileAlreadyLoaded(file_name):
                MainView(cmd_package=self.package, cmd_case = file_name).show()
        else:
            self.loadFile(file_name)


    def loadingAborted(self, msg, fn):
        """Show a message window dialog.
        Delete the case if it is already loaded, but non conformal.

        @param msg text to display in the popup window
        @param fn name of the file pf parameters
        """
        msg += self.tr("\n\nThe loading of %s is aborted." % fn)
        title = self.tr("File of parameters reading error")

        QMessageBox.critical(self, title, msg)

        if hasattr(self, 'case'):
            delattr(self, 'case')


    def loadFile(self, file_name=None):
        """
        Private method

        Load an existing file.

        @type fname: C{str}
        @param fname: file name to load
        """
        file_name = os.path.abspath(str(file_name))
        fn = os.path.basename(file_name)
        log.debug("loadFile -> %s" % file_name)

        # XML syntax checker

        msg =  XMLengine.xmlChecker(file_name)
        if msg:
            self.loadingAborted(msg, fn)
            return

        # Instantiate a new case

        try:
            self.case = XMLengine.Case(package=self.package, file_name=file_name, studymanager=True)
        except:
            msg = self.tr("This file is not in accordance with XML specifications.")
            self.loadingAborted(msg, fn)
            return

        # Cleaning the '\n' and '\t' from file_name (except in formula)
        self.case.xmlCleanAllBlank(self.case.xmlRootNode())

        # we consider we are in calculation mode when we open an xml file
        msg = self.initCase()
        if msg:
            self.loadingAborted(msg, fn)
            return

        # All checks are fine, wan can continue...

        self.addRecentFile(fn)
        self.Browser.configureTree(self.case)
        self.dockWidgetBrowserDisplay(True)

        # Update the case and the StudyIdBar
        self.case['xmlfile'] = file_name
        title = fn + " - " + self.tr(self.package.code_name) + self.tr(" STUDYMANAGER GUI") \
                   + " - " + self.package.version
        self.setWindowTitle(title)

        msg = self.tr("Loaded: %s" % fn)
        self.statusbar.showMessage(msg, 2000)

        self.case['saved'] = "yes"


    def fileOpen(self):
        """
        public slot

        open an existing file
        """
        msg = self.tr("Opening an existing case.")
        self.statusbar.showMessage(msg, 2000)

        title = self.tr("Open existing file.")

        if hasattr(self, 'case') and os.path.isdir(self.case['data_path']):
            path = self.case['data_path']
        else:
            path = os.getcwd()
            dataPath = os.path.join(path, "..", "DATA")
            if os.path.isdir(dataPath): path = dataPath

        filetypes = self.tr(self.package.code_name) + self.tr(" GUI files (*.xml);;""All Files (*)")

        file_name, _selfilter = getopenfilename(self, title, path, filetypes)

        if not file_name:
            msg = self.tr("Loading aborted")
            self.statusbar.showMessage(msg, 2000)
            file_name = None
            return
        else:
            file_name = str(file_name)
            log.debug("fileOpen -> %s" % file_name)

        if hasattr(self, 'case'):
            if not self.fileAlreadyLoaded(file_name):
                MainView(cmd_package=self.package, cmd_case = file_name).show()
        else:
            self.loadFile(file_name)

        self.statusbar.clearMessage()


    def updateStudyId(self):
        """
        private method

        update the XML file name
        """
        file_name = XMLengine._encode(self.case['xmlfile'])


    def fileSave(self):
        """
        public slot

        save the current case
        """
        log.debug("fileSave()")

        if not hasattr(self, 'case'):
            return

        file_name = self.case['xmlfile']
        log.debug("fileSave(): %s" % file_name)
        if not file_name:
            self.fileSaveAs()
            return

        log.debug("fileSave(): %s" % os.path.dirname(file_name))
        log.debug("fileSave(): %s" % os.access(os.path.dirname(file_name), os.W_OK))
        if not os.access(os.path.dirname(file_name), os.W_OK):
            title = self.tr("Save error")
            msg   = self.tr("Failed to write %s " % file_name)
            QMessageBox.critical(self, title, msg)
            msg = self.tr("Saving aborted")
            self.statusbar.showMessage(msg, 2000)
            return

        self.updateStudyId()
        self.case.xmlSaveDocument()

        log.debug("fileSave(): ok")

        msg = self.tr("%s saved" % file_name)
        self.statusbar.showMessage(msg, 2000)


    def fileSaveAs(self):
        """
        public slot

        save the current case with a new name
        """
        log.debug("fileSaveAs()")

        if hasattr(self,'case'):
            filetypes = self.tr(self.package.code_name) + self.tr(" GUI files (*.xml);;""All Files (*)")
            fname, _selfilter = getsavefilename(self,
                                                self.tr("Save File As"),
                                                self.case['data_path'],
                                                filetypes)

            if fname:
                f = str(fname)
                self.case['xmlfile'] = f
                self.addRecentFile(f)
                self.fileSave()
                self.updateStudyId()
                self.case.xmlSaveDocument()
                title = os.path.basename(self.case['xmlfile']) + " - " + self.tr(self.package.code_name) + self.tr(" GUI") \
                     + " - " + self.package.version
                self.setWindowTitle(title)
            else:
                msg = self.tr("Saving aborted")
                self.statusbar.showMessage(msg, 2000)


    def displayManual(self, pkg, manual, reader = None):
        """
        private method

        open a manual
        """
        argv_info = ['--guide']
        argv_info.append(manual)
        cs_info.main(argv_info, pkg)


    def displayNewPage(self, index):
        """
        private slot

        display a new page when the Browser send the order

        @type index: C{QModelIndex}
        @param index: index of the item in the C{QTreeView} clicked in the browser
        """
        # stop if the entry is a folder or a file

        if self.Browser.isFolder(): return

        # warning and stop if is no case
        if not hasattr(self, 'case'):
            log.debug("displayNewPage(): no attr. 'case', return ")

            msg = self.tr("You have to create a new case or load\n"\
                          "an existing case before selecting an item")
            w = QMessageBox(self)
            w.information(self,
                          self.tr("Warning"),
                          msg,
                          self.tr('OK'))
            return

        self.page = self.Browser.display(self,
                                         self.case,
                                         self.statusbar,
                                         self.Browser)

        if self.page is not None:
            self.scrollArea.setWidget(self.page)

        else:
            log.debug("displayNewPage() self.page == None")
            raise


    def displayAbout(self):
        """
        public slot

        the About dialog window shows:
         - title
         - version
         - contact
        """
        msg = self.package.code_name + "\n"                 +\
              "version " + self.package.version + "\n\n"    +\
              "For information about this application "     +\
              "please contact:\n\n"                         +\
              self.package.bugreport + "\n\n"               +\
              "Please visit our site:\n"                    +\
              self.package.url
        QMessageBox.about(self, self.package.name + ' study manager', msg)


    def displayLicence(self):
        """
        public slot

        GNU GPL license dialog window
        """
        QMessageBox.about(self, self.package.code_name + ' study manager', "see COPYING file") # TODO


    def displayConfig(self):
        """
        public slot

        configuration information window
        """
        QMessageBox.about(self, self.package.code_name + ' study manager', "see config.py") # TODO


    def setColor(self):
        """
        public slot

        choose GUI color
        """
        c = self.palette().color(QPalette.Window)
        color = QColorDialog.getColor(c, self)
        if color.isValid():
            app = QCoreApplication.instance()
            if not self.palette_default:
                self.palette_default = QPalette().resolve(self.palette())
            app.setPalette(QPalette(color))
            settings = QSettings()
            settings.setValue("MainWindow/Color",
                              self.palette().color(QPalette.Window).name())


    def setFontSize(self):
        """
        public slot

        choose GUI font
        """
        font, ok = QFontDialog.getFont(self)
        log.debug("setFont -> %s" % ok)
        if ok:
            if not self.font_default:
                self.font_default = self.font()
            self.setFont(font)
            app = QCoreApplication.instance()
            app.setFont(font)
            settings = QSettings()
            settings.setValue("MainWindow/Font",
                              self.font().toString())


    def restoreStyleDefaults(self):
        """
        public slot

        Restore default style.
        """

        reply = QMessageBox.question(self, "Restore defaults",
                                     "Restore default color and font ?",
                                     QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            app = QCoreApplication.instance()
            if self.palette_default:
                app.setPalette(self.palette_default)
            if self.font_default:
                print(self.font_default)
                print(self.font())
                self.setFont(self.font_default)
                app.setFont(self.font_default)
            settings = QSettings()
            settings.remove("MainWindow/Color")
            settings.remove("MainWindow/Font")


    def tr(self, text):
        """
        private method

        translation

        @param text: text to translate
        @return: translated text
        """
        return text

#-------------------------------------------------------------------------------
# Main Window for Code_Saturne
#-------------------------------------------------------------------------------

class MainViewSaturne(QMainWindow, Ui_MainForm, MainView):

    def __new__(cls, cmd_package = None, cmd_case = ""):
        return super(MainViewSaturne, cls). __new__(cls, cmd_package, cmd_case)


    def __init__(self,
                 cmd_package      = None,
                 cmd_case         = ""):
        """
        Initializes a Main Window for a new document:
          1. finish the Main Window layout
          2. connection betwenn signal and slot
          3. Ctrl+C signal handler
          4. create some instance variables
          5. restore system settings

        @type cmd_case:
        @param cmd_case:
        """
        QMainWindow.__init__(self)
        Ui_MainForm.__init__(self)

        self.setupUi(self)

        # create some instance variables

        self.cmd_case    = cmd_case
        self.package     = cmd_package

        self.XML_DOC_VERSION = XML_DOC_VERSION

        self.Browser = BrowserView()
        self.ui_initialize()

        self.displayCSManualAction.triggered.connect(self.displayCSManual)
        self.displayCSTutorialAction.triggered.connect(self.displayCSTutorial)
        self.displayCSTheoryAction.triggered.connect(self.displayCSTheory)
        self.displayCSSmgrAction.triggered.connect(self.displayCSSmgr)
        self.displayCSRefcardAction.triggered.connect(self.displayCSRefcard)
        self.displayCSDoxygenAction.triggered.connect(self.displayCSDoxygen)

        docdir = self.package.get_dir('docdir')
        if os.path.isdir(docdir):
            liste = os.listdir(docdir)
        else:
            liste = []

        if 'user.pdf' not in liste:
            self.displayCSManualAction.setEnabled(False)
        if 'theory.pdf' not in liste:
            self.displayCSTheoryAction.setEnabled(False)
        if 'studymanager.pdf' not in liste:
            self.displayCSSmgrAction.setEnabled(False)
        if 'refcard.pdf' not in liste:
            self.displayCSRefcardAction.setEnabled(False)
        if 'doxygen' not in liste:
            self.displayCSDoxygenAction.setEnabled(False)
        self.displayNCManualAction.setVisible(False)


    def initCase(self):
        """
        Initializes the new case with default xml nodes.
        If previous case, just check if all mandatory nodes exist.
        """
        XMLinit(self.case).initialize()


    def displayCSManual(self):
        """
        public slot

        open the user manual
        """
        self.displayManual(self.package, 'user')


    def displayCSTutorial(self):
        """
        public slot

        open the tutorial for Code_Saturne
        """
        msg = "See " + self.package.url + " web site for tutorials."
        QMessageBox.about(self, self.package.name + ' study manager', msg)


    def displayCSTheory(self):
        """
        public slot

        open the theory and programmer's guide
        """
        self.displayManual(self.package, 'theory')

    def displayCSSmgr(self):
        """
        public slot

        open the studymanager guide
        """
        self.displayManual(self.package, 'studymanager')

    def displayCSRefcard(self):
        """
        public slot

        open the quick reference card for Code_Saturne
        """
        self.displayManual(self.package, 'refcard')


    def displayCSDoxygen(self):
        """
        public slot

        open the quick doxygen for Code_Saturne
        """
        self.displayManual(self.package, 'Doxygen')


#-------------------------------------------------------------------------------

def isAlive(qobj):
    """
    return True if the object qobj exist

    @param qobj: the name of the attribute
    @return: C{True} or C{False}
    """
    import sip
    try:
        sip.unwrapinstance(qobj)
    except RuntimeError:
        return False
    return True

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
