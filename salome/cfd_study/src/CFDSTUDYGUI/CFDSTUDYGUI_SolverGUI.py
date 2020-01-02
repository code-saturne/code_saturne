# -*- coding: utf-8 -*-

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
Solver GUI
==========

The two solvers I{Code_Saturne} and C{NEPTUNE_CFD} have their own GUI. The
purpose of the class C{CFDSTUDYGUI_SolverGUI} is to display the solver GUI of
the selected code in the SALOME workspace.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from CFDSTUDYGUI_Commons import CFD_Code, CFD_Saturne, CFD_Neptune, sgPyQt, sg
from CFDSTUDYGUI_Commons import LoggingMgr
import CFDSTUDYGUI_DataModel
from CFDSTUDYGUI_Management import CFDGUI_Management
from code_saturne import cs_info

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_SolverGUI")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

mw = None

_c_CFDGUI = CFDGUI_Management()

#-------------------------------------------------------------------------------
# Function definitions
#-------------------------------------------------------------------------------


def findObjectBrowserDockWindow():
    dsk = sgPyQt.getDesktop()
    ldock = []
    if dsk != None:
        ldock = dsk.findChildren(QDockWidget)
    objectBrowserDockWindow = None
    if ldock != []:
        for i in ldock:
            if 'Object Browser' in str(i.windowTitle()):
                objectBrowserDockWindow = i
    return objectBrowserDockWindow


def tabifyCfdGui():
    """
    tabify DockWidgets which contains CFD study CASE QMainview :
    CFDSTUDY Main Window
    """
    log.debug("tabifyCfdGui")
    dsk = sgPyQt.getDesktop()
    ldockMainWin = []

    ldockMainWin = _c_CFDGUI.getDockListe()

    objectBrowserDockWindow = findObjectBrowserDockWindow()
    if len(ldockMainWin) >= 1:
        dsk.splitDockWidget(objectBrowserDockWindow,ldockMainWin[0],Qt.Horizontal)
        dsk.tabifyDockWidget(objectBrowserDockWindow,ldockMainWin[0])
    for i in range(1, len(ldockMainWin)):
        dsk.tabifyDockWidget(ldockMainWin[0], ldockMainWin[i])
        dsk.tabifyDockWidget(objectBrowserDockWindow, ldockMainWin[i])


def updateObjectBrowser():
    """
    force le regroupement en onglets des QTreeView apres updateObjBrowser
    """
    sg.updateObjBrowser()
    tabifyCfdGui()


def findDockWindow(xmlName, caseName, studyCFDName):
    """
    Find if the dockwindow corresponding to this xmlcase is already opened
    """
    log.debug("findDockWindow")
    bool_findDockWindow = False

    if _c_CFDGUI != None:
        bool_findDockWindow = _c_CFDGUI.findElem(xmlName, caseName, studyCFDName)

    return bool_findDockWindow


#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_SolverGUI(QObject):
    """
    Auxilliary class for interaction with solvers GUI
    """
    def __init__(self):
        log.debug("CFDSTUDY_SolverGUI.__init__: ")
        QObject.__init__(self, None)
        self._CurrentWindow = None
        self.dockMainWin = None

    def ExecGUI(self, WorkSpace, sobjXML, aCase, Args=''):
        """
        Executes GUI for solver relatively CFDCode
        """
        log.debug("CFDSTUDY_SolverGUI.ExecGUI: ")
        mw = None
        if sobjXML != None:
            #searching
            aTitle = sobjXML.GetName()
            if aCase != None:
                if findDockWindow(aTitle, aCase.GetName(), aCase.GetFather().GetName()):
                    fileN = str(aCase.GetFather().GetName() + "." + aCase.GetName()) + '.' + str(aTitle)
                    mess = "Case file " + fileN + " is already opened"
                    QMessageBox.warning(None, "Warning: ", mess)
                    return
        else:
            aTitle = "unnamed"
            if aCase != None:
                if findDockWindow(aTitle, aCase.GetName(), aCase.GetFather().GetName()):
                    mess = "A new case is already opened"
                    QMessageBox.warning(None, "Warning: ",mess)
                    return
        if aCase != None:
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "^DATA$")
            if not len(aChildList)== 1:
                # no DATA folder
                mess = "DATA directory is not present in the case"
                QMessageBox.warning(None, "Warning: ", mess)
                return None

            aStartPath = CFDSTUDYGUI_DataModel._GetPath(aChildList[0])
            if aStartPath != None and aStartPath != '':
                os.chdir(aStartPath)
        mw = self.launchGUI(WorkSpace, aCase, sobjXML, Args)
        if mw != None:
            self._CurrentWindow = mw

        return mw


    def isActive(self):
        if _c_CFDGUI.getDocks() == {}:
            self._CurrentWindow = None
        if self._CurrentWindow != None:
            return True
        else:
            return False


    def okToContinue(self):
        log.debug("okToContinue")
        if self._CurrentWindow != None and self._CurrentWindow.okToContinue():
            if self._CurrentWindow.case['probes']:
                self._CurrentWindow.case['probes'].removeActors()
            return True
        else:
            return False


    def SaveXmlFile(self):
        log.debug("SaveXmlFile")
        xml_file = None
        if self._CurrentWindow != None:
            self._CurrentWindow.fileSave()
            xml_file = self._CurrentWindow.case['xmlfile']
        return xml_file

    def SaveAsXmlFile(self):
        """
        First: get the xmlfile name with the case (whose path is stored into the MainView Object)
        then save as into tne new xml file (the new name is stored into the case of the MainView Object instead of the old one)
        return old_xml_file,new_xml_file
        """
        old_xml_file = None
        xml_file = None
        if self._CurrentWindow != None:
            old_xml_file = self._CurrentWindow.case['xmlfile']
            self._CurrentWindow.fileSaveAs()
            xml_file = self._CurrentWindow.case['xmlfile']
            if old_xml_file == "":
                old_xml_file = None
            if xml_file == "":
                xml_file = None

        return old_xml_file, xml_file


    def getDockTitleName(self, xml_file):
        """
        Build the Dock Title Name STUDY.CASE.file.xml with the entire file Name path
        """
        lnames = xml_file.split("/")
        if len(lnames) < 4:
            return None
        xmlname   = lnames[-1]
        casename  = lnames[-3]
        studyname = lnames[-4]
        return '.'.join([studyname, casename, xmlname])


    def getDockTitleNameFromOB(self, studyname, casename, xmlname):
        return '.'.join([studyname, casename, xmlname])


    def onUndo(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotUndo()


    def onRedo(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotRedo()

    def onOTStudyMode(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotOpenTurnsMode()

    def onOpenShell(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.openXterm()


    def onDisplayCase(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayCase()


    def onEditSRCFiles(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.fileEditorOpen()


    def onCheckSRCFiles(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.testUserFilesCompilation()


    def onViewLogFiles(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.fileViewerOpen()


    def onLaunchSolver(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.runOrSubmit()


    def onLaunchOT(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.runOTMode()

    def onHelpAbout(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayAbout()


    def onSaturneHelpLicense(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayLicence()


    def onSaturneHelpManual(self):
        from code_saturne.cs_package import package
        argv_info = ['--guide', 'user']
        cs_info.main(argv_info, package())


    def onSaturneHelpTutorial(self):
        from code_saturne.cs_package import package
        msg = "See http://code-saturne.org web site for tutorials."
        QMessageBox.about(self._CurrentWindow, 'code_saturne Interface', msg)


    def onSaturneHelpKernel(self):
        from code_saturne.cs_package import package
        argv_info = ['--guide', 'theory']
        cs_info.main(argv_info, package())


    def onSaturneHelpRefcard(self):
        from code_saturne.cs_package import package
        argv_info = ['--guide', 'refcard']
        cs_info.main(argv_info, package())


    def onSaturneHelpDoxygen(self):
        from code_saturne.cs_package import package
        argv_info = ['--guide', 'Doxygen']
        cs_info.main(argv_info, package())


    def onNeptuneHelpManual(self):
        from neptune_cfd.nc_package import package
        argv_info = ['--guide', 'user']
        cs_info.main(argv_info, package())


    def onNeptuneHelpTutorial(self):
        from neptune_cfd.nc_package import package
        argv_info = ['--guide', 'tutorial']
        cs_info.main(argv_info, package())


    def onNeptuneHelpKernel(self):
        from neptune_cfd.nc_package import package
        argv_info = ['--guide', 'theory']
        cs_info.main(argv_info, package())


    def onNeptuneHelpDoxygen(self):
        from neptune_cfd.nc_package import package
        argv_info = ['--guide', 'Doxygen']
        cs_info.main(argv_info, package())


    def setWindowTitle_CFD(self,mw,aCase,baseTitleName):
        if aCase != None:
            fatherName = aCase.GetFather().GetName()
            aTitle = str(fatherName + "." + aCase.GetName()) + '.' + str(baseTitleName)
            if mw != None:
                mw.setWindowTitle(aTitle)
        return aTitle


    def launchGUI(self, WorkSpace, aCase, sobjXML, Args):
        """
        mw.dockWidgetBrowser is the Browser of the CFD MainView
        """
        log.debug("launchGUI")
        from code_saturne.cs_gui import process_cmd_line
        from code_saturne.Base.MainView import MainView
        if CFD_Code() == CFD_Saturne:
            from code_saturne.cs_package import package
        elif CFD_Code() == CFD_Neptune:
            from neptune_cfd.nc_package import package

        if sobjXML == None:
            Title = "unnamed"
        else:
            Title = sobjXML.GetName()

        self.Workspace = WorkSpace
        pkg = package()
        case, splash = process_cmd_line(Args)
        try:
            mw = MainView(pkg, case, aCase)
        except:
            mess = "Error in Opening CFD GUI"
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, QMessageBox.NoButton)
            return None

        # Put the standard panel of the MainView inside a QDockWidget
        # in the SALOME Desktop
        aTitle = self.setWindowTitle_CFD(mw, aCase, Title)
        dsk = sgPyQt.getDesktop()
#####
        objectBrowserDockWindow = findObjectBrowserDockWindow()

        self.mainWin = QMainWindow()
        self.mainWin.setWindowTitle(aTitle)
        self.mainWin.setCentralWidget(mw.centralwidget)
        self.mainWin.addDockWidget(Qt.LeftDockWidgetArea,mw.dockWidgetBrowser)
#####
        self.dockMainWin = QDockWidget(aTitle)
        self.dockMainWin.setWidget(self.mainWin)
##

        dsk.addDockWidget(Qt.LeftDockWidgetArea,self.dockMainWin)
        self.dockMainWin.setVisible(True)
        self.dockMainWin.show()
        self.dockMainWin.raise_()

        objectBrowserDockWindow.visibilityChanged["bool"].connect(self.resizeObjectBrowserDock)

        #Add Dock windows are managed by CFDGUI_Management class
        aStudyCFD = aCase.GetFather()
        aCaseCFD  = aCase
        xmlFileName = str(Title)
        _c_CFDGUI.set_d_CfdCases(self.dockMainWin, mw, aStudyCFD, aCaseCFD, xmlFileName, sobjXML)
        dockMain = _c_CFDGUI.getDockWithCFDNames(aStudyCFD.GetName(), aCaseCFD.GetName(), xmlFileName)
        if dockMain != None:
            dockMain.visibilityChanged["bool"].connect(self.resizeMainWindowDock)

        updateObjectBrowser()
        return mw


    def resizeObjBrowserDock(self):
        """
        called by closeStudy in CFDSTUDYGUI.py because the Object Browser size is stored into ~/.config/salome/SalomeApprc.xxx file xxx is the salome version
        """
        log.debug("resizeObjectBrowserDock")
        dsk = sgPyQt.getDesktop()
        if dsk != None:
            objectBrowserDockWindow = findObjectBrowserDockWindow()
            if objectBrowserDockWindow != None:
                dsk.resizeDocks({objectBrowserDockWindow}, {300},Qt.Horizontal)


    def resizeMainWindowDock(self,visible):
        """
        visible referred to Object Browser dock widget
        """
        log.debug("resizeMainWindowDock")
        dsk = sgPyQt.getDesktop()
        if dsk != None:
            dock = self.sender()
            if visible:
                dsk.resizeDocks({dock}, {900},Qt.Horizontal)

    def resizeObjectBrowserDock(self,visible):
        """
        visible referred to Object Browser dock widget
        """
        log.debug("resizeObjectBrowserDock")
        dsk = sgPyQt.getDesktop()
        if dsk != None:
            if visible:
                objectBrowserDockWindow = findObjectBrowserDockWindow()
                dsk.resizeDocks({objectBrowserDockWindow}, {300},Qt.Horizontal)

            else:
                if self.dockMainWin != None :
                    dsk.resizeDocks({self.dockMainWin}, {900},Qt.Horizontal)

    def hideDocks(self):
        _c_CFDGUI.hideDocks()
        ob = sgPyQt.getObjectBrowser()
        # Clear the current selection in the SALOME object browser, which does not match with the shown dock window
        if ob != None:
            ob.clearSelection()


    def showDocks(self):
        _c_CFDGUI.showDocks()
        ob = sgPyQt.getObjectBrowser()
        # Clear the current selection in the SALOME object browser, which does not match with the shown dock window
        if ob != None:
            ob.clearSelection()


    def disconnectDockWindows(self):
        """
        Hide the dock windows of CFDSTUDY GUI, when activating another Salome module
        We can have one or several of them with the right click on the main menu bar of
        Salome
        """
        if _c_CFDGUI != None:
            if _c_CFDGUI.d_CfdCases != []:
                self.hideDocks()


    def connectDockWindows(self):
        """
        Show all the dock windows of CFDSTUDY GUI, when activating Salome CFDSTUDY module
        """
        if _c_CFDGUI != None:
            if _c_CFDGUI.d_CfdCases != []:
                self.showDocks()
                tabifyCfdGui()


    def getStudyCaseXmlNames(self, mw):
        if _c_CFDGUI != None:
            studyCFDName, caseName, xmlName  = _c_CFDGUI.getStudyCaseXmlNames(mw)
        return studyCFDName, caseName, xmlName


    def getCase(self, mw):
        if _c_CFDGUI != None:
            case  = _c_CFDGUI.getCase(mw)
        return case


    def removeDockWindowfromStudyAndCaseNames(self, studyCFDName, caseName):
        """
        Close the CFD_study_dock_windows if opened from close Study popup menu into the object browser
        """
        log.debug("removeDockWindowfromStudyAndCaseNames -> %s %s" % (studyCFDName, caseName))
        dsk = sgPyQt.getDesktop()
        if _c_CFDGUI != None:
            _c_CFDGUI.delDockfromStudyAndCaseNames(dsk, studyCFDName, caseName)


    def removeDockWindow(self, studyCFDName, caseName, xmlName):
        """
        Close the CFD_study_dock_windows from remove  popup menu in object browser
        """
        log.debug("removeDockWindow -> %s %s %s" % (studyCFDName, caseName, xmlName))
        dsk = sgPyQt.getDesktop()
        if _c_CFDGUI != None:
            _c_CFDGUI.delDock(dsk, studyCFDName, caseName, xmlName)
#-------------------------------------------------------------------------------
