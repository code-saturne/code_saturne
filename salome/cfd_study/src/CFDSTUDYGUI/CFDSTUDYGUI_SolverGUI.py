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
Solver GUI
==========

The two solvers I{Code_Saturne} and C{NEPTUNE_CFD} have their own GUI. The
purpose of the class C{CFDSTUDYGUI_SolverGUI} is to display the solver GUI of
the selected code in the SALOME workspace.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, string, logging

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

from CFDSTUDYGUI_Commons import CFD_Code, CFD_Saturne, CFD_Neptune, sgPyQt
from CFDSTUDYGUI_Commons import LoggingMgr
import CFDSTUDYGUI_DataModel
from CFDSTUDYGUI_Management import CFDGUI_Management
import cs_info

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

def getObjectBrowserDock():
    dock = None
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    for dock in dsk.findChildren(QDockWidget):
        dockTitle = str(dock.windowTitle())
        if (dockTitle == 'Object Browser'):
            return dock


def tabObjectBrowser():
    """
    tabify DockWidgets which contains QTreeView:
    Object Browser and CFDSTUDY tree
    """
    dsk = sgPyQt.getDesktop()
    ldock = dsk.findChildren(QDockWidget)
    ldocktree = []
    for i in ldock:
        lo = i.findChildren(QTreeView)
        if len(lo):
            ldocktree.append(i)
    for i in range(1, len(ldocktree)):
        dsk.tabifyDockWidget(ldocktree[0], ldocktree[i])


def updateObjectBrowser():
    """
    force le regroupement en onglets des QTreeView apres updateObjBrowser
    """
    studyId = sgPyQt.getStudyId()
    sgPyQt.updateObjBrowser(studyId, 1)
    tabObjectBrowser()


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
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI.getDocks(studyId) == {}:
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
        lnames = string.split(xml_file, "/")
        if len(lnames) < 4:
            return None
        xmlname   = lnames[-1]
        casename  = lnames[-3]
        studyname = lnames[-4]
        return string.join([studyname, casename, xmlname], ".")


    def getDockTitleNameFromOB(self, studyname, casename, xmlname):
        return string.join([studyname, casename, xmlname], ".")


    def onUndo(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotUndo()


    def onRedo(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotRedo()

    def onPreproMode(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotPreproMode()

    def onCalculationMode(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotCalculationMode()

    def onOTStudyMode(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.slotOpenTurnsMode()

    def onOpenShell(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.openXterm()


    def onDisplayCase(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayCase()


    def onHelpAbout(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayAbout()


    def onSaturneHelpLicense(self):
        if self._CurrentWindow != None:
            self._CurrentWindow.displayLicence()


    def onSaturneHelpManual(self):
        from cs_package import package
        argv_info = ['--guide', 'user']
        cs_info.main(argv_info, package())


    def onSaturneHelpTutorial(self):
        from cs_package import package
        msg = "See http://code-saturne.org web site for tutorials."
        QMessageBox.about(self._CurrentWindow, 'code_saturne Interface', msg)


    def onSaturneHelpKernel(self):
        from cs_package import package
        argv_info = ['--guide', 'theory']
        cs_info.main(argv_info, package())


    def onSaturneHelpRefcard(self):
        from cs_package import package
        argv_info = ['--guide', 'refcard']
        cs_info.main(argv_info, package())


    def onSaturneHelpDoxygen(self):
        from cs_package import package
        argv_info = ['--guide', 'Doxygen']
        cs_info.main(argv_info, package())


    def onNeptuneHelpManual(self):
        from nc_package import package
        argv_info = ['--guide', 'user']
        cs_info.main(argv_info, package())


    def onNeptuneHelpTutorial(self):
        from nc_package import package
        argv_info = ['--guide', 'tutorial']
        cs_info.main(argv_info, package())


    def onNeptuneHelpKernel(self):
        from nc_package import package
        argv_info = ['--guide', 'theory']
        cs_info.main(argv_info, package())


    def onNeptuneHelpDoxygen(self):
        from nc_package import package
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
        from cs_gui import process_cmd_line
        if CFD_Code() == CFD_Saturne:
            from cs_package import package
            from code_saturne.Base.MainView import MainView
        elif CFD_Code() == CFD_Neptune:
            from nc_package import package
            from neptune_cfd.core.MainView import MainView

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
        dock = QDockWidget(aTitle)

        dock.setWidget(mw.frame)
        dock.setMinimumWidth(520)
        dsk.addDockWidget(Qt.RightDockWidgetArea, dock)
        dock.setVisible(True)
        dock.show()

        # Put the QTreeView of the MainView which is already inside a QDockWidget
        # in the SALOME Desktop
        BrowserTitle = aTitle  + " Browser"
        mw.dockWidgetBrowser.setWindowTitle(BrowserTitle)
        dsk.addDockWidget(Qt.LeftDockWidgetArea, mw.dockWidgetBrowser)

        mw.dockWidgetBrowser.setVisible(True)
        mw.dockWidgetBrowser.show()
        mw.dockWidgetBrowser.raise_()
        dock.raise_()

        #Add Dock windows are managed by CFDGUI_Management class
        studyId = sgPyQt.getStudyId()
        aStudyCFD = aCase.GetFather()
        aCaseCFD  = aCase
        xmlFileName = str(Title)
        _c_CFDGUI.set_d_CfdCases(studyId, dock, mw.dockWidgetBrowser, mw, aStudyCFD, aCaseCFD, xmlFileName, sobjXML)
        dock.visibilityChanged["bool"].connect(self.setdockWindowBrowserActivated)
        mw.dockWidgetBrowser.visibilityChanged["bool"].connect(self.setdockWindowActivated)
        dock.toggleViewAction().toggled["bool"].connect(self.setdockWB)
        mw.dockWidgetBrowser.toggleViewAction().toggled["bool"].connect(self.setdock)
        _c_CFDGUI.tabifyDockWindows(dsk, studyId)
        self.showDockWindows(studyId, xmlFileName, aCaseCFD.GetName(), aStudyCFD.GetName())
        updateObjectBrowser()

        return mw


    def setdockWB(self, istoggled):
        """
        istoggled referred to CFD Window dock widget
        """
        log.debug("setdockWB")
        studyId = sgPyQt.getStudyId()
        dock = self.sender().parent()
        log.debug("setdockWB -> %s" % (dock,))

        if _c_CFDGUI != None:
            dockWB = _c_CFDGUI.getdockWB(studyId, dock)
            if dockWB != None:
                dockWB.setVisible(istoggled) #dock.isVisible())
                if istoggled:
                    dock.show()
                    dock.raise_()
                    dockWB.show()
                    dockWB.raise_()
                mw = _c_CFDGUI.getMW(studyId, dock)
                self._CurrentWindow = mw
                mw.activateWindow()
                log.debug("setdockWB -> mw = %s" % (mw,))
        else:
            self._CurrentWindow = None


    def setdock(self, istoggled):
        """
        istoggled referred to Window browser dock widget
        """
        log.debug("setdock")
        studyId = sgPyQt.getStudyId()
        dockWB = self.sender().parent()
        log.debug("setdock -> %s" % (dockWB,))

        if _c_CFDGUI != None:
            dock = _c_CFDGUI.getdock(studyId, dockWB)
            if dock != None:
                dock.setVisible(istoggled) #dockWB.isVisible())
                if istoggled:
                    dock.show()
                    dock.raise_()
                    dockWB.show()
                    dockWB.raise_()
                mw = _c_CFDGUI.getMW(studyId, dock)
                self._CurrentWindow = mw
                mw.activateWindow()
                log.debug("setdock -> mw = %s" % (mw,))
        else:
            self._CurrentWindow = None


    def setdockWindowBrowserActivated(self, visible):
        """
        mw is the Main CFD window allocated by MainView code
        When we click on a cfd study window tab, the cfd study window appears and the associated CFD window browser raises too
        """

        studyId = sgPyQt.getStudyId()
        dock = self.sender()
        log.debug("setdockWindowBrowserActivated -> %s" % (dock,))

        if not visible:
            return
        #if dock.isActiveWindow() == False:
            #return
        if _c_CFDGUI != None:
            dockWB = _c_CFDGUI.getdockWB(studyId, dock)
            if dockWB != None:
                dockWB.activateWindow()
                dockWB.show()
                dockWB.raise_()
                mw = _c_CFDGUI.getMW(studyId, dock)
                self._CurrentWindow = mw
                mw.activateWindow()
                log.debug("setdockWindowBrowserActivated -> mw = %s" % (mw,))
                ob = sgPyQt.getObjectBrowser()
                # Clear the current selection in the SALOME object browser, which does not match with the shown dock window
                ob.clearSelection()
        else:
            self._CurrentWindow = None


    def setdockWindowActivated(self, visible):
        """
        mv is the Main CFD window allocated by MainView code
        When we click on a  CFD window browser tab, the CFD window browser appears and the associated cfd study window raises too
        """
        studyId = sgPyQt.getStudyId()
        dockWB = self.sender()
        log.debug("setdockWindowActivated -> %s" % (dockWB,))

        if not visible:
            return
        #if dockWB.isActiveWindow() == False:
            #return
        if _c_CFDGUI != None:
            dock = _c_CFDGUI.getdock(studyId, dockWB)
            if dock != None:
                dock.activateWindow()
                dock.show()
                dock.raise_()
                mw = _c_CFDGUI.getMW(studyId, dock)
                self._CurrentWindow = mw
                mw.activateWindow()
                log.debug("setdockWindowActivated -> mw = %s" % (mw,))
                ob = sgPyQt.getObjectBrowser()
                # effacer la selection en cours
                ob.clearSelection()
        else:
            self._CurrentWindow = None


    def disconnectDockWindows(self):
        """
        Hide the dock windows of CFDSTUDY GUI, when activating another Salome module
        We can have one or several of them with the right click on the main menu bar of
        Salome
        """
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
          _c_CFDGUI.hideDocks(studyId)


    def connectDockWindows(self):
        """
        Show all the dock windows of CFDSTUDY GUI, when activating Salome CFDSTUDY module
        """
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
            _c_CFDGUI.showDocks(studyId)
            _c_CFDGUI.tabifyDockWindows(sgPyQt.getDesktop(),studyId)


    def showDockWindows(self, studyId, xmlName, caseName, studyCFDName):
        """
        Find if the dockwindow corresponding to this xmlcase is already opened
        """
        if _c_CFDGUI != None:
            _c_CFDGUI.showDockWindows(studyId, xmlName, caseName, studyCFDName)


    def getStudyCaseXmlNames(self, mw):
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
            studyCFDName, caseName, xmlName  = _c_CFDGUI.getStudyCaseXmlNames(studyId, mw)
        return studyCFDName, caseName, xmlName


    def getCase(self, mw):
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
            case  = _c_CFDGUI.getCase(studyId,mw)
        return case


    def removeDockWindowfromStudyAndCaseNames(self, studyCFDName, caseName):
        """
        Close the CFD_study_dock_windows if opened from close Study popup menu into the object browser
        """
        log.debug("removeDockWindow -> %s %s" % (studyCFDName, caseName))
        dsk = sgPyQt.getDesktop()
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
           _c_CFDGUI.delDockfromStudyAndCaseNames(dsk, studyId, studyCFDName, caseName)


    def removeDockWindow(self, studyCFDName, caseName, xmlName):
        """
        Close the CFD_study_dock_windows from remove  popup menu in object browser
        """
        log.debug("removeDockWindow -> %s %s %s" % (studyCFDName, caseName, xmlName))
        dsk = sgPyQt.getDesktop()
        studyId = sgPyQt.getStudyId()
        if _c_CFDGUI != None:
           _c_CFDGUI.delDock(dsk, studyId, studyCFDName, caseName, xmlName)

#-------------------------------------------------------------------------------
