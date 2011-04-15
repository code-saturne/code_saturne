# -*- coding: utf-8 -*-
#============================================================================
#
#     This file is part of CFDSTUDY the plug-in for Salome
#     of Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     CFDSTUDY is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     CFDSTUDY is distributed in the hope that it will be
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

import os, sys
import string

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtGui import QApplication, QMainWindow, QDockWidget, QTreeView, QMessageBox
from PyQt4.QtCore import Qt, QObject, QEvent, SIGNAL, SLOT

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from CFDSTUDYGUI_Commons import CFD_Code, Trace, CFD_Saturne, CFD_Neptune, sgPyQt
from CFDSTUDYGUI_Commons import LogModeOn, LogModeOff, LoggingMgr, LoggingAgent
import CFDSTUDYGUI_DataModel

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

_d_DockWindows = {}
_d_DockWindowsBrowser = {}
_d_DockWindowsRuncase = {}

_selectedMainViewCase = []
mw=None

_LoggingMgr = LoggingMgr()

#-------------------------------------------------------------------------------
# Function definitions
#-------------------------------------------------------------------------------

def getDockWindowsLists():
    """
    """
    studyId = sgPyQt.getStudyId()
    if studyId in _d_DockWindows.keys() and studyId in _d_DockWindowsBrowser.keys():
        return True #_d_DockWindows[studyId],_d_DockWindowsBrowser[studyId]
    else :
        return False


def tabObjectBrowser():
    """
    Regroupe en onglets les DockWidgets contenant des QTreeView:
    Object Browser et arbres CFDSTUDY
    """
    dsk = sgPyQt.getDesktop()
    ldock=dsk.findChildren(QDockWidget)
    ldocktree=[]
    for i in ldock:
        lo=i.findChildren(QTreeView)
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


def updateDockWindows() :
    """
    force le regroupement en onglets des fenetres d'etudes CFD
    """
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()

    if not getDockWindowsLists():
        return
    if len(_d_DockWindows[studyId]) > 1:
        for i in range(1,len(_d_DockWindows[studyId])):
            dsk.tabifyDockWidget(_d_DockWindows[studyId][0], _d_DockWindows[studyId][i])


def updateDockWindowsBrowser() :
    """
    force le regroupement en onglets des fenetres Browser d'etudes CFD
    """
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()

    if not getDockWindowsLists():
        return
    if len(_d_DockWindowsBrowser[studyId]) > 1 :
        for i in range(1,len(_d_DockWindowsBrowser[studyId])):
            dsk.tabifyDockWidget(_d_DockWindowsBrowser[studyId][0], _d_DockWindowsBrowser[studyId][i])


def findDockWindow(xmlName,caseName,studyCFDName) :
    """
    find if the dockwindow corresponding to this xmlcase is already opened
    """
    import string
    name = string.join([studyCFDName,caseName,xmlName],".")
    bool_findDockWindow = False
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    ldockWindows =dsk.findChildren(QDockWidget)
    for dock in ldockWindows:
        dockTitle = dock.windowTitle()
        if studyId not in _d_DockWindowsBrowser.keys():
            _d_DockWindowsBrowser[studyId] = []
        if (str(dockTitle) == 'Object Browser') and (dock not in _d_DockWindowsBrowser[studyId]):
            _d_DockWindowsBrowser[studyId].append(dock)
        if studyId not in _d_DockWindows.keys():
            _d_DockWindows[studyId] = []

        if _d_DockWindows[studyId] != [] :
            for dock in _d_DockWindows[studyId] :
                if str(name) in str(dock.windowTitle()) :
                    bool_findDockWindow = True
                    dock.show()
                    dock.raise_()

    return bool_findDockWindow


def removeDockWindow(caseName) :
    """
    A TESTER
    remove the CFD_study_dock_windows from remove  popup menu in object browser
    """
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    ldockWindows =dsk.findChildren(QDockWidget)
    for dock in ldockWindows:
        dockTitle = dock.windowTitle()
        if string.split(str(dockTitle),".")[0] == str(caseName):
            if "Browser" not in str(dockTitle) :
                _d_DockWindows[studyId].remove(dock)
                dsk.removeDockWidget(dock)
                dock.setParent(None)
            else :
                _d_DockWindowsBrowser[studyId].remove(dock)
                dsk.removeDockWidget(dock)
                dock.setParent(None)
        updateObjectBrowser()


#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_SolverGUI(QObject):
    """
    Auxilliary class for interaction with solvers GUI
    """
    def __init__(self):
        if Trace() : print "CFDSTUDY_SolverGUI.__init__ : "
        QObject.__init__(self, None)
        self._WindowsMap = {}
        self._CurrentWindow = None
        _d_DockWindows = {}
        _d_DockWindowsBrowser = {}

    def ExecGUI(self, WorkSpace, aTitle, aCase, Args=''):
        """
        Executes GUI for solver relatively CFDCode
        """
        if Trace() : print "CFDSTUDY_SolverGUI.ExecGUI : "
        mw = None
        if aTitle != None and aTitle != "":
            #searching in existing tabs
            win = self._findWindow(aTitle,aCase)
            if win != None:
                #need for activation
                QApplication.postEvent(win, QEvent(QEvent.FocusIn))
                return win
        else:
            aTitle = "unnamed"
            if findDockWindow(aTitle,aCase.GetName(),aCase.GetFather().GetName()) :
                mess = "A case is not finished to be set"
                QMessageBox.warning(None, "Warning : ",mess)
                return

        if aCase != None:
            if CFD_Code() == CFD_Saturne:
                # object of DATA folder
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "^DATA$")
                if not len(aChildList)== 1:
                    # no DATA folder
                    if Trace(): print "CFDSTUDYGUI_SolverGUI.ExecGUI :There are not data folder in selected by user case"
                    return None
                aStartPath = CFDSTUDYGUI_DataModel._GetPath(aChildList[0])
            elif CFD_Code() == CFD_Neptune:
                aStartPath = CFDSTUDYGUI_DataModel._GetPath(aCase)
            if aStartPath != None and aStartPath != '':
                os.chdir(aStartPath)
        if CFD_Code() == CFD_Saturne:
            mw = self._ExecICS(WorkSpace, aCase, aTitle, Args)
        elif CFD_Code() == CFD_Neptune:
            mw = self._ExecIPB(WorkSpace, aTitle, Args)
        if mw != None:
            self._WindowsMap[mw] = aCase
            self._CurrentWindow = mw
            _selectedMainViewCase.append(mw)
        return mw


    def eventFilter(self, anObject, anEvent):
        if Trace() : print "CFDSTUDY_SolverGUI.eventFilter : "
        if anEvent.type() == QEvent.Close:
            del self._WindowsMap[anObject]
        elif anEvent.type() == QEvent.Hide:
            if self._CurrentWindow == anObject:
                self._CurrentWindow = None
        elif anEvent.type() == QEvent.Show:
            self._CurrentWindow = anObject
        return False


    def isActive(self):
        return self._CurrentWindow != None


    def onSaveXmlFile(self):
        if Trace() : print "CFDSTUDY_SolverGUI.onSaveXmlFile : "
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.fileSave()

    def onSaveAsXmlFile(self):
        """
        """
        old_xml_file = None
        xml_file = None
        if  len(_selectedMainViewCase) != 0:
            _sMainViewCase = self._CurrentWindow
            if CFD_Code() == CFD_Saturne:
                old_xml_file = _sMainViewCase.case['xmlfile']
                _sMainViewCase.fileSaveAs()
                xml_file = _sMainViewCase.case['xmlfile']
            if xml_file != None and xml_file != "" and xml_file != old_xml_file:
                #need update object browser
                case = self._WindowsMap[_sMainViewCase]
                study = CFDSTUDYGUI_DataModel.GetStudyByObj(case)

                obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(study, xml_file)
                if obj:
                    #changes Title of the tab
                    new_title = os.path.basename(xml_file)
                    _sMainViewCase.setWindowTitle(new_title)

                    #updates Object Browser
                    CFDSTUDYGUI_DataModel._RebuildTreeRecursively(obj)
                    #if os.path.exists(old_xml_file) and os.path.exists(xml_file) :
                    if os.path.exists(xml_file) :
                        self.replaceDockTitleName(xml_file,old_xml_file,case,study)

                    updateObjectBrowser()


    def getDockTitleName(self,xml_file) :
        """
        Build the Dock Title Name STUDY.CASE.file.xml with the entire file Name
        """
        lnames = string.split(xml_file,"/")
        if len(lnames) < 4 : return None
        xmlname   = lnames[-1]
        casename  = lnames[-3]
        studyname = lnames[-4]
        return string.join([studyname,casename,xmlname],".")


    def replaceDockTitleName(self,new_xml_file,old_xml_file,case,study) :
        """
        replace dock title name in the title of the dock widget and update
        the _d_DockWindows and _d_DockWindowsBrowser dictionary
        """
        OldDockTitleName = self.getDockTitleName(old_xml_file)
        NewDockTitleName = self.getDockTitleName(new_xml_file)

        if NewDockTitleName == None :
            mess = "File : "+xml_file+ \
                   " is not stored into a CFDSTUDY directory structure like .../STUDY_CFD/CASE/DATA/filename.xml"
            QMessageBox.warning(None, "File Error : ",mess)
            return
        if OldDockTitleName == None :
            studyname = study.GetName()
            casename  = case.GetName()
            xmlname = "unnamed"
            OldDockTitleName = string.join([studyname,casename,xmlname],".")
        if NewDockTitleName != None :
            studyId = sgPyQt.getStudyId()
            if studyId in _d_DockWindows.keys():
                for dock in _d_DockWindows[studyId]:
                    if str(OldDockTitleName) in str(dock.windowTitle()) :
                        dock.setWindowTitle(str(NewDockTitleName))
                        dock.show()
                        dock.raise_()
            if studyId in _d_DockWindowsBrowser.keys():
                for dock in _d_DockWindowsBrowser[studyId]:
                        if str(OldDockTitleName) in str(dock.windowTitle()) :
                        dock.setWindowTitle(string.join([str(NewDockTitleName),"Browser"]))
                        dock.show()
                        dock.raise_()
        return


    def onOpenShell(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onOpenShell : "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.openXterm()


    def onDisplayCase(self):
        if Trace() : print "CFDSTUDY_SolverGUI.onDisplayCase : "
        _LoggingMgr.start(sys)
        if self._CurrentWindow != None:

            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCase()
        _LoggingMgr.finish(sys)

    def onHelpAbout(self):
        if Trace() : print "CFDSTUDY_SolverGUI.onHelpAbout : "
        if self._CurrentWindow != None:
           if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayAbout()

#-----------------------------------------------------------------------------

    def onSaturneReloadModule(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneReloadModule: "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.reload_modules()
        return

    def onSaturneReloadPage(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneReloadPage: "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.reload_page()
        return

    def onSaturneHelpLicense(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneHelpLicense: "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayLicence()
        return

    def onSaturneHelpCS(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneHelpcs: "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCSManual()
        return

    def onSaturneHelpSD(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneHelpSD: "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayECSManual()
        return


    def onSaturneHelpCS_Kernel(self):
        """
        """
        if Trace() : print "CFDSTUDY_SolverGUI.onSaturneHelpCS_Kernel : "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCSKernel()

        return


    def onSaturneHelpCS_Infos(self):
        """
        """
        if Trace() : print "CFDSTUDYGUI_SolverGUI.onSaturneHelpCS_INFOS : "
        if self._CurrentWindow != None :
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayECSInfos()

        return


    def onNeptuneWinBrowser(self, flag):
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Neptune:
                self._CurrentWindow.browserDockDisplay(flag)


    def onNeptuneWinIdenty(self,flag):
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Neptune:
                self._CurrentWindow.identityDockDisplay(flag)


    def _ExecIPB(self, WorkSpace, Title, Args):
        """
        A developper
        """
        pass


    def _ExecICS(self, WorkSpace, aCase, Title, Args):
        """
        """
        from Base.MainView import MainView
        from cs_gui import process_cmd_line

        self.Workspace = WorkSpace
        case, splash, batch_window, batch_file, tree_window, read_only = process_cmd_line(Args)

        mw = MainView(case, batch_window, batch_file, tree_window, read_only, aCase)
        fatherName = aCase.GetFather().GetName()

        aTitle = str(fatherName + "." + aCase.GetName()) + '.' + str(Title)
        mw.setWindowTitle(aTitle)

        dsk = sgPyQt.getDesktop()
        dock = QDockWidget(aTitle)

        dock.setWidget(mw.frame)
        dock.setMinimumWidth(520)
        dsk.addDockWidget(Qt.RightDockWidgetArea, dock)

        studyId = sgPyQt.getStudyId()

        if studyId not in _d_DockWindows.keys():
            _d_DockWindows[studyId] = []
        _d_DockWindows[studyId].append(dock)

        dock.setVisible(True)
        dock.show()
        updateDockWindows()

        BrowserTitle = aTitle  + " Browser"
        mw.dockWidgetBrowser.setWindowTitle(BrowserTitle)
        dsk.addDockWidget(Qt.LeftDockWidgetArea,mw.dockWidgetBrowser)

        if studyId not in _d_DockWindowsBrowser.keys():
            _d_DockWindowsBrowser[studyId] = []
        _d_DockWindowsBrowser[studyId].append(mw.dockWidgetBrowser)

        mw.dockWidgetBrowser.setVisible(True)
        mw.dockWidgetBrowser.show()
        mw.dockWidgetBrowser.raise_()
        dock.raise_()

        self.connect(dock, SIGNAL("visibilityChanged(bool)"), self.setdockWindowBrowserActivated)
        self.connect(mw.dockWidgetBrowser, SIGNAL("visibilityChanged(bool)"),self.setdockWindowActivated)

        updateDockWindowsBrowser()
        updateObjectBrowser()
        return mw


    def _findWindow(self, aTitle, aCase):

        if aTitle == '' or aCase == None:
            return None

        for win in self._WindowsMap.keys():
            if aTitle == str(win.windowTitle().toLatin1()) and \
                       aCase.GetID() == self._WindowsMap[win].GetID():
                return win

        return None


    def setdockWindowBrowserActivated(self,visible) :

        if not visible : return
        dock = self.sender()
        if dock.isActiveWindow() == False : return
        titledock = str(dock.windowTitle())
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        for i in range(len(_d_DockWindowsBrowser[studyId])):
            if titledock in str(_d_DockWindowsBrowser[studyId][i].windowTitle()):
                _d_DockWindowsBrowser[studyId][i].activateWindow()
                _d_DockWindowsBrowser[studyId][i].setVisible(True)
                _d_DockWindowsBrowser[studyId][i].show()
                _d_DockWindowsBrowser[studyId][i].raise_()
                self.activateCurrentWindow(titledock)


    def setdockWindowActivated(self,visible) :

        if not visible : return
        dock = self.sender()
        if dock.isActiveWindow() == False : return
        title = str(dock.windowTitle())
        titledock,br = string.split(title," Browser")
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        for i in range(len(_d_DockWindows[studyId])):
            if str(_d_DockWindows[studyId][i].windowTitle()) in str(title):
                _d_DockWindows[studyId][i].activateWindow()
                _d_DockWindows[studyId][i].setVisible(True)
                _d_DockWindows[studyId][i].show()
                _d_DockWindows[studyId][i].raise_()
                self.activateCurrentWindow(titledock)


    def activateCurrentWindow(self,title) :
        """
        """
        for mw in self._WindowsMap.keys() :
            casename = mw.case['salome'].GetName()
            xmlfile = mw.case['xmlfile']
            if xmlfile == "" or xmlfile == None :
                xmlfileName = "unnamed"
            else :
                xmlfileName = os.path.basename(xmlfile)
            fatherCaseName = mw.case['salome'].GetFather().GetName()
            if title == string.join([fatherCaseName,casename,xmlfileName],".") :
                self._CurrentWindow = mw
                if mw not in _selectedMainViewCase :
                    _selectedMainViewCase.append(mw)
                mw.activateWindow()


    def disconnectDockWindows(self) :
        """
        Hide the dock windows of CFDSTUDY GUI, when activating another Salome Component
        We can have one or several of them with the right click on the main menu bar of
        Salome
        """
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        if studyId not in _d_DockWindows.keys():
            return

        if len(_d_DockWindows[studyId]) != 0 :
            for dock in _d_DockWindows[studyId]:
                dock.hide()

        if len(_d_DockWindowsBrowser[studyId]) != 0 :
            for dock in _d_DockWindowsBrowser[studyId]:
                if dock.windowTitle() != 'Object Browser':
                    dock.hide()
        if studyId not in _d_DockWindowsRuncase.keys():
            return
        if len(_d_DockWindowsRuncase[studyId]) != 0:
            for dock in _d_DockWindowsRuncase[studyId]:
                dock.hide()


    def connectDockWindows(self) :
        """
        Show all the dock windows of CFDSTUDY GUI, when activating another Salome Component
        visualise les dock windows  lors d'un changement de composant
        """
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        if studyId not in _d_DockWindows.keys():
            return

        if len(_d_DockWindows[studyId]) != 0:
            for dock in _d_DockWindows[studyId]:
                dock.show()
                dock.setVisible(True)

        if len(_d_DockWindowsBrowser[studyId]) != 0:
            for dock in _d_DockWindowsBrowser[studyId]:
                dock.show()
                dock.setVisible(True)

        if studyId not in _d_DockWindowsRuncase.keys():
            return
        if len(_d_DockWindowsRuncase[studyId]) != 0:
            for dock in _d_DockWindowsRuncase[studyId]:
                dock.show()
                dock.setVisible(True)
