# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_SolverGUI")
log.setLevel(logging.DEBUG)
#log.setLevel(logging.NOTSET)

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
    else:
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


def updateDockWindows():
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


def updateDockWindowsBrowser():
    """
    force le regroupement en onglets des fenetres Browser d'etudes CFD
    """
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()

    if not getDockWindowsLists():
        return
    if len(_d_DockWindowsBrowser[studyId]) > 1:
        for i in range(1,len(_d_DockWindowsBrowser[studyId])):
            dsk.tabifyDockWidget(_d_DockWindowsBrowser[studyId][0], _d_DockWindowsBrowser[studyId][i])


def findDockWindow(xmlName, caseName, studyCFDName):
    """
    Find if the dockwindow corresponding to this xmlcase is already opened
    """
    bool_findDockWindow = False
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    for dock in dsk.findChildren(QDockWidget):
        dockTitle = str(dock.windowTitle())
        if studyId not in _d_DockWindowsBrowser.keys():
            _d_DockWindowsBrowser[studyId] = []
        if (dockTitle == 'Object Browser') and (dock not in _d_DockWindowsBrowser[studyId]):
            _d_DockWindowsBrowser[studyId].append(dock)
        if studyId not in _d_DockWindows.keys():
            _d_DockWindows[studyId] = []

        if _d_DockWindows[studyId] != []:
            for dock in _d_DockWindows[studyId]:
                if string.rstrip(string.join([studyCFDName, caseName, xmlName], ".")) == string.rstrip(dockTitle) :
                    bool_findDockWindow = True
                    dock.show()
                    dock.raise_()

    return bool_findDockWindow

def update_selectedMainViewCase_list(studyCFDName, caseName, xmlName) :
    """
    """
    lind = []
    ind = 0
    for win in _selectedMainViewCase :
        xmlfile = os.path.basename(win.case['xmlfile'])
        boo = xmlfile.rstrip() == xmlName and win.salome.GetName()==caseName and win.salome.GetFather().GetName() == studyCFDName
        if boo :
            lind.append(ind)
        ind = ind + 1
    if len(lind) > 0 :
        for i in lind :
            del _selectedMainViewCase[i]

def removeDockWindow(studyCFDName, caseName, xmlName=""):
    """
    Close the CFD_study_dock_windows from remove  popup menu in object browser
    """
    log.debug("removeDockWindow -> caseName = %s" % caseName)
    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    for dock in dsk.findChildren(QDockWidget):
        dockTitle = str(dock.windowTitle())
        log.debug("removeDockWindow -> dockTitle = %s" % dockTitle)
        stringName = string.rstrip(string.join([studyCFDName, caseName, xmlName], "."))
        stringNameB = stringName + " Browser"
        if stringName  == string.rstrip(dockTitle) or stringNameB == string.rstrip(dockTitle):
            log.debug("removeDockWindow -> widget to close = %s" % dockTitle)
            if "Browser" not in str(dockTitle):
                _d_DockWindows[studyId].remove(dock)
            else:
                _d_DockWindowsBrowser[studyId].remove(dock)
            dsk.removeDockWidget(dock)
            dock.setParent(None)
            dock.close()
        updateObjectBrowser()

#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_SolverGUI(QObject):
    """
    Auxilliary class for interaction with solvers GUI
    """
    def __init__(self):
        if Trace(): print "CFDSTUDY_SolverGUI.__init__: "
        QObject.__init__(self, None)
        self._WindowsMap = {}
        self._CurrentWindow = None
        _d_DockWindows = {}
        _d_DockWindowsBrowser = {}


    def ExecGUI(self, WorkSpace, sobjXML, aCase, Args=''):
        """
        Executes GUI for solver relatively CFDCode
        """
        if Trace(): print "CFDSTUDY_SolverGUI.ExecGUI: "
        mw = None
        if sobjXML != None :
            #searching in existing tabs
            aTitle = sobjXML.GetName()
            win = self._findWindow(aTitle,aCase)
            if win != None:
                #need for activation
                QApplication.postEvent(win, QEvent(QEvent.FocusIn))
                return win
        else:
            aTitle = "unnamed"
            if findDockWindow(aTitle,aCase.GetName(),aCase.GetFather().GetName()):
                mess = "A case is not finished to be set"
                QMessageBox.warning(None, "Warning: ",mess)
                return

        if aCase != None:
            if CFD_Code() == CFD_Saturne:
                # object of DATA folder
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "^DATA$")
                if not len(aChildList)== 1:
                    # no DATA folder
                    if Trace(): print "CFDSTUDYGUI_SolverGUI.ExecGUI:There are not data folder in selected by user case"
                    return None
                aStartPath = CFDSTUDYGUI_DataModel._GetPath(aChildList[0])
            elif CFD_Code() == CFD_Neptune:
                aStartPath = CFDSTUDYGUI_DataModel._GetPath(aCase)
            if aStartPath != None and aStartPath != '':
                os.chdir(aStartPath)
        if CFD_Code() == CFD_Saturne:
            mw = self._ExecICS(WorkSpace, aCase, sobjXML, Args)
        elif CFD_Code() == CFD_Neptune:
            mw = self._ExecIPB(WorkSpace, aTitle, Args)
        if mw != None:
            self._WindowsMap[mw] = aCase
            self._CurrentWindow = mw
            _selectedMainViewCase.append(mw)
        return mw


    def eventFilter(self, anObject, anEvent):
        log.debug("eventFilter")
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
        log.debug("onSaveXmlFile")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                #print "self._CurrentWindow.case['xmlfile'] = ",self._CurrentWindow.case['xmlfile']
                if self._CurrentWindow.case['xmlfile'] != "":
                    self._CurrentWindow.fileSave()
                else:
                    self.SaveAsXmlFile()

    def SaveAsXmlFile(self):
        """
        First : get the xmlfile name with the case (whose path is stored into the MainView Object)
        then save as into tne new xml file (the new name is stored into the case of the MainView Object instead of the old one)
        return old_xml_file,new_xml_file
        """
        old_xml_file = None
        xml_file = None
        NewSObj = None
        OldSobj = None
        if  len(_selectedMainViewCase) != 0:
            _sMainViewCase = self._CurrentWindow
            if CFD_Code() == CFD_Saturne:
                old_xml_file = _sMainViewCase.case['xmlfile']
                _sMainViewCase.fileSaveAs()
                xml_file = _sMainViewCase.case['xmlfile']
                if old_xml_file == "" :
                    old_xml_file = None


        return old_xml_file,xml_file


    def getDockTitleName(self,xml_file):
        """
        Build the Dock Title Name STUDY.CASE.file.xml with the entire file Name
        """
        lnames = string.split(xml_file,"/")
        if len(lnames) < 4: return None
        xmlname   = lnames[-1]
        casename  = lnames[-3]
        studyname = lnames[-4]
        return string.join([studyname,casename,xmlname],".")

    def getDockTitleNameFromOB(self,studyname,casename,xmlname) :
        return string.join([studyname,casename,xmlname],".")

    def replaceDockTitleName(self,new_xml_file,old_xml_file,case,study):
        """
        replace dock title name in the title of the dock widget and update
        the _d_DockWindows and _d_DockWindowsBrowser dictionary
        """
        OldDockTitleName = self.getDockTitleName(old_xml_file)
        NewDockTitleName = self.getDockTitleName(new_xml_file)

        if NewDockTitleName == None:
            mess = "File: "+xml_file+ \
                   " is not stored into a CFDSTUDY directory structure like .../STUDY_CFD/CASE/DATA/filename.xml"
            QMessageBox.warning(None, "File Error: ",mess)
            return
        if OldDockTitleName == None:
            studyname = study.GetName()
            casename  = case.GetName()
            xmlname = "unnamed"
            OldDockTitleName = string.join([studyname,casename,xmlname],".")
        if NewDockTitleName != None:
            studyId = sgPyQt.getStudyId()
            if studyId in _d_DockWindows.keys():
                for dock in _d_DockWindows[studyId]:
                    if str(OldDockTitleName) in str(dock.windowTitle()):
                        dock.setWindowTitle(str(NewDockTitleName))
                        dock.show()
                        dock.raise_()
            if studyId in _d_DockWindowsBrowser.keys():
                for dock in _d_DockWindowsBrowser[studyId]:
                    if str(OldDockTitleName) in str(dock.windowTitle()):
                        dock.setWindowTitle(string.join([str(NewDockTitleName),"Browser"]))
                        dock.show()
                        dock.raise_()
        return


    def onOpenShell(self):
        """
        """
        log.debug("onOpenShell")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.openXterm()


    def onDisplayCase(self):
        log.debug("onDisplayCase")
        _LoggingMgr.start(sys)
        if self._CurrentWindow != None:

            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCase()
        _LoggingMgr.finish(sys)

    def onHelpAbout(self):
        log.debug("onHelpAbout")
        if self._CurrentWindow != None:
           if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayAbout()

#-----------------------------------------------------------------------------

    def onSaturneReloadModule(self):
        """
        """
        log.debug("onSaturneReloadModule")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.reload_modules()
        return

    def onSaturneReloadPage(self):
        """
        """
        log.debug("CFDSTUDY_SolverGUI.onSaturneReloadPage")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.reload_page()
        return

    def onSaturneHelpLicense(self):
        """
        """
        log.debug("onSaturneHelpLicense")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayLicence()
        return

    def onSaturneHelpCS(self):
        """
        """
        log.debug("onSaturneHelpcs")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCSManual()
        return

    def onSaturneHelpSD(self):
        """
        """
        log.debug("onSaturneHelpSD")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayECSManual()
        return


    def onSaturneHelpCS_Kernel(self):
        """
        """
        log.debug("onSaturneHelpCS_Kernel")
        if self._CurrentWindow != None:
            if CFD_Code() == CFD_Saturne:
                self._CurrentWindow.displayCSKernel()

        return


    def onSaturneHelpCS_Infos(self):
        """
        """
        log.debug("onSaturneHelpCS_INFOS")
        if self._CurrentWindow != None:
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


    def setWindowTitle_CFD(self,mw,aCase,baseTitleName) :
        """
        """
        if aCase != None :
            fatherName = aCase.GetFather().GetName()
            aTitle = str(fatherName + "." + aCase.GetName()) + '.' + str(baseTitleName)
            if mw != None :
                mw.setWindowTitle(aTitle)
        return aTitle


    def _ExecICS(self, WorkSpace, aCase, sobjXML, Args):
        """
        """
        log.debug("_ExecICS")
        from cs_gui import process_cmd_line
        from cs_package import package
        from Base.MainView import MainView
        if sobjXML == None :
            Title = "unnamed"
        else :
            Title = sobjXML.GetName()
        self.Workspace = WorkSpace
        pkg = package()
        case, splash, batch_window, batch_file, tree_window, read_only = process_cmd_line(Args)
        mw = MainView(pkg, case, batch_window, batch_file, tree_window, read_only, aCase)

        aTitle = self.setWindowTitle_CFD(mw,aCase,Title)
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
        ##dock.installEventFilter(self)
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
        """
        """
        log.debug("_findWindow")
        if aTitle == '' or aCase == None:
            return None

        for win in self._WindowsMap.keys():
            if aTitle == str(win.windowTitle().toLatin1()) and \
                       aCase.GetID() == self._WindowsMap[win].GetID():
                return win

        return None


    def setdockWindowBrowserActivated(self,visible):
        """
        """
        if not visible: return
        dock = self.sender()
        if dock.isActiveWindow() == False: return
        titledock = str(dock.windowTitle())
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        titledockBr = string.rstrip(titledock) + " Browser"
        for i in range(len(_d_DockWindowsBrowser[studyId])):

            if string.rstrip(titledockBr) == string.rstrip(str(_d_DockWindowsBrowser[studyId][i].windowTitle())):
                _d_DockWindowsBrowser[studyId][i].activateWindow()
                _d_DockWindowsBrowser[studyId][i].setVisible(True)
                _d_DockWindowsBrowser[studyId][i].show()
                _d_DockWindowsBrowser[studyId][i].raise_()
                self.activateCurrentWindow(titledockBr)


    def setdockWindowActivated(self,visible):
        """
        """
        if not visible: return
        dock = self.sender()
        if dock.isActiveWindow() == False: return
        title = str(dock.windowTitle())
        titledock,br = string.split(title," Browser")
        studyId = sgPyQt.getStudyId()
        if studyId not in _d_DockWindowsBrowser.keys():
            return
        for i in range(len(_d_DockWindows[studyId])):
            if string.rstrip(str(_d_DockWindows[studyId][i].windowTitle())) == string.rstrip(str(titledock)):
                _d_DockWindows[studyId][i].activateWindow()
                _d_DockWindows[studyId][i].setVisible(True)
                _d_DockWindows[studyId][i].show()
                _d_DockWindows[studyId][i].raise_()
                self.activateCurrentWindow(titledock)


    def activateCurrentWindow(self,title):
        """
        """
        for mw in self._WindowsMap.keys():

            casename = mw.salome.GetName()
            xmlfile = os.path.basename(mw.case['xmlfile'])
            xmlfile = xmlfile.rstrip()
            if xmlfile == "" or xmlfile == None:
                xmlfileName = "unnamed"
            else:
                xmlfileName = xmlfile
            fatherCaseName = mw.salome.GetFather().GetName()
            if string.rstrip(title) == string.join([fatherCaseName,casename,xmlfileName],"."):
                self._CurrentWindow = mw
                if mw not in _selectedMainViewCase:
                    _selectedMainViewCase.append(mw)
                mw.activateWindow()


    def disconnectDockWindows(self):
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

        if len(_d_DockWindows[studyId]) != 0:
            for dock in _d_DockWindows[studyId]:
                dock.hide()

        if len(_d_DockWindowsBrowser[studyId]) != 0:
            for dock in _d_DockWindowsBrowser[studyId]:
                if dock.windowTitle() != 'Object Browser':
                    dock.hide()
        if studyId not in _d_DockWindowsRuncase.keys():
            return
        if len(_d_DockWindowsRuncase[studyId]) != 0:
            for dock in _d_DockWindowsRuncase[studyId]:
                dock.hide()


    def connectDockWindows(self):
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


    def update_WindowsMap_dict(self,aCase,aXmlFileName) :
        l_win = []
        for win in self._WindowsMap.keys() :
            xmlfile = os.path.basename(win.case['xmlfile'])
            boo = self._WindowsMap[win].GetName() == aCase.GetName() and xmlfile.rstrip() == aXmlFileName and self._WindowsMap[win].GetFather().GetName() == aCase.GetFather().GetName()
            if boo :
                l_win.append(win)
        if len(l_win) != 0 :
            for w in l_win :
                del self._WindowsMap[w]

