# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
Dialog Collector
================

This file gathers the C{QDialog} definitions of the CFD_STUDY module.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, re, shutil, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QMessageBox, QDialog, QPushButton, QToolTip, QColor, \
                        QRadioButton, QTableWidget, QTableWidgetItem
from PyQt4.QtCore import Qt, QStringList, QString, SIGNAL

#For Testing
from PyQt4.QtCore import QTranslator

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from ui_InfoDialog            import Ui_InfoDialog
from ui_SetTreeLocationDialog import Ui_SetTreeLocationDialog
from ui_RunCaseDialog         import Ui_RunCaseDialog
from ui_ECSConversionDialog   import Ui_ECSConversionDialog
from ui_CopyDialog            import Ui_CopyDialog
from ui_GUIActivationDialog   import Ui_GUIActivationDialog
import CFDSTUDYGUI_DataModel

from CFDSTUDYGUI_Commons import _SetCFDCode, CFD_Code, sgPyQt
from CFDSTUDYGUI_Commons import CFD_Saturne, CFD_Neptune, CheckCFD_CodeEnv

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_DialogCollector")
log.setLevel(logging.DEBUG)
#log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Dialog definitions
#-------------------------------------------------------------------------------

class InfoDialog(QtGui.QDialog, Ui_InfoDialog):
    """
    Dialog informations about solver installation.
    """
    def __init__(self, parent = None):
        """
        """
        QtGui.QDialog.__init__(self, parent)
        Ui_InfoDialog.__init__(self)

        self.setupUi(self)


class InfoDialogHandler(InfoDialog):
    """
    """
    def __init__(self, parent = None):
        """
        """
        InfoDialog.__init__(self,parent)

        self.status = 1 #initial access status

        aBtn = self.findChild(QtGui.QPushButton, "OKButton")
        if aBtn != None:
            aBtn.setText(self.tr("DLG_OK_BUTTON_TEXT"))

        codeBG = self.findChild(QtGui.QButtonGroup, "CodeBG")
        if codeBG != None:
            codeBG.setTitle(self.tr("INFO_DLG_CFDCODE_TITLE"))

        self.SaturneRB = self.findChild(QtGui.QRadioButton, "SaturneRB")
        if self.SaturneRB != None:
            self.SaturneRB.setText(self.tr("INFO_DLG_SATURNE_TEXT"))
            self.SaturneRB.setChecked(True)

        self.NeptuneRB = self.findChild(QtGui.QRadioButton, "NeptuneRB")
        if self.NeptuneRB != None:
            self.NeptuneRB.setText(self.tr("INFO_DLG_NEPTUNE_TEXT"))

        self.setWindowTitle(self.tr("INFO_DLG_CAPTION"))


    def accept(self):
        iok, mess = CheckCFD_CodeEnv(CFD_Code())
        if iok :
            if mess != "" :
                Error = "Error : "+ self.tr("CFDSTUDY_INVALID_ENV")
                QMessageBox.critical(ActionHandler.dskAgent().workspace(),
                                 Error, mess, QMessageBox.Ok, 0)
            else :
                InfoDialog.accept(self)
                #block other code
                self.setCode(CFD_Code(), True)
        else:
            Error = "Error : " + self.tr("INFO_DLG_INVALID_ENV")
            QMessageBox.critical(self, Error, mess, QMessageBox.Ok, 0)


    def setCode(self, code, isDisableOther):
        if code == CFD_Saturne:
            self.SaturneRB.setEnabled(True)
            self.SaturneRB.setChecked(True)
            #self.NeptuneRB.setEnabled(not isDisableOther)
            from cs_package import package
        elif code == CFD_Neptune:
            self.NeptuneRB.setEnabled(True)
            self.NeptuneRB.setChecked(True)
            #self.SaturneRB.setEnabled(not isDisableOther)
            from nc_package import package
        else:
            raise DialogError, "Invalid CFD_Code in InfoDialog class"
        pkg = package()
        self.labelVersionValue.setText(pkg.version)
        self.labelPrefixValue.setText(pkg.dirs['prefix'][1])
        _SetCFDCode(code)


    def onCodeChanged(self, currenBtnId):
        codeBG = self.findChild(QtGui.QButtonGroup, "CodeBG")
        if codeBG != None:
            if codeBG.selected() == self.SaturneRB:
                _SetCFDCode(CFD_Saturne)
                from cs_package import package
            if codeBG.selected() == self.NeptuneRB:
                _SetCFDCode(CFD_Neptune)
                from nc_package import package
            pkg = package()
            self.labelVersionValue.setText(pkg.version)
            self.labelPrefixValue.setText(pkg.dirs['prefix'][1])


#-----------------------------------------------------------------------------------------------------------

class SetTreeLocationDialog(QtGui.QDialog, Ui_SetTreeLocationDialog):
    """
    Tree Location Dialog informations about environment variables
    """
    def __init__(self, parent = None):
        """
        """
        QtGui.QDialog.__init__(self, parent)
        Ui_SetTreeLocationDialog.__init__(self)

        self.setupUi(self)


class SetTreeLocationDialogHandler(SetTreeLocationDialog):
    """
    Load the CFD study location. If the name of the CFD study
    does not exists, the corresponding folder is created.
    """
    def __init__(self, parent = None):
        """
        Constructor. Initialize text and label of the QDialog.
        """
        SetTreeLocationDialog.__init__(self, parent)
        #AB add case mode
        self.isCaseMode = False

        aBtn = self.findChild(QtGui.QPushButton,"OKButton")
        if not aBtn == None:
            aBtn.setText(self.tr("DLG_OK_BUTTON_TEXT"))

        aBtn = self.findChild(QtGui.QPushButton,"CancelButton")
        if not aBtn == None:
            aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        self.setWindowTitle(self.tr("LOCATION_DLG_CAPTION"))

        aLabel = self.findChild(QtGui.QLabel,"NameLabel")
        if not aLabel == None:
            aLabel.setText(self.tr("LOCATION_DLG_STUDY_NAME"))

        aLabel = self.findChild(QtGui.QLabel,"CaseLabel")
        if not aLabel == None:
            aLabel.setText(self.tr("LOCATION_DLG_CASE_NAME"))

        aGB = self.findChild(QtGui.QGroupBox, "CaseGroupBox")
        if not aGB == None:
            aGB.setTitle(self.tr("LOCATION_DLG_ADD_CASE"))

        aLabel = self.findChild(QtGui.QLabel,"StudyDirLabel")
        if not aLabel == None:
            aLabel.setText(self.tr("LOCATION_DLG_STUDY_DIR_LABEL"))

        aLE = self.findChild(QtGui.QLineEdit,"StudyDirName")
        if aLE != None:
            aLE.clear()

        # Installing validator on case name
        aLE = self.findChild(QtGui.QLineEdit,"CaseLineEdit")
        if aLE != None:
            aLE.clear()
        aLE = self.findChild(QtGui.QLineEdit,"StudyLineEdit")
        if aLE != None:
            aLE.clear()

        self.StudyPath = ''
        self.StudyName = ''
        self.adjustSize()


    def onBrowsePath(self):
        """
        Call into ui_SetTreeLocationDialog.py from setTreeLocationDialog.ui built with qtdesigner
        """
        aLE = self.findChild(QtGui.QLineEdit,"StudyDirName")
        if aLE != None:
            new_path = aLE.text()

            new_path = sgPyQt.getExistingDirectory(self, new_path, str(self.tr("SET_STUDY_LOCATION_BROWSE_CAPTION").toLatin1()))

            if not new_path or new_path == "":
                return
        new_path = os.path.abspath(str(new_path))
        if os.path.exists(os.path.join(new_path , 'MAILLAGE')) or os.path.exists(os.path.join(new_path, 'MESH')):
            new_path, self.StudyName = os.path.split(new_path)
        aLE.setText(new_path)
        self.findChild(QtGui.QLineEdit, "StudyLineEdit").setText(self.StudyName)


    def setCaseMode(self, flag):
        """
        modify the Dialog look:
        flag == False -> study dans cases creation
        flag == True -> only cases creation
        """
        self.isCaseMode = flag
        if self.isCaseMode == True:
            self.findChild(QtGui.QPushButton,"BrowseButton").hide()
            self.findChild(QtGui.QLineEdit,"StudyLineEdit").hide()
            self.findChild(QtGui.QLineEdit,"StudyDirName").hide()
            self.findChild(QtGui.QLabel,"NameLabel").hide()
            self.findChild(QtGui.QLabel,"StudyDirLabel").hide()
            self.findChild(QtGui.QGroupBox,"CaseGroupBox").setCheckable(False)
        else:
            self.findChild(QtGui.QPushButton,"BrowseButton").show()
            self.findChild(QtGui.QLineEdit,"StudyLineEdit").show()
            self.findChild(QtGui.QLineEdit,"StudyDirName").show()
            self.findChild(QtGui.QLabel,"NameLabel").show()
            self.findChild(QtGui.QLabel,"StudyDirLabel").show()

            self.findChild(QtGui.QGroupBox,"CaseGroupBox").setCheckable(True)
            self.findChild(QtGui.QGroupBox,"CaseGroupBox").setChecked(True)

        self.adjustSize();


    def accept(self):
        aDirLE = self.findChild(QtGui.QLineEdit,"StudyDirName")
        aNameLE = self.findChild(QtGui.QLineEdit,"StudyLineEdit")
        aCaseGB = self.findChild(QtGui.QGroupBox,"CaseGroupBox")
        aCaseLE = self.findChild(QtGui.QLineEdit,"CaseLineEdit")

        if aDirLE == None  or aNameLE == None \
               or aCaseGB == None or aCaseLE == None:
            raise DialogError, "Can't find control widget!"

        if self.isCaseMode == False:
            # check study directory
            aStudyDir = str(aDirLE.text().toLatin1())

            # create from study dir + study name

            aStudyDirName = str(aNameLE.text().toUpper().toLatin1())

            self.StudyPath = os.path.join(aStudyDir, aStudyDirName)
            self.StudyName = aStudyDirName

            studyObj = CFDSTUDYGUI_DataModel.FindStudyByPath(self.StudyPath)
            if studyObj != None:
                mess = self.tr("LOCATION_DLG_ERROR_OPEN_MESS")
                mess = mess.arg(self.StudyPath)
                QMessageBox.critical( self, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)
                return False

            if os.path.isdir(self.StudyPath) and self.StudyName != '':
                mess = self.tr("LOCATION_DLG_WARN_MESS")
                mess = mess.arg(self.StudyPath)
                if QMessageBox.information(self, "Information", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                    return False
            elif os.path.isdir(aStudyDir) and self.StudyName != '':
                mess = self.tr("LOCATION_DLG_WARN_MESS1")
                mess = mess.arg(aStudyDir)
                if QMessageBox.warning(self, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                    return False
            else:
                mess = self.tr("LOCATION_DLG_ERROR_MESS")
                mess = mess.arg(aStudyDir)
                QMessageBox.critical(self, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)
                return False

        # ckeck case name
        if aCaseGB.isChecked() or self.isCaseMode == True:
            if aCaseLE.text() == "":
                QMessageBox.critical(self, "Error", "Case name is empty!", QMessageBox.Ok, QMessageBox.NoButton)
                return False
            self.CaseNames = str(aCaseLE.text().toLatin1())
        else:
            self.CaseNames = ""

        SetTreeLocationDialog.accept(self)

#----------------------------------------------------------------------------------------------------------------------

class RunCaseDialog(QtGui.QDialog, Ui_RunCaseDialog):
    """
    Dialog informations about environment variables
    """
    def __init__(self, parent = None):
        """
        """
        QtGui.QDialog.__init__(self, parent)
        Ui_RunCaseDialog.__init__(self)

        self.setupUi(self)

class RunCaseDialogHandler(RunCaseDialog):
    """
    """
    def __init__(self, parent = None):
        """
        """
        RunCaseDialog.__init__(self,parent)
        self.aRunBtn = self.findChild(QtGui.QPushButton,"RunCaseBtn")
        self.aRunBtn.setText(self.tr("RUNCASE_DLG_RUN_BUTTON_TEXT"))

        aBtn = self.findChild(QtGui.QPushButton,"CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        aBtnGroup = self.findChild(QtGui.QGroupBox,"ModeButtonGroup")
        aBtnGroup.setTitle(self.tr("RUNCASE_DLG_MODE_TITLE"))

        aBtn = self.CompileModeBtn
        aBtn.setText(self.tr("RUNCASE_DLG_COMPILE_MODE_BTN_TEXT"))

        aBtn = self.RunModeBtn
        aBtn.setText(self.tr("RUNCASE_DLG_RUN_MODE_BTN_TEXT"))
        aBtn.setChecked(True)

        aGroupBox = self.findChild(QtGui.QGroupBox,"MainGroupBox")
        aGroupBox.setTitle(self.tr("RUNCASE_DLG_MAIN_GROUP_TITLE"))

        aLbl = self.findChild(QtGui.QLabel,"StudyLabel")
        aLbl.setText(self.tr("RUNCASE_DLG_STUDY_LABEL_TEXT"))

        aLbl = self.findChild(QtGui.QLabel,"CaseLabel")
        aLbl.setText(self.tr("RUNCASE_DLG_CASE_LABEL_TEXT"))

        self.setWindowTitle(self.tr("RUNCASE_DLG_CAPTION"))

        self.StudyCB = self.findChild(QtGui.QComboBox,"StudyCB")
        self.connect(self.StudyCB, SIGNAL("activated(int)"), self.slotUpdateStudy)

        self.CaseCB = self.findChild(QtGui.QComboBox,"CaseCB")
        self.connect(self.CaseCB, SIGNAL("activated(int)"), self.slotUpdateCase)

        self.adjustSize()

        self.InitCase = None


    def fillData(self):
        self.StudyCB.clear()
        aStudyList = CFDSTUDYGUI_DataModel.GetStudyNameList()
        if len(aStudyList) == 0:
            self.__SetInvalidMode(True)
            return
        self.StudyCB.clear()
        self.StudyCB.addItem(self.CurrentStudy.GetName())
        aCaseList = CFDSTUDYGUI_DataModel.GetCaseNameList(self.CurrentStudy)
        self.CaseCB.clear()

        for i in aCaseList:
            self.CaseCB.addItem(i)

    def slotUpdateStudy(self, newStudyIndex):

        aStudyList = CFDSTUDYGUI_DataModel.GetStudyList()

        if newStudyIndex > len(aStudyList)-1:
            return
        #obtain study object
        aStudyObj = aStudyList[newStudyIndex]
        #init case-combobox by available cases
        aCaseList = CFDSTUDYGUI_DataModel.GetCaseNameList(aStudyObj)

        self.CaseCB.clear()

        for item in aCaseList:
            self.CaseCB.addItem(item)


    def slotUpdateCase(self, newCaseIndex):
        """
        """
        aStudyList = CFDSTUDYGUI_DataModel.GetStudyList()
        aStudyIndex = self.StudyCB.currentIndex()
        if aStudyIndex > len(aStudyList)-1:
            print "Error: not correct index of study"
            return
        #obtain study object
        aStudyObj = aStudyList[aStudyIndex]

        # current case
        aCaseList =  CFDSTUDYGUI_DataModel.GetCaseList(aStudyObj)

        if newCaseIndex > len(aCaseList)-1:
            print "Error: not correct index of case"
            return
        #obtain study object
        aCaseObj = aCaseList[newCaseIndex]

        # object of DATA folder
        aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCaseObj, "DATA")
        if not len(aChildList) == 1:
            return
        aDataObj =  aChildList[0]

        aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, "THCH")
        if not len(aChildList) == 1:
            return
        aTHCHObj =  aChildList[0]

        #check for SaturneGUI or NeptuneGUI file in case:
        aRunBtn = self.findChild(QtGui.QRadioButton,"RunModeBtn")

        flag = CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCaseObj)
        aRunBtn.setEnabled(flag)
        if not flag:
            aCompBtn = self.findChild(QtGui.QRadioButton,"RunModeBtn")
            aCompBtn.setChecked(True)


    def __SetEnableAllUnderCase(self, param):
        self.aRunBtn.setEnabled(param)


    def show(self):
        self.fillData()
        RunCaseDialog.show(self)


    def setInitData(self, aCaseObj):
        self.CurrentCase = aCaseObj

    def setCurrentCase(self, aCase):
        self.CurrentCase = aCase

    # use before fill data
    def setCurrentStudy(self, theStudy):
        self.CurrentStudy = theStudy
#----------------------------------------------------------------------------------------------------------------------

class ECSConversionDialog(QtGui.QDialog,Ui_ECSConversionDialog):
    """
    Tree Location Dialog informations about environment variables
    """
    def __init__(self, parent = None):
        """
        """
        QtGui.QDialog.__init__(self, parent)
        Ui_ECSConversionDialog.__init__(self)

        self.setupUi(self)

class ECSConversionDialogHandler(ECSConversionDialog):
    """
    """
    def __init__(self, parent = None):
        """
        Constructor. Initialize text and label of the QDialog.
        """
        ECSConversionDialog.__init__(self, parent)

        self.ConvertBtn = self.findChild(QtGui.QPushButton,"ConvertBtn")
        self.ConvertBtn.setText(self.tr("ECSCONVERT_DLG_CONVERT_BUTTON"))

        aBtn = self.findChild(QtGui.QPushButton,"CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        aLabel = self.findChild(QtGui.QLabel,"CaseLabel")
        aLabel.hide()

        self.CaseCB = self.findChild(QtGui.QComboBox,"CaseCB")
        self.CaseCB.hide()

        self.ResultLabel = self.findChild(QtGui.QLabel,"ResultLabel")
        self.ResultLabel.setText(self.tr("ECSCONVERT_DLG_RESULT_LABEL_TEXT"))

        self.ResultName = self.findChild(QtGui.QLineEdit,"ResultNameLE")

        self.connect(self.ResultName, SIGNAL('textChanged(const QString&)'), self.slotResNameChanged)

        self.setWindowTitle(self.tr("ECSCONVERT_DLG_CAPTION"))

        self.adjustSize()


    def GetCaseName(self):
        return self.CaseCB.currentText()


    def show(self):
        aStudy = CFDSTUDYGUI_DataModel.GetFirstStudy()

        self.CaseCB.clear()

        aCaseList = CFDSTUDYGUI_DataModel.GetCaseNameList(aStudy)
        if len(aCaseList) == 0:
            self.CaseCB.setEnabled(False)
            self.ConvertBtn.setEnabled(False)
        else:
            self.CaseCB.setEnabled(True)
            self.ConvertBtn.setEnabled(True)
            for i in aCaseList:
                self.CaseCB.insertItem(i)

        ECSConversionDialog.show(self)


    def resultFileName(self):
        aResName = str(self.ResultName.text().toLatin1())
        if re.match(".*\.med$", aResName) and len(aResName) > 4:
            aResName = aResName[0:len(aResName)-4]

        return aResName


    def setResultFileName(self, NewName):
        if NewName == None:
            self.ResultName.clear()
        else:
            self.ResultName.setText(NewName)

        self.slotResNameChanged()


    def slotResNameChanged(self):
        self.ConvertBtn.setEnabled(str(self.ResultName.text().toLatin1())!= '')

#----------------------------------------------------------------------------------------------------------------------

class CopyDialog(QtGui.QDialog, Ui_CopyDialog):
    """
    Dialog informations about environment variables
    """
    def __init__(self, parent=None):
        """
        """
        QtGui.QDialog.__init__(self, parent)
        Ui_CopyDialog.__init__(self)

        self.setupUi(self)


class CopyDialogHandler(CopyDialog):
    """
    """
    def __init__(self, parent=None):
        CopyDialog.__init__(self, parent)
        self.CopyBtn = self.findChild(QtGui.QPushButton, "CopyBtn")
        self.CopyBtn.setText(self.tr("COPY_DLG_COPY_BUTTON"))

        aBtn = self.findChild(QtGui.QPushButton, "CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        aLabel = self.findChild(QtGui.QLabel, "SourceCaseLabel")
        aLabel.setText(self.tr("Case"))

        aLabel = self.findChild(QtGui.QLabel, "SourceFileLabel")
        aLabel.setText(self.tr("File"))

        aLabel = self.findChild(QtGui.QLabel, "DestCaseLabel")
        aLabel.setText(self.tr("DATA directory"))

        aLabel = self.findChild(QtGui.QLabel, "DestFilelabel")
        aLabel.setText(self.tr("New name"))

        self.SourceFileName = self.findChild(QtGui.QLabel,      "SourceFileName")
        self.SourceCaseName = self.findChild(QtGui.QLabel,      "SourceCaseName")
        self.DestDirLE      = self.findChild(QtGui.QLineEdit,   "DataDirectoryLineEdit")
        self.DestDirPB      = self.findChild(QtGui.QPushButton, "DataDirectoryPushButton")
        self.DestFileLE     = self.findChild(QtGui.QLineEdit,   "NewNameLineEdit")

        self.setWindowTitle(self.tr("COPY_DLG_CAPTION"))
        self.connect(self.DestDirPB, SIGNAL("clicked()"), self.onBrowsePath)


    def onBrowsePath(self):
        new_path = self.DestDirLE.text()
        new_path = sgPyQt.getExistingDirectory(self, new_path, str(self.tr("DATA directory").toLatin1()))

        if not new_path or new_path == "":
            return
        self.DestDirLE.setText(os.path.abspath(str(new_path)))



    def show(self):
        #aStudyList = CFDSTUDYGUI_DataModel.GetStudyList()
        #aCaseList  = []
        #for s in aStudyList:
            #aCaseList += CFDSTUDYGUI_DataModel.GetCaseNameList(s)

        #self.DestDirLE.clear()
        #self.DestFileLE.clear()
        #if self.CopyBtn.isEnabled():
            #if len(aCaseList) == 0:
                #self.DestDirLE.setEnabled(False)
                #self.DestFileLE.setEnabled(False)
                #self.CopyBtn.setEnabled(False)
            #else:
                #self.DestDirLE.setEnabled(True)
                #self.DestFileLE.setEnabled(True)
                #self.CopyBtn.setEnabled(True)

        CopyDialog.exec_(self)


    def setCurrentObject(self, sobj):
        self.Object = sobj
        aCase  = CFDSTUDYGUI_DataModel.GetCase(sobj)
        if not sobj or not aCase:
            CopyDialog.reject(self)
        else:
            c = aCase.GetName()
            p = CFDSTUDYGUI_DataModel._GetPath(aCase)
            self.SourceFileName.setText(sobj.GetName())
            self.SourceCaseName.setText(c)
            self.DestFileLE.setText(sobj.GetName())
            self.DestDirLE.setText(os.path.join(str(p), "DATA"))


    def accept(self):
        aDestDirName    = str(self.DestDirLE.text().toLatin1())
        aDestFileName   = str(self.DestFileLE.text().toLatin1())
        aDestFilePath = os.path.join(aDestDirName, aDestFileName)

        if os.path.exists(aDestFilePath) and os.path.isfile(aDestFilePath):
            QMessageBox.critical(self, self.tr("COPY_DLG_EXISTS_ERROR_CAPTION"), self.tr("COPY_DLG_EXISTS_ERROR_TEXT"), 1, 0)
            return False

        aSourceFilePath = CFDSTUDYGUI_DataModel._GetPath(self.Object)
        shutil.copyfile(aSourceFilePath, aDestFilePath)
        CopyDialog.accept(self)


    def destCaseName(self):
        return str(self.DestDirLE.text().toLatin1())


#class CopyDialog(QtGui.QDialog, Ui_CopyDialog):
    #"""
    #Dialog informations about environment variables
    #"""
    #def __init__(self, parent=None):
        #"""
        #"""
        #QtGui.QDialog.__init__(self, parent)
        #Ui_CopyDialog.__init__(self)

        #self.setupUi(self)


#class CopyDialogHandler(CopyDialog):
    #"""
    #"""
    #def __init__(self, parent = None):
        #CopyDialog.__init__(self, parent)
        #self.CopyBtn = self.findChild(QtGui.QPushButton, "CopyBtn")
        #self.CopyBtn.setText(self.tr("COPY_DLG_COPY_BUTTON"))

        #aBtn = self.findChild(QtGui.QPushButton, "CancelBtn")
        #aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        #aLabel = self.findChild(QtGui.QLabel, "DestCaseLabel")
        #aLabel.setText(self.tr("COPY_DLG_DEST_CASE_LABEL"))

        #aLabel = self.findChild(QtGui.QLabel, "SourceCaseLabel")
        #aLabel.setText(self.tr("COPY_DLG_SOURCE_CASE_LABEL"))

        #aLabel = self.findChild(QtGui.QLabel, "FileNameLabel")
        #aLabel.setText(self.tr("COPY_DLG_FILE_NAME_LABEL"))

        #self.FileName = self.findChild(QtGui.QLabel, "FileName")
        #self.SourceCaseName = self.findChild(QtGui.QLabel, "SourceCaseName")
        #self.DestDirLE = self.findChild(QtGui.QComboBox, "DestDirLE")

        #self.setWindowTitle(self.tr("COPY_DLG_CAPTION"))


    #def show(self):
        #aStudy = CFDSTUDYGUI_DataModel.GetFirstStudy()

        #self.DestDirLE.clear()
        #if self.CopyBtn.isEnabled():
            #aCaseList = CFDSTUDYGUI_DataModel.GetCaseNameList(aStudy)

            #if len(aCaseList) == 0:
                #self.DestDirLE.setEnabled(False)
                #self.CopyBtn.setEnabled(False)
            #else:
                #self.DestDirLE.setEnabled(True)
                #self.CopyBtn.setEnabled(True)

                #sourceCase = str(self.SourceCaseName.text().toLatin1())

                #for i in aCaseList:
                    #if not i == sourceCase:
                        #self.DestDirLE.addItem(i)

        #CopyDialog.exec_(self)


    #def setCurrentObject(self, sobj):
        #self.Object = sobj
        #aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
        #if not sobj == None and not aCase == None:
            #self.FileName.setText(sobj.GetName())
            #self.SourceCaseName.setText(aCase.GetName())

            #self.DestDirLE.setEnabled(True)
            #self.CopyBtn.setEnabled(True)

        #else:
            #self.FileName.setText("")
            #self.SourceCaseName.SetText("")
            #self.DestDirLE.setEnabled(False)
            #self.CopyBtn.setEnabled(False)


    #def accept(self):
        #aSourceFilePath = CFDSTUDYGUI_DataModel._GetPath(self.Object)
        #aSourceCase = CFDSTUDYGUI_DataModel.GetCase(self.Object)
        #aSourceCasePath = CFDSTUDYGUI_DataModel._GetPath(aSourceCase)

        #aSourceCaseName = str(self.SourceCaseName.text().toLatin1())
        #aDestCaseName = str(self.DestDirLE.currentText().toLatin1())

        ##check for existing of file in destinate CASE
        #aDestCasePath = str(QString(aSourceCasePath).left(len(aSourceCasePath)-len(aSourceCaseName)).toLatin1())
        #aDestCasePath += aDestCaseName

        #aDestFilePath = aDestCasePath + str(QString(aSourceFilePath).right(len(aSourceFilePath)-len(aSourceCasePath)).toLatin1())

        #if os.path.exists(aDestFilePath) and os.path.isfile(aDestFilePath):
            #QMessageBox.critical(self, self.tr("COPY_DLG_EXISTS_ERROR_CAPTION"), self.tr("COPY_DLG_EXISTS_ERROR_TEXT"), 1, 0)
            #return False

        #aCmd = "cp " + aSourceFilePath + " " + aDestFilePath

        #status = os.system(aCmd)
        #if not status == 0:
            #QMessageBox.critical(self, self.tr("COPY_DLG_COPY_ERROR_CAPTION"), self.tr("COPY_DLG_COPY_ERROR_TEXT"), 1, 0)
            #return False

        #CopyDialog.accept(self)


    #def destCaseName(self):
        #return str(self.DestDirLE.currentText().toLatin1())


#----------------------------------------------------------------------------------------------------------------------


class GUIActivationDialog(QtGui.QDialog,Ui_GUIActivationDialog):
    """
    Set environment variables Dialog informations
    """
    def __init__(self, parent = None):
        """
        """
        QtGui.QDialog.__init__(self)
        Ui_GUIActivationDialog.__init__(self)
        self.setupUi(self)


class GUIActivationDialogHandler(GUIActivationDialog):
    """
    """
    def __init__(self, parent = None):

        self.xmlfile = ""
        GUIActivationDialog.__init__(self, parent)

        self.ActivateBtn = self.findChild(QtGui.QPushButton,"ActivateBtn")
        self.ActivateBtn.setText(self.tr("GUIACTIVATE_DLG_ACTIVATE_BTN"))

        aBtn = self.findChild(QtGui.QPushButton,"CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        self.CaseCB = self.findChild(QtGui.QComboBox,"CaseCB")

        aGB = self.findChild(QtGui.QGroupBox,"OptionsGB")
        aGB.setTitle(self.tr("GUIACTIVATE_DLG_OPTIONS_TITLE"))

        self.FileCheckBox = self.findChild(QtGui.QCheckBox,"FileCheckBox")
        self.FileCB = self.findChild(QtGui.QComboBox,"FileCB")

        self.connect(self.FileCheckBox, SIGNAL("clicked()"), self.slotUseXMLfileChange)

        self.connect(self.CaseCB, SIGNAL("activated(int)"), self.slotUpdateData)

        self.adjustSize()


    def fillData(self, xmlFileName):
        """
        """
        self.disconnect(self.CaseCB, SIGNAL("activated(int)"), self.slotUpdateData)
        self.CaseCB.clear()
        self.xmlfile = xmlFileName

        for i in CFDSTUDYGUI_DataModel.GetCaseNameList(self.CurrentStudy):
            if self.CurrentCase and self.CurrentCase.GetName() == i:
                self.CaseCB.addItem(i)

        if self.xmlfile == "" :
            self.slotUpdateData()
        else:
            self.slotUpdateData()

        if CFD_Code() == CFD_Saturne:
            self.setWindowTitle(self.tr("ICSACTIVATE_DLG_CAPTION"))
            self.CaseLabel.setTitle(self.tr("ICSACTIVATE_DLG_CASE_LABEL"))
            #self.FileCheckBox.setChecked(False)
        elif CFD_Code() == CFD_Neptune:
            self.setWindowTitle(self.tr("IPBACTIVATE_DLG_CAPTION"))
            self.CaseLabel.setTitle(self.tr("IPBACTIVATE_DLG_CASE_LABEL"))


    def slotUpdateData(self):
        self.ActivateBtn.setEnabled(True)
        self.FileCB.clear()

        if self.CaseCB.currentText() == None:
            self.FileCheckBox.setEnabled(False)
            self.FileCB.setEnabled(False)
        else:
            self.FileCheckBox.setEnabled(True)
            self.FileCB.clear()

            # current case
            aCaseName = str(self.CaseCB.currentText().toLatin1())
            if aCaseName == "":
                self.ActivateBtn.setEnabled(False)
                return

            aCase = None
            aCases =  CFDSTUDYGUI_DataModel.GetCaseList(self.CurrentStudy)
            for c in aCases:
                if c.GetName() == aCaseName:
                    aCase = c
                    break

            if aCase == None:
                self.ActivateBtn.setEnabled(False)
                return

            # object of DATA folder
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "DATA")
            if not len(aChildList) == 1:
                self.ActivateBtn.setEnabled(False)
                return
            aDataObj =  aChildList[0]

            #fill File combo-box
            if self.xmlfile == "" :
                aFileList = CFDSTUDYGUI_DataModel.ScanChildNames(aDataObj, ".*\.xml$")
                if len(aFileList) == 0:
                    self.FileCheckBox.setEnabled(False);

                for i in aFileList:
                    self.FileCB.addItem(i)

                self.FileCB.setEnabled(self.FileCheckBox.isChecked())
            else :
                self.FileCB.addItem(self.xmlfile)
                self.FileCB.setEnabled(self.FileCheckBox.isChecked())

            #check for activation file SaturneGUI or NeptuneGUI
            if not CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCase):
                #Warning message
                if CFD_Code() == CFD_Saturne:
                    mess = self.tr("ICSACTIVATE_DLG_BAD_CASE_MESS")
                elif CFD_Code() == CFD_Neptune:
                    mess = self.tr("IPBACTIVATE_DLG_BAD_CASE_MESS")
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                self.ActivateBtn.setEnabled(False)


    # use before fill data
    def setCurrentCase(self, aCase):
        self.CurrentCase = aCase


    # use before fill data
    def setCurrentStudy(self, theStudy):
        self.CurrentStudy = theStudy


    # use after fill data
    def setCurrentXMLfile(self, aFile):
        if self.CurrentCase != None:
            self.FileCheckBox.setChecked(True)
            self.FileCB.setEnabled(True)
            #self.FileCB.setEditText(aFile)
            self.xmlfile = aFile


    def currentCaseName(self):
        if self.CaseCB.currentText() == 0:
            return None

        return str(self.CaseCB.currentText().toLatin1())


    def currentXMLfile(self):
        if not self.FileCB.isEnabled():
            return None
        return self.FileCB.currentText().toLatin1()


    def isUseXmlFile(self):
        return self.FileCheckBox.isChecked()


    def slotUseXMLfileChange(self):
        self.FileCB.setEnabled(self.FileCheckBox.isChecked())



class CFDSTUDYGUI_DialogCollector:
    def __init__(self):
        self.SetTreeLocationDialog = SetTreeLocationDialogHandler()
        self.InfoDialog = InfoDialogHandler()
        self.RunCaseDialog = RunCaseDialogHandler()
        self.ECSConversionDialog = ECSConversionDialogHandler()
        self.CopyDialog = CopyDialogHandler()
        self.GUIActivationDialog = GUIActivationDialogHandler()


#def CFD_Code():
#    return CFD_Saturne


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    lang = "en"
    qm = ["en","fr"]
    if len(sys.argv) == 2:
        lang = sys.argv[1]

    QM_FILE = "CFDSTUDY_msg_"+lang+".qm"
    fi = QFileInfo(QM_FILE)

    if not fi.exists():
        QMessageBox.warning(None, "File error",
            "Cannot find translation for language: "+lang)
    else:
        translator = QTranslator()
        translator.load(QM_FILE)
        app.installTranslator(translator)
    w = InfoDialogHandler()

    w.show()
    sys.exit(app.exec_())
