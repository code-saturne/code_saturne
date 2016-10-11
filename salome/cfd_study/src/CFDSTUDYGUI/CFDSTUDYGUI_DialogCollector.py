# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from ui_InfoDialog            import Ui_InfoDialog
from ui_SetTreeLocationDialog import Ui_SetTreeLocationDialog
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
log.setLevel(logging.NOTSET)
#-------------------------------------------------------------------------------
# Dialog definitions
#-------------------------------------------------------------------------------

class InfoDialog(QDialog, Ui_InfoDialog):
    """
    Dialog informations about solver installation.
    """
    def __init__(self, parent = None):
        """
        """
        QDialog.__init__(self, parent)
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

        aBtn = self.findChild(QPushButton, "OKButton")
        aBtn.setText(self.tr("DLG_OK_BUTTON_TEXT"))

        self.setWindowTitle(self.tr("INFO_DLG_CAPTION"))


    def accept(self):
        iok, mess = CheckCFD_CodeEnv(CFD_Code())
        if iok:
            if mess != "" :
                Error = "Error : "+ self.tr("CFDSTUDY_INVALID_ENV")
                QMessageBox.critical(ActionHandler.dskAgent().workspace(),
                                 Error, mess, QMessageBox.Ok, 0)
            else :
                InfoDialog.accept(self)
        else:
            Error = "Error : " + self.tr("INFO_DLG_INVALID_ENV")
            QMessageBox.critical(self, Error, mess, QMessageBox.Ok, 0)


    def setCode(self, env_saturne, env_neptune):
        if env_neptune:
            code = CFD_Neptune
            from nc_package import package

        elif env_saturne:
            code = CFD_Saturne
            from cs_package import package

        else:
            raise DialogError, "Invalid CFD_Code in InfoDialog class"

        pkg = package()
        self.labelVersionValue.setText(pkg.version)
        self.labelPrefixValue.setText(pkg.get_dir('prefix'))
        self.labelCodeValue.setText(pkg.name)
        _SetCFDCode(code)


    def update(self, code):
        if code == CFD_Saturne:
            from cs_package import package
        if code == CFD_Neptune:
            from nc_package import package

        pkg = package()
        self.labelVersionValue.setText(pkg.version)
        self.labelPrefixValue.setText(pkg.get_dir('prefix'))
        self.labelCodeValue.setText(pkg.name)



#-----------------------------------------------------------------------------------------------------------

class SetTreeLocationDialog(QDialog, Ui_SetTreeLocationDialog):
    """
    Tree Location Dialog informations about environment variables
    """
    def __init__(self, parent = None):
        """
        """
        QDialog.__init__(self, parent)
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
        self.isCaseMode = False

        aBtn = self.findChild(QPushButton,"OKButton")
        if not aBtn == None:
            aBtn.setText(self.tr("DLG_OK_BUTTON_TEXT"))

        aBtn = self.findChild(QPushButton,"CancelButton")
        if not aBtn == None:
            aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        self.setWindowTitle(self.tr("LOCATION_DLG_CAPTION"))

        aLabel = self.findChild(QLabel,"CaseLabel")
        if not aLabel == None:
            aLabel.setText(self.tr("LOCATION_DLG_CASE_NAME"))

        aLabel = self.findChild(QLabel,"StudyDirLabel")
        if not aLabel == None:
            aLabel.setText(self.tr("LOCATION_DLG_STUDY_DIR_LABEL"))

        aLE = self.findChild(QLineEdit,"StudyDirName")
        if aLE != None:
            aLE.clear()

        # Installing validator on case name
        aLE = self.findChild(QLineEdit,"CaseLineEdit")
        if aLE != None:
            aLE.clear()
        aLE = self.findChild(QLineEdit,"StudyLineEdit")
        if aLE != None:
            aLE.clear()

        self.StudyPath = ''
        self.StudyName = ''
        self.CaseRefName = ''
        neptune_status, mess2 = CheckCFD_CodeEnv(CFD_Neptune)
        if not neptune_status:
            self.findChild(QRadioButton,"radioButtonNeptune").setEnabled(False)

        self.Create = self.findChild(QCheckBox,"checkBoxCreate")
        self.Create.clicked.connect(self.slotCreateCase)
        self.Load = self.findChild(QCheckBox,"checkBoxLoad")
        self.Load.clicked.connect(self.slotLoadCase)
        self.CopyFrom = self.findChild(QCheckBox,"checkBoxCopyFrom")
        self.CopyFrom.clicked.connect(self.slotCopyFrom)

        # Define option when openning
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(True)
        self.findChild(QRadioButton,"radioButtonSaturne").setChecked(True)
        self.findChild(QRadioButton,"radioButtonNeptune").setChecked(False)
        self.findChild(QGroupBox,"CaseGroupBox").setEnabled(False)
        self.findChild(QCheckBox, "checkBoxMesh").setChecked(True)
        self.findChild(QCheckBox, "checkBoxPOST").setChecked(True)
        self.findChild(QPushButton, "BrowseButtonCopy").setEnabled(False)

        self.setCaseMode(self.isCaseMode)


    def slotCreateCase(self):
        self.Load.setChecked(False)
        self.findChild(QGroupBox,"CaseGroupBox").setEnabled(True)
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(True)
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(False)


    def slotLoadCase(self):
        self.Create.setChecked(False)
        self.findChild(QGroupBox,"CaseGroupBox").setEnabled(False)
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(True)


    def slotCopyFrom(self):
        if self.CopyFrom.isChecked():
            self.findChild(QPushButton, "BrowseButtonCopy").setEnabled(True)
        else:
            self.findChild(QPushButton, "BrowseButtonCopy").setEnabled(False)


    def onBrowsePath(self):
        """
        Call into ui_SetTreeLocationDialog.py from setTreeLocationDialog.ui built with qtdesigner
        """
        aLE = self.findChild(QLineEdit,"StudyDirName")
        if aLE != None:
            new_path = aLE.text()

            new_path = sgPyQt.getExistingDirectory(self, new_path, str(self.tr("SET_STUDY_LOCATION_BROWSE_CAPTION")))

            if not new_path or new_path == "":
                return
        new_path = os.path.abspath(str(new_path))
        if os.path.exists(os.path.join(new_path, 'MESH')):
            new_path, self.StudyName = os.path.split(new_path)
        aLE.setText(new_path)
        self.findChild(QLineEdit, "StudyLineEdit").setText(self.StudyName)


    def onBrowsePathCopy(self):
        """
        Call into ui_SetTreeLocationDialog.py from setTreeLocationDialog.ui built with qtdesigner
        for option --copy-from
        """
        aLE = self.findChild(QLineEdit,"StudyDirName")
        if aLE != None:
            new_path = aLE.text()

            new_path = sgPyQt.getExistingDirectory(self, new_path, str(self.tr("SET_CASE_LOCATION_BROWSE_CAPTION")))

            if not new_path or new_path == "":
                return
        # TODO check if it is a case directory
        self.CaseRefName = os.path.abspath(str(new_path))


    def setCaseMode(self, flag):
        """
        modify the Dialog look:
        flag == False -> study dans cases creation
        flag == True -> only cases creation
        """
        self.isCaseMode = flag
        if self.isCaseMode == True:
            self.findChild(QCheckBox,"checkBoxCreate").setChecked(True)
            self.findChild(QCheckBox,"checkBoxCreate").setCheckable(False)
            self.findChild(QCheckBox,"checkBoxLoad").setChecked(False)
            self.findChild(QCheckBox,"checkBoxLoad").setCheckable(False)
            self.findChild(QGroupBox,"CaseGroupBox").setEnabled(True)

            self.findChild(QGroupBox,"groupBox").hide()
            self.findChild(QPushButton,"BrowseButton").hide()
            self.findChild(QLineEdit,"StudyDirName").hide()
            self.findChild(QLabel,"StudyDirLabel").hide()
        else:
            self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
            self.findChild(QCheckBox,"checkBoxCreate").setCheckable(True)
            self.findChild(QCheckBox,"checkBoxLoad").setChecked(True)
            self.findChild(QCheckBox,"checkBoxLoad").setCheckable(True)

            self.findChild(QGroupBox,"groupBox").show()
            self.findChild(QPushButton,"BrowseButton").show()
            self.findChild(QLineEdit,"StudyDirName").show()
            self.findChild(QLabel,"StudyDirLabel").show()

        self.adjustSize();


    def accept(self):
        aDirLE       = self.findChild(QLineEdit,"StudyDirName")
        aNameLE      = self.findChild(QLineEdit,"StudyLineEdit")
        aCaseLE      = self.findChild(QLineEdit,"CaseLineEdit")
        CreateOption = self.findChild(QCheckBox,"checkBoxCreate")
        Neptune      = self.findChild(QRadioButton,"radioButtonNeptune")
        Saturne      = self.findChild(QRadioButton,"radioButtonSaturne")
        MeshOpt      = self.findChild(QCheckBox, "checkBoxMesh")
        PostOpt      = self.findChild(QCheckBox, "checkBoxPOST")
        CopyFromOpt  = self.findChild(QCheckBox, "checkBoxCopyFrom")

        if aDirLE == None or aNameLE == None:
            raise DialogError, "Can't find control widget!"

        if self.isCaseMode == False:
            # check study directory
            aStudyDir = str(aDirLE.text())

            # create from study dir + study name
            if aNameLE.text() != aNameLE.text():
                raise DialogError, "Names must not contain special characters."

            aStudyDirName = str(aNameLE.text())

            self.StudyPath = os.path.join(aStudyDir, aStudyDirName)
            self.StudyName = aStudyDirName

            studyObj = CFDSTUDYGUI_DataModel.FindStudyByPath(self.StudyPath)
            if studyObj != None:
                mess = str(self.tr("LOCATION_DLG_ERROR_OPEN_MESS"))
                mess = mess%(self.StudyPath)
                QMessageBox.critical( self, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)
                return False

            if os.path.isdir(self.StudyPath) and self.StudyName != '':
                mess = str(self.tr("LOCATION_DLG_WARN_MESS"))
                mess = mess%(self.StudyPath)
                if QMessageBox.information(self, "Information", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                    return False
            elif os.path.isdir(aStudyDir) and self.StudyName != '':
                mess = str(self.tr("LOCATION_DLG_WARN_MESS1"))
                mess = mess%(aStudyDir)
                if QMessageBox.warning(self, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                    return False
            else:
                mess = str(self.tr("LOCATION_DLG_ERROR_MESS"))
                mess = mess%(aStudyDir)
                QMessageBox.critical(self, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)
                return False

        # ckeck case name
        if CreateOption.isChecked() or self.isCaseMode == True:
            if aCaseLE.text() == "":
                QMessageBox.critical(self, "Error", "Case name is empty!", QMessageBox.Ok, QMessageBox.NoButton)
                return False
            self.CaseNames = str(aCaseLE.text())
            self.CreateOption = True
            if Neptune.isChecked():
                self.code = CFD_Neptune
            else:
                self.code = CFD_Saturne
            self.meshOption = MeshOpt.isChecked()
            self.postOption = PostOpt.isChecked()
            self.CopyFromOption = CopyFromOpt.isChecked()
        else:
            self.CaseNames = ""
            self.CreateOption = False
            self.code = None
            self.meshOption = None
            self.postOption = None
            self.CopyFromOption = None

        SetTreeLocationDialog.accept(self)


#----------------------------------------------------------------------------------------------------------------------

class ECSConversionDialog(QDialog,Ui_ECSConversionDialog):
    """
    Tree Location Dialog informations about environment variables
    """
    def __init__(self, parent = None):
        """
        """
        QDialog.__init__(self, parent)
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

        self.ConvertBtn = self.findChild(QPushButton,"ConvertBtn")
        self.ConvertBtn.setText(self.tr("ECSCONVERT_DLG_CONVERT_BUTTON"))

        aBtn = self.findChild(QPushButton,"CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        aLabel = self.findChild(QLabel,"CaseLabel")
        aLabel.hide()

        self.CaseCB = self.findChild(QComboBox,"CaseCB")
        self.CaseCB.hide()

        self.ResultLabel = self.findChild(QLabel,"ResultLabel")
        self.ResultLabel.setText(self.tr("ECSCONVERT_DLG_RESULT_LABEL_TEXT"))

        self.ResultName = self.findChild(QLineEdit,"ResultNameLE")
        self.ResultName.textChanged.connect(self.slotResNameChanged)
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
        aResName = str(self.ResultName.text())
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
        self.ConvertBtn.setEnabled(str(self.ResultName.text())!= '')

#----------------------------------------------------------------------------------------------------------------------

class CopyDialog(QDialog, Ui_CopyDialog):
    """
    Dialog informations about environment variables
    """
    def __init__(self, parent=None):
        """
        """
        QDialog.__init__(self, parent)
        Ui_CopyDialog.__init__(self)

        self.setupUi(self)


class CopyDialogHandler(CopyDialog):
    """
    """
    def __init__(self, parent=None):
        CopyDialog.__init__(self, parent)
        self.CopyBtn = self.findChild(QPushButton, "CopyBtn")
        self.CopyBtn.setText(self.tr("COPY_DLG_COPY_BUTTON"))

        aBtn = self.findChild(QPushButton, "CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        aLabel = self.findChild(QLabel, "SourceCaseLabel")
        aLabel.setText(self.tr("Case"))

        aLabel = self.findChild(QLabel, "SourceFileLabel")
        aLabel.setText(self.tr("File"))

        aLabel = self.findChild(QLabel, "DestCaseLabel")
        aLabel.setText(self.tr("DATA directory"))

        aLabel = self.findChild(QLabel, "DestFilelabel")
        aLabel.setText(self.tr("New name"))

        self.SourceFileName = self.findChild(QLabel,      "SourceFileName")
        self.SourceCaseName = self.findChild(QLabel,      "SourceCaseName")
        self.DestDirLE      = self.findChild(QLineEdit,   "DataDirectoryLineEdit")
        self.DestDirPB      = self.findChild(QPushButton, "DataDirectoryPushButton")
        self.DestFileLE     = self.findChild(QLineEdit,   "NewNameLineEdit")

        self.setWindowTitle(self.tr("COPY_DLG_CAPTION"))
        self.DestDirPB.clicked.connect(self.onBrowsePath)

    def onBrowsePath(self):
        new_path = self.DestDirLE.text()

        new_path = sgPyQt.getExistingDirectory(self, new_path, str(self.tr("DATA directory")))

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
        aDestDirName    = str(self.DestDirLE.text())
        aDestFileName   = str(self.DestFileLE.text())
        aDestFilePath = os.path.join(aDestDirName, aDestFileName)

        if os.path.exists(aDestFilePath) and os.path.isfile(aDestFilePath):
            QMessageBox.critical(self, self.tr("COPY_DLG_EXISTS_ERROR_CAPTION"), self.tr("COPY_DLG_EXISTS_ERROR_TEXT"), 1, 0)
            return False

        aSourceFilePath = CFDSTUDYGUI_DataModel._GetPath(self.Object)
        shutil.copyfile(aSourceFilePath, aDestFilePath)
        CopyDialog.accept(self)


    def destCaseName(self):
        return str(self.DestDirLE.text())


#----------------------------------------------------------------------------------------------------------------------


class GUIActivationDialog(QDialog,Ui_GUIActivationDialog):
    """
    Set environment variables Dialog informations
    """
    def __init__(self, parent = None):
        """
        """
        QDialog.__init__(self)
        Ui_GUIActivationDialog.__init__(self)
        self.setupUi(self)


class GUIActivationDialogHandler(GUIActivationDialog):
    """
    """
    def __init__(self, parent = None):
        log.debug("__init__")
        self.xmlfile = ""
        GUIActivationDialog.__init__(self, parent)

        self.ActivateBtn = self.findChild(QPushButton,"ActivateBtn")
        self.ActivateBtn.setText(self.tr("GUIACTIVATE_DLG_ACTIVATE_BTN"))

        aBtn = self.findChild(QPushButton,"CancelBtn")
        aBtn.setText(self.tr("DLG_CANCEL_BUTTON_TEXT"))

        self.CaseCB = self.findChild(QComboBox,"CaseCB")

        aGB = self.findChild(QGroupBox,"OptionsGB")
        aGB.setTitle(self.tr("GUIACTIVATE_DLG_OPTIONS_TITLE"))

        self.FileCheckBox = self.findChild(QCheckBox,"FileCheckBox")
        self.FileCB = self.findChild(QComboBox,"FileCB")

        self.FileCheckBox.clicked.connect(self.slotUseXMLfileChange)
        self.CaseCB.activated["int"].connect(self.slotUpdateData)

        self.adjustSize()


    def fillData(self, xmlFileName):
        """
        """
        log.debug("fillData")
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
        log.debug("slotUpdateData")
        self.ActivateBtn.setEnabled(True)
        self.FileCB.clear()

        if self.CaseCB.currentText() == None:
            self.FileCheckBox.setEnabled(False)
            self.FileCB.setEnabled(False)
        else:
            self.FileCheckBox.setEnabled(True)
            self.FileCB.clear()

            # current case
            aCaseName = str(self.CaseCB.currentText())
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

        return str(self.CaseCB.currentText())


    def currentXMLfile(self):
        if not self.FileCB.isEnabled():
            return None
        return self.FileCB.currentText()


    def isUseXmlFile(self):
        return self.FileCheckBox.isChecked()


    def slotUseXMLfileChange(self):
        self.FileCB.setEnabled(self.FileCheckBox.isChecked())



class CFDSTUDYGUI_DialogCollector:
    def __init__(self):
        self.SetTreeLocationDialog = SetTreeLocationDialogHandler()
        self.InfoDialog = InfoDialogHandler()
        self.ECSConversionDialog = ECSConversionDialogHandler()
        self.CopyDialog = CopyDialogHandler()
        self.GUIActivationDialog = GUIActivationDialogHandler()

if __name__ == "__main__":
    app = QApplication(sys.argv)
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
