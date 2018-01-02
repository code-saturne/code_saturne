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
# ObjectTR is a convenient object for traduction purpose

ObjectTR = QObject()
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
import CFDSTUDYGUI_Commons
from CFDSTUDYGUI_Commons import _SetCFDCode, CFD_Code, sgPyQt
from CFDSTUDYGUI_Commons import CFD_Saturne, CFD_Neptune, CheckCFD_CodeEnv
from CFDSTUDYGUI_Message import cfdstudyMess
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
                mess = cfdstudyMess.trMessage(self.tr("CFDSTUDY_INVALID_ENV"),[]) + mess
                cfdstudyMess.criticalMessage(mess)
            else :
                InfoDialog.accept(self)
        else:
            mess = cfdstudyMess.trMessage(self.tr("INFO_DLG_INVALID_ENV"),[]) + mess
            cfdstudyMess.criticalMessage(mess)


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
        self.StudyPath = ''
        self.StudyName = ''
        self.CaseRefName = ''
        neptune_status = False
        cs_root_dir = os.getenv("CS_ROOT_DIR")
        if cs_root_dir:
            if "neptune_cfd" in os.listdir(os.path.join(cs_root_dir, "bin")):
                neptune_status = True
        if not neptune_status:
            self.findChild(QRadioButton,"radioButtonNeptune").setEnabled(False)

        self.findChild(QCheckBox,"checkBoxCreate").clicked.connect(self.slotCreateStudy)
        self.findChild(QCheckBox,"checkBoxLoad").clicked.connect(self.slotLoadStudy)
        self.findChild(QCheckBox,"checkBoxCopyFrom").clicked.connect(self.slotCopyFrom)
        self.findChild(QWidget,"copyFromCase_widget").hide()
        # Define option when openning
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
        self.findChild(QRadioButton,"radioButtonSaturne").setChecked(True)
        self.findChild(QRadioButton,"radioButtonNeptune").setChecked(False)
        self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
        self.findChild(QCheckBox,"checkBoxCouplingSaturneSyrthes").setChecked(False)
        self.findChild(QWidget,"SyrthesGroupBox").setVisible(False)
        self.findChild(QLineEdit, "NprocsLineEdit").setText("1")
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(False)
        self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(False)
        self.findChild(QLineEdit, "CaseLineEdit").setEnabled(False)
        self.CaseNames              = ""
        self.CreateOption           = False
        self.code                   = None
        self.CopyFromOption         = False
        self.CouplingSaturneSyrthes = False
        self.SyrthesCase            = ""
        self.adjustSize()

    def reinit(self) :
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(False)
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
        self.findChild(QRadioButton,"radioButtonSaturne").setChecked(True)
        self.findChild(QRadioButton,"radioButtonNeptune").setChecked(False)
        self.findChild(QWidget,"CaseWidget").setEnabled(False)
        self.findChild(QWidget,"CaseWidget").setVisible(True)
        self.findChild(QLineEdit,"StudyDirName").setText("")
        self.findChild(QLineEdit, "StudyLineEdit").setEnabled(False)
        self.findChild(QLineEdit, "StudyLineEdit").setText("")
        self.findChild(QLineEdit, "CaseLineEdit").setText("")
        self.findChild(QLineEdit, "CaseLineEdit").setEnabled(False)
        self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
        self.findChild(QCheckBox,"checkBoxCouplingSaturneSyrthes").setChecked(False)
        self.findChild(QLineEdit, "NprocsLineEdit").setText("1")
        self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(False)
        self.findChild(QWidget,"SyrthesGroupBox").setVisible(False)
        self.findChild(QWidget,"copyFromCase_widget").setEnabled(False)
        self.CaseNames              = ""
        self.CreateOption           = False
        self.code                   = None
        self.CopyFromOption         = False
        self.CouplingSaturneSyrthes = False
        self.StudyPath = ''
        self.StudyName = ''
        self.CaseRefName = ''
        self.adjustSize()


    def on_checkBoxCreate_pressed(self) :
        self.reinit()
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(False)
        self.findChild(QLineEdit,"StudyDirName").setText("")
        self.findChild(QLineEdit, "StudyLineEdit").setEnabled(True)
        self.findChild(QLineEdit, "StudyLineEdit").setText("")
        self.findChild(QLineEdit, "CaseLineEdit").setText("")
        self.findChild(QLineEdit, "copyFromCase").setText("")
        self.findChild(QWidget,"StudyGB").setEnabled(True)
        self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)


    def slotCreateStudy(self):
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(True)
        new_path = QFileDialog.getExistingDirectory(None, self.tr("Select a directory location to create a new CFD Study"))
        if str(new_path) == "" :
            self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
            self.reinit()
            return
        self.findChild(QLineEdit,"StudyDirName").setText(str(new_path))

        if self.findChild(QLineEdit,"StudyDirName").text() == "" :
            self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
            self.findChild(QLineEdit, "StudyLineEdit").setEnabled(False)

    def on_checkBoxLoad_pressed(self) :
        self.findChild(QCheckBox,"checkBoxCreate").setChecked(False)
        self.findChild(QWidget,"StudyGB").setEnabled(False)
        self.findChild(QLineEdit,"StudyDirName").setText("")
        self.findChild(QLineEdit, "StudyLineEdit").setText("")
        self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
        self.findChild(QCheckBox,"checkBoxCouplingSaturneSyrthes").setChecked(False)
        self.findChild(QLineEdit, "copyFromCase").setText("")

    def slotLoadStudy(self):
        self.findChild(QCheckBox,"checkBoxLoad").setChecked(True)
        self.adjustSize()
        new_path = QFileDialog.getExistingDirectory(None, self.tr("Select an existing CFD Study"))
        if str(new_path) == "" :
            self.reinit()
            return
        new_path, self.StudyName = os.path.split(str(new_path))
        self.findChild(QLineEdit,"StudyDirName").setText(new_path)
        self.findChild(QLineEdit, "StudyLineEdit").setText(self.StudyName)
        self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(True)


    def on_StudyLineEdit_textChanged(self,stringValue):
        if str(stringValue) != "" :
            self.findChild(QLineEdit, "CaseLineEdit").setEnabled(True)
            self.findChild(QWidget,"CaseWidget").setEnabled(True)
            self.findChild(QCheckBox,"checkBoxCopyFrom").setEnabled(True)
            self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(True)
        else:
            self.findChild(QLineEdit, "CaseLineEdit").setEnabled(False)
            self.findChild(QWidget,"CaseWidget").setEnabled(False)
            self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(False)


    def on_checkBoxCopyFrom_pressed(self):
        if self.findChild(QCheckBox,"checkBoxCopyFrom").isChecked() and self.findChild(QLineEdit, "copyFromCase").text() != "":
            self.findChild(QLineEdit,"copyFromCase").setText("")
            return

#####################################################################################
    def slotCopyFrom(self):
        """
        Call into ui_SetTreeLocationDialog.py from setTreeLocationDialog.ui built with qtdesigner
        for option --copy-from
        """

        if not self.findChild(QCheckBox,"checkBoxCopyFrom").isChecked() :
            self.findChild(QWidget,"copyFromCase_widget").hide()
            return

        CopyFromCasePath = ""
        if self.findChild(QLineEdit,"StudyDirName").text() != "":
            CopyFromCasePath = QFileDialog.getExistingDirectory(None, ObjectTR.tr("SET_CASE_LOCATION_BROWSE_CAPTION"))

            if CopyFromCasePath == None or str(CopyFromCasePath) == "":
                self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
                return

        self.CaseRefName = os.path.abspath(str(CopyFromCasePath))
        # check if it is a case directory
        if not self.isCfdCaseDir(self.CaseRefName):
            mess = cfdstudyMess.trMessage(self.tr("CASE_DLG_ERROR_MESS"),[self.CaseRefName])
            cfdstudyMess.aboutMessage(mess)
            self.findChild(QWidget,"copyFromCase_widget").hide()
            self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
            return
        CaseRefDATAPath = os.path.join(self.CaseRefName,"DATA")
        if "NeptuneGUI" in os.listdir(CaseRefDATAPath) and self.findChild(QRadioButton,"radioButtonSaturne").isChecked():
            mess = cfdstudyMess.trMessage(self.tr("CASE_COPYFROM_NEPTUNE_DLG_ERROR_MESS"),[])
            cfdstudyMess.aboutMessage(mess)
            self.findChild(QWidget,"copyFromCase_widget").hide()
            self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
            return
        if "SaturneGUI" in os.listdir(CaseRefDATAPath) and self.findChild(QRadioButton,"radioButtonNeptune").isChecked():
            mess = cfdstudyMess.trMessage(self.tr("CASE_COPYFROM_SATURNE_DLG_ERROR_MESS"),[])
            cfdstudyMess.aboutMessage(mess)
            self.findChild(QWidget,"copyFromCase_widget").hide()
            self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
            return

        self.findChild(QWidget,"copyFromCase_widget").show()
        self.findChild(QLineEdit,"copyFromCase").setText(self.CaseRefName)
        self.findChild(QDialogButtonBox,"buttonBox").button(QDialogButtonBox.Ok).setEnabled(True)


    def on_checkBoxCouplingSaturneSyrthes_clicked(self):
        if self.findChild(QCheckBox,"checkBoxCouplingSaturneSyrthes").isChecked() and self.findChild(QLineEdit, "syrthesCase").text() != "":
            self.findChild(QLineEdit,"syrthesCase").setText("")


    def isCfdCaseDir(self,CfdCaseRefDir) :
        listDir_required = ["DATA","SRC"]
        boolDir = True
        for i in listDir_required :
            if i not in os.listdir(CfdCaseRefDir) :
                boolDir = False
                return boolDir
        return boolDir

#####################################################################################
    def setCaseMode(self):
        """
        Called for adding cases
        """
        self.setWindowTitle(self.tr("Add new case to study"))
        self.findChild(QCheckBox,"checkBoxLoad").hide()
        self.findChild(QCheckBox,"checkBoxCreate").hide()
        self.findChild(QWidget,"CaseWidget").setEnabled(True)
        self.findChild(QCheckBox,"checkBoxCopyFrom").setEnabled(False)
        self.findChild(QCheckBox,"checkBoxCopyFrom").setChecked(False)
        self.findChild(QGroupBox,"groupBox").show()
        self.findChild(QGroupBox,"groupBox").setEnabled(True)
        self.findChild(QWidget,"StudyGB").setEnabled(False)
        self.findChild(QLineEdit,"CaseLineEdit").setEnabled(True)
        self.findChild(QRadioButton,"radioButtonSaturne").setEnabled(True)
        self.findChild(QRadioButton,"radioButtonNeptune").setEnabled(True)
        self.findChild(QCheckBox,"checkBoxCouplingSaturneSyrthes").hide()
        self.findChild(QLineEdit,"StudyDirName").show()
        self.findChild(QLabel,"StudyDirLabel").show()
        self.adjustSize()



    def accept(self):
        aDirLE                      = self.findChild(QLineEdit,"StudyDirName")
        aNameLE                     = self.findChild(QLineEdit,"StudyLineEdit")
        aCaseLE                     = self.findChild(QLineEdit,"CaseLineEdit")
        CreateOption                = self.findChild(QCheckBox,"checkBoxCreate")
        Neptune                     = self.findChild(QRadioButton,"radioButtonNeptune")
        Saturne                     = self.findChild(QRadioButton,"radioButtonSaturne")
        self.CaseNames              = str(self.findChild(QLineEdit,"CaseLineEdit").text())
        self.CopyFromOption         = self.findChild(QCheckBox, "checkBoxCopyFrom").isChecked()
        self.CouplingSaturneSyrthes = self.findChild(QCheckBox, "checkBoxCouplingSaturneSyrthes").isChecked()
        self.Nprocs                 = str(self.findChild(QLineEdit,"NprocsLineEdit").text())
        if  aNameLE.text() == "" :
            mess = cfdstudyMess.trMessage(self.tr("LOCATION_DLG_ERROR_MESS"),[])
            cfdstudyMess.criticalMessage(mess)
            return False

        # check study directory
        aStudyDir = str(aDirLE.text())
        # Load from study dir + study name
        if aNameLE.text() != aNameLE.text():
            raise DialogError, "Names must not contain special characters."

        aStudyDirName = str(aNameLE.text())
        self.StudyPath = os.path.join(aStudyDir, aStudyDirName)
        self.StudyName = aStudyDirName

        if Neptune.isChecked():
            self.code = "NEPTUNE_CFD"
        else:
            self.code = "Code_Saturne"

        if self.checkBoxLoad.isChecked():
            studyObj = CFDSTUDYGUI_DataModel.FindStudyByPath(self.StudyPath)
            if studyObj != None:
                mess = cfdstudyMess.trMessage(self.tr("LOCATION_DLG_ERROR_OPEN_MESS"),[self.StudyPath])
                cfdstudyMess.aboutMessage(mess)
                self.reinit()
                return False
            if self.StudyName == '':
                mess = cfdstudyMess.trMessage(self.tr("LOCATION_DLG_ERROR_MESS"),[])
                cfdstudyMess.criticalMessage(mess)
                self.reinit()
                return False
            if not CFDSTUDYGUI_Commons.isaSaturneSyrthesCouplingStudy(self.StudyPath):
                if not CFDSTUDYGUI_Commons.isaCFDStudy(self.StudyPath):
                    mess = cfdstudyMess.trMessage(self.tr("NOT_A_STUDY_DIRECTORY"),[self.StudyPath,"CFD","SYRTHES"])
                    cfdstudyMess.criticalMessage(mess)
                    self.reinit()
                    return False

        # ckeck case name
        if self.checkBoxCreate.isChecked() :
            if self.StudyName == '':
                mess = cfdstudyMess.trMessage(self.tr("LOCATION_DLG_ERROR_MESS"),[aStudyDir])
                cfdstudyMess.criticalMessage(mess)
                return False
            self.CaseNames = str(aCaseLE.text())
            self.CreateOption = True
            if self.CouplingSaturneSyrthes :
                self.SyrthesCase = str(self.findChild(QLineEdit,"syrthesCase").text())
                if self.SyrthesCase == "":
                    mess = cfdstudyMess.trMessage(self.tr("EMPTY_SYRTHES_CASENAME_MESS"),[])
                    cfdstudyMess.criticalMessage(mess)
                    return False
        SetTreeLocationDialog.accept(self)

    def reject(self):
        self.reinit()
        SetTreeLocationDialog.reject(self)

#----------------------------------------------------------------------------------------------

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
        new_path = QFileDialog.getExistingDirectory(None, self.tr("DATA directory"),new_path)
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
            mess = cfdstudyMess.trMessage(self.tr("COPY_DLG_EXISTS_ERROR_TEXT"),[])
            cfdstudyMess.criticalMessage(mess)
            return False

        aSourceFilePath = CFDSTUDYGUI_DataModel._GetPath(self.Object)
        shutil.copyfile(aSourceFilePath, aDestFilePath)

        CopyDialog.accept(self)


    def destCaseName(self):
        return str(self.findChild(QLineEdit,   "DataDirectoryLineEdit").text())


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
                    mess = cfdstudyMess.trMessage(self.tr("ICSACTIVATE_DLG_BAD_CASE_MESS"),[])
                elif CFD_Code() == CFD_Neptune:
                    mess = cfdstudyMess.trMessage(self.tr("IPBACTIVATE_DLG_BAD_CASE_MESS"),[])
                cfdstudyMess.warningMessage(mess)
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
    import sys
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
#    w = InfoDialogHandler()
    w = SetTreeLocationDialogHandler()
    w.show()
    sys.exit(app.exec_())
