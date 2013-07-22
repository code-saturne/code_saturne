# -*- coding: iso-8859-1 -*-
#
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
This module contains the following classes and function:
- BatchRunningUserFilesDialogView
- BatchRunningAdvancedOptionsDialogView
- BatchRunningView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys
import string, types
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from BatchRunningForm import Ui_BatchRunningForm
from BatchRunningUserFilesDialogForm import Ui_BatchRunningUserFilesDialogForm
from BatchRunningAdvancedOptionsDialogForm import Ui_BatchRunningAdvancedOptionsDialogForm

from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, IntValidator, RegExpValidator, setGreenColor
from BatchRunningModel import RuncaseModel, BatchRunningModel
from Pages.LocalizationModel import LocalizationModel, Zone

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BatchRunningView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Popup window class: Data and results user files
#-------------------------------------------------------------------------------


class BatchRunningUserFilesDialogView(QDialog, Ui_BatchRunningUserFilesDialogForm):
    """
    Class for data and results user files
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningUserFilesDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("User files"))

        self.default = default
        self.result  = {}
        for key in self.default.keys():
            self.result[key] = []

        # Models
        rows = 0
        columns = 1
        self.modelData = QListModel(rows, columns)
        self.modelResu = QListModel(rows, columns)

        # associated with views.
        self.viewData.setModel(self.modelData)
        self.viewResu.setModel(self.modelResu)

        self.viewData.setItemDelegate(DataDelegate(self, self.default['data_path']))
        self.viewResu.setItemDelegate(ResuDelegate(self))

        # Connections
        self.connect(self.buttonNewData, SIGNAL("clicked()"), self.slotAddData)
        self.connect(self.buttonAddData, SIGNAL("clicked()"), self.slotNewData)
        self.connect(self.buttonAddResu, SIGNAL("clicked()"), self.slotNewResu)
        self.connect(self.buttonDeleteData, SIGNAL("clicked()"), self.slotDeleteData)
        self.connect(self.buttonDeleteResu, SIGNAL("clicked()"), self.slotDeleteResu)

        # Previous values
        for item in self.default['data']:
            self.setFileData(item)
        for item in self.default['results']:
            self.setFileResu(item)


    def setFileData(self, item):
        # Verify that the input is not already in the QListView
        indexList = self.modelData.search(QString(item))

        if indexList:
            title = self.tr("Warning")
            msg   = self.tr("%s is already in the list." % str(item))
            QMessageBox.warning(self, title, msg)
        else:
            std_item = QStandardItem(QString(item))
            self.modelData.appendRow(std_item)


    def setFileResu(self, item):
        # Verify that the input is not already in the QListView
        indexList = self.modelResu.search(QString(item))

        if indexList:
            title = self.tr("Warning")
            msg   = self.tr("%s is already in the list." % str(item))
            QMessageBox.warning(self, title, msg)
        else:
            std_item = QStandardItem(QString(item))
            self.modelResu.appendRow(std_item)


    @pyqtSignature("")
    def slotAddData(self):
        """
        Add data users files input in entries in the good list.
        """
        title = self.tr("Select user data files.")
        filetypes = self.tr("User data files (*);;""All Files (*)")
        list = QFileDialog.getOpenFileNames(self,
                                            title,
                                            self.default['data_path'],
                                            filetypes)
        for item in list:
            self.setFileData(os.path.basename(str(item)))


    @pyqtSignature("")
    def slotNewData(self):
        std_item = QStandardItem(QString(""))
        self.modelData.appendRow(std_item)
        index = self.modelData.indexFromItem(std_item)
        self.viewData.edit(index)


    @pyqtSignature("")
    def slotNewResu(self):
        std_item = QStandardItem(QString(""))
        self.modelResu.appendRow(std_item)
        index = self.modelResu.indexFromItem(std_item)
        self.viewResu.edit(index)


    @pyqtSignature("")
    def slotDeleteData(self):
        """
        Delete the selection from the listbox (one by one).
        """
        index = self.viewData.currentIndex()
        if index.isValid():
            self.modelData.removeRow(index.row())


    @pyqtSignature("")
    def slotDeleteResu(self):
        """
        Delete the selection from the listbox (one by one).
        """
        index = self.viewResu.currentIndex()
        if index.isValid():
            self.modelResu.removeRow(index.row())


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        column = 0

        for row in range(self.modelData.rowCount()):
            index = self.modelData.index(row, column, QModelIndex())
            qstring = index.data(Qt.DisplayRole).toString()
            self.result['data'].append(str(qstring))

        for row in range(self.modelResu.rowCount()):
            index = self.modelResu.index(row, column, QModelIndex())
            qstring = index.data(Qt.DisplayRole).toString()
            self.result['results'].append(str(qstring))

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        self.result  = self.default.copy()
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text


class QListModel(QStandardItemModel):
    def __init__(self, row,  column,  parent=None):
        super(QListModel, self).__init__(row, column, parent)

    def search(self, item):
        result = []
        column = 0
        for row in range(self.rowCount()):
            index = self.index(row, column, QModelIndex())
            qstring = index.data(Qt.DisplayRole).toString()
            if item == qstring:
                result.append(index)

        return result


class DataDelegate(QItemDelegate):
    def __init__(self, parent=None, path=None):
        super(DataDelegate, self).__init__(parent)
        self.path = path
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        vd = RegExpValidator(editor, QRegExp("[_A-Za-z0-9\-\*\!\?\.]*"))
        editor.setValidator(vd)
        editor.setFrame(False)
        self.connect(editor, SIGNAL("returnPressed()"), self.commitAndCloseEditor)
        editor.setCursorPosition(0)
        return editor


    def commitAndCloseEditor(self):
        editor = self.sender()
        if isinstance(editor, QLineEdit):
            self.emit(SIGNAL("commitData(QWidget*)"), editor)
            self.emit(SIGNAL("closeEditor(QWidget*)"), editor)


    def setEditorData(self, editor, index):
        text = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(text)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        item = editor.text()

        if model.search(item):
            model.removeRow(index.row())
            title = self.tr("Warning")
            msg   = self.tr("%s is already in the list." % str(item))
            QMessageBox.warning(self.parent, title, msg)
            return

        path = self.path + "/" + str(item)
        if not os.path.isfile(path) and not os.path.islink(path):
            model.removeRow(index.row())
            title = self.tr("Information")
            msg   = self.tr("%s is not in the data directory:\n\n%s"
                            "\n\nCheck location of this file.\n"
                            "(Note: wildcards are authorized)" % (str(item), self.path))
            QMessageBox.information(self.parent, title, msg)

        model.setData(index, QVariant(item), Qt.DisplayRole)


class ResuDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ResuDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        vd = RegExpValidator(editor, QRegExp("[_A-Za-z0-9\-\*\!\?\.]*"))
        editor.setValidator(vd)
        editor.setFrame(False)
        self.connect(editor, SIGNAL("returnPressed()"), self.commitAndCloseEditor)
        editor.setCursorPosition(0)
        return editor


    def commitAndCloseEditor(self):
        editor = self.sender()
        if isinstance(editor, QLineEdit):
            self.emit(SIGNAL("commitData(QWidget*)"), editor)
            self.emit(SIGNAL("closeEditor(QWidget*)"), editor)


    def setEditorData(self, editor, index):
        text = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(text)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        item = editor.text()

        if model.search(item):
            model.removeRow(index.row())
            title = self.tr("Warning")
            msg   = self.tr("%s is already in the list." % str(item))
            QMessageBox.warning(self.parent, title, msg)
        else:
            model.setData(index, QVariant(item), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Popup advanced options
#-------------------------------------------------------------------------------


class BatchRunningAdvancedOptionsDialogView(QDialog, Ui_BatchRunningAdvancedOptionsDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningAdvancedOptionsDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result  = self.default.copy()

        self.lineEdit.setReadOnly(True)

        # Combo models
        self.modelCSTMPPREFIX  = ComboModel(self.comboBoxCSTMPPREFIX, 2, 1)
        self.modelExecPrepro   = ComboModel(self.comboBox, 2, 1)
        self.modelExecPartit   = ComboModel(self.comboBox_2, 2, 1)
        self.modelExecKernel   = ComboModel(self.comboBox_3, 2, 1)
        self.modelArg_cs_verif = ComboModel(self.comboBox_5, 2, 1)
        self.modelCSOUT1       = ComboModel(self.comboBox_6, 2, 1)
        self.modelCSOUT2       = ComboModel(self.comboBox_7, 3, 1)

        # Combo items
        self.modelCSTMPPREFIX.addItem(self.tr("automatic"), 'automatic')
        self.modelCSTMPPREFIX.addItem(self.tr("prescribed"), 'prescribed')

        self.modelExecPrepro.addItem(self.tr("Run the preprocessor"), 'yes')
        self.modelExecPrepro.addItem(self.tr("Use existing DATA/preprocessor_output"), 'no')

        self.modelExecPartit.addItem(self.tr("Run the partioner"), 'yes')
        self.modelExecPartit.addItem(self.tr("Use existing domain_number_<p> file in DATA/PARTITION_OUTPUT/\n"\
                                             "if present, unoptimized partition otherwise"), 'no')

        self.modelExecKernel.addItem(self.tr("Setup data and run the calculation"), 'yes')
        self.modelExecKernel.addItem(self.tr("Do not setup data and run the calculation"), 'no')

        self.modelArg_cs_verif.addItem(self.tr("Off"), 'standard')
        self.modelArg_cs_verif.addItem(self.tr("Mesh quality criteria"), 'mesh_quality')

        self.modelCSOUT1.addItem(self.tr("to standard output"), 'standard')
        self.modelCSOUT1.addItem(self.tr("to listing"), 'listing')

        self.modelCSOUT2.addItem(self.tr("no output"), 'shunte')
        self.modelCSOUT2.addItem(self.tr("to standard output"), 'standard')
        self.modelCSOUT2.addItem(self.tr("to listing_n<N>"), 'listing')

        # connections
        self.connect(self.comboBoxCSTMPPREFIX, SIGNAL("activated(const QString&)"), self.slotCSTMPPREFIX)
        self.connect(self.toolButton, SIGNAL("clicked()"), self.slotSearchDirectory)
        self.connect(self.comboBox, SIGNAL("activated(const QString&)"), self.slotExePrepro)
        self.connect(self.comboBox_2, SIGNAL("activated(const QString&)"), self.slotExePartit)
        self.connect(self.comboBox_3, SIGNAL("activated(const QString&)"), self.slotExeKernel)
        self.connect(self.toolButton_2, SIGNAL("clicked()"), self.slotSearchFile)
        self.connect(self.lineEdit_2, SIGNAL("textChanged(const QString &)"), self.slotPartitionList)
        self.connect(self.lineEdit_3, SIGNAL("textChanged(const QString &)"), self.slotValgrind)
        self.connect(self.lineEdit_4, SIGNAL("textChanged(const QString &)"), self.slotCs_lib_add)
        self.connect(self.comboBox_5, SIGNAL("activated(const QString&)"), self.slotArgCsVerif)
        self.connect(self.comboBox_6, SIGNAL("activated(const QString&)"), self.slotArgCsOutput)
        self.connect(self.comboBox_7, SIGNAL("activated(const QString&)"), self.slotArgCsOutput)

        # Previous values
        self.dir_name = self.default['CS_TMP_PREFIX']
        self.lineEdit.setText(QString(self.dir_name))
        if self.dir_name == "":
            self.lineEdit.setEnabled(False)
            self.toolButton.setEnabled(False)
            self.modelCSTMPPREFIX.setItem(str_model='automatic')
        else:
            self.lineEdit.setEnabled(True)
            self.toolButton.setEnabled(True)
            self.modelCSTMPPREFIX.setItem(str_model='prescribed')

        self.exe_prepro = self.default['EXEC_PREPROCESS']
        self.modelExecPrepro.setItem(str_model=self.exe_prepro)

        self.exe_partit = self.default['EXEC_PARTITION']
        self.modelExecPartit.setItem(str_model=self.exe_partit)

        self.exe_kernel = self.default['EXEC_KERNEL']
        self.modelExecKernel.setItem(str_model=self.exe_kernel)

        self.partition_list = self.default['PARTITION_LIST']
        self.lineEdit_2.setText(QString(self.partition_list))

        self.valgrind = self.default['VALGRIND']
        self.lineEdit_3.setText(QString(self.valgrind))

        self.cs_lib_add = self.default['CS_LIB_ADD']
        self.lineEdit_4.setText(QString(self.cs_lib_add))

        self.setArgCsVerif()
        self.setArgCsOutput()


    @pyqtSignature("const QString &")
    def slotCSTMPPREFIX(self, text):
        """
        Select mode for CS_TMP_PREFIX.
        """
        if self.modelCSTMPPREFIX.dicoV2M[str(text)] == 'prescribed':
            self.dir_name = self.default['CS_TMP_PREFIX']
            self.lineEdit.setEnabled(True)
            self.toolButton.setEnabled(True)
            setGreenColor(self.toolButton, True)
        else:
            self.dir_name = ""
            self.lineEdit.setEnabled(False)
            self.toolButton.setEnabled(False)
            setGreenColor(self.toolButton, False)
        self.lineEdit.setText(QString(self.dir_name))


    @pyqtSignature("const QString &")
    def slotPartitionList(self, text):
        """
        Input for Partitioner.
        """
        self.partition_list = str(text)


    @pyqtSignature("const QString &")
    def slotValgrind(self, text):
        """
        Input for Valgrind.
        """
        self.valgrind = str(text)


    @pyqtSignature("const QString &")
    def slotCs_lib_add(self, text):
        """
        Input for external libraries.
        """
        self.cs_lib_add = str(text)


    @pyqtSignature("const QString &")
    def slotExePrepro(self, text):
        """
        Preprocessor execution mode option.
        """
        self.exe_prepro = self.modelExecPrepro.dicoV2M[str(text)]


    @pyqtSignature("const QString &")
    def slotExePartit(self, text):
        """
        Partitioner execution mode option.
        """
        self.exe_partit = self.modelExecPartit.dicoV2M[str(text)]


    @pyqtSignature("const QString &")
    def slotExeKernel(self, text):
        """
        Kernel execution mode option.
        """
        self.exe_kernel = self.modelExecKernel.dicoV2M[str(text)]


    def setArgCsVerif(self):
        """
        Put ARG_CS_VERIF option from "lance" file.
        """
        if self.default['ARG_CS_VERIF'] == '':
            self.arg_cs_verif = 'standard'
            self.val_verif = ""
        if self.default['ARG_CS_VERIF'] == '--quality' or self.default['ARG_CS_VERIF'] == '-q':
            self.arg_cs_verif = 'mesh_quality'
            self.val_verif = '--quality'
        self.modelArg_cs_verif.setItem(str_model=self.arg_cs_verif)


    @pyqtSignature("const QString &")
    def slotArgCsVerif(self, text):
        """
        Input ARG_CS_VERIF option.
        """
        self.val_verif = ''
        self.arg_cs_verif = self.modelArg_cs_verif.dicoV2M[str(text)]
        arg_verif = self.arg_cs_verif

        if arg_verif == 'standard'     : self.val_verif = ''
        if arg_verif == 'mesh_quality' : self.val_verif = '--quality'


    def setArgCsOutput(self):
        """
        Put ARG_CS_OUTPUT options from 'lancer' file.
        """
        self.val_output = self.default['ARG_CS_OUTPUT']
        if self.default['ARG_CS_OUTPUT'] == '':
            self.modelCSOUT1.setItem(str_model='listing')
            self.modelCSOUT2.setItem(str_model='shunte')
        else:
            list = self.default['ARG_CS_OUTPUT'].split()
            l1 = 0
            l2 = 0
            for n in range(len(list)):
                if list[n] == '--log':
                    l1 = 1
                    if list[n+1] == '0': self.modelCSOUT1.setItem(str_model='standard')
                    if list[n+1] == '1': self.modelCSOUT1.setItem(str_model='listing')
                if list[n] == '--logp':
                    l2 = 1
                    if list[n+1] == '0': self.modelCSOUT2.setItem(str_model='standard')
                    if list[n+1] == '1': self.modelCSOUT2.setItem(str_model='listing')
                    if list[n+1] == '-1': self.modelCSOUT2.setItem(str_model='shunte')
            if l1 == 0: self.modelCSOUT1.setItem(str_model='listing')
            if l2 == 0: self.modelCSOUT2.setItem(str_model='shunte')


    @pyqtSignature("const QString &")
    def slotArgCsOutput(self, text):
        """
        Input ARG_CS_OUTPUT options.
        """
        self.val_output =''
        out1 = ''
        out2 = ''
        arg_out1 = self.modelCSOUT1.dicoV2M[str(self.comboBox_6.currentText())]
        arg_out2 = self.modelCSOUT2.dicoV2M[str(self.comboBox_7.currentText())]
        if arg_out1 == 'listing': out1 = ''
        if arg_out1 == 'standard': out1 = '--log 0'
        if arg_out2 == 'shunte': out2 = ''
        if arg_out2 == 'standard': out2 = '--logp 0'
        if arg_out2 == 'listing': out2 = '--logp 1'
        self.val_output = out1 + ' ' + out2


    @pyqtSignature("")
    def slotSearchDirectory(self):
        """
        Choice temporary directory for batch
        """
        title    = self.tr("Select directory")
        default  = os.getcwd()
        options  = QFileDialog.ShowDirsOnly # | QFileDialog.DontResolveSymlinks
        dir_name = QFileDialog.getExistingDirectory(self, title, default, options)

        dir = str(dir_name)
        if dir:
            self.dir_name = dir
            setGreenColor(self.toolButton, False)
        self.lineEdit.setText(QString(self.dir_name))

        return self.dir_name


    @pyqtSignature("")
    def slotSearchFile(self):
        """
        Choice temporary directory for batch
        """
        file_name = ""

        title = self.tr("Select file for use VALGRIND option")
        path  = os.getcwd()
        filetypes = self.tr("All Files (*)")
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes)
        file_name = str(file_name)

        # TO CHECK ...
        if file_name:
            self.valgrind = str(self.lineEdit_3.text())
            if not self.valgrind:
                new = file_name + " --tool=memcheck"
            else:
                new = ""
                for i in string.split(self.valgrind):
                    if i == string.split(self.valgrind)[0]:
                        i = file_name
                        new = new + i + ' '
            self.valgrind = new
            self.lineEdit_3.setText(QString(self.valgrind))


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        self.result['CS_TMP_PREFIX']   = self.dir_name
        self.result['EXEC_PREPROCESS'] = self.exe_prepro
        self.result['EXEC_PARTITION']  = self.exe_partit
        self.result['EXEC_KERNEL']     = self.exe_kernel
        self.result['PARTITION_LIST']  = self.partition_list
        self.result['VALGRIND']        = self.valgrind
        self.result['CS_LIB_ADD']      = self.cs_lib_add
        self.result['ARG_CS_VERIF']    = self.val_verif
        self.result['ARG_CS_OUTPUT']   = self.val_output

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class BatchRunningView(QWidget, Ui_BatchRunningForm):
    """
    This class is devoted to the script selection.
    When a new script is selected, The old data frame is deleted and
    a new apropriate frame is open.
    If the batch script file name is known, informations are displayed
    in the apropiate widget.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BatchRunningForm.__init__(self)
        self.setupUi(self)

        self.case = case
        RuncaseModel(self.case)

        # Validators and connections

        if self.case['batch_type']:

            validatorSimpleName = RegExpValidator(self.lineEditJobName,
                                                  QRegExp("[_A-Za-z0-9]*"))
            self.lineEditJobName.setValidator(validatorSimpleName)
            self.lineEditJobGroup.setValidator(validatorSimpleName)
            self.pushButtonRunSubmit.setText("Submit job")

            self.connect(self.lineEditJobName, SIGNAL("textChanged(const QString &)"),
                         self.slotJobName)
            self.connect(self.spinBoxNodes, SIGNAL("valueChanged(int)"),
                         self.slotJobNodes)
            self.connect(self.spinBoxPpn, SIGNAL("valueChanged(int)"),
                         self.slotJobPpn)
            self.connect(self.spinBoxProcs, SIGNAL("valueChanged(int)"),
                         self.slotJobProcs)
            self.connect(self.spinBoxDays, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxHours, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxMinutes, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxSeconds, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.comboBoxClass, SIGNAL("activated(const QString&)"),
                         self.slotClass)
            self.connect(self.lineEditJobGroup, SIGNAL("textChanged(const QString &)"),
                         self.slotJobGroup)

        else:

            self.pushButtonRunSubmit.setText("Start calculation")

            self.connect(self.spinBoxNProcs, SIGNAL("valueChanged(int)"), self.slotParallelComputing)

        self.connect(self.toolButtonSearchBatch, SIGNAL("clicked()"), self.slotSearchBatchScriptFile)
        self.connect(self.toolButtonFiles, SIGNAL("clicked()"), self.slotUserFiles)
        self.connect(self.toolButtonAdvanced, SIGNAL("clicked()"), self.slotAdvancedOptions)
        self.connect(self.pushButtonRunSubmit, SIGNAL("clicked()"), self.slotBatchRunning)

        # initialize Widgets

        name = self.case['batch_type']
        if name:
            self.labelBatchName.setText(QString(name))

        self.hideBatchInfo()

        self.class_list = None

        if self.case['scripts_path']:
            if "runcase" in os.listdir(self.case['scripts_path']):
                self.case['batchScript'] = "runcase"
                self.displayBatchInfo()
                setGreenColor(self.toolButtonSearchBatch, False)


    @pyqtSignature("int")
    def slotParallelComputing(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.mdl.dicoValues['NUMBER_OF_PROCESSORS'] = v
        self.mdl.updateBatchScriptFile('NUMBER_OF_PROCESSORS')


    @pyqtSignature("const QString &")
    def slotJobName(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobName.validator().state == QValidator.Acceptable:
            self.mdl.dicoValues['job_name'] = str(v)
            self.mdl.updateBatchScriptFile('job_name')


    @pyqtSignature("int")
    def slotJobNodes(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.mdl.dicoValues['job_nodes'] = str(self.spinBoxNodes.text())
        self.mdl.updateBatchScriptFile('job_nodes')


    @pyqtSignature("int")
    def slotJobPpn(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.mdl.dicoValues['job_ppn']  = str(self.spinBoxPpn.text())
        self.mdl.updateBatchScriptFile('job_ppn')


    @pyqtSignature("int")
    def slotJobProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.mdl.dicoValues['job_procs']  = str(self.spinBoxProcs.text())
        self.mdl.updateBatchScriptFile('job_procs')


    @pyqtSignature("")
    def slotJobWallTime(self):

        h_cput = self.spinBoxDays.value()*24 + self.spinBoxHours.value()
        m_cput = self.spinBoxMinutes.value()
        s_cput = self.spinBoxSeconds.value()
        self.mdl.dicoValues['job_walltime'] = h_cput*3600 + m_cput*60 + s_cput
        self.mdl.updateBatchScriptFile('job_walltime')


    @pyqtSignature("")
    def slotClass(self):

        self.mdl.dicoValues['job_class'] = str(self.comboBoxClass.currentText())
        if len(self.mdl.dicoValues['job_class']) > 0:
            self.mdl.updateBatchScriptFile('job_class')


    @pyqtSignature("const QString &")
    def slotJobGroup(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobName.validator().state == QValidator.Acceptable:
            self.mdl.dicoValues['job_group'] = str(v)
            self.mdl.updateBatchScriptFile('job_group')


    @pyqtSignature("")
    def slotUserFiles(self):
        """
        Input 'USER_INPUT_FILES' and 'USER_OUTPUT_FILES'
        """
        default = {}
        default['data_path'] = self.case['data_path']
        default['data']      = string.split(self.mdl.dicoValues['USER_INPUT_FILES'])
        default['results']   = string.split(self.mdl.dicoValues['USER_OUTPUT_FILES'])
        log.debug("slotUserFiles -> %s" % str(default))

        dialog = BatchRunningUserFilesDialogView(self, default)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotUserFiles -> %s" % str(result))
            self.mdl.dicoValues['USER_INPUT_FILES']   = string.join(result['data'])
            self.mdl.dicoValues['USER_OUTPUT_FILES'] = string.join(result['results'])
            self.mdl.updateBatchScriptFile('USER_INPUT_FILES')
            self.mdl.updateBatchScriptFile('USER_OUTPUT_FILES')


    @pyqtSignature("")
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        default = {}
        list = ['CS_TMP_PREFIX', 'EXEC_PREPROCESS', 'EXEC_PARTITION', 'EXEC_KERNEL',
                'PARTITION_LIST', 'VALGRIND', 'CS_LIB_ADD',
                'ARG_CS_VERIF', 'ARG_CS_OUTPUT', ]
        for option in list:
            default[option] = self.mdl.dicoValues[option]
        log.debug("slotAdvancedOptions result = %s "%str(default))

        dialog = BatchRunningAdvancedOptionsDialogView(self, default)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions result = %s "%str(result))
            for option in list:
                self.mdl.dicoValues[option] = result[option]
                self.mdl.updateBatchScriptFile(option)


    @pyqtSignature("")
    def slotBatchRunning(self):
        """
        Launch Code_Saturne batch running.
        """
        # Is the file saved?

        if self.case['new'] == "yes" or self.case.isModified():

            title = self.tr("Warning")
            msg   = self.tr("The current case must be saved before "\
                            "running a Code_Saturne batch job.")
            QMessageBox.information(self, title, msg)
            return

        # Test 2: have we a mesh?

        if not GuiParam.matisse :
            node_ecs = self.case.xmlGetNode('solution_domain')
            if not node_ecs.xmlGetNode('meshes_list'):
                if not node_ecs.xmlGetNode('meshes_list').xmlGetNodeList('mesh'):
                    title = self.tr("Warning")
                    msg   = self.tr("You have to select a mesh.\n\n")
                    QMessageBox.information(self, title, msg)
                    return

        # Test 3: have we a trouble with the mesh generation?

        if GuiParam.matisse :
            import Pages.Matisse as Matisse
            if not Matisse.MatisseMeshRunning(self.case).ok :
                title = self.tr("Warning")
                msg   = self.tr("Mesh generation error.\nSee the file 'listsim'")
                QMessageBox.information(self, title, msg)
                return

        # Test 4: verify if boundary definition exists

        if not GuiParam.matisse:
            bd = LocalizationModel('BoundaryZone', self.case)
            if not bd.getZones():
                if self.case['no_boundary_conditions'] == False:
                    title = self.tr("Warning")
                    msg   = self.tr("No boundary definition declared.\n\n")
                    QMessageBox.warning(self, title, msg)
                    self.case['no_boundary_conditions'] = True

        # Build command line

        batch_type = self.case['batch_type']

        script = os.path.join(self.case['scripts_path'], self.case['batchScript'])
        batch1 = self.case['batchScript']

        if not batch_type:
            batch1 = os.path.join(self.case['scripts_path'], "batch")
            batch2 = batch1 + '~'
            try:
                os.rename(batch1, batch2)
            except:
                pass
            cmd = 'nice nohup ' + script + ' | tee ' + batch1 + ' &'
        elif batch_type[0:3] == 'CCC':
            cmd = 'msub ' + batch1
        elif batch_type[0:5] == 'LOADL':
            cmd = 'llsubmit ' + batch1
        elif batch_type[0:3] == 'LSF':
            cmd = 'bsub < ' + batch1
        elif batch_type[0:3] == 'PBS' or batch_type[0:3] == 'SGE':
            cmd = 'qsub ' + batch1
        elif batch_type[0:5] == 'SLURM':
            cmd = 'sbatch ' + batch1
        else:
            pass

        if self.case['salome']:
            import  SalomeHandler
            SalomeHandler.runSolver(self.case, script)
        else:
            saved_path = os.getcwd()
            try:
                os.chdir(self.case['scripts_path'])
                os.system(cmd)
            except Exception:
                pass
            os.chdir(saved_path)


    def getCommandOutput(self, cmd):
        """
        Run a command and return it's standard output.
        """
        import subprocess
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        lines = []
        while True:
            l = p.stdout.readline()
            lines.append(l.strip())
            if len(l) == 0 and p.poll() != None:
                break
        output = p.communicate()

        if p.returncode == 0:
            return lines


    def getClassList(self):
        """
        Layout of the second part of this page.
        """

        self.class_list = []

        try:

            if self.case['batch_type'][0:3] == 'CCC':
                output = self.getCommandOutput('class')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:5] == 'LOADL':
                output = self.getCommandOutput('llclass')
                ignore = True
                for l in output:
                    if l[0:3] == '---':
                        ignore = not ignore
                    elif ignore == False:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'LSF':
                output = self.getCommandOutput('bqueues')
                ignore = True
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'PBS':
                output = self.getCommandOutput('qstat -q')
                ignore = True
                for l in output:
                    if l[0:3] == '---':
                        ignore = not ignore
                    elif ignore == False:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'SGE':
                output = self.getCommandOutput('qconf -sc')
                for l in output:
                    if l[0:1] != '#':
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:5] == 'SLURM':
                output = self.getCommandOutput('sinfo -s')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        name = l.split(' ')[0]
                        if name[-1:] == '*':
                            name = name[:-1]
                        self.class_list.append(name)

        except Exception:
            pass


    def hideBatchInfo(self):
        """
        hide all batch info before displaying a selected subset
        """

        self.groupBoxJob.hide()

        self.labelJobName.hide()
        self.lineEditJobName.hide()
        self.labelNodes.hide()
        self.spinBoxNodes.hide()
        self.labelPpn.hide()
        self.spinBoxPpn.hide()
        self.labelProcs.hide()
        self.spinBoxProcs.hide()
        self.labelClass.hide()
        self.labelWTime.hide()
        self.spinBoxDays.hide()
        self.labelDays.hide()
        self.spinBoxHours.hide()
        self.labelHours.hide()
        self.spinBoxMinutes.hide()
        self.labelMinutes.hide()
        self.spinBoxSeconds.hide()
        self.labelSeconds.hide()
        self.comboBoxClass.hide()
        self.labelJobGroup.hide()
        self.lineEditJobGroup.hide()

        self.labelNProcs.hide()
        self.spinBoxNProcs.hide()


    def displayBatchInfo(self):
        """
        Layout of the second part of this page.
        """
        if hasattr(self, 'mdl'):
            del self.mdl
        self.mdl = BatchRunningModel(self.case)
        self.mdl.readBatchScriptFile()

        self.labelBatchName.show()

        name = self.case['batchScript']
        self.labelBatchName.setText(QString(name))

        self.job_name  = self.mdl.dicoValues['job_name']
        self.job_nodes = self.mdl.dicoValues['job_nodes']
        self.job_ppn  = self.mdl.dicoValues['job_ppn']
        self.job_procs = self.mdl.dicoValues['job_procs']
        self.job_walltime = self.mdl.dicoValues['job_walltime']
        self.job_class  = self.mdl.dicoValues['job_class']
        self.job_group  = self.mdl.dicoValues['job_group']

        if self.job_name != None:
            self.labelJobName.show()
            self.lineEditJobName.setText(QString(self.job_name))
            self.lineEditJobName.show()

        if self.job_nodes != None:
            self.labelNodes.show()
            self.spinBoxNodes.setValue(int(self.job_nodes))
            self.spinBoxNodes.show()

        if self.job_ppn != None:
            self.labelPpn.show()
            self.spinBoxPpn.setValue(int(self.job_ppn))
            self.spinBoxPpn.show()

        if self.job_procs != None:
            self.labelProcs.show()
            self.spinBoxProcs.setValue(int(self.job_procs))
            self.spinBoxProcs.show()

        if self.job_walltime != None:
            seconds = self.job_walltime
            minutes = seconds / 60
            hours = minutes / 60
            days = hours / 24
            seconds = seconds % 60
            minutes = minutes % 60
            hours = hours % 24
            self.spinBoxDays.setValue(days)
            self.spinBoxHours.setValue(hours)
            self.spinBoxMinutes.setValue(minutes)
            self.spinBoxSeconds.setValue(seconds)
            self.labelWTime.show()
            self.spinBoxDays.show()
            self.labelDays.show()
            self.spinBoxHours.show()
            self.labelHours.show()
            self.spinBoxMinutes.show()
            self.labelMinutes.show()
            self.spinBoxSeconds.show()
            self.labelSeconds.show()

        if self.job_class != None:

            # Only one pass here
            if self.class_list == None:
                self.getClassList()
                if len(self.class_list) > 0:
                    for c in self.class_list:
                        self.comboBoxClass.addItem(self.tr(c), QVariant(c))
                else:
                    c = self.job_class
                    self.comboBoxClass.addItem(self.tr(c), QVariant(c))

            # All passes
            try:
                index = self.class_list.index(self.job_class)
                self.comboBoxClass.setCurrentIndex(index)
            except Exception:
                if len(self.class_list) > 0:
                    self.job_class = self.class_list[0]
            self.labelClass.show()
            self.comboBoxClass.show()

        if self.job_group != None:
            self.labelJobGroup.show()
            self.lineEditJobGroup.setText(QString(self.job_group))
            self.lineEditJobGroup.show()

        # Show Job management box

        if self.case['batch_type']:
            if self.case['batch_type'][0:5] == 'LOADL':
                self.groupBoxJob.setTitle("Load Leveler job parameters")
            elif self.case['batch_type'][0:3] == 'LSF':
                self.groupBoxJob.setTitle("LSF job parameters")
            elif self.case['batch_type'][0:3] == 'PBS':
                self.groupBoxJob.setTitle("PBS job parameters")
            elif self.case['batch_type'][0:3] == 'SGE':
                self.groupBoxJob.setTitle("Sun Grid Engine job parameters")
            elif self.case['batch_type'][0:5] == 'SLURM':
                self.groupBoxJob.setTitle("SLURM job parameters")
            else:
                self.groupBoxJob.setTitle("Batch job parameters")
            self.groupBoxJob.show()
            self.labelNProcs.hide()
            self.spinBoxNProcs.hide()
        else:
            nb_proc_str = self.mdl.dicoValues['NUMBER_OF_PROCESSORS']
            if nb_proc_str == "":
                nb_proc_str = "1"
            self.spinBoxNProcs.setValue(int(nb_proc_str))
            self.labelNProcs.show()
            self.spinBoxNProcs.show()

        self.mdl.updateBatchScriptFile()


    @pyqtSignature("const QString &")
    def slotBatchCalculation(self, text):
        """
        1) Look if the batch script file name is allready known
        2) Display the apropriate frame
        """
        log.debug("slotBatchCalculation -> %s" % str(text))
        self.groupBoxBatch.hide()

        self.labelLauncher.setEnabled(True)

        if self.case['batchScript']:
            self.displayBatchInfo()
            setGreenColor(self.toolButtonSearchBatch, False)
        else:
            self.labelBatchName.hide()
            setGreenColor(self.toolButtonSearchBatch, True)


    @pyqtSignature("")
    def slotSearchBatchScriptFile(self):
        """
        Open a FileDialog in order to search the batch script file
        in the system file.
        """
        file_name = ""
        if self.case['scripts_path'] and os.path.isdir(self.case['scripts_path']):
            path = self.case['scripts_path']
        else:
            path = os.getcwd()
        title = self.tr("Select the batch script")
        filetypes = self.tr("All Files (*)")
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes)
        file_name = str(file_name)

        if file_name:

            launcher = os.path.basename(file_name)
            setGreenColor(self.toolButtonSearchBatch, False)

            if self.case['scripts_path'] == os.path.dirname(file_name):
                self.case['batchScript'] = os.path.basename(file_name)
                self.displayBatchInfo()

            else:
                title = self.tr("Warning")
                msg   = self.tr("The new batch script file is not in scripts "\
                                "directory given in the 'Identity and paths' "\
                                "section.\n\n" + \
                                "Verify the existence and location of these files, "\
                                "and the 'Identity and Pathes' section")
                QMessageBox.warning(self, title, msg)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
