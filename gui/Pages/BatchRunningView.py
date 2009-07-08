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
- BatchRunningPBSJobManagementDialogView
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
from BatchRunningPBSJobManagementDialogForm import Ui_BatchRunningPBSJobManagementDialogForm
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
        title = self.tr("Search user data files.")
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
# Popup window class: Cluster job management
#-------------------------------------------------------------------------------


class BatchRunningPBSJobManagementDialogView(QDialog, Ui_BatchRunningPBSJobManagementDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningPBSJobManagementDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("PBS job management"))

        self.default = default
        self.result  = self.default.copy()

        # Validators
        validatorJobName  = RegExpValidator(rx, self.lineEditJobName, QRegExp("[_A-Za-z0-9]*"))
        validatorNodes    = IntValidator(self.lineEditNodes, min=1)
        validatorCPUNodes = IntValidator(self.lineEditCPUNodes, min=1, max=2)
        validatorHours    = IntValidator(self.lineEditHours, min=0, max=999)
        validatorMinutes  = IntValidator(self.lineEditMinutes, min=0, max=59)
        validatorSeconds  = IntValidator(self.lineEditSeconds, min=0, max=59)
        validatorMemory   = IntValidator(self.lineEditMemory, min=1, max=9999999)

        self.lineEditJobName.setValidator(validatorJobName)
        self.lineEditNodes.setValidator(validatorNodes)
        self.lineEditCPUNodes.setValidator(validatorCPUNodes)

        self.lineEditHours.setValidator(validatorHours)
        self.lineEditMinutes.setValidator(validatorMinutes)
        self.lineEditSeconds.setValidator(validatorSeconds)
        self.lineEditMemory.setValidator(validatorMemory)

        # Previous values
        self.job_name     = self.default['PBS_JOB_NAME']
        self.cluster_node = self.default['PBS_nodes']
        self.cluster_ppn  = self.default['PBS_ppn']
        self.job_mem      = self.default['PBS_mem']
        L = string.split(self.default['PBS_walltime'], ":")
        self.h_cput = L[0]
        self.m_cput = L[1]
        self.s_cput = L[2]

        self.lineEditJobName.setText(QString(self.job_name))
        self.lineEditNodes.setText(QString(str(self.cluster_node)))

        self.lineEditCPUNodes.setText(QString(str(self.cluster_ppn)))
        self.lineEditHours.setText(QString(str(self.h_cput)))
        self.lineEditMinutes.setText(QString(str(self.m_cput)))
        self.lineEditSeconds.setText(QString(str(self.s_cput)))
        self.lineEditMemory.setText(QString(str(self.job_mem)))


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        self.result['PBS_JOB_NAME'] = str(self.lineEditJobName.text())
        if self.lineEditNodes.validator().state == QValidator.Acceptable:
            self.result['PBS_nodes'] = str(self.lineEditNodes.text())
        if self.lineEditCPUNodes.validator().state == QValidator.Acceptable:
            self.result['PBS_ppn']  = str(self.lineEditCPUNodes.text())
        if self.lineEditCPUNodes.validator().state == QValidator.Acceptable and \
            self.lineEditMinutes.validator().state == QValidator.Acceptable and \
            self.lineEditSeconds.validator().state == QValidator.Acceptable:
            h_cput = str(self.lineEditHours.text())
            m_cput = str(self.lineEditMinutes.text())
            s_cput = str(self.lineEditSeconds.text())
            self.result['PBS_walltime'] = h_cput + ":" + m_cput + ":" + s_cput
        if self.lineEditMemory.validator().state == QValidator.Acceptable:
            self.result['PBS_mem']  = str(self.lineEditMemory.text())

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

        # Combo models
        self.modelExecPrepro   = ComboModel(self.comboBox, 2, 1)
        self.modelExecPartit   = ComboModel(self.comboBox_2, 2, 1)
        self.modelExecKernel   = ComboModel(self.comboBox_3, 2, 1)
        self.modelArg_cs_verif = ComboModel(self.comboBox_5, 2, 1)
        self.modelCSOUT1       = ComboModel(self.comboBox_6, 2, 1)
        self.modelCSOUT2       = ComboModel(self.comboBox_7, 3, 1)

        # Combo items
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
        self.exe_name = self.default['CS_TMP_PREFIX']
        self.lineEdit.setText(QString(self.exe_name))

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

        self.dir_name = self.default['CS_TMP_PREFIX']
        self.lineEdit.setEnabled(False)


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
        self.dir_name = ''

        title    = self.tr("Select directory")
        default  = os.getcwd()
        options  = QFileDialog.ShowDirsOnly # | QFileDialog.DontResolveSymlinks
        dir_name = QFileDialog.getExistingDirectory(self, title, default, options)

        self.dir_name = str(dir_name)
        if self.dir_name:
            self.exe_name = self.dir_name
        else:
            self.exe_name = ""
        self.lineEdit.setText(QString(self.exe_name))

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
    This class is devoted to the Computer selection.
    When a new computer is selected, The old data frame is deleted and
    a new apropriate frame is open.
    If the batch script file name is known, informations are display
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

        # Combo model

        self.modelArchi = ComboModel(self.comboBoxArchi, 3, 1)
        self.modelArchi.addItem(self.tr("Workstation"),                   'station')
        self.modelArchi.addItem(self.tr("Cluster with PBS queue system"), 'pbs')
        self.modelArchi.addItem(self.tr("Cluster with LSF queue system"), 'lsf')
        self.modelArchi.addItem(self.tr("Cluster with SGE queue system"), 'sge')
        self.modelArchi.disableItem(str_model='lsf')
        self.modelArchi.disableItem(str_model='sge')

        # Connections

        self.connect(self.comboBoxArchi, SIGNAL("activated(const QString &)"), self.slotBatchCalculation)
        self.connect(self.toolButtonSearch, SIGNAL("clicked()"), self.slotSearchBatchScriptFile)
        self.connect(self.spinBoxProcs, SIGNAL("valueChanged(int)"), self.slotParallelComputing)
        self.connect(self.toolButtonFiles, SIGNAL("clicked()"), self.slotUserFiles)
        self.connect(self.toolButtonAdvanced, SIGNAL("clicked()"), self.slotAdvancedOptions)
        self.connect(self.pushButtonRun, SIGNAL("clicked()"), self.slotBatchRunning)

        # initialize Widgets

        if self.case['computer'] == "":
            key = 'station'
            self.case['computer'] = key
            self.modelArchi.setItem(str_model=key)
            self.slotBatchCalculation(self.tr('Workstation'))
            if self.case['scripts_path']:
                if "runcase" in os.listdir(self.case['scripts_path']):
                    if key in self.case['batchScript']:
                        self.case['batchScript'][key] = "runcase"
                        self.displayBatchScriptInfo()
                        setGreenColor(self.toolButtonSearch, False)
        else:
            self.modelArchi.setItem(str_model=self.case['computer'])
            self.slotBatchCalculation(self.modelArchi.dicoM2V[self.case['computer']])


    @pyqtSignature("int")
    def slotParallelComputing(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.mdl.dicoValues['NUMBER_OF_PROCESSORS'] = v
        self.mdl.updateBatchScriptFile('NUMBER_OF_PROCESSORS')


    @pyqtSignature("")
    def pbsJobManagement(self):
        """
        Get PBS card informations.
        """
        default = {}
        list = ['PBS_JOB_NAME', 'PBS_nodes', 'PBS_ppn', 'PBS_walltime', 'PBS_mem']
        for opt in list:
            default[opt] = self.mdl.dicoValues[opt]
        log.debug("pbsJobManagement -> %s" % str(default))

        dialog = BatchRunningPBSJobManagementDialogView(self, default)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("pbsJobManagement -> %s" % str(result))
            for option in list:
                self.mdl.dicoValues[option] = result[option]

        self.mdl.updateBatchScriptFile()

    @pyqtSignature("")
    def lsfJobManagement(self):
        """
        Get LSF card informations
        """
        pass
##        default = {}
##        default['job_name'] = self.mdl.dicoValues['PBS_JOB_NAME']
##        default['NQS_cput'] =
##        default['NQS_cpuT'] =
##        default['NQS_mem'] =
##        windowTitle = self.tr("CaThy")
##        dialog = LSFJobManagementDialog(self.myPage, title=t.CATHY, default=default)
##
##        self.job_name = dialog.result['job_name']
##        self.NQS_cpult = dialog.result['NQS_cput']
##        self.NQS_cpulT = dialog.result['NQS_cpuT'] + dialog.result['NQS_cput']
##        self.job_memory = dialog.result['NQS_mem']
##
##        self.mdl.updateBatchScriptFile()


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
        # Test 1: is the file saved?

        if self.case['new'] == "yes" or self.case.isModified():

            title = self.tr("Warning")
            msg   = self.tr("The current case must be saved before "\
                            "Code_Saturne batch running.")
            QMessageBox.information(self, title, msg)
            return
#        if self.case.saved() == "no":
#            self.case.xmlSaveDocument()
#            self.mdl.dicoValues['PARAM'] = os.path.basename(self.case['xmlfile'])
#            self.mdl.updateBatchScriptFile('PARAM')

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

##        if not GuiParam.matisse:
##            from DefineBoundaryRegionsModel import DefBCModel
##            groupList = DefBCModel(self.case).getListLabel()
##            if not groupList:
##                if self.case['no_boundary_conditions'] == False:
##                    title = self.tr("Warning")
##                    msg   = self.tr("No boundary definition declared.\n\n")
##                    QMessageBox.warning(self, title, msg)
##                    self.case['no_boundary_conditions'] = True

        if not GuiParam.matisse:
            bd = LocalizationModel('BoundaryZone', self.case)
            if not bd.getZones():
                if self.case['no_boundary_conditions'] == False:
                    title = self.tr("Warning")
                    msg   = self.tr("No boundary definition declared.\n\n")
                    QMessageBox.warning(self, title, msg)
                    self.case['no_boundary_conditions'] = True

        # Command line building

        key = self.case['computer']

        script = self.case['scripts_path'] + "/" + self.case['batchScript'][key]
        batch1 = self.case['scripts_path'] + "/" + "batch"
        batch2 = batch1 + '~'

        if key == 'station':
            try:
                os.rename(batch1, batch2)
            except:
                pass
            cmd = 'nice nohup ' + script + ' | tee ' + batch1 + ' &'
# FIXME: Work in progress
##            dialog = CalculationDialog(self.master, title=t.BATCH, stbar=self.stbar,
##                                       script, batch1)
##            cmd = 'nice nohup ' + script + 'dialog' + ' &'
        elif key == 'pbs':
            #cmd = 'qsub ' + script + ' ' + self.case['batchScript'][key] + ' &'
            cmd = 'qsub ' + script
        elif key == 'lsf':
            cmd = 'bsub ' + script + ' ' + self.case['batchScript'][key] + ' &'
        elif key == 'sge':
            pass

        if self.case['salome']:
            from SalomeHandler import runSolver
            cmd = ['nice', 'nohup', script]
            runSolver(self.case, cmd, self.mdl, batch1)
        else:
            os.system(cmd)


    def displayBatchScriptInfo(self):
        """
        Layout of the second part of this page.
        """
        self.groupBoxBatch.show()

        self.labelJob.hide()
        self.toolButtonJob.hide()

        if hasattr(self, 'mdl'):
            del self.mdl
        self.mdl = BatchRunningModel(self.case)
        self.mdl.readBatchScriptFile()

        self.labelFilename.show()
        self.computer = self.modelArchi.dicoV2M[str(self.comboBoxArchi.currentText())]
        name = self.case['batchScript'][self.computer]
        self.labelFilename.setText(QString(name))

        if self.case['computer'] == 'station':
            if GuiParam.matisse :
                self.labelFiles.show()
                self.toolButtonFiles.show()
                self.labelProcs.hide()
                self.spinBoxProcs.hide()
            else:
                self.labelFiles.show()
                self.toolButtonFiles.show()
                self.labelProcs.show()
                self.spinBoxProcs.show()
        else:
            self.labelProcs.hide()
            self.spinBoxProcs.hide()
            self.labelJob.show()
            self.toolButtonJob.show()
            if self.case['computer'] == "pbs":
                self.connect(self.toolButtonJob, SIGNAL("clicked()"), self.pbsJobManagement)
            if self.case['computer'] == "lsf":
                self.connect(self.toolButtonJob, SIGNAL("clicked()"), self.lsfJobManagement)
            if self.case['computer'] == "sge":
                pass

#            self.lineBatch1.show()
#            self.labelAdvanced.show()
#            self.toolButtonAdvanced.show()

        dico = self.mdl.dicoValues

        if dico['NUMBER_OF_PROCESSORS'] == "":
            dico['NUMBER_OF_PROCESSORS'] = "1"
        if self.case['computer'] == 'station':
            self.spinBoxProcs.setValue(int(dico['NUMBER_OF_PROCESSORS']))

        self.mdl.updateBatchScriptFile()


    @pyqtSignature("const QString &")
    def slotBatchCalculation(self, text):
        """
        1) Look if the batch script file name is allready known
        for the current computer
        2) Display the apropriate frame for the selected computer
        """
        log.debug("slotBatchCalculation -> %s" % str(text))
        self.groupBoxBatch.hide()

        key = self.modelArchi.dicoV2M[str(text)]
        self.case['computer'] = key
        log.debug("slotBatchCalculation -> %s" % key)

        self.labelLauncher.setEnabled(True)

        if key in self.case['batchScript'] and self.case['batchScript'][key]:
            self.displayBatchScriptInfo()
            setGreenColor(self.toolButtonSearch, False)
        else:
            self.labelFilename.hide()
            setGreenColor(self.toolButtonSearch, True)


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
        title = self.tr("Search the batch script")
        filetypes = self.tr("All Files (*)")
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes)
        file_name = str(file_name)

        if file_name:
            launcher = os.path.basename(file_name)
            setGreenColor(self.toolButtonSearch, False)

            if self.case['scripts_path'] == os.path.dirname(file_name):

                self.computer = self.modelArchi.dicoV2M[str(self.comboBoxArchi.currentText())]
                key = self.computer
                if key in self.case['batchScript']:
                    self.case['batchScript'][key] = launcher
                else:
                    print "Warning: slotSearchBatchScriptFile\n Error with key:", key
                self.displayBatchScriptInfo()
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
