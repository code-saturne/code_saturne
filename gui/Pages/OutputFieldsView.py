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
This module defines the 'Output fields' page.

This module contains the following classes :
- ProbesValidator
- ProbesDelegate
- NameDelegate
- StandardItemModelGlobalVariables
- OutputFieldsView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.QtPage import ComboModel, RegExpValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from OutputFields import Ui_OutputFields
from code_saturne.model.OutputFieldsModel import *

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputFieldsView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbesValidator(QRegExpValidator):
    """
    Validator for real data.
    """
    def __init__(self, parent, xml_model):
        """
        Initialization for validator
        """
        regExp = QRegExp("^[0-9 ]*$")
        super(ProbesValidator, self).__init__(regExp, parent)
        self.parent = parent
        self.mdl = xml_model
        self.state = QValidator.Invalid


    def validate(self, stri, pos):
        """
        Validation method.

        QValidator.Invalid       0  The string is clearly invalid.
        QValidator.Intermediate  1  The string is a plausible intermediate value during editing.
        QValidator.Acceptable    2  The string is acceptable as a final result; i.e. it is valid.
        """
        state = QRegExpValidator.validate(self, stri, pos)[0]

        valid = True
        for probe in str(stri).split():
            if probe not in self.mdl.getVariableProbeList():
                valid = False

        if state == QValidator.Acceptable:
            if not valid:
                state = QValidator.Intermediate

        palette = self.parent.palette()

        if state == QValidator.Intermediate:
            palette.setColor(QPalette.Text, QColor("red"))
            self.parent.setPalette(palette)
        else:
            palette.setColor(QPalette.Text, QColor("black"))
            self.parent.setPalette(palette)

        self.state = state

        return (state, stri, pos)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbesDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None, xml_model=None):
        super(ProbesDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = ProbesValidator(editor, self.mdl)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = editor.text()
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for the name
#-------------------------------------------------------------------------------

class NameDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_plabel = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_label = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_plabel = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_plabel = str(editor.text())

            if new_plabel in model.mdl.getVariableLabelsList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getVariableLabelsList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, new_plabel, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandardItemModelGlobalVariables class
#-------------------------------------------------------------------------------

class StandardItemModelGlobalVariables(QStandardItemModel):

    def __init__(self, parent, fieldId, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Name"),
                         self.tr("Listing"),
                         self.tr("Writer"),
                         self.tr("Probes")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self.currentid = fieldId

        self._data  = []
        self.mdl    = mdl


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return None

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return data
            else:
                return None

        elif role == Qt.CheckStateRole:
            data = self._data[index.row()][index.column()]
            if index.column() == 1 or index.column() == 2 :
                if data == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        row = index.row()
        if index.column() == 1 or index.column() == 2 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()

        # Name
        if col == 0:
            oldlabel = self._data[row][col]
            new_plabel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_plabel
            self.mdl.setVariableLabel(self.currentid, new_plabel, oldlabel)
        # Listing
        elif col == 1:
            state = from_qvariant(value, int)
            label = self._data[row][0]
            if state == Qt.Unchecked:
                self._data[row][col] = "off"
                self.mdl.setListingStatus(self.currentid, label, "off")
            else:
                self._data[row][col] = "on"
                self.mdl.setListingStatus(self.currentid, label, "on")
        # Writer
        elif col == 2:
            state = from_qvariant(value, int)
            label = self._data[row][0]
            if state == Qt.Unchecked:
                self._data[row][col] = "off"
                self.mdl.setPostProcessingStatus(self.currentid, label, "off")
            else:
                self._data[row][col] = "on"
                self.mdl.setPostProcessingStatus(self.currentid, label, "on")
        # Probe
        elif col == 3:
            label = self._data[row][0]
            lst = from_qvariant(value, to_text_string)
            self._data[row][col] = lst
            self.mdl.setProbesList(self.currentid, label, lst)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, fieldId, var):
        """
        Load a variable in the model.
        """
        row = self.rowCount()

        label     = self.mdl.getVariableLabel(fieldId, var)
        listing   = self.mdl.getListingStatus(fieldId, var)
        writer    = self.mdl.getPostProcessingStatus(fieldId, var)
        listProbe = self.mdl.getProbesList(fieldId, var)
        probes = string.join(listProbe," ")

        variable = [label, listing, writer, probes]

        self._data.append(variable)
        self.setRowCount(row+1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def updateField(self, field):
        """
        update field when selection change
        """
        self.currentid = field


#-------------------------------------------------------------------------------
# NonCondensableView class
#-------------------------------------------------------------------------------

class OutputFieldsView(QWidget, Ui_OutputFields):
    """
    Output Fields layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_OutputFields.__init__(self)
        self.setupUi(self)

        self.case   = case
        self.case.undoStopGlobal()
        self.mdl    = OutputFieldsModel(self.case)

        # Combo box models
        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for fieldId in self.mdl.getFieldIdList() :
            label = self.mdl.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)

        self.currentid = -1
        if len(self.mdl.getFieldIdList()) > 0 :
            self.currentid = self.mdl.getFieldIdList()[0]
            self.modelField.setItem(str_model = self.currentid)

        self.tableModelGlobalVariables = StandardItemModelGlobalVariables(self, "none", self.mdl)
        self.tableViewGlobalVariables.setModel(self.tableModelGlobalVariables)
        self.tableViewGlobalVariables.resizeColumnsToContents()
        self.tableViewGlobalVariables.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewGlobalVariables.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewGlobalVariables.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewGlobalVariables.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewGlobalVariables.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewGlobalVariables.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewGlobalVariables.setSelectionMode(QAbstractItemView.SingleSelection)

        self.tableModelFieldsVariables = StandardItemModelGlobalVariables(self, self.currentid, self.mdl)
        self.tableViewFieldsVariables.setModel(self.tableModelFieldsVariables)
        self.tableViewFieldsVariables.resizeColumnsToContents()
        self.tableViewFieldsVariables.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewFieldsVariables.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFieldsVariables.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFieldsVariables.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFieldsVariables.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewFieldsVariables.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewFieldsVariables.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateNameG   = NameDelegate(self.tableViewGlobalVariables)
        delegateNameF   = NameDelegate(self.tableViewFieldsVariables)
        probesDelegateG = ProbesDelegate(self.tableViewGlobalVariables, self.mdl)
        probesDelegateF = ProbesDelegate(self.tableViewFieldsVariables, self.mdl)

        self.tableViewGlobalVariables.setItemDelegateForColumn(0, delegateNameG)
        self.tableViewGlobalVariables.setItemDelegateForColumn(3, probesDelegateG)
        self.tableViewFieldsVariables.setItemDelegateForColumn(0, delegateNameF)
        self.tableViewFieldsVariables.setItemDelegateForColumn(3, probesDelegateF)

        # Connect signals to slots
        self.tableModelGlobalVariables.dataChanged.connect(self.dataChangedGlobalVariables)
        self.tableModelFieldsVariables.dataChanged.connect(self.dataChangedFieldsVariables)
        self.comboBoxField.activated[str].connect(self.slotField)

        for var in self.mdl.getGlobalVariables() :
            self.tableModelGlobalVariables.newItem("none", var)

        if (len(self.mdl.getFieldIdList())> 0) :
            self.groupBoxField.show()
            self.initializeVariables(self.currentid)
        else :
            self.groupBoxField.hide()

        self.case.undoStartGlobal()


    def dataChangedGlobalVariables(self, topLeft, bottomRight) :
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewGlobalVariables.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewGlobalVariables.resizeColumnToContents(col)


    def dataChangedFieldsVariables(self, topLeft, bottomRight) :
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewFieldsVariables.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewFieldsVariables.resizeColumnToContents(col)


    @pyqtSlot(str)
    def slotField(self, text):
        """
        INPUT label for choice of field
        """
        self.currentid = self.modelField.dicoV2M[text]
        self.initializeVariables(self.currentid)


    def initializeVariables(self, FieldId) :
        """
        Update variable and property in table
        """
        while self.tableModelFieldsVariables.rowCount():
            self.tableModelFieldsVariables.deleteItem(0)

        for var in self.mdl.getFieldVariables(self.currentid) :
            self.tableModelFieldsVariables.newItem(self.currentid, var)

        self.tableModelFieldsVariables.updateField(self.currentid)
