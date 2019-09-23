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
This module defines the 'Main fields' page.

This module contains the following classes:
- LabelDelegate
- TypeDelegate
- FieldDelegate
- ValueDelegate
- StandardItemModelNonCondensable
- NoncondensableView
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

from code_saturne.Base.QtPage import ComboModel, RegExpValidator, DoubleValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Pages.Noncondensable import Ui_NonCondensable
from code_saturne.model.NonCondensableModel import *
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NonCondensableView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Combo box delegate for the nature of non condensable
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
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

            if new_plabel in model.mdl.getNonCondensableLabelList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getNonCondensableLabelList()
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


class TypeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(TypeDelegate, self).__init__(parent)
        self.parent   = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 5, 1)
        self.modelCombo.addItem(self.tr("H2"), 'H2')
        self.modelCombo.addItem(self.tr("N2"), 'N2')
        self.modelCombo.addItem(self.tr("HE"), 'HE')
        self.modelCombo.addItem(self.tr("O2"), 'O2')
        self.modelCombo.addItem(self.tr("Air"), 'Air')

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_model=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("TypeDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the field selection
#-------------------------------------------------------------------------------

class FieldDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, thermo):
        super(FieldDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl    = mdl
        self.thermo = thermo


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        for fieldId in self.mdl.getGasPhaseList() :
            if self.thermo.getMaterials(fieldId) == "Water" :
                label     = self.mdl.getLabel(fieldId)
                self.modelCombo.addItem(self.tr(label), label)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_model=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("FieldDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Line edit delegate for the value
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor, min=0.)
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = index.model().data(index, Qt.DisplayRole)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandardItemModelNonCondensable class
#-------------------------------------------------------------------------------

class StandardItemModelNonCondensable(QStandardItemModel):

    def __init__(self, parent, mdl, thermo):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Label"),
                         self.tr("Type"),
                         self.tr("Field label"),
                         self.tr("Molar mass"),
                         self.tr("Cobin 1"),
                         self.tr("Cobin 2")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl
        self.thermo = thermo


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

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        row = index.row()
        if index.column() > 3 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else :
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
        name = "mass_fraction_non_condensable_gas_" + str(row + 1)

        # Label
        if col == 0:
            new_label = from_qvariant(value, to_text_string)
            self._data[row][col] = new_label
            self.mdl.setNonCondLabel(name, new_label)
        # Type
        elif col == 1:
            new_type = from_qvariant(value, to_text_string)
            self._data[row][col] = new_type
            self.mdl.setNonCondType(name, new_type)
            self.updateItem(index)
        # field Id
        elif col == 2:
            new_field = from_qvariant(value, to_text_string)
            self._data[row][col] = new_field
            self.mdl.setNonCondFieldId(name, new_field)
        # mass molar
        elif col == 3:
            coeff  = from_qvariant(value, float)
            self._data[row][col] = coeff
            self.mdl.setNonCondMassMol(name, coeff)
        # cobin 1
        elif col == 4:
            coeff  = from_qvariant(value, float)
            self._data[row][col] = coeff
            self.mdl.setNonCondCobin1(name, coeff)
        # cobin 2
        elif col == 5:
            coeff  = from_qvariant(value, float)
            self._data[row][col] = coeff
            self.mdl.setNonCondCobin2(name, coeff)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_noncond=None):
        """
        Add/load a field in the model.
        """
        row = self.rowCount()

        name = ""
        if existing_noncond == None :
           name = self.mdl.addNonCondensable()
        else :
           name = existing_noncond
        label      = self.mdl.getNonCondLabel(name)
        fieldId    = self.mdl.getNonCondFieldId(name)
        labelfield = self.mdl.getLabel(fieldId)
        type       = self.mdl.getNonCondType(name)
        massmol    = self.mdl.getNonCondMassMol(name)
        cobin1     = self.mdl.getNonCondCobin1(name)
        cobin2     = self.mdl.getNonCondCobin2(name)

        nonCond = [label, type, labelfield, massmol, cobin1, cobin2]

        self._data.append(nonCond)
        self.setRowCount(row+1)


    def updateItem(self, index):
        """
        update item
        """
        row = index.row()
        name = "mass_fraction_non_condensable_gas_" + str(row + 1)
        self._data[row][3] = self.mdl.getNonCondMassMol(name)
        self._data[row][4] = self.mdl.getNonCondCobin1(name)
        self._data[row][5] = self.mdl.getNonCondCobin2(name)

        self.flags(index)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        self.mdl.deleteNonCondensable(row + 1)
        row = self.rowCount()
        self.setRowCount(row-1)


#-------------------------------------------------------------------------------
# NonCondensableView class
#-------------------------------------------------------------------------------

class NonCondensableView(QWidget, Ui_NonCondensable):
    """
    Non condensable layout.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NonCondensable.__init__(self)
        self.setupUi(self)

        self.case   = case
        self.case.undoStopGlobal()
        self.mdl    = NonCondensableModel(self.case)
        self.thermo = ThermodynamicsModel(self.case)
        self.browser = tree

        self.tableModelNonCondensable = StandardItemModelNonCondensable(self, self.mdl, self.thermo)
        self.tableViewNonCondensable.setModel(self.tableModelNonCondensable)
        self.tableViewNonCondensable.resizeColumnsToContents()
        self.tableViewNonCondensable.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewNonCondensable.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewNonCondensable.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewNonCondensable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewNonCondensable.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewNonCondensable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewNonCondensable.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateLabel  = LabelDelegate(self.tableViewNonCondensable)
        delegateType   = TypeDelegate(self.tableViewNonCondensable)
        delegateField  = FieldDelegate(self.tableViewNonCondensable, self.mdl, self.thermo)
        delegateMasMol = ValueDelegate(self.tableViewNonCondensable)
        delegateCobin1 = ValueDelegate(self.tableViewNonCondensable)
        delegateCobin2 = ValueDelegate(self.tableViewNonCondensable)
        self.tableViewNonCondensable.setItemDelegateForColumn(0, delegateLabel)
        self.tableViewNonCondensable.setItemDelegateForColumn(1, delegateType)
        self.tableViewNonCondensable.setItemDelegateForColumn(2, delegateField)
        self.tableViewNonCondensable.setItemDelegateForColumn(3, delegateMasMol)
        self.tableViewNonCondensable.setItemDelegateForColumn(4, delegateCobin1)
        self.tableViewNonCondensable.setItemDelegateForColumn(5, delegateCobin2)

        # Connect signals to slots
        self.pushButtonAdd.clicked.connect(self.slotAddNonCondensable)
        self.pushButtonDelete.clicked.connect(self.slotDeleteNonCondensable)
        self.tableModelNonCondensable.dataChanged.connect(self.dataChanged)

        for noncond in self.mdl.getNonCondensableNameList():
            self.tableModelNonCondensable.newItem(noncond)

        self.case.undoStartGlobal()


    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewNonCondensable.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewNonCondensable.resizeColumnToContents(col)


    @pyqtSlot()
    def slotAddNonCondensable(self):
        """
        Add a non condensable
        """
        self.tableViewNonCondensable.clearSelection()
        self.tableModelNonCondensable.newItem()
        self.browser.configureTree(self.case)


    @pyqtSlot()
    def slotDeleteNonCondensable(self):
        """
        Delete a non condensable from the list (one by one).
        """
        row = self.tableViewNonCondensable.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteNonCondensable -> %s" % row)
            self.tableModelNonCondensable.deleteItem(row)
        self.browser.configureTree(self.case)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
