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
This module defines the 'Species transport' page.

This module contains the following classes:
- LabelDelegate
- FieldDelegate
- StandardItemModelUserScalar
- SpeciesView
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator, RegExpValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from Species import Ui_Species
from SpeciesModel import SpeciesModel
from code_saturne.Base.Common import LABEL_LENGTH_MAX

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SpeciesView")
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

            if new_plabel in model.mdl.getScalarLabelList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getScalarLabelList()
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

            model.setData(index, to_qvariant(new_plabel), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the field selection
#-------------------------------------------------------------------------------

class FieldDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl):
        super(FieldDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl    = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        for fieldId in self.mdl.getFieldIdList() :
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
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# StandardItemModelUserScalar class
#-------------------------------------------------------------------------------

class StandardItemModelUserScalar(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Label"),
                         self.tr("Carrier Field")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(data)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()
        name = row + 1

        # Label
        if col == 0:
            new_label = from_qvariant(value, to_text_string)
            self._data[row][col] = new_label
            self.mdl.setScalarLabel(name, new_label)
        # field Id
        elif col == 1:
            new_field= from_qvariant(value, to_text_string)
            self._data[row][col] = new_field
            id = self.mdl.getFieldId(self._data[row][col])
            self.mdl.setScalarFieldId(name, id)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_scalar=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        name = ""
        if existing_scalar == None :
            name = self.mdl.addScalar(row + 1)
        else :
            name = existing_scalar

        fieldId    = self.mdl.getScalarFieldId(row + 1)
        labelfield = self.mdl.getLabel(fieldId)
        label      = self.mdl.getScalarLabel(row + 1)

        scalar = [label, labelfield]

        self._data.append(scalar)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        self.mdl.deleteScalar(row + 1)
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# SpeciesView class
#-------------------------------------------------------------------------------

class SpeciesView(QWidget, Ui_Species):
    """
    Species creation layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_Species.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = SpeciesModel(self.case)

        self.tableModelScalar = StandardItemModelUserScalar(self, self.mdl)
        self.tableViewPassifScalar.setModel(self.tableModelScalar)
        self.tableViewPassifScalar.resizeColumnsToContents()
        self.tableViewPassifScalar.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewPassifScalar.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewPassifScalar.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewPassifScalar.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewPassifScalar.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewPassifScalar.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewPassifScalar.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateLabel  = LabelDelegate(self.tableViewPassifScalar)
        delegateField  = FieldDelegate(self.tableViewPassifScalar, self.mdl)

        self.tableViewPassifScalar.setItemDelegateForColumn(0, delegateLabel)
        self.tableViewPassifScalar.setItemDelegateForColumn(1, delegateField)

        # Validators
        validatorDiffusionCoef   = DoubleValidator(self.lineEditDiffusionCoef, min = 0.0)
        validatorSchmidt = DoubleValidator(self.lineEditSchmidt, min = 0.0)
        validatorMin     = DoubleValidator(self.lineEditMinValue, min = 0.0)
        validatorMax     = DoubleValidator(self.lineEditMaxValue, min = 0.0)

        self.lineEditDiffusionCoef.setValidator(validatorDiffusionCoef)
        self.lineEditSchmidt.setValidator(validatorSchmidt)
        self.lineEditMinValue.setValidator(validatorMin)
        self.lineEditMaxValue.setValidator(validatorMax)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddScalar)
        self.pushButtonDelete.clicked.connect(self.slotDeleteScalar)
        self.tableModelScalar.dataChanged.connect(self.dataChanged)
        self.lineEditDiffusionCoef.textChanged[str].connect(self.slotDiffusionCoef)
        self.lineEditSchmidt.textChanged[str].connect(self.slotSchmidt)
        self.lineEditMinValue.textChanged[str].connect(self.slotMinValue)
        self.lineEditMaxValue.textChanged[str].connect(self.slotMaxValue)
        self.checkBoxTimeDepend.clicked.connect(self.slotTimeDepend)
        self.checkBoxDiffusion.clicked.connect(self.slotDiffusion)
        self.tableViewPassifScalar.clicked.connect(self.slotChangeSelection)

        # hide properties by default
        self.groupBoxScalarProperties.hide()

        # load values
        self.currentid = -1

        for scalar in self.mdl.getScalarNameList():
            self.tableModelScalar.newItem(scalar)

        self.case.undoStartGlobal()


    def slotChangeSelection(self, text=None):
        """
        detect change selection to update constant properties
        """
        row = self.tableViewPassifScalar.currentIndex().row()
        self.update(row)


    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewPassifScalar.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewPassifScalar.resizeColumnToContents(col)

        row = self.tableViewPassifScalar.currentIndex().row()
        self.update(row)


    def update(self, row):
        """
        show groupBoxScalarProperties if necessary
        """
        self.currentid = row + 1
        self.groupBoxScalarProperties.show()
        self.lineEditDiffusionCoef.setText(str(self.mdl.getDiffusionCoef(self.currentid)))
        self.lineEditSchmidt.setText(str(self.mdl.getSchmidt(self.currentid)))
        self.lineEditMinValue.setText(str(self.mdl.getMinValue(self.currentid)))
        self.lineEditMaxValue.setText(str(self.mdl.getMaxValue(self.currentid)))
        isTimeDep  = self.mdl.getTimeDependStatus(self.currentid) == "on"
        self.checkBoxTimeDepend.setChecked(isTimeDep)
        isDiffus  = self.mdl.getDiffusionStatus(self.currentid) == "on"
        self.checkBoxDiffusion.setChecked(isDiffus)


    @pyqtSlot()
    def slotAddScalar(self):
        """
        Add a scalar
        """
        self.tableViewPassifScalar.clearSelection()
        self.tableModelScalar.newItem()


    @pyqtSlot()
    def slotDeleteScalar(self):
        """
        Delete the a scalar from the list (one by one).
        """
        row = self.tableViewPassifScalar.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteScalar -> %s" % row)
            self.tableModelScalar.deleteItem(row)
        self.groupBoxScalarProperties.hide()


    @pyqtSlot(str)
    def slotDiffusionCoef(self, text):
        """
        Update the diffusion coefficient
        """
        if self.lineEditDiffusionCoef.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setDiffusionCoef(self.currentid, value)


    @pyqtSlot(str)
    def slotSchmidt(self, text):
        """
        Update the schmidt
        """
        if self.lineEditSchmidt.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setSchmidt(self.currentid, value)


    @pyqtSlot(str)
    def slotMinValue(self, text):
        """
        Update the minimum value
        """
        if self.lineEditMinValue.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMinValue(self.currentid, value)


    @pyqtSlot(str)
    def slotMaxValue(self, text):
        """
        Update the maximum value
        """
        if self.lineEditMaxValue.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMaxValue(self.currentid, value)


    @pyqtSlot(bool)
    def slotTimeDepend(self, checked):
        """
        check box for time depend
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setTimeDependStatus(self.currentid, status)


    @pyqtSlot(bool)
    def slotDiffusion(self, checked):
        """
        check box for diffusion
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setDiffusionStatus(self.currentid, status)

