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
This module defines the conjugate heat transfer view data management.

This module contains the following classes and function:
- SyrthesVerbosityDelegate
- ProjectionAxisDelegate
- SelectionCriteriaDelegate
- StandardItemModelSyrthes
- ConjugateHeatTransferView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

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

from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.ConjugateHeatTransferForm import Ui_ConjugateHeatTransferForm
from code_saturne.Pages.ConjugateHeatTransferModel import ConjugateHeatTransferModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ConjugateHeatTransferView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QLineEdit delegate for validation of Syrthes verbosity or visualization
#-------------------------------------------------------------------------------

class SyrthesVerbosityDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(SyrthesVerbosityDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QComboBox delegate for Axis Projection in Conjugate Heat Transfer table
#-------------------------------------------------------------------------------

class ProjectionAxisDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent = None):
        super(ProjectionAxisDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem("off")
        editor.addItem("X")
        editor.addItem("Y")
        editor.addItem("Z")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataSyrthes[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for location
#-------------------------------------------------------------------------------

class SelectionCriteriaDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(SelectionCriteriaDelegate, self).__init__(parent)
        self.parent = parent
        self.__model = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "" :
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelSyrthes(QStandardItemModel):

    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)
        self.setColumnCount(5)
        self.headers = [self.tr("Instance name"),
                        self.tr("Verbosity"),
                        self.tr("Visualization"),
                        self.tr("Projection Axis"),
                        self.tr("Selection criteria")]
        self.tooltip = [self.tr("Name of coupled instance"),
                        self.tr("Verbosity level"),
                        self.tr("Visualization output level (0 for none)"),
                        self.tr("Projection axis to match 2D Solid domain"),
                        self.tr("Selection criteria for coupled boundary faces")]
        self.setColumnCount(len(self.headers))
        self.dataSyrthes = []
        self.__model = model


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.ToolTipRole:
            return to_qvariant(self.tooltip[index.column()])
        if role == Qt.DisplayRole:
            return to_qvariant(self.dataSyrthes[index.row()][index.column()])
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
            return

        row = index.row()
        if index.column() in (0, 3, 4):
            self.dataSyrthes[row][index.column()] = str(from_qvariant(value, to_text_string))
        else:
            self.dataSyrthes[row][index.column()] = from_qvariant(value, int)

        num = row + 1
        self.__model.setSyrthesInstanceName(num, self.dataSyrthes[row][0])
        self.__model.setSyrthesVerbosity(num, self.dataSyrthes[row][1])
        self.__model.setSyrthesVisualization(num, self.dataSyrthes[row][2])
        self.__model.setSyrthesProjectionAxis(num, self.dataSyrthes[row][3])
        self.__model.setSelectionCriteria(num, self.dataSyrthes[row][4])

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def addItem(self, syrthes_name,
                verbosity, visualization, proj_axis, location):
        """
        Add a row in the table.
        """
        self.dataSyrthes.append([syrthes_name,
                                 verbosity, visualization, proj_axis, location])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataSyrthes[row]
        row = self.rowCount()
        self.setRowCount(row-1)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ConjugateHeatTransferView(QWidget, Ui_ConjugateHeatTransferForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ConjugateHeatTransferForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.case.undoStopGlobal()

        self.__model = ConjugateHeatTransferModel(self.case)

        # Models
        self.modelSyrthes = StandardItemModelSyrthes(self.__model)
        self.tableViewSyrthes.setModel(self.modelSyrthes)

        if QT_API == "PYQT4":
            self.tableViewSyrthes.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setResizeMode(4, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSyrthes.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)

        delegateSyrthesVerbosity = SyrthesVerbosityDelegate(self.tableViewSyrthes)
        self.tableViewSyrthes.setItemDelegateForColumn(1, delegateSyrthesVerbosity)
        self.tableViewSyrthes.setItemDelegateForColumn(2, delegateSyrthesVerbosity)
        delegateProjectionAxis = ProjectionAxisDelegate(self.tableViewSyrthes)
        self.tableViewSyrthes.setItemDelegateForColumn(3, delegateProjectionAxis)
        delegateSelectionCriteria = SelectionCriteriaDelegate(self.tableViewSyrthes, self.__model)
        self.tableViewSyrthes.setItemDelegateForColumn(4, delegateSelectionCriteria)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddSyrthes)
        self.pushButtonDelete.clicked.connect(self.slotDeleteSyrthes)

        # Insert list of Syrthes couplings for view
        for c in self.__model.getSyrthesCouplingList():
            [syrthes_name, verbosity, visualization, proj_axis, location] = c
            self.modelSyrthes.addItem(syrthes_name,
                                      verbosity, visualization, proj_axis, location)

        if len(self.__model.getSyrthesCouplingList()) < 2:
            self.tableViewSyrthes.hideColumn(0)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddSyrthes(self):
        """
        Set in view label and variables to see on profile
        """
        syrthes_name  = self.__model.defaultValues()['syrthes_name']
        verbosity     = self.__model.defaultValues()['verbosity']
        visualization = self.__model.defaultValues()['visualization']
        proj_axis     = self.__model.defaultValues()['projection_axis']
        location      = self.__model.defaultValues()['selection_criteria']
        num = self.__model.addSyrthesCoupling(syrthes_name,
                                              verbosity, visualization,
                                              proj_axis, location)
        self.modelSyrthes.addItem(syrthes_name, verbosity, visualization,
                                  proj_axis, location)
        if len(self.__model.getSyrthesCouplingList()) > 1:
            self.tableViewSyrthes.showColumn(0)


    @pyqtSlot()
    def slotDeleteSyrthes(self):
        """
        Delete the profile from the list (one by one).
        """
        row = self.tableViewSyrthes.currentIndex().row()
        log.debug("slotDeleteSyrthes -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing coupling")
            QMessageBox.information(self, title, msg)
        else:
            self.modelSyrthes.deleteRow(row)
            self.__model.deleteSyrthesCoupling(row+1)
            if len(self.__model.getSyrthesCouplingList()) < 2:
                self.tableViewSyrthes.hideColumn(0)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
