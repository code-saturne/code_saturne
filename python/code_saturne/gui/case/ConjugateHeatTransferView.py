# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.ConjugateHeatTransferForm import Ui_ConjugateHeatTransferForm
from code_saturne.model.ConjugateHeatTransferModel import ConjugateHeatTransferModel
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ConjugateHeatTransferView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelSyrthes(QStandardItemModel):

    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)
        self.setColumnCount(1)
        self.headers = [self.tr("Instance name"),
                        self.tr("Boundary zones")
                        ]
        self.tooltip = [self.tr("Name of coupled instance"),
                        self.tr("Boundary zones where the selected instance applies")
                        ]
        self.setColumnCount(len(self.headers))
        self.dataSyrthes = []
        self.__model = model


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]
        if role == Qt.DisplayRole:
            return self.dataSyrthes[index.row()][index.column()]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return

        row = index.row()
        self.dataSyrthes[row][index.column()] = str(from_qvariant(value, to_text_string))
        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True

    def addItem(self, row):
        """
        Add a row in the table.
        """
        self.dataSyrthes.append(row)
        row = self.rowCount()
        self.setRowCount(row+1)

    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataSyrthes[row]
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ConjugateHeatTransferView(QWidget, Ui_ConjugateHeatTransferForm):
    """
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ConjugateHeatTransferForm.__init__(self)
        self.setupUi(self)

    def setup(self, case):

        self.case = case

        self.case.undoStopGlobal()

        self.__model = ConjugateHeatTransferModel(self.case)

        # Models
        self.modelSyrthes = StandardItemModelSyrthes(self.__model)
        self.tableViewSyrthes.setModel(self.modelSyrthes)

        if QT_API == "PYQT4":
            self.tableViewSyrthes.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setResizeMode(1, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSyrthes.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewSyrthes.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

        self.modelProjectionAxis = ComboModel(self.comboBoxProjectionAxis, 1, 1)
        self.modelProjectionAxis.addItem("off")
        self.modelProjectionAxis.addItem("x")
        self.modelProjectionAxis.addItem("y")
        self.modelProjectionAxis.addItem("z")

        # connections
        self.comboBoxProjectionAxis.currentTextChanged[str].connect(self.slotProjectionAxis)
        self.lineEditVerbosity.textChanged[str].connect(self.slotVerbosity)
        self.lineEditVisualization.textChanged[str].connect(self.slotVisualization)
        self.lineEditTolerance.editingFinished.connect(self.slotTolerance)
        self.lineEditVerbosity.setValidator(IntValidator(self.lineEditVerbosity))
        self.lineEditVisualization.setValidator(IntValidator(self.lineEditVisualization))

        _tolValidator = DoubleValidator(self.lineEditTolerance, min=0.)
        _tolValidator.setExclusiveMin()
        self.lineEditTolerance.setValidator(_tolValidator)

        self.initializeParameters()

        self.case.undoStartGlobal()

    def initializeParameters(self):
        projection_axis = self.__model.getSyrthesProjectionAxis()
        if projection_axis not in ["", None]:
            self.comboBoxProjectionAxis.setCurrentText(projection_axis)
        else:
            self.comboBoxProjectionAxis.setCurrentText("off")
            self.slotProjectionAxis("off")

        verbosity = self.__model.getSyrthesVerbosity()
        if verbosity not in ["", None]:
            self.lineEditVerbosity.setText(verbosity)
        else:
            self.lineEditVerbosity.setText("0")

        visualization = self.__model.getSyrthesVisualization()
        if visualization not in ["", None]:
            self.lineEditVisualization.setText(visualization)
        else:
            self.lineEditVisualization.setText("1")

        tolerance = self.__model.getSyrthesTolerance()
        self.lineEditTolerance.setText(tolerance)

        for syrthes_name, boundary_labels in self.__model.getSyrthesCouplingList():
            self.modelSyrthes.addItem([syrthes_name, ", ".join(boundary_labels)])

    @pyqtSlot(str)
    def slotProjectionAxis(self, value):
        self.__model.setSyrthesProjectionAxis(value)
        return

    @pyqtSlot(str)
    def slotVerbosity(self, value):
        self.__model.setSyrthesVerbosity(value)
        pass

    @pyqtSlot(str)
    def slotVisualization(self, value):
        self.__model.setSyrthesVisualization(value)
        pass

    @pyqtSlot()
    def slotTolerance(self):
        """
        Input tolerance value.
        """
        if self.lineEditTolerance.validator().state == QValidator.Acceptable:
            text = self.lineEditTolerance.text()
            self.__model.setSyrthesTolerance(text)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
