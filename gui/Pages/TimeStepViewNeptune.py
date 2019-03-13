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
This module defines the 'Time step control' page.

This module contains the following classes :
- ValueDelegate
- StandardItemModelCourantFourier
- TimeStepView
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

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, IntValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from TimeStep import Ui_TimeStep
from code_saturne.model.TimeStepModelNeptune import TimeStepModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TimeStepView")
log.setLevel(GuiParam.DEBUG)


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
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
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
# StandardItemModelCourantFourier class
#-------------------------------------------------------------------------------

class StandardItemModelCourantFourier(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Field label"),
                         self.tr("Maximum Courant"),
                         self.tr("Maximum Fourier")]

        self.setColumnCount(len(self.headers))

        self.tooltip = []

        self._data = []
        self.mdl = mdl


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
        if index.column() == 0 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()
        FieldId = row + 1

        # Maximum courant
        if col == 1:
            new_pvalue = from_qvariant(value, float)
            self._data[row][col] = new_pvalue
            self.mdl.setMaxCourant(FieldId, new_pvalue)

        # Maximum fourier
        elif col == 2:
            new_pvalue = from_qvariant(value, float)
            self._data[row][col] = new_pvalue
            self.mdl.setMaxFourier(FieldId,new_pvalue)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, fieldId):
        """
        Add/load a field in the model.
        """
        row = self.rowCount()

        label   = self.mdl.getLabel(fieldId)
        courant = self.mdl.getMaxCourant(fieldId)
        fourier = self.mdl.getMaxFourier(fieldId)

        field = [label, courant, fourier]

        self._data.append(field)
        self.setRowCount(row+1)


#-------------------------------------------------------------------------------
# TimeStepView class
#-------------------------------------------------------------------------------

class TimeStepView(QWidget, Ui_TimeStep):
    """
    Time step control layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_TimeStep.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = TimeStepModel(self.case)

        # Combo box models

        self.modelTimeStepOption = ComboModel(self.comboBoxTimeStepOption, 3, 1)
        self.modelTimeStepOption.addItem(self.tr("Constant"), "constant")
        self.modelTimeStepOption.addItem(self.tr("Adaptive"), "uniform")
        self.modelTimeStepOption.addItem(self.tr("Steady"), "steady")

        self.modelTimeStop = ComboModel(self.comboBoxTimeStopType, 2, 1)
        self.modelTimeStop.addItem(self.tr("Number of time steps"), "iteration")
        self.modelTimeStop.addItem(self.tr("Time analysis (s)"), "time")

        # Validators

        validatorRefT   = DoubleValidator(self.lineEditReferenceTimeStep, min = 0.0)
        validatorNumT   = IntValidator(self.lineEditNumberTimeStep, min = 0)
        validatorAnaT   = DoubleValidator(self.lineEditTimeAnalysis)
        validatorDtMin  = DoubleValidator(self.lineEditDtMin, min = 0.0)
        validatorDtMax  = DoubleValidator(self.lineEditDtMax, min = 0.0)
        validatorIncMax = DoubleValidator(self.lineEditDtIncreasingMax, min = 0.0)
        validatorDecMax = DoubleValidator(self.lineEditDtDecreasingMax, min = 0.0)

        validatorRefT.setExclusiveMin(True)
        validatorNumT.setExclusiveMin(True)
        validatorDtMin.setExclusiveMin(True)
        validatorDtMax.setExclusiveMin(True)
        validatorIncMax.setExclusiveMin(True)
        validatorDecMax.setExclusiveMin(True)

        self.lineEditReferenceTimeStep.setValidator(validatorRefT)
        self.lineEditNumberTimeStep.setValidator(validatorNumT)
        self.lineEditTimeAnalysis.setValidator(validatorAnaT)
        self.lineEditDtMin.setValidator(validatorDtMin)
        self.lineEditDtMax.setValidator(validatorDtMax)
        self.lineEditDtIncreasingMax.setValidator(validatorIncMax)
        self.lineEditDtDecreasingMax.setValidator(validatorDecMax)

        # tableViewCourantFourier

        self.tableModelCourantFourier = StandardItemModelCourantFourier(self.mdl)
        self.tableViewCourantFourier.setModel(self.tableModelCourantFourier)
        self.tableViewCourantFourier.resizeColumnsToContents()
        self.tableViewCourantFourier.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewCourantFourier.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCourantFourier.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCourantFourier.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewCourantFourier.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCourantFourier.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCourantFourier.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewCourantFourier.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewCourantFourier.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateMaxFourier = ValueDelegate(self.tableViewCourantFourier)
        delegateMaxCourant = ValueDelegate(self.tableViewCourantFourier)

        self.tableViewCourantFourier.setItemDelegateForColumn(1, delegateMaxFourier)
        self.tableViewCourantFourier.setItemDelegateForColumn(2, delegateMaxCourant)

        # Connect signals to slots
        self.comboBoxTimeStepOption.activated[str].connect(self.slotTimeStepOption)
        self.comboBoxTimeStopType.activated[str].connect(self.slotTimeStop)
        self.lineEditReferenceTimeStep.textChanged[str].connect(self.slotReferenceTimeStep)
        self.lineEditNumberTimeStep.textChanged[str].connect(self.slotNumberTimeStep)
        self.lineEditTimeAnalysis.textChanged[str].connect(self.slotTimeAnalysis)
        self.lineEditDtMin.textChanged[str].connect(self.slotDtMin)
        self.lineEditDtMax.textChanged[str].connect(self.slotDtMax)
        self.lineEditDtIncreasingMax.textChanged[str].connect(self.slotDtIncreasingMax)
        self.lineEditDtDecreasingMax.textChanged[str].connect(self.slotDtDecreasingMax)
        self.tableModelCourantFourier.dataChanged.connect(self.dataChanged)

        # Initialize widget
        self.initializeVariables()

        self.case.undoStartGlobal()


    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewCourantFourier.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewCourantFourier.resizeColumnToContents(col)


    @pyqtSlot(str)
    def slotTimeStepOption(self, text):
        """
        INPUT time option
        """
        model = self.modelTimeStepOption.dicoV2M[text]
        self.mdl.setTimePassingChoice(model)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotTimeStop(self, text):
        """
        INPUT time stop option
        """
        model = self.modelTimeStop.dicoV2M[text]
        self.mdl.setTimeStopChoice(model)

        if model == "time" :
            self.lineEditNumberTimeStep.hide()
            self.lineEditTimeAnalysis.show()
            value = self.mdl.getMaximumTime()
            self.lineEditTimeAnalysis.setText(str(value))
        else :
            self.lineEditNumberTimeStep.show()
            value = self.mdl.getTimeStepsNumber()
            self.lineEditNumberTimeStep.setText(str(value))
            self.lineEditTimeAnalysis.hide()


    @pyqtSlot(str)
    def slotReferenceTimeStep(self, var):
        """
        """
        if self.lineEditReferenceTimeStep.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setTimeStep(value)


    @pyqtSlot(str)
    def slotNumberTimeStep(self, var):
        """
        """
        if self.lineEditNumberTimeStep.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, int)
            self.mdl.setTimeStepsNumber(value)


    @pyqtSlot(str)
    def slotTimeAnalysis(self, var):
        """
        """
        if self.lineEditTimeAnalysis.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMaximumTime(value)


    @pyqtSlot(str)
    def slotDtMin(self, var):
        """
        """
        if self.lineEditDtMin.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMinDtDt0Variation(value)


    @pyqtSlot(str)
    def slotDtMax(self, var):
        """
        """
        if self.lineEditDtMax.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMaxDtDt0Variation(value)


    @pyqtSlot(str)
    def slotDtIncreasingMax(self, var):
        """
        """
        if self.lineEditDtIncreasingMax.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMaxDtVariationIncreasing(value)


    @pyqtSlot(str)
    def slotDtDecreasingMax(self, var):
        """
        """
        if self.lineEditDtDecreasingMax.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMaxDtVariationDecreasing(value)


    def initializeVariables(self):
        """
        Initialize when a timee step option is choosen
        """
        value = self.mdl.getTimeStep()
        self.lineEditReferenceTimeStep.setText(str(value))

        model = self.mdl.getTimePassingChoice()
        self.modelTimeStepOption.setItem(str_model = model)

        if model == "constant" or model == "steady":
            self.groupBoxAdvancedParameters.hide()
            self.groupBoxCourantFourierParameters.hide()
        else :
            self.groupBoxAdvancedParameters.show()
            self.groupBoxCourantFourierParameters.show()

            value = self.mdl.getMinDtDt0Variation()
            self.lineEditDtMin.setText(str(value))

            value = self.mdl.getMaxDtDt0Variation()
            self.lineEditDtMax.setText(str(value))

            value = self.mdl.getMaxDtVariationIncreasing()
            self.lineEditDtIncreasingMax.setText(str(value))

            value = self.mdl.getMaxDtVariationDecreasing()
            self.lineEditDtDecreasingMax.setText(str(value))

            if (self.tableModelCourantFourier.rowCount() <= 0) :
                for fieldId in self.mdl.getFieldIdList():
                    self.tableModelCourantFourier.newItem(fieldId)

        model = self.mdl.getTimeStopChoice()
        self.modelTimeStop.setItem(str_model = model)

        if model == "time" :
            self.lineEditNumberTimeStep.hide()
            self.lineEditTimeAnalysis.show()
            value = self.mdl.getMaximumTime()
            self.lineEditTimeAnalysis.setText(str(value))
        else :
            self.lineEditNumberTimeStep.show()
            value = self.mdl.getTimeStepsNumber()
            self.lineEditNumberTimeStep.setText(str(value))
            self.lineEditTimeAnalysis.hide()

