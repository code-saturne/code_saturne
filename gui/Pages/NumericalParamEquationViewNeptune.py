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
This module defines the 'Equation parameters' page.

This module contains the following classes and function:
- SchemeDelegate
- SolverDelegate
- PrecisionDelegate
- IterationDelegate
- StandardItemModelScheme
- StandardItemModelSolver
- NumericalParamEquationView
"""

#-------------------------------------------------------------------------------
# Library modules import
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, IntValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from NumericalParamEquationNeptune import Ui_NumericalParamEquation
from code_saturne.Pages.NumericalParamEquationModelNeptune import NumericalParamEquatModel
from code_saturne.Pages.GlobalNumericalParametersModel import GlobalNumericalParametersModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NumericalParamEquationView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Combo box delegate for the scheme
#-------------------------------------------------------------------------------

class SchemeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, dicoM2V, dicoV2M):
        super(SchemeDelegate, self).__init__(parent)
        self.parent   = parent
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 3, 1)
        self.modelCombo.addItem(self.tr(self.dicoM2V["centered"]), 'centered')
        self.modelCombo.addItem(self.tr(self.dicoM2V["upwind"]), 'upwind')
        self.modelCombo.addItem(self.tr(self.dicoM2V["solu"]), 'solu')

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("SchemeDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the solver
#-------------------------------------------------------------------------------

class SolverDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, dicoM2V, dicoV2M):
        super(SolverDelegate, self).__init__(parent)
        self.parent   = parent
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 12, 1)
        self.modelCombo.addItem(self.tr(self.dicoM2V["automatic"]), 'automatic')
        self.modelCombo.addItem(self.tr(self.dicoM2V["jacobi"]), 'jacobi')
        self.modelCombo.addItem(self.tr(self.dicoM2V["pcg"]), 'pcg')
        self.modelCombo.addItem(self.tr(self.dicoM2V["cgstab"]), 'cgstab')
        self.modelCombo.addItem(self.tr(self.dicoM2V["jacobi_saturne"]), 'jacobi_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["pcg_saturne"]), 'pcg_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["bicgstab_saturne"]), 'bicgstab_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["bicgstab2_saturne"]), 'bicgstab2_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["gmres_saturne"]), 'gmres_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["gauss_seidel_saturne"]), 'gauss_seidel_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["sym_gauss_seidel_saturne"]), 'sym_gauss_seidel_saturne')
        self.modelCombo.addItem(self.tr(self.dicoM2V["pcr3_saturne"]), 'pcr3_saturne')

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("SolverDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Line edit delegate for the value
#-------------------------------------------------------------------------------

class PrecisionDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(PrecisionDelegate, self).__init__(parent)
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
                    model.setData(idx, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for the value
#-------------------------------------------------------------------------------

class IterationDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(IterationDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = IntValidator(editor, min=1)
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
                    model.setData(idx, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Scheme class
#-------------------------------------------------------------------------------

class StandardItemModelScheme(QStandardItemModel):

    def __init__(self, NPE, dicoM2V, dicoV2M):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = NPE
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

        self.headers = [self.tr("Name"),
                        self.tr("Scheme"),
                        self.tr("Slope\nTest")]

        self.setColumnCount(len(self.headers))

        self._data  = []


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

        elif role == Qt.CheckStateRole:
            data = self._data[index.row()][index.column()]
            if index.column() == 2:
                if data == 'on':
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        elif index.column() == 0 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 2:
            data = self._data[index.row()][1]
            if data == 'upwind':
                return Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
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
        var = self._data[row][0]

        # Scheme
        if col == 1:
            new_scheme = from_qvariant(value, to_text_string)
            self._data[row][col] = new_scheme
            self.mdl.setSchemeModel(var, self.dicoV2M[new_scheme])
            self._data[row][2] = str(self.mdl.getSlopeTestStatus(var))

        # Slope Test
        elif col == 2:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self._data[row][col] = "off"
                self.mdl.setSlopeTestStatus(var,"off")
            else:
                self._data[row][col] = "on"
                self.mdl.setSlopeTestStatus(var,"on")

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, variable):
        """
        Add/load a variable in the model.
        """
        row = self.rowCount()

        scheme = self.dicoM2V[self.mdl.getSchemeModel(variable)]
        slope  = self.mdl.getSlopeTestStatus(variable)

        var = [variable, scheme, slope]

        self._data.append(var)
        self.setRowCount(row+1)


#-------------------------------------------------------------------------------
# Solver class
#-------------------------------------------------------------------------------

class StandardItemModelSolver(QStandardItemModel):

    def __init__(self, NPE, dicoM2V, dicoV2M):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = NPE
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

        self.headers = [self.tr("Name"),
                        self.tr("Solver"),
                        self.tr("Solver\nprecision"),
                        self.tr("Maximum\niteration")]

        self.setColumnCount(len(self.headers))

        self._data  = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(str(data))
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        if index.column() == 0 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
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
        var = self._data[row][0]

        # Solver
        if col == 1:
            new_solver = from_qvariant(value, to_text_string)
            self._data[row][col] = new_solver
            self.mdl.setSolverModel(var, self.dicoV2M[new_solver])

        # Precision solver
        elif col == 2:
            coeff  = from_qvariant(value, float)
            self._data[row][col] = coeff
            self.mdl.setSolverPrecision(var, coeff)

        # Maximum iteration
        elif col == 3:
            coeff  = from_qvariant(value, int)
            self._data[row][col] = coeff
            self.mdl.setMaximumIteration(var, coeff)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, variable):
        """
        Add/load a variable in the model.
        """
        row = self.rowCount()

        solver = self.dicoM2V[self.mdl.getSolverModel(variable)]
        prec   = self.mdl.getSolverPrecision(variable)
        iter   = self.mdl.getMaximumIteration(variable)

        var = [variable, solver, prec, iter]

        self._data.append(var)
        self.setRowCount(row+1)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class NumericalParamEquationView(QWidget, Ui_NumericalParamEquation):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_NumericalParamEquation.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = NumericalParamEquatModel(self.case)

        # dico
        self.dicoM2V = {"automatic"                              : 'Automatic',
                        "jacobi"                                 : 'Jacobi (NEPTUNE_CFD)',
                        "pcg"                                    : 'PCG (NEPTUNE_CFD)',
                        "cgstab"                                 : 'CGSTAB (NEPTUNE_CFD)',
                        "jacobi_saturne"                         : 'Jacobi (Code_Saturne)',
                        "pcg_saturne"                            : 'PCG (Code_Saturne)',
                        "bicgstab_saturne"                       : 'BICGSTAB (Code_Saturne)',
                        "bicgstab2_saturne"                      : 'BICGSTAB-2 (Code_Saturne)',
                        "gmres_saturne"                          : 'GMRES (Code_Saturne)',
                        "gauss_seidel_saturne"                   : 'Gauss-Seidel (Code_Saturne)',
                        "sym_gauss_seidel_saturne"               : 'Sym. Gauss-Seidel (Code_Saturne)',
                        "pcr3_saturne"                           : 'PCR3 (Code_Saturne)',
                        "centered"                               : 'Centered',
                        "upwind"                                 : 'Upwind',
                        "solu"                                   : 'SOLU'}
        self.dicoV2M = {"Automatic"                                    : 'automatic',
                        "Jacobi (NEPTUNE_CFD)"                         : 'jacobi',
                        "PCG (NEPTUNE_CFD)"                            : 'pcg',
                        "CGSTAB (NEPTUNE_CFD)"                         : 'cgstab',
                        "Jacobi (Code_Saturne)"                        : 'jacobi_saturne',
                        "PCG (Code_Saturne)"                           : 'pcg_saturne',
                        "BICGSTAB (Code_Saturne)"                      : 'bicgstab_saturne',
                        "BICGSTAB-2 (Code_Saturne)"                    : 'bicgstab2_saturne',
                        "GMRES (Code_Saturne)"                         : 'gmres_saturne',
                        "Gauss-Seidel (Code_Saturne)"                  : 'gauss_seidel_saturne',
                        "Sym. Gauss-Seidel (Code_Saturne)"             : 'sym_gauss_seidel_saturne',
                        "PCR3 (Code_Saturne)"                          : 'pcr3_saturne',
                        "Centered"                                     : 'centered',
                        "Upwind"                                       : 'upwind',
                        "SOLU"                                         : 'solu'}

        # Scheme
        self.modelScheme = StandardItemModelScheme(self.mdl, self.dicoM2V, self.dicoV2M)
        self.tableViewScheme.setModel(self.modelScheme)
        self.tableViewScheme.setAlternatingRowColors(True)
        self.tableViewScheme.resizeColumnToContents(0)
        self.tableViewScheme.resizeRowsToContents()
        self.tableViewScheme.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewScheme.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewScheme.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewScheme.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewScheme.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateScheme   = SchemeDelegate(self.tableViewScheme, self.dicoM2V, self.dicoV2M)

        self.tableViewScheme.setItemDelegateForColumn(1, delegateScheme)

        # Solveur
        self.modelSolver = StandardItemModelSolver(self.mdl, self.dicoM2V, self.dicoV2M)
        self.tableViewSolver.setModel(self.modelSolver)
        self.tableViewSolver.setAlternatingRowColors(True)
        self.tableViewSolver.resizeColumnToContents(0)
        self.tableViewSolver.resizeRowsToContents()
        self.tableViewSolver.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewSolver.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewSolver.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewSolver.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSolver.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateSolver    = SolverDelegate(self.tableViewSolver, self.dicoM2V, self.dicoV2M)
        delegatePrecision = PrecisionDelegate(self.tableViewSolver)
        delegateIteration = IterationDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(1, delegateSolver)
        self.tableViewSolver.setItemDelegateForColumn(2, delegatePrecision)
        self.tableViewSolver.setItemDelegateForColumn(3, delegateIteration)
        if QT_API == "PYQT4":
            self.tableViewSolver.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSolver.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # Connect signals to slots
        self.modelSolver.dataChanged.connect(self.dataChangedSolver)
        self.modelScheme.dataChanged.connect(self.dataChangedScheme)
        self.tabWidgetScheme.currentChanged[int].connect(self.slotchanged)

        # load variables
        for var in self.mdl.getVariableList() :
            self.modelSolver.newItem(var)
            self.modelScheme.newItem(var)
        self.tableViewSolver.resizeRowsToContents()

        self.tabWidgetScheme.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


    def dataChangedSolver(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewSolver.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewSolver.resizeColumnToContents(col)


    def dataChangedScheme(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewScheme.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewScheme.resizeColumnToContents(col)


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index

