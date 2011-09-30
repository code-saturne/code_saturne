# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
This module contains the following classes:
- NatureScalarDelegate
- DoubleValueDelegate
- StandardItemModelScalars
- ScalarsBoundariesView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from BoundaryConditionsScalarsForm import Ui_BoundaryConditionsScalarsForm

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsScalarsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Combo box delegate for nature of the scalar
#-------------------------------------------------------------------------------
# TODO: Use a  model to enable/disable items in the combo

class NatureScalarDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None):
        super(NatureScalarDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for item in index.model().dico.values():
            editor.addItem(QString(item))
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"dirichlet": 0, "neumann": 1, "exchange_coefficient": 2}
        row = index.row()
        col = index.column()
        string = index.model().getItem(row)[col]
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = str(comboBox.currentText())

        if value == self.tr("Prescribed value"):
            d = "dirichlet"
        elif value == self.tr("Prescribed flux"):
            d = "neumann"
        elif value == self.tr("Exchange coefficient"):
            d = "exchange_coefficient"

        model.setData(index, QVariant(d), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Delegate for a double value QTableView
#-------------------------------------------------------------------------------

class DoubleValueDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(DoubleValueDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value, ok = editor.text().toDouble()
        if editor.validator().state == QValidator.Acceptable:
            model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class to display scalars values in a QTableView
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):

    def __init__(self, case, boundary):
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Scalar Name"),
                        self.tr("Type"),
                        self.tr("Value"),
                        self.tr("Exchange\nCoefficient")]
        self.setColumnCount(len(self.headers))

        self.case      = case
        self.nature    = boundary.getNature()
        self.boundary  = boundary
        self.sca_model = DefineUserScalarsModel(self.case)

        self.dico = {}
        self.dico["dirichlet"] = self.tr("Prescribed value")
        self.dico["neumann"]   = self.tr("Prescribed flux")
        if self.nature == "wall":
            self.dico["exchange_coefficient"] = self.tr("Exchange coefficient")

        self._data      = []
        self._disabled  = []

        self._populateModel()


    def _populateModel(self):
        for s_label in self.sca_model.getScalarLabelsList():
            log.debug("_initData for scalar label %s " % s_label)
            row = self.rowCount()
            self._newItem(s_label, row)
            self.setRowCount(row+1)


    def _newItem(self, s_label, row):
        line = ["", "", "", ""]
        line[0] = s_label

        if self.nature == 'inlet':
            line[1] = 'dirichlet'
            line[2] = self.boundary.getScalar(s_label)
            if (row,1) not in self._disabled: self._disabled.append((row,1))
            if (row,3) not in self._disabled: self._disabled.append((row,3))

        elif self.nature == 'wall':
            choice = self.boundary.getScalarChoice(s_label)
            line[1] = choice
            if choice == 'dirichlet':
                line[2] = self.boundary.getScalarImposedValue(s_label)
                if (row,3) not in self._disabled: self._disabled.append((row,3))
            elif choice == 'neumann':
                line[2] = self.boundary.getScalarImposedFlux(s_label)
                if (row,3) not in self._disabled: self._disabled.append((row,3))
            elif choice == 'exchange_coefficient':
                line[2] = self.boundary.getScalarImposedValue(s_label)
                line[3] = self.boundary.getScalarExchangeCoefficient(s_label)
                if (row,3) in self._disabled: self._disabled.remove((row,3))

        elif self.nature == 'outlet':
            choice = self.boundary.getScalarChoice(s_label)
            line[1] = choice
            if (row,3) not in self._disabled: self._disabled.append((row,3))
            if choice == 'dirichlet':
                line[2] = self.boundary.getScalar(s_label)
            elif choice == 'neumann':
                line[2] = self.boundary.getScalar(s_label)

        self._data.append(line)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        row = index.row()
        col = index.column()

        if role == Qt.DisplayRole:
            if col == 1:
                if self._data[row][col] in self.dico:
                    return QVariant(self.dico[self._data[row][col]])
                else:
                    return QVariant()
            else:
                return QVariant(self._data[row][col])

        if role == Qt.ToolTipRole:
            if col == 0 :
                return QVariant(self.tr("Code_Saturne keyword: NOMVAR"))

        if role == Qt.StatusTipRole:
            if col == 0:
                return QVariant(self.tr("Scalar name"))

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        elif (index.row(), index.column()) in self._disabled:
            return Qt.ItemIsEnabled

        else:
            if index.column() == 0:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()
        s_label = self._data[row][0]

        # type in dirichlet, neumann or exchange_coefficient
        if col == 1:
            type = str(value.toString())
            self._data[row][col] = type
            self.boundary.setScalarChoice(s_label, type)

            if type == 'dirichlet':
                if self.nature == 'wall':
                    self._data[row][2] = self.boundary.getScalarImposedValue(s_label)
                elif self.nature == 'outlet':
                           self._data[row][2] = self.boundary.getScalar(s_label)
                if (row,3) not in self._disabled: self._disabled.append((row,3))
                self._data[row][3] = ""

            elif type == 'neumann':
                if self.nature == 'wall':
                    self._data[row][2] = self.boundary.getScalarImposedFlux(s_label)
                elif self.nature == 'outlet':
                    self._data[row][2] = self.boundary.getScalar(s_label)
                if (row,3) not in self._disabled: self._disabled.append((row,3))
                self._data[row][3] = ""

            elif type == 'exchange_coefficient':
                if (row,3) in self._disabled: self._disabled.remove((row,3))
                self._data[row][2] = self.boundary.getScalarImposedValue(s_label)
                self._data[row][3] = self.boundary.getScalarExchangeCoefficient(s_label)

        # value(s) associated to the choice
        if col == 2:
            choice = self._data[row][1]
            val, ok = value.toDouble()
            self._data[row][col] = val

            if self.nature == 'inlet':
                self.boundary.setScalar(s_label, val)
            elif self.nature == 'wall':
                if choice == 'dirichlet' or choice == 'exchange_coefficient':
                    self.boundary.setScalarImposedValue(s_label, val)
                elif choice == 'neumann':
                    self.boundary.setScalarImposedFlux(s_label, val)
            elif self.nature == 'outlet':
                self.boundary.setScalar(s_label, val)

        # exchange_coefficient
        if col == 3:
            choice = self._data[row][1]
            coeff, ok = value.toDouble()
            self._data[row][col] = coeff

            if self.nature == 'wall':
                if choice == 'exchange_coefficient':
                    self.boundary.setScalarExchangeCoefficient(s_label, coeff)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def insertItem(self, label, codeNumber, var_nature, local):
        """
        Insert an element in the table view.
        """
        line = [label, codeNumber, var_nature, local]
        self._data.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self._data[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsScalarsView(QWidget, Ui_BoundaryConditionsScalarsForm):
    """
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsScalarsForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None


    def __setBoundary(self, boundary):
        """
        Set the current boundary
        """
        self.__boundary = boundary

        if hasattr(self, "modelScalars"):
            del self.modelScalars

        # Model and QTableView for Scalars
        self.modelScalars = StandardItemModelScalars(self.__case, self.__boundary)
        self.tableViewScalars.setModel(self.modelScalars)
        self.tableViewScalars.resizeColumnsToContents()
        self.tableViewScalars.resizeRowsToContents()
        self.tableViewScalars.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        # Delegates
        delegateNatureScalar = NatureScalarDelegate(self.tableViewScalars)
        delegateDouble = DoubleValueDelegate(self.tableViewScalars)
        self.tableViewScalars.setItemDelegateForColumn(1, delegateNatureScalar)
        self.tableViewScalars.setItemDelegateForColumn(2, delegateDouble)
        self.tableViewScalars.setItemDelegateForColumn(3, delegateDouble)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        if DefineUserScalarsModel(self.__case).getScalarLabelsList():
            self.__setBoundary(boundary)
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
