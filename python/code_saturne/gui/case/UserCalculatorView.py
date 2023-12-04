# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This module defines the 'User define calculator formulae' page.

This module contains the following classes:
- NameDelegate
- LocationDelegate
- StandardItemModelCalculator
- UserCalculatorView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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

from code_saturne.gui.base.QtPage import ComboModel, RegExpValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.UserCalculatorForm import Ui_UserCalculator
from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.case.QMegEditorView import QMegEditorView

from code_saturne.model.UserCalculatorModel import UserCalculatorModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("UserCalculator")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# lineEdit delegate for the name
#-------------------------------------------------------------------------------

class NameDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_pname = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_name = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_pname = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_pname = str(editor.text())

            if new_pname in model.mdl.getFunctionsNamesList():
                default = {}
                default['name']  = self.old_pname
                default['list']   = model.mdl.getFunctionsNamesList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_pname = result['name']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_pname = self.old_pname

            model.setData(index, new_pname, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Location delegate
#-------------------------------------------------------------------------------

class LocationDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(LocationDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 4, 1)
        self.modelCombo.addItem(self.tr("cells"), "cells")
        self.modelCombo.addItem(self.tr("boundary faces"), "boundary")
        # For the moment deactivate internal faces or vertices
#        self.modelCombo.addItem(self.tr("interior faces"), "internal")
#        self.modelCombo.addItem(self.tr("vertices"), "vertices")

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
        log.debug("LocationDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Dimension delegate
#-------------------------------------------------------------------------------

class FuncDimensionDelegate(QItemDelegate):
    """
    Use of a combo box in the table for the user array dimension.
    """

    def __init__(self, parent):
        super(FuncDimensionDelegate, self).__init__(parent)
        self.parent = parent

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 5, 1)

        self.modelCombo.addItem(self.tr("1"), "1")
        self.modelCombo.addItem(self.tr("2"), "2")
        self.modelCombo.addItem(self.tr("3"), "3")
        self.modelCombo.addItem(self.tr("6"), "6")
        self.modelCombo.addItem(self.tr("9"), "9")

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
        log.debug("FuncDimensionDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandardItemModelCalculator class
#-------------------------------------------------------------------------------

class StandardItemModelCalculator(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Name"),
                        self.tr("Location"),
                        self.tr("Dimension")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

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

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
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

        name = self._data[row][0]
        new_val = from_qvariant(value, to_text_string)

        # Name
        if col == 0:
            self.mdl.setFunctionName(name, new_val)
            self._data[row][col] = new_val
        # location
        elif col == 1:
            self.mdl.setFunctionLocation(name, new_val)
            self._data[row][col] = new_val
        # dimension
        elif col == 2:
            self.mdl.setFunctionDim(name, new_val)
            self._data[row][col] = new_val

        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, name=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        if name is None:
            name = "Calc_{}".format(str(row + 1))

        self.mdl.addFunction(name)
        location = self.mdl.getFunctionLocation(name)
        dim      = self.mdl.getFunctionDim(name)

        scalar = [name, location, dim]

        self._data.append(scalar)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        name = self._data[row][0]
        del self._data[row]
        self.mdl.deleteFunction(name)
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# UserCalculatorView class
#-------------------------------------------------------------------------------

class UserCalculatorView(QWidget, Ui_UserCalculator):
    """
    UserCalculatorView layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_UserCalculator.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.mdl = UserCalculatorModel(self.case)

        self.tableModelCalc = StandardItemModelCalculator(self, self.mdl)
        self.tableViewUserCalculator.setModel(self.tableModelCalc)
        if QT_API == "PYQT4":
            self.tableViewUserCalculator.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
            self.tableViewUserCalculator.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
            self.tableViewUserCalculator.horizontalHeader().setResizeMode(2,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewUserCalculator.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
            self.tableViewUserCalculator.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)
            self.tableViewUserCalculator.horizontalHeader().setSectionResizeMode(2,QHeaderView.Stretch)
        self.tableViewUserCalculator.resizeColumnsToContents()
        self.tableViewUserCalculator.resizeRowsToContents()
        self.tableViewUserCalculator.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewUserCalculator.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateName     = NameDelegate(self.tableViewUserCalculator)
        delegateLocation = LocationDelegate(self.tableViewUserCalculator)
        delegateDim      = FuncDimensionDelegate(self.tableViewUserCalculator)

        self.tableViewUserCalculator.setItemDelegateForColumn(0, delegateName)
        self.tableViewUserCalculator.setItemDelegateForColumn(1, delegateLocation)
        self.tableViewUserCalculator.setItemDelegateForColumn(2, delegateDim)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddCalc)
        self.pushButtonDelete.clicked.connect(self.slotDeleteCalc)

        self.pushButtonFormula.clicked.connect(self.slotFormulaCalculator)

        self.tableViewUserCalculator.clicked[QModelIndex].connect(self.__slotSelectCalc)

        # load values
        for user in self.mdl.getFunctionsNamesList():
            self.tableModelCalc.newItem(user)

        self.current_name = None
        self._updatePostView()

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddCalc(self):
        """
        Add a scalar
        """
        self.tableViewUserCalculator.clearSelection()
        self.tableModelCalc.newItem()


    @pyqtSlot()
    def slotDeleteCalc(self):
        """
        Delete the a scalar from the list (one by one).
        """
        row = self.tableViewUserCalculator.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteScalar -> %s" % row)
            self.tableModelCalc.deleteItem(row)


    @pyqtSlot()
    def slotFormulaCalculator(self):
        """
        Modify postprocessing formula
        """

        exp, req, sca, sym = self.mdl.getFunctionFormulaComponents(self.current_name)

        exa = "{} = 1.;".format(self.current_name)

        print(exp)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'pca',
                                zone_name     = 'None',
                                variable_name = self.current_name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaCalculator -> {}".format(str(result)))
            self.mdl.setFunctionFormula(self.current_name, str(result))
            self.pushButtonFormula.setToolTip(result)
            self.pushButtonFormula.setStyleSheet("background-color: green")


    @pyqtSlot("QModelIndex")
    def __slotSelectCalc(self, index):
        """
        Select the user array in the QTable
        """

        _id = self.tableViewUserCalculator.currentIndex();

        row_id = _id.row()
        self.current_name = self.tableModelCalc.getData(_id)[0]

        self._updatePostView()


    def _updatePostView(self):
        """
        Update view
        """

        _is_selected = self.current_name != None
        self.labelFormula.setVisible(_is_selected)
        self.pushButtonFormula.setVisible(_is_selected)

