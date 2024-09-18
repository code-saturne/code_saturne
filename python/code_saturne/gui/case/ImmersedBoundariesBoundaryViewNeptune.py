# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the Immersed boundaries view data management.

This module contains the following classes and function:
 - ImmersedBoundariesBoundaryViewNeptune
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging, os

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
from code_saturne.gui.case.ImmersedBoundariesBoundaryFormNeptune import Ui_ImmersedBoundariesBoundaryFormNeptune
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesBoundaryViewNeptune")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QComboBox delegate for the moving type : Set motion or computed from fluid forces
#-------------------------------------------------------------------------------

class BoundaryConditionNatureTypeDelegate(QItemDelegate):
    """
    Use of a combobox to set the fsi moving type
    """

    def __init__(self, parent, mdl):
        super(BoundaryConditionNatureTypeDelegate, self).__init__(parent)
        self.parent  = parent
        self.mdl     = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ["Wall"]:
            editor.addItem(itm)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().data[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemBoundary(QStandardItemModel):

    def __init__(self, model, case):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Object name"),
                        self.tr("Nature")]
        self.setColumnCount(len(self.headers))

        self.data = []
        self.__model = model
        self.case = case

    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        # Tooltips
        #if role == Qt.ToolTipRole:
        #    return self.tooltip[col]

        # Display
        if role == Qt.DisplayRole:
            return self.data[row][col]

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        if col == 0:
            return Qt.ItemIsSelectable

        elif (col == 1):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            #return
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        num = row + 1

        if col == 1:
            self.data[row][col] = str(from_qvariant(value, to_text_string))
            self.__model.setObjectBoundaryConditionNature(num, self.data[row][col])

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getData(self, index):
        row = index.row()
        return self.data[row]

    def addItem(self, object_name, object_bc_nature):
        """
        Add a row in the table.
        """
        self.data.append([object_name, object_bc_nature])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.data[row]
        row = self.rowCount()
        self.setRowCount(row-1)

    def getItem(self, row):
        return self.data[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesBoundaryViewNeptune(QWidget, Ui_ImmersedBoundariesBoundaryFormNeptune):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesBoundaryFormNeptune.__init__(self)
        self.setupUi(self)

        self.case = case
        self.stbar = stbar
        self.ibm = ImmersedBoundariesModel(self.case)
        self.current_obj = None

        # Models
        self.model_vol = StandardItemBoundary(self.ibm, self.case)
        self.tableViewIBMBoundaryzone.setModel(self.model_vol)

        for obj in range(1,self.ibm.getNumberOfObjects()+1):
            self.model_vol.addItem(self.ibm.getObjectName(obj),
                                   self.ibm.getObjectBoundaryConditionNature(obj))

        if QT_API == "PYQT4":
            self.tableViewIBMBoundaryzone.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewIBMBoundaryzone.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            #self.tableViewIBMBoundaryzone.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewIBMBoundaryzone.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewIBMBoundaryzone.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            #self.tableViewIBMBoundaryzone.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        delegateBCType = BoundaryConditionNatureTypeDelegate(self.tableViewIBMBoundaryzone, self.ibm)
        self.tableViewIBMBoundaryzone.setItemDelegateForColumn(1, delegateBCType)

        self.model_vol.dataChanged.connect(self.dataChanged)

        self.tableViewIBMBoundaryzone.clicked[QModelIndex].connect(self.slotChangedSelection)

        self.updatePageView()

        self.case.undoStartGlobal()


    @pyqtSlot("QModelIndex")
    def slotChangedSelection(self, index):
        """
        detect change in selection and update view
        """
        row = self.tableViewIBMBoundaryzone.currentIndex().row()
        self.current_obj = row + 1
        self.updatePageView()

    def dataChanged(self, topLeft, bottomRight):
        self.updatePageView()

    def updatePageView(self):
        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            self.groupBoxLocalization.hide()
            return

        self.groupBoxLocalization.show()

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
