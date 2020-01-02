# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import RegExpValidator, IntValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.BalanceForm import Ui_BalanceForm
from code_saturne.Pages.FacesSelectionView import StandardItemModelFaces

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BalanceView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Line edit delegate for selection
#-------------------------------------------------------------------------------

class LineEditDelegateSelector(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = editor.text()
        model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for index
#-------------------------------------------------------------------------------

class LineEditDelegateIndex(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


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
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelPressureDrop(QStandardItemModel):

    def __init__(self, mdl=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.mdl = mdl

        self.headers = [self.tr("Zone id"),
                        self.tr("Selection criteria")]

        self.tooltip = [self.tr("Zone id"),
                        self.tr("Selection criteria string")]

        self.setColumnCount(len(self.headers))

        self._data = []


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            if index.column() in (0, 1):
                return self._data[row][col]

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        if col == 0:
            new_code = from_qvariant(value, int)
            self._data[row][col] = new_code
        elif col == 1:
            criteria = from_qvariant(value, to_text_string)
            self._data[row][col] = criteria
            self.mdl.setPressureDropCriteria(row, criteria)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def addItem(self, idx):
        """
        Add an item in the QListView.
        """
        row = self.rowCount()
        crit = self.mdl.getPressureDropCriteria(idx)
        pres = [idx, crit]
        self._data.append(pres)
        self.setRowCount(row+1)


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)
        for id in range(0, len(self.mdl.getPressureDropList())):
            self._data[id][0] = id


#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelScalarBalance(QStandardItemModel):

    def __init__(self, mdl=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.mdl = mdl

        self.headers = [self.tr("Zone id"),
                        self.tr("variables"),
                        self.tr("Selection criteria")]

        self.tooltip = [self.tr("Zone id"),
                        self.tr("variables"),
                        self.tr("Selection criteria string")]

        self.setColumnCount(len(self.headers))

        self._data = []


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()

            if index.column() in (0, 1, 2):
                return self._data[row][col]

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        if col == 2:
            criteria = from_qvariant(value, to_text_string)
            self._data[row][col] = criteria
            self.mdl.setScalarBalanceCriteria(row, criteria)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def addItem(self, idx):
        """
        Add an item in the QListView.
        """
        row = self.rowCount()
        lst = self.mdl.getVariable(idx)
        crit = self.mdl.getScalarBalanceCriteria(idx)
        pres = [idx, " ; ".join(lst), crit]
        self._data.append(pres)
        self.setRowCount(row+1)


    def replaceItem(self, idx, lst):
        """
        Replace a row in the table.
        """
        self._data[idx][1] = lst


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)
        for id in range(0, len(self.mdl.getScalarBalanceList())):
            self._data[id][0] = id


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------
class BalanceView(QWidget, Ui_BalanceForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)
        Ui_BalanceForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
            from code_saturne.model.BalanceModel import BalanceModel
            self.mdl = BalanceModel(self.case)
        elif self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            from code_saturne.model.BalanceModelNeptune import BalanceModelNeptune
            self.mdl = BalanceModelNeptune(self.case)

        # tableView Pressure Drop
        self.pressureModel = StandardItemModelPressureDrop(self.mdl)
        self.tableViewPressureDrop.setModel(self.pressureModel)
        self.tableViewPressureDrop.setAlternatingRowColors(True)
        self.tableViewPressureDrop.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewPressureDrop.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateIdxPres = LineEditDelegateIndex(self.tableViewPressureDrop)
        self.tableViewPressureDrop.setItemDelegateForColumn(0, delegateIdxPres)

        delegateSelector = LineEditDelegateSelector(self.tableViewPressureDrop)
        self.tableViewPressureDrop.setItemDelegateForColumn(1, delegateSelector)

        self.tableViewPressureDrop.resizeColumnsToContents()
        self.tableViewPressureDrop.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewPressureDrop.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewPressureDrop.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)

        # tableView scalar balance
        self.modelScalarBalance = StandardItemModelScalarBalance(self.mdl)
        self.tableViewScalarBalance.setModel(self.modelScalarBalance)
        self.tableViewScalarBalance.resizeColumnToContents(0)
        self.tableViewScalarBalance.setAlternatingRowColors(True)
        self.tableViewScalarBalance.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewScalarBalance.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateIdxScalar = LineEditDelegateIndex(self.tableViewScalarBalance)
        self.tableViewScalarBalance.setItemDelegateForColumn(0, delegateIdxScalar)

        delegateVariable = LineEditDelegateSelector(self.tableViewScalarBalance)
        self.tableViewScalarBalance.setItemDelegateForColumn(1, delegateVariable)

        delegateSelector = LineEditDelegateSelector(self.tableViewScalarBalance)
        self.tableViewScalarBalance.setItemDelegateForColumn(1, delegateSelector)

        self.tableViewScalarBalance.resizeColumnsToContents()
        self.tableViewScalarBalance.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewScalarBalance.horizontalHeader().setResizeMode(2,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewScalarBalance.horizontalHeader().setSectionResizeMode(2,QHeaderView.Stretch)

        # QListView layout
        self.gridlayout1 = QGridLayout(self.widgetDrag)
        self.gridlayout1.setContentsMargins(0, 0, 0, 0)
        self.DragList = QListView(self.widgetDrag)
        self.gridlayout1.addWidget(self.DragList,0,0,1,1)

        self.gridlayout2 = QGridLayout(self.widgetDrop)
        self.gridlayout2.setContentsMargins(0, 0, 0, 0)
        self.DropList = QListView(self.widgetDrop)
        self.gridlayout2.addWidget(self.DropList,0,0,1,1)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.DragList.setModel(self.modelDrag)
        self.DropList.setModel(self.modelDrop)
        self.DragList.setAlternatingRowColors(True)
        self.DragList.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.DropList.setAlternatingRowColors(True)
        self.DropList.setEditTriggers(QAbstractItemView.NoEditTriggers)

        # Connections
        self.pushButtonAddPressureDrop.clicked.connect(self.slotAddPressureDrop)
        self.pushButtonDeletePressureDrop.clicked.connect(self.slotDeletePressureDrop)
        self.tableViewScalarBalance.pressed[QModelIndex].connect(self.slotSelectScalarBalance)
        self.pushButtonAddScalarBalance.clicked.connect(self.slotAddScalarBalance)
        self.pushButtonDeleteScalarBalance.clicked.connect(self.slotDeleteScalarBalance)
        self.pushButtonAddVar.clicked.connect(self.slotAddVarProfile)
        self.pushButtonSuppressVar.clicked.connect(self.slotDeleteVarProfile)

        if self.mdl.getPressureDropList() != None:
            for i in range(len(self.mdl.getPressureDropList())):
                self.pressureModel.addItem(i)

        if self.mdl.getScalarBalanceList() != None:
            for i in range(len(self.mdl.getScalarBalanceList())):
                self.modelScalarBalance.addItem(i)

        #update list of variables, properties, scalars ...
        liste_label = []
        for label in self.mdl.getScalarVariables():
            liste_label.append(label)
        self.modelDrag.setStringList(sorted(liste_label, key=str.lower))

        self.__eraseEntries()

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddPressureDrop(self):
        """
        Add Pressure Drop
        """
        self.mdl.addPressureDrop()
        self.pressureModel.addItem(len(self.mdl.getPressureDropList()) -1)
        self.tableViewPressureDrop.clearSelection()


    @pyqtSlot()
    def slotDeletePressureDrop(self):
        """
        Delete the selected Pressure Drop from the list
        """
        idx = self.tableViewPressureDrop.currentIndex().row()
        self.mdl.delPressureDrop(idx)
        self.pressureModel.delItem(idx)
        self.tableViewPressureDrop.clearSelection()


    @pyqtSlot()
    def slotAddScalarBalance(self):
        """
        Set in view label and variables to see on scalar balance
        """
        self.mdl.addScalarBalance()
        self.modelScalarBalance.addItem(len(self.mdl.getScalarBalanceList()) -1)
        self.__eraseEntries()


    @pyqtSlot()
    def slotDeleteScalarBalance(self):
        """
        Delete the scalar balance from the list (one by one).
        """
        row = self.tableViewScalarBalance.currentIndex().row()
        log.debug("slotDeleteScalarBalance -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing scalar balance")
            QMessageBox.information(self, title, msg)
        else:
            self.mdl.deleteScalarBalance(row)
            self.modelScalarBalance.delItem(row)
            self.__eraseEntries()


    @pyqtSlot("QModelIndex")
    def slotSelectScalarBalance(self, index):
        """
        Return the selected item from the list.
        """
        self.groupBoxScalarBalance.show()

        row = index.row()
        log.debug("slotSelectScalarBalance -> %s" % (row,))

        liste = self.mdl.getVariable(row)

        self.modelDrop.setStringList([])
        liste = [str(s) for s in liste]

        self.modelDrop.setStringList(liste)


    @pyqtSlot()
    def slotAddVarProfile(self):
        """
        Add a new var from list to profile
        """
        row = self.tableViewScalarBalance.currentIndex().row()
        if (self.DragList.currentIndex().row() >=0) :
            liste = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.DragList.currentIndex().row()]
            if var not in liste :
                liste.append(var)
            liste = [str(s) for s in liste]
            self.modelDrop.setStringList(liste)
            self.mdl.setVariable(row, liste)

            row = self.tableViewScalarBalance.currentIndex().row()
            liste = self.mdl.getVariable(row)
            self.modelScalarBalance.replaceItem(row, " ; ".join(liste))


    @pyqtSlot()
    def slotDeleteVarProfile(self):
        """
        Supress a var from profile
        """
        row = self.tableViewScalarBalance.currentIndex().row()
        self.modelDrop.removeRows(self.DropList.currentIndex().row(), 1)
        liste = self.modelDrop.stringList()
        liste = [str(s) for s in liste]
        self.mdl.setVariable(row, liste)

        row = self.tableViewScalarBalance.currentIndex().row()
        liste = self.mdl.getVariable(row)
        self.modelScalarBalance.replaceItem(row, " ; ".join(liste))


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.groupBoxScalarBalance.hide()
        self.tableViewScalarBalance.clearSelection()


    def tr(self, text):
        """
        Translation
        """
        return text



#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    BalanceView = BalanceView(app)
    BalanceView.show()
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
