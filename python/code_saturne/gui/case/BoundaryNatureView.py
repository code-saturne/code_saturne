# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------

"""
This module contains the following classes:
- LabelDelegate
- LocalizationSelectorDelegate
- StandardItemModelLocalization
- LocalizationView
"""

# -------------------------------------------------------------------------------
# Standard modules
# -------------------------------------------------------------------------------

import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.base.QtPage import IntValidator, RegExpValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.BoundaryNatureForm import Ui_BoundaryNatureForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryNatureView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# QComboBox delegate for the boundary nature
# -------------------------------------------------------------------------------

class BoundaryNatureDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """

    def __init__(self, parent, dicoM2V):
        super(BoundaryNatureDelegate, self).__init__(parent)
        self.parent = parent
        self.dicoM2V = dicoM2V

        self.dicoV2M = {}
        for k, v in list(self.dicoM2V.items()):
            self.dicoV2M[v] = k

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for k in list(self.dicoV2M.keys()):
            editor.addItem(k)
        editor.installEventFilter(self)
        editor.setMinimumWidth(120)
        return editor

    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        str_model = index.model().getData(row, col)
        idx = list(self.dicoM2V.keys()).index(str_model)
        comboBox.setCurrentIndex(idx)

    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.dicoV2M[txt]
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


# -------------------------------------------------------------------------------
# StandarItemModel class
# -------------------------------------------------------------------------------

class StandardItemModelLocalization(QStandardItemModel):
    def __init__(self, mdl, zoneType, dicoM2V, tree=None, case=None):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"),
                        self.tr("Nature")]
        self.setColumnCount(len(self.headers))

        self.mdl = mdl
        self.zoneType = zoneType
        self.dicoM2V = dicoM2V
        self.browser = tree
        self.case = case

        self._data = []
        self._disable = []

    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()

            if col == 0:
                return self._data[row][col]

            elif col == 1:
                key = self._data[row][col]
                return self.dicoM2V[key]

        return None

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0:
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None

    def setData(self, index, value, role):
        col = index.column()

        if col == 1:
            row = index.row()

            nature = str(from_qvariant(value, to_text_string))
            self._data[row][1] = nature

            self.mdl.setNature(self._data[row][0], nature)

        self.dataChanged.emit(index, index)
        self.browser.configureTree(self.case)
        return True

    def addItem(self, zone):
        """
        Add an element in the table view.
        """
        line = [zone.getLabel(),
                zone.getNature()]
        self._data.append(line)
        row = self.rowCount()
        self.setRowCount(row + 1)

        # If row is disabled, it appears in light gray
        # self._disable.append((row, 0))
        self.browser.configureTree(self.case)
        return zone

    def getItem(self, row):
        return self._data[row]

    def getData(self, row, column):
        return self._data[row][column]


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class BoundaryNatureView(QWidget, Ui_BoundaryNatureForm):
    """
    Main class
    """

    def __init__(self, parent, case, tree=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryNatureForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        dicoM2V = Zone("BoundaryZone", case=self.case).getModel2ViewDictionary()
        zoneType = "BoundaryZone"
        self.zoneType = "BoundaryZone"

        self.mdl = LocalizationModel("BoundaryZone", case)
        self.case['dump_python'].append([self.mdl.__module__, "BoundaryZone", ()])

        self.browser = tree

        # Delegates
        delegateNature = BoundaryNatureDelegate(self.tableView, dicoM2V)

        # Model for table View
        self.modelLocalization = StandardItemModelLocalization(self.mdl,
                                                               zoneType,
                                                               dicoM2V,
                                                               tree,
                                                               case)
        self.tableView.setModel(self.modelLocalization)
        self.tableView.setItemDelegateForColumn(1, delegateNature)
        last_section = 1

        # Populate QTableView model
        for zone in self.mdl.getZones():
            self.modelLocalization.addItem(zone)

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(last_section, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(last_section, QHeaderView.Stretch)

        # Connections
        self.modelLocalization.dataChanged.connect(self.dataChanged)
        self.tableView.clicked.connect(self.slotChangeSelection)

        # Context menu
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested[QPoint].connect(self.slotContextMenu)

        self.case.undoStartGlobal()

    def slotChangeSelection(self):
        """
        """
        current = self.tableView.currentIndex()

    @pyqtSlot()
    def slotContextMenu(self):
        """
        Public slot

        Create the popup menu of the Browser
        """
        fileMenu = QMenu(self.tableView)

        self.actionInlet = QAction(self.tr("Select all inlets"), self.tableView)
        self.actionInlet.triggered.connect(self.slotSelectBoundaries)
        fileMenu.addAction(self.actionInlet)

        self.actionOutlet = QAction(self.tr("Select all outlets"), self.tableView)
        self.actionOutlet.triggered.connect(self.slotSelectBoundaries)
        fileMenu.addAction(self.actionOutlet)

        self.actionWall = QAction(self.tr("Select all walls"), self.tableView)
        self.actionWall.triggered.connect(self.slotSelectBoundaries)
        fileMenu.addAction(self.actionWall)

        self.actionSymmetry = QAction(self.tr("Select all symmetries"), self.tableView)
        self.actionSymmetry.triggered.connect(self.slotSelectBoundaries)
        fileMenu.addAction(self.actionSymmetry)

        fileMenu.popup(QCursor().pos())
        fileMenu.show()

    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row() + 1):
            self.tableView.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column() + 1):
            self.tableView.resizeColumnToContents(col)

    @pyqtSlot()
    def slotSelectBoundaries(self):
        """
        Public slot.

        Warning: works only if the selection mode of the view is set to MultiSelection.
        """
        previous_selecion_mode = self.tableView.selectionMode()
        self.tableView.setSelectionMode(QAbstractItemView.MultiSelection)
        self.tableView.clearSelection()

        if self.sender() == self.actionInlet:
            select = "inlet"
        elif self.sender() == self.actionOutlet:
            select = "outlet"
        elif self.sender() == self.actionWall:
            select = "wall"
        elif self.sender() == self.actionSymmetry:
            select = "symmetry"

        for row in range(self.modelLocalization.rowCount()):
            [label, nature] = self.modelLocalization.getItem(row)
            if nature == select:
                self.tableView.selectRow(row)

        self.tableView.setSelectionMode(previous_selecion_mode)


# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
