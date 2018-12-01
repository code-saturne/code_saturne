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
This module defines the 'Boundary regions definition' page.

This module contains the following classes:
- LocalizationView
- BoundaryConditionsView
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

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.LocalizationModel import Zone
from code_saturne.Pages.LocalizationView import *
from code_saturne.Pages.LocalizationModelNeptune import LocalizationModel
from code_saturne.Pages.LocalizationForm import Ui_LocalizationForm

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LocalizationView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelLocalization(QStandardItemModel):
    def __init__(self, mdl, zoneType, dicoM2V, tree = None, case = None):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"),
                        self.tr("Zone"),
                        self.tr("Nature"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))

        self.mdl      = mdl
        self.zoneType = zoneType
        self.dicoM2V  = dicoM2V
        self.browser = tree
        self.case = case

        self._data = []
        self._disable = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()

            if col in [0, 1, 3]:
                return to_qvariant(self._data[row][col])

            elif col == 2:
                key = self._data[row][col]
                return to_qvariant(self.dicoM2V[key])

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if (index.row(), index.column()) in self._disable:
            return Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        [old_label, old_code, old_nature, old_local] = self._data[row]

        old_zone = Zone(self.zoneType,
                        case         = self.case,
                        label        = old_label,
                        codeNumber   = old_code,
                        localization = old_local,
                        nature       = old_nature)

        new_label  = old_label
        new_code   = old_code
        new_nature = old_nature
        new_local  = old_local

        if col == 0:
            new_label = from_qvariant(value, to_text_string)
            self._data[row][col] = new_label

        elif col == 1:
            new_code = from_qvariant(value, int)
            self._data[row][col] = new_code

        elif col == 2:
            new_nature = str(from_qvariant(value, to_text_string))
            self._data[row][col] = new_nature

        elif col == 3:
            new_local = str(from_qvariant(value, to_text_string))
            self._data[row][col] = new_local

        new_zone = Zone(self.zoneType,
                        case         = self.case,
                        label        = new_label,
                        codeNumber   = new_code,
                        localization = new_local,
                        nature       = new_nature)

        self.mdl.replaceZone(old_zone, new_zone)

        self.dataChanged.emit(index, index)
        self.browser.configureTree(self.case)
        return True


    def addItem(self, zone=None):
        """
        Add an element in the table view.
        """
        if not zone:
            zone = self.mdl.addZone(Zone(self.zoneType, case = self.case))

        line = [zone.getLabel(),
                zone.getCodeNumber(),
                zone.getNature(),
                zone.getLocalization()]
        self._data.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)

        # Warning: the Volume region 'all_cells' is mandatory, and can not be removed.
        if self.zoneType == "VolumicZone":
            if zone.getLabel() == "all_cells":
                for c in [0, 3]:
                    self._disable.append((row, c))
            # self._disable.append((row, 2))
        self._disable.append((row, 1))
        self.browser.configureTree(self.case)
        return zone

    def getItem(self, row):
        return self._data[row]


    def deleteItem(self, irow):
        del self._data[irow]
        nb_rows = self.rowCount()
        self.setRowCount(nb_rows-1)


    def getData(self, row, column):
        return self._data[row][column]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class LocalizationView(QWidget, Ui_LocalizationForm):
    """
    Main class
    """
    def __init__(self, zoneType, parent, case, dicoM2V, tree = None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LocalizationForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.zoneType = zoneType

        self.mdl = LocalizationModel(zoneType, case)
        self.case['dump_python'].append([self.mdl.__module__, zoneType, ()])

        self.browser = tree

        # Model for table View
        self.modelLocalization = StandardItemModelLocalization(self.mdl, zoneType, dicoM2V, tree, case)
        self.tableView.setModel(self.modelLocalization)

        # Delegates
        delegateLabel = LabelDelegate(self.tableView, self.mdl)
        delegateCode  = CodeNumberDelegate(self.tableView, self.mdl)
        delegateLocal = LocalizationSelectorDelegate(self.tableView, self.mdl)

        self.tableView.setItemDelegateForColumn(0, delegateLabel)
        self.tableView.setItemDelegateForColumn(1, delegateCode)
        self.tableView.setItemDelegateForColumn(3, delegateLocal)

        # Populate QTableView model
        for zone in self.mdl.getZones():
            self.modelLocalization.addItem(zone)

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(3, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)

        # Connections
        self.pushButtonNew.clicked.connect(self.slotAddZone)
        self.pushButtonDelete.clicked.connect(self.slotDeleteZone)
        self.toolButtonCreation.clicked.connect(self.slotAddFromPrePro)
        self.modelLocalization.dataChanged.connect(self.dataChanged)

        self.pushButtonSalome.hide()
        if case['salome']:
            self.pushButtonSalome.show()
            self.pushButtonSalome.clicked.connect(self.slotAddFromSalome)

        # Context menu
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested[QPoint].connect(self.slotContextMenu)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddZone(self):
        """
        Insert a new item in the table view.
        """
        self.modelLocalization.addItem()


    @pyqtSlot()
    def slotDeleteZone(self):
        """
        Private Slot.
        Warning: the Volume region 'all_cells' is mandatory, therefore it can not be deleted.
        """
        lst = []
        for index in self.tableView.selectionModel().selectedRows():
            row = index.row()
            lst.append(row)

        lst.sort()
        lst.reverse()

        for row in lst:
            [label, codeNumber, nature, localization] = self.modelLocalization.getItem(row)
            if not (label == "all_cells" and self.zoneType == 'VolumicZone'):
                self.mdl.deleteZone(label)
                self.modelLocalization.deleteItem(row)


    @pyqtSlot()
    def slotAddFromPrePro(self):
        """
        Research a preprocessor log to pick colors or groups of cells or faces.
        """
        entity = 'faces'

        file_name = preprocessorFile(self, self.case['resu_path'])

        if file_name:
            for loc in Informations(file_name, entity).getLocalizations():
                if loc not in self.mdl.getLocalizationsZonesList():
                    zone = Zone(self.zoneType, case = self.case, localization = loc)
                    self.mdl.addZone(zone)
                    self.modelLocalization.addItem(zone)


    @pyqtSlot()
    def slotContextMenu(self):
        """
        Public slot

        Create the popup menu of the Browser
        """
        fileMenu = QMenu(self.tableView)

        actionMerge = QAction(self.tr("Merge selected zones"), self.tableView)
        actionMerge.triggered.connect(self.slotMerge)
        fileMenu.addAction(actionMerge)

        if self.zoneType == 'BoundaryZone':
            fileMenu.addSeparator()

            self.actionInlet = QAction(self.tr("Select all inlets"), self.tableView)
            self.actionInlet.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionInlet)

            self.actionOutlet = QAction(self.tr("Select all outlets"), self.tableView)
            self.actionOutlet.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionOutlet)

            self.actionWall = QAction(self.tr("Select all walls"), self.tableView)
            self.actionWall.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionWall)

            self.actionSymmetry = QAction(self.tr("Select all symmetries"), self.tableView)
            self.actionSymmetry.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionSymmetry)

        fileMenu.popup(QCursor().pos())
        fileMenu.show()


    @pyqtSlot()
    def slotMerge(self):
        """
        public slot
        """
        list = []
        for index in self.tableView.selectionModel().selectedRows():
            list.append(index.row())

        row = list.pop(0)
        [label, code, nature, new_localization] = self.modelLocalization.getItem(row)

        new_zone = Zone(self.zoneType,
                        case = self.case,
                        label        = label,
                        codeNumber   = code,
                        localization = new_localization,
                        nature       = nature)

        for row in list:
            [label, code, nature, localization] = self.modelLocalization.getItem(row)
            if "all[]" not in string.split(new_localization, " "):
                new_localization += " or " + localization
            if localization == "all[]":
                new_localization = "all[]"

        new_zone.setLocalization(new_localization)
        self.slotDeleteZone()
        self.mdl.addZone(new_zone)
        self.modelLocalization.addItem(new_zone)


    @pyqtSlot()
    def slotAddFromSalome(self):
        """
        When GUI is embeded in the Salome desktop. Add selection criteria from
        graphical selection in the VTK viwver, or in the ObjectBrowser.
        """
        if self.case['salome']:
            from code_saturne.Pages.SalomeHandler import BoundaryGroup, VolumeGroup

            log.debug("slotAddFromSalome: zoneType -> %s" % self.zoneType)
            loc = BoundaryGroup()

            log.debug("slotAddFromSalome: selection criteria -> %s" % loc)
            if loc not in self.mdl.getLocalizationsZonesList():
                zone = Zone(self.zoneType, case = self.case, localization = loc)
                self.mdl.addZone(zone)
                self.modelLocalization.addItem(zone)


    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableView.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableView.resizeColumnToContents(col)


    def tr(self, text):
        """
        Translation.
        """
        return text


#-------------------------------------------------------------------------------
# Define Boundary regions class
#-------------------------------------------------------------------------------

class BoundaryLocalizationView(LocalizationView):
    """
    Define boundary regions class.
    """
    def __init__(self, parent=None, case=None, tree = None):
        """
        Constructor.
        """
        self.case = case
        dicoM2V = Zone("BoundaryZone", case = self.case).getModel2ViewDictionary()
        self.browser = tree

        LocalizationView.__init__(self, "BoundaryZone", parent, self.case, dicoM2V, tree)

        title = self.tr("Boundary regions definition")
        self.groupBoxLocalization.setTitle(title)

        # Delegates
        delegateNature = BoundaryNatureDelegate(self.tableView, dicoM2V)
        self.tableView.setItemDelegateForColumn(2, delegateNature)


    @pyqtSlot()
    def slotSelectBoudaries(self):
        """
        Public slot.

        Warning: works only if the selection mode of the view is set to MultiSelection.
        """
        previous_selecion_mode = self.tableView.selectionMode()
        self.tableView.setSelectionMode(QAbstractItemView.MultiSelection)
        self.tableView.clearSelection()

        if self.sender() == self.actionInlet:
            select  = "inlet"
        elif self.sender() == self.actionOutlet:
            select  = "outlet"
        elif self.sender() == self.actionWall:
            select  = "wall"
        elif self.sender() == self.actionSymmetry:
            select  = "symmetry"

        for row in range(self.modelLocalization.rowCount()):
            [label, code, nature, localization] = self.modelLocalization.getItem(row)
            if nature == select:
                self.tableView.selectRow(row)

        self.tableView.setSelectionMode(previous_selecion_mode)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
