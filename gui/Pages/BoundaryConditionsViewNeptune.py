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
This module defines the 'Main fields boundary conditions' page.

This module contains the following classes:
- StandardItemModelBoundaries
- StandardItemModelMainFields
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
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, to_qvariant
from code_saturne.model.LocalizationModelNeptune import LocalizationModel
from code_saturne.model.LocalizationModel import Zone
from BoundaryConditionsNeptune import Ui_BoundaryConditions
from code_saturne.model.BoundaryNeptune import *
from code_saturne.model.BoundaryConditionsModelNeptune import *
from code_saturne.model.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# StandarItemModel class to display boundaries in a QTreeView
#-------------------------------------------------------------------------------

class StandardItemModelBoundaries(QStandardItemModel):
    def __init__(self):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"),
                        self.tr("Zone"),
                        self.tr("Nature"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))
        self.dataBoundary = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            return to_qvariant(self.dataBoundary[index.row()][index.column()])
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        self.dataChanged.emit(index, index)
        return True


    def insertItem(self, label, codeNumber, var_nature, local):
        line = [label, codeNumber, var_nature, local]
        self.dataBoundary.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataBoundary[row]


#-------------------------------------------------------------------------------
# StandardItemModelMainFields class
#-------------------------------------------------------------------------------

class StandardItemModelMainFields(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Field\nlabel"),
                         self.tr("Phase of\nfield"),
                         self.tr("Interfacial criterion")]

        self.setColumnCount(len(self.headers))

        self.tooltip = []

        self._data = []
        self.mdl = mdl


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

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, fieldId):
        """
        Add/load a field in the model.
        """
        row = self.rowCount()

        label        = self.mdl.getLabel(fieldId)
        nature       = self.mdl.getFieldNature(fieldId)
        criterion    = self.mdl.getCriterion(fieldId)

        field = [label, nature, criterion]

        self._data.append(field)
        self.setRowCount(row+1)


#-------------------------------------------------------------------------------
# BoundaryConditionsView class
#-------------------------------------------------------------------------------

class BoundaryConditionsView(QWidget, Ui_BoundaryConditions):
    """
    Boundary conditions layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditions.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = BoundaryConditionsModel(self.case)

        # Model and QTreeView for Boundaries

        self.__modelBoundaries = StandardItemModelBoundaries()
        self.treeViewBoundaries.setModel(self.__modelBoundaries)

        # Fill the model with the boundary zone
        list = ('wall', 'inlet', 'outlet')

        d = LocalizationModel('BoundaryZone', self.case)
        for zone in d.getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            codeNumber = zone.getCodeNumber()
            local = zone.getLocalization()
            if nature in list:
                self.__modelBoundaries.insertItem(label, codeNumber, nature, local)


        # Main fields definition
        self.tableModelFields = StandardItemModelMainFields(self.mdl)
        self.tableViewFields.setModel(self.tableModelFields)
        self.tableViewFields.resizeColumnsToContents()
        self.tableViewFields.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewFields.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFields.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFields.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFields.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFields.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFields.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewFields.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewFields.setSelectionMode(QAbstractItemView.SingleSelection)

        for fieldId in self.mdl.getFieldIdList():
            self.tableModelFields.newItem(fieldId)

        self.__currentField = -1
        self.__nature = ""
        self.__label = ""

        # Connect signals to slots
        self.treeViewBoundaries.clicked[QModelIndex].connect(self.__slotSelectBoundary)
        self.tableViewFields.clicked[QModelIndex].connect(self.__slotSelectField)

        self.__hideAllWidgets()
        self.tableViewFields.hide()

        self.case.undoStartGlobal()


    @pyqtSlot("QModelIndex")
    def __slotSelectField(self, index):
        """
        Select a Field in the QTable
        """

        self.__currentField = self.tableViewFields.currentIndex().row() + 1
        log.debug("slotSelectField current field %s" %str(self.__currentField))

        self.__hideAllWidgets()
        boundary = Boundary(self.__nature, self.__label, self.case, self.__currentField)

        if self.__nature == 'wall':
            self.__selectWallBoundary(boundary)
            self.tableViewFields.hide()
            if len(MainFieldsModel(self.case).getSolidFieldIdList()) > 0:
                self.tableViewFields.show()
        elif self.__nature == 'inlet':
            self.__selectInletBoundary(boundary)
            self.tableViewFields.show()
        elif self.__nature == 'outlet':
            self.__selectOutletBoundary(boundary)
            self.tableViewFields.show()


    @pyqtSlot("QModelIndex")
    def __slotSelectBoundary(self, index):
        """
        Select a boundary in the QTreeView.
        """
        label, codeNumber, nature, local = self.__modelBoundaries.getItem(index.row())
        log.debug("slotSelectBoundary label %s (%s)" % (label, nature))

        self.__nature = nature
        self.__label = label

        self.__hideAllWidgets()

        if (self.__nature == 'inlet' or self.__nature == 'outlet') :
            boundary = Boundary(self.__nature, label, self.case, self.__currentField)
            self.tableViewFields.show()

            if self.__currentField > 0 :
                if self.__nature == 'inlet':
                    self.__selectInletBoundary(boundary)
                elif self.__nature == 'outlet':
                    self.__selectOutletBoundary(boundary)
        elif self.__nature == 'wall':
            boundary = Boundary(self.__nature, label, self.case, self.__currentField)

            self.tableViewFields.hide()
            if len(MainFieldsModel(self.case).getSolidFieldIdList()) > 0:
                self.tableViewFields.show()
            self.__selectWallBoundary(boundary)


    def __selectInletBoundary(self, boundary):
        """
        Shows widgets for inlet.
        """
        self.PressureWidget.hideWidget()
        self.VelocityWidget.setup(self.case, self.__currentField)
        self.VelocityWidget.showWidget(boundary)
        self.TurbulenceWidget.setup(self.case, self.__currentField)
        self.TurbulenceWidget.showWidget(boundary)
        self.EnergyWidget.setup(self.case, self.__currentField)
        self.EnergyWidget.showWidget(boundary)
        self.FractionWidget.setup(self.case, self.__currentField)
        self.FractionWidget.showWidget(boundary)
        self.NonCondensableWidget.setup(self.case, self.__currentField)
        self.NonCondensableWidget.showWidget(boundary)
        self.InterfacialAreaWidget.setup(self.case, self.__currentField)
        self.InterfacialAreaWidget.showWidget(boundary)
        self.ScalarWidget.setup(self.case, self.__currentField)
        self.ScalarWidget.showWidget(boundary)
        self.WallWidget.hideWidget()


    def __selectWallBoundary(self, boundary):
        """
        Shows widgets for wall.
        """
        self.PressureWidget.hideWidget()
        self.VelocityWidget.hideWidget()
        self.TurbulenceWidget.hideWidget()
        self.EnergyWidget.setup(self.case, "none")
        self.EnergyWidget.showWidget(boundary)
        self.FractionWidget.hideWidget()
        self.NonCondensableWidget.hideWidget()
        self.InterfacialAreaWidget.hideWidget()
        self.ScalarWidget.hideWidget()
        if self.__currentField > 0 :
            self.WallWidget.setup(self.case, self.__currentField)
            self.WallWidget.showWidget(boundary)


    def __selectOutletBoundary(self, boundary):
        """
        Shows widgets for wall.
        """
        self.PressureWidget.setup(self.case)
        self.PressureWidget.showWidget(boundary)
        self.VelocityWidget.hideWidget()
        self.TurbulenceWidget.hideWidget()
        self.EnergyWidget.setup(self.case, self.__currentField)
        self.EnergyWidget.showWidget(boundary)
        self.FractionWidget.setup(self.case, self.__currentField)
        self.FractionWidget.showWidget(boundary)
        self.NonCondensableWidget.setup(self.case, self.__currentField)
        self.NonCondensableWidget.showWidget(boundary)
        self.InterfacialAreaWidget.hideWidget()
        self.ScalarWidget.setup(self.case, self.__currentField)
        self.ScalarWidget.showWidget(boundary)
        self.WallWidget.hideWidget()


    def __hideAllWidgets(self):
        """
        Hides all promoted QWidgets, for all nature.
        """
        self.PressureWidget.hideWidget()
        self.VelocityWidget.hideWidget()
        self.TurbulenceWidget.hideWidget()
        self.EnergyWidget.hideWidget()
        self.FractionWidget.hideWidget()
        self.NonCondensableWidget.hideWidget()
        self.InterfacialAreaWidget.hideWidget()
        self.ScalarWidget.hideWidget()
        self.WallWidget.hideWidget()


    def tr(self, text):
        """
        Translation
        """
        return text
