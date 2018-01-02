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
This module contains the following classes:
- StandardItemModelBoundaries
- StandardItemModelScalars
- RadiativeBoundariesView
"""

#-------------------------------------------------------------------------------
# Standard modules
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

from code_saturne.Pages.BoundaryConditionsWallRadiativeTransferForm import \
     Ui_BoundaryConditionsWallRadiativeTransferForm
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsWallRadiativeTransferView")

#-------------------------------------------------------------------------------
# StandarItemModel class to display scalars properties
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):
    def __init__(self, bdModel):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Wall radiative\ncaracteristics"),
                        self.tr("Value"),
                        self.tr("Unit")]
        self.setColumnCount(len(self.headers))
        self.bdModel = bdModel
        self.lst = self.getListVariablesForCondition()
        self.setRowCount(len(self.lst))
        log.debug("StandardItemModelScalars.__init__  lst = %s " % str(self.lst))

        self.dataScalars = {}
        self.dataScalars["EPSP"]  = bdModel.getEmissivity()
        self.dataScalars["XLAMP"] = bdModel.getThermalConductivity()
        self.dataScalars["EPAP"]  = bdModel.getThickness()
        self.dataScalars["TEXTP"] = bdModel.getExternalTemperatureProfile()
        self.dataScalars["TINTP"] = bdModel.getInternalTemperatureProfile()
        self.dataScalars["FLUX"]  = bdModel.getFlux()


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            if index.column() == 0:
                return to_qvariant(self.lst[index.row()][1])
            elif index.column() == 1:
                key = self.lst[index.row()][3]
                return to_qvariant(self.dataScalars[key])
            elif index.column() == 2:
                return to_qvariant(self.lst[index.row()][2])
        if role == Qt.ToolTipRole:
            kword = self.lst[index.row()][3]
            return to_qvariant(self.tr("Code_Saturne keyword: " + kword))
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        if index.column() == 1:
            row = index.row()
            key = self.lst[row][3]
            tag = self.lst[row][4]
            val = from_qvariant(value, float)
            self.bdModel.setValRay(val, tag)
            self.dataScalars[key] = val
        self.dataChanged.emit(index, index)
        return True


    def deleteAll(self):
        self.dataScalars = []
        self.setRowCount(0)


    def insertItem(self):
        self.dataScalars.append()
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataScalars[row]


    def getListVariablesForCondition(self):
        """
        Get list of variables for condition choosed
        """
        cond = self.bdModel.getRadiativeChoice()

        if cond == 'itpimp':
            lst = [(0, self.tr("Emissivite"), '',  'EPSP',  'emissivity'),
                   (1, self.tr("Initial temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
        if cond == 'ipgrno':
            lst = [(0, self.tr("Emissivity"), '',  'EPSP',  'emissivity'),
                   (1, self.tr("Conductivity"), 'W/m/K', 'XLAMP', 'wall_thermal_conductivity'),
                   (2, self.tr("Thickness"), 'm', 'EPAP' , 'thickness'),
                   (3, self.tr("Profile of external temperature"), 'K', 'TEXTP', 'external_temperature_profile'),
                   (4, self.tr("Profile of internal temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
        if cond == 'ifgrno':
            lst = [(0, self.tr("Emissivity"),'', 'EPSP', 'emissivity'),
                   (1, self.tr("Flux of conduction"), 'W/m2', 'FLUX',  'flux'),
                   (2, self.tr("Inital temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
        return lst

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsWallRadiativeTransferView(QWidget,
                                                  Ui_BoundaryConditionsWallRadiativeTransferForm):
    """
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsWallRadiativeTransferForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        # Create the Page layout.

        # Combo
        self.modelRadiative = ComboModel(self.comboBoxRadiative,3,1)
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and profile of fixed internal temperature"), 'itpimp')
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and profile of fixed external temperature"), 'ipgrno')
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and flux of fixed conduction"), 'ifgrno')

        # Connections
        self.comboBoxRadiative.activated[str].connect(self.slotRadiativeChoice)

        self.__case.undoStartGlobal()


    def showWidget(self, b):
        """
        Show the widget
        """
        if ThermalRadiationModel(self.__case).getRadiativeModel() != "off":
            label = b.getLabel()
            self.__boundary = Boundary('radiative_wall', label, self.__case)
            choice = self.__boundary.getRadiativeChoice()
            self.modelRadiative.setItem(str_model=choice)

            if hasattr(self, "modelScalars"):
                del self.modelScalars
            self.modelScalars = StandardItemModelScalars(self.__boundary)
            self.tableViewScalars.setModel(self.modelScalars)

            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all the widget
        """
        self.hide()


    @pyqtSlot(str)
    def slotRadiativeChoice(self, text):
        cond = self.modelRadiative.dicoV2M[str(text)]
        log.debug("slotRadiativeChoice cond = %s "%cond)
        self.__boundary.setRadiativeChoice(cond)
        self.modelScalars.deleteAll()
        self.modelScalars = StandardItemModelScalars(self.__boundary)
        self.tableViewScalars.setModel(self.modelScalars)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
