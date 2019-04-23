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
This module contains the following classes:
- StandardItemModelBoundaries
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
from code_saturne.model.ThermalRadiationModel import ThermalRadiationModel

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.Base.QtPage import from_qvariant
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary

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
            return None
        if role == Qt.DisplayRole:
            if index.column() == 0:
                return self.lst[index.row()][1]
            elif index.column() == 1:
                key = self.lst[index.row()][3]
                return self.dataScalars[key]
            elif index.column() == 2:
                return self.lst[index.row()][2]
        return None


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


    def getListVariablesForCondition(self):
        """
        Get list of variables for condition choosed
        """
        cond = self.bdModel.getRadiativeChoice()

        if cond == 'itpimp':
            lst = [(0, self.tr("Emissivity"), '',  'EPSP',  'emissivity'),
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

        validatorEmissivity = DoubleValidator(self.lineEditEmissivity, min=0.0)
        self.lineEditEmissivity.setValidator(validatorEmissivity)

        validatorConductivity = DoubleValidator(self.lineEditConductivity, min=0.0)
        self.lineEditConductivity.setValidator(validatorConductivity)

        validatorThickness = DoubleValidator(self.lineEditThickness, min=0.0)
        self.lineEditThickness.setValidator(validatorThickness)

        validatorExtTemperature = DoubleValidator(self.lineEditExtTemperature, min=0.0)
        self.lineEditExtTemperature.setValidator(validatorExtTemperature)

        validatorIntTemperature = DoubleValidator(self.lineEditIntTemperature, min=0.0)
        self.lineEditIntTemperature.setValidator(validatorIntTemperature)

        validatorConductionFlux = DoubleValidator(self.lineEditConductionFlux, min=0.0)
        self.lineEditConductionFlux.setValidator(validatorConductionFlux)


    def __updateView__(self):
        cond = self.__boundary.getRadiativeChoice()

        #self.labelEmissivity.show()
        #self.lineEditEmissivity.show()
        self.lineEditEmissivity.setText(str(self.__boundary.getEmissivity()))

        #self.labelIntTemperature.hide()
        #self.lineEditIntTemperature.hide()
        #self.labelIntTemperatureUnit.hide()
        self.lineEditIntTemperature.setText(str(self.__boundary.getInternalTemperatureProfile()))

        self.labelConductivity.hide()
        self.lineEditConductivity.hide()
        self.labelConductivityUnit.hide()
        self.lineEditConductivity.setText(str(self.__boundary.getThermalConductivity()))

        self.labelThickness.hide()
        self.lineEditThickness.hide()
        self.labelThicknessUnit.hide()
        self.lineEditThickness.setText(str(self.__boundary.getThickness()))

        self.labelExtTemperature.hide()
        self.lineEditExtTemperature.hide()
        self.labelExtTemperatureUnit.hide()
        self.lineEditExtTemperature.setText(str(self.__boundary.getExternalTemperatureProfile()))

        self.labelConductionFlux.hide()
        self.lineEditConductionFlux.hide()
        self.labelConductionFluxUnit.hide()
        self.lineEditConductionFlux.setText(str(self.__boundary.getFlux()))

        if cond == 'ipgrno':
            self.labelConductivity.show()
            self.lineEditConductivity.show()
            self.labelConductivityUnit.show()

            self.labelThickness.show()
            self.lineEditThickness.show()
            self.labelThicknessUnit.show()

            self.labelExtTemperature.show()
            self.lineEditExtTemperature.show()
            self.labelExtTemperatureUnit.show()
        elif cond == 'ifgrno':
            self.labelConductionFlux.show()
            self.lineEditConductionFlux.show()
            self.labelConductionFluxUnit.show()


    def setup(self, case):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None

        self.case.undoStopGlobal()

        # Create the Page layout.

        # Combo
        self.modelRadiative = ComboModel(self.comboBoxRadiative,3,1)
        self.modelRadiative.addItem(self.tr("Fixed interior temperature"), 'itpimp')
        self.modelRadiative.addItem(self.tr("Fixed exterior temperature"), 'ipgrno')
        self.modelRadiative.addItem(self.tr("Fixed conduction flux"), 'ifgrno')

        # Connections
        self.comboBoxRadiative.activated[str].connect(self.slotRadiativeChoice)

        self.lineEditEmissivity.textChanged[str].connect(self.slotEmissivity)
        self.lineEditConductivity.textChanged[str].connect(self.slotConductivity)
        self.lineEditThickness.textChanged[str].connect(self.slotThickness)
        self.lineEditExtTemperature.textChanged[str].connect(self.slotExtTemperature)
        self.lineEditIntTemperature.textChanged[str].connect(self.slotIntTemperature)
        self.lineEditConductionFlux.textChanged[str].connect(self.slotConductionFlux)

        self.case.undoStartGlobal()


    def showWidget(self, b):
        """
        Show the widget
        """
        if ThermalRadiationModel(self.case).getRadiativeModel() != "off":
            label = b.getLabel()
            self.__boundary = Boundary('radiative_wall', label, self.case)
            choice = self.__boundary.getRadiativeChoice()
            self.modelRadiative.setItem(str_model=choice)
            self.__updateView__()
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
        self.__updateView__()


    @pyqtSlot(str)
    def slotEmissivity(self, text):
        """
        """
        if self.lineEditEmissivity.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setEmissivity(c)


    @pyqtSlot(str)
    def slotConductivity(self, text):
        """
        """
        if self.lineEditConductivity.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setThermalConductivity(c)


    @pyqtSlot(str)
    def slotThickness(self, text):
        """
        """
        if self.lineEditThickness.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setThickness(c)


    @pyqtSlot(str)
    def slotExtTemperature(self, text):
        """
        """
        if self.lineEditExtTemperature.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setExternalTemperatureProfile(c)


    @pyqtSlot(str)
    def slotIntTemperature(self, text):
        """
        """
        if self.lineEditIntTemperature.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setInternalTemperatureProfile(c)


    @pyqtSlot(str)
    def slotConductionFlux(self, text):
        """
        """
        if self.lineEditConductionFlux.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.__boundary.setFlux(c)


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
