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
This module defines the 'Nucleate boiling' page.

This module contains the following classes:
- NucleateBoilingView
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
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from NucleateBoiling import Ui_NucleateBoiling
from code_saturne.model.NucleateBoilingModel import NucleateBoilingModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NucleateBoilingView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
#  class InterfacialForces
#-------------------------------------------------------------------------------

class NucleateBoilingView(QWidget, Ui_NucleateBoiling):
    """
    Nucleate boiling model layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NucleateBoiling.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = NucleateBoilingModel(self.case)

        self.modelHeatTransferModel = ComboModel(self.comboBoxHeatTransferModel, 2, 1)
        self.modelHeatTransferModel.addItem(self.tr("Extended Kurul-Podowski model"),"extended_kurul-podowski")
        self.modelHeatTransferModel.addItem(self.tr("Standard Kurul-Podowski model"),"standard_kurul-podowski")

        self.modelWallFunctionModel = ComboModel(self.comboBoxWallFunctionModel, 3, 1)
        self.modelWallFunctionModel.addItem(self.tr("standard (single phase wall function)"), "standard")
        self.modelWallFunctionModel.addItem(self.tr("Koncar Tiselj-NED 2008"), "koncar")
        self.modelWallFunctionModel.addItem(self.tr("Mimouni et al-NED 2008"), "mimouni")

        self.modelYPlus = ComboModel(self.comboBoxYPlus, 3, 1)
        self.modelYPlus.addItem(self.tr("Boundary cell center"), "center")
        self.modelYPlus.addItem(self.tr("Y+ = "), "Yplus_value")
        self.modelYPlus.addItem(self.tr("Nucleate bubble diameter"), "diameter")

        # Validators

        validatorYplus = DoubleValidator(self.lineEditYPlus, min = 0.0)
        validatorRad   = DoubleValidator(self.lineEditMaxRadius, min = 0.0)
        validatorDiam  = DoubleValidator(self.lineEditMaxDiam, min = 0.0)
        validatorSat   = DoubleValidator(self.lineEditMaxOverSaturation, min = 0.0)
        validatorLam   = DoubleValidator(self.lineEditThermalConductivity, min = 0.0)
        validatorRho   = DoubleValidator(self.lineEditDensity, min = 0.0)
        validatorCp    = DoubleValidator(self.lineEditSpecificHeat, min = 0.0)
        validatorTh    = DoubleValidator(self.lineEditThickness, min = 0.0)

        validatorYplus.setExclusiveMin(True)
        validatorRad.setExclusiveMin(True)
        validatorDiam.setExclusiveMin(True)
        validatorSat.setExclusiveMin(True)
        validatorLam.setExclusiveMin(True)
        validatorRho.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorTh.setExclusiveMin(True)

        self.lineEditYPlus.setValidator(validatorYplus)
        self.lineEditMaxRadius.setValidator(validatorRad)
        self.lineEditMaxDiam.setValidator(validatorDiam)
        self.lineEditMaxOverSaturation.setValidator(validatorSat)
        self.lineEditThermalConductivity.setValidator(validatorLam)
        self.lineEditDensity.setValidator(validatorRho)
        self.lineEditSpecificHeat.setValidator(validatorCp)
        self.lineEditThickness.setValidator(validatorTh)


        # Connect signals to slots
        self.comboBoxHeatTransferModel.activated[str].connect(self.slotHeatTransferModel)
        self.comboBoxWallFunctionModel.activated[str].connect(self.slotWallFunctionModel)
        self.comboBoxYPlus.activated[str].connect(self.slotYPlus)
        self.checkBoxThickness.clicked.connect(self.slotThickness)
        self.lineEditYPlus.textChanged[str].connect(self.slotYPlusValue)
        self.lineEditMaxRadius.textChanged[str].connect(self.slotMaxRadius)
        self.lineEditMaxDiam.textChanged[str].connect(self.slotMaxDiam)
        self.lineEditMaxOverSaturation.textChanged[str].connect(self.slotMaxOverSaturation)
        self.lineEditThermalConductivity.textChanged[str].connect(self.slotThermalConductivity)
        self.lineEditDensity.textChanged[str].connect(self.slotDensity)
        self.lineEditSpecificHeat.textChanged[str].connect(self.slotSpecificHeat)
        self.lineEditThickness.textChanged[str].connect(self.slotThicknessValue)

        # load values
        isYPlus = self.mdl.getYPlusModel()
        self.modelYPlus.setItem(str_model=isYPlus)

        if isYPlus == "Yplus_value" :
           self.lineEditYPlus.show()
           self.lineEditYPlus.setText(str(self.mdl.getYPlusValue()))
        else :
           self.lineEditYPlus.hide()

        self.lineEditMaxRadius.setText(str(self.mdl.getMaxRadius()))
        self.lineEditMaxDiam.setText(str(self.mdl.getMaxDiameter()))
        self.lineEditMaxOverSaturation.setText(str(self.mdl.getMaxOverSaturation()))
        self.lineEditThermalConductivity.setText(str(self.mdl.getThermalConductivity()))
        self.lineEditDensity.setText(str(self.mdl.getDensity()))
        self.lineEditSpecificHeat.setText(str(self.mdl.getSpecificHeat()))

        model = self.mdl.getHeatTransferModel()
        self.modelHeatTransferModel.setItem(str_model=model)
        if model == "standard_kurul-podowski" :
            self.labelMaxRadius.setEnabled(0)
            self.lineEditMaxRadius.setEnabled(0)
            self.labelMaxRadiusUnit.setEnabled(0)
            self.labelMaxDiam.setEnabled(0)
            self.lineEditMaxDiam.setEnabled(0)
            self.labelMaxDiamUnit.setEnabled(0)
        else :
            self.labelMaxRadius.setEnabled(1)
            self.lineEditMaxRadius.setEnabled(1)
            self.labelMaxRadiusUnit.setEnabled(1)
            self.labelMaxDiam.setEnabled(1)
            self.lineEditMaxDiam.setEnabled(1)
            self.labelMaxDiamUnit.setEnabled(1)

        isThickness = self.mdl.getThicknessStatus() == "on"
        self.checkBoxThickness.setChecked(isThickness)

        if isThickness :
            self.lineEditThickness.show()
            self.labelThickness1.show()
            self.labelThickness2.show()
            self.lineEditThickness.setText(str(self.mdl.getThicknessValue()))
        else :
            self.lineEditThickness.hide()
            self.labelThickness1.hide()
            self.labelThickness2.hide()

        model = self.mdl.getWallFunctionModel()
        self.modelWallFunctionModel.setItem(str_model=model)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotHeatTransferModel(self, text):
        """
        configure standard or extend kurul-podowski model
        """
        value = self.modelHeatTransferModel.dicoV2M[text]
        log.debug("slotHeatTransferModel -> %s" % value)
        self.mdl.setHeatTransferModel(value)
        if value == "standard_kurul-podowski" :
            self.labelMaxRadius.setEnabled(0)
            self.lineEditMaxRadius.setEnabled(0)
            self.labelMaxRadiusUnit.setEnabled(0)
            self.labelMaxDiam.setEnabled(0)
            self.lineEditMaxDiam.setEnabled(0)
            self.labelMaxDiamUnit.setEnabled(0)
        else :
            self.labelMaxRadius.setEnabled(1)
            self.lineEditMaxRadius.setEnabled(1)
            self.labelMaxRadiusUnit.setEnabled(1)
            self.labelMaxDiam.setEnabled(1)
            self.lineEditMaxDiam.setEnabled(1)
            self.labelMaxDiamUnit.setEnabled(1)
            self.lineEditMaxRadius.setText(str(self.mdl.getMaxRadius()))
            self.lineEditMaxDiam.setText(str(self.mdl.getMaxDiameter()))


    @pyqtSlot(str)
    def slotWallFunctionModel(self, text):
        """
        configure wall function model
        """
        value = self.modelWallFunctionModel.dicoV2M[text]
        log.debug("slotWallFunctionModel -> %s" % value)
        self.mdl.setWallFunctionModel(value)


    @pyqtSlot(str)
    def slotYPlus(self, text):
        """
        configure Y Plus model
        """
        value = self.modelYPlus.dicoV2M[text]
        log.debug("slotYPlus -> %s" % value)
        self.mdl.setYPlusModel(value)

        if value == "Yplus_value" :
           self.lineEditYPlus.show()
           self.lineEditYPlus.setText(str(self.mdl.getYPlusValue()))
        else :
           self.lineEditYPlus.hide()


    @pyqtSlot(str)
    def slotYPlusValue(self, text):
        """
        Update the Yplus value
        """
        if self.lineEditYPlus.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setYPlusValue(value)


    @pyqtSlot(str)
    def slotMaxRadius(self, text):
        """
        Update the max radius
        """
        if self.lineEditMaxRadius.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMaxRadius(value)


    @pyqtSlot(str)
    def slotMaxDiam(self, text):
        """
        Update the max diameter
        """
        if self.lineEditMaxDiam.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMaxDiameter(value)


    @pyqtSlot(str)
    def slotMaxOverSaturation(self, text):
        """
        Update the maximum oversaturation temperature
        """
        if self.lineEditMaxOverSaturation.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMaxOverSaturation(value)


    @pyqtSlot(str)
    def slotThermalConductivity(self, text):
        """
        Update the thermal conductivity
        """
        if self.lineEditThermalConductivity.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setThermalConductivity(value)


    @pyqtSlot(str)
    def slotDensity(self, text):
        """
        Update the density
        """
        if self.lineEditDensity.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setDensity(value)


    @pyqtSlot(str)
    def slotSpecificHeat(self, text):
        """
        Update the specific heat
        """
        if self.lineEditSpecificHeat.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setSpecificHeat(value)


    @pyqtSlot(bool)
    def slotThickness(self, checked):
        """
        check box for Y Plus
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setThicknessStatus(status)

        if status == 'on' :
            self.lineEditThickness.show()
            self.labelThickness1.show()
            self.labelThickness2.show()
            self.lineEditThickness.setText(str(self.mdl.getThicknessValue()))
        else :
            self.lineEditThickness.hide()
            self.labelThickness1.hide()
            self.labelThickness2.hide()


    @pyqtSlot(str)
    def slotThicknessValue(self, text):
        """
        Update the thickness value
        """
        if self.lineEditThickness.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setThicknessValue(value)


