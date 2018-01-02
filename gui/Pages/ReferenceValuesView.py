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
This module defines the values of reference.

This module contains the following classes and function:
- ReferenceValuesView
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from code_saturne.Pages.ReferenceValuesForm import Ui_ReferenceValuesForm
from code_saturne.Pages.ReferenceValuesModel import ReferenceValuesModel
from code_saturne.Pages.GasCombustionModel import GasCombustionModel
from code_saturne.Pages.CompressibleModel import CompressibleModel
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.GroundwaterModel import GroundwaterModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ReferenceValuesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ReferenceValuesView(QWidget, Ui_ReferenceValuesForm):
    """
    Class to open Reference Pressure Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ReferenceValuesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = ReferenceValuesModel(self.case)

        # Combo models
        self.modelLength = ComboModel(self.comboBoxLength,2,1)
        self.modelLength.addItem(self.tr("Automatic"), 'automatic')
        self.modelLength.addItem(self.tr("Prescribed"), 'prescribed')
        self.comboBoxLength.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        # Connections
        self.lineEditP0.textChanged[str].connect(self.slotPressure)
        self.lineEditV0.textChanged[str].connect(self.slotVelocity)
        self.comboBoxLength.activated[str].connect(self.slotLengthChoice)
        self.lineEditL0.textChanged[str].connect(self.slotLength)
        self.lineEditT0.textChanged[str].connect(self.slotTemperature)
        self.lineEditOxydant.textChanged[str].connect(self.slotTempOxydant)
        self.lineEditFuel.textChanged[str].connect(self.slotTempFuel)
        self.lineEditMassMolar.textChanged[str].connect(self.slotMassemol)

        # Validators

        validatorP0 = DoubleValidator(self.lineEditP0, min=0.0)
        self.lineEditP0.setValidator(validatorP0)

        validatorV0 = DoubleValidator(self.lineEditV0, min=0.0)
        self.lineEditV0.setValidator(validatorV0)

        validatorL0 = DoubleValidator(self.lineEditL0, min=0.0)
        self.lineEditL0.setValidator(validatorL0)

        validatorT0 = DoubleValidator(self.lineEditT0,  min=0.0)
        validatorT0.setExclusiveMin(True)
        self.lineEditT0.setValidator(validatorT0)

        validatorOxydant = DoubleValidator(self.lineEditOxydant,  min=0.0)
        validatorOxydant.setExclusiveMin(True)
        self.lineEditOxydant.setValidator(validatorOxydant)

        validatorFuel = DoubleValidator(self.lineEditFuel,  min=0.0)
        validatorFuel.setExclusiveMin(True)
        self.lineEditFuel.setValidator(validatorFuel)

        validatorMM = DoubleValidator(self.lineEditMassMolar, min=0.0)
        validatorMM.setExclusiveMin(True)
        self.lineEditMassMolar.setValidator(validatorMM)

        # Display

        model = self.mdl.getParticularPhysical()

        self.groupBoxMassMolar.hide()

        if model == "atmo":
            self.labelInfoT0.hide()
        elif model == "comp" or model == "coal":
            self.groupBoxMassMolar.show()
        elif model == "off":
            thmodel = ThermalScalarModel(self.case).getThermalScalarModel()
            if thmodel == "enthalpy":
                self.labelT0.setText("enthalpy")
                self.labelUnitT0.setText("J/kg")
                self.groupBoxTemperature.setTitle("Reference enthalpy")
            elif thmodel == "temperature_celsius":
                self.labelUnitT0.setText("C")

            if FluidCharacteristicsModel(self.case).getMaterials() != "user_material":
                self.labelInfoT0.hide()
        else:
            self.groupBoxTemperature.hide()

        gas_comb = GasCombustionModel(self.case).getGasCombustionModel()
        if gas_comb == 'd3p':
            self.groupBoxTempd3p.show()
            t_oxy  = self.mdl.getTempOxydant()
            t_fuel = self.mdl.getTempFuel()
            self.lineEditOxydant.setText(str(t_oxy))
            self.lineEditFuel.setText(str(t_fuel))
        else:
            self.groupBoxTempd3p.hide()

        # Initialization

        darc = GroundwaterModel(self.case).getGroundwaterModel()
        if darc != 'off':
            self.groupBoxPressure.hide()
        else:
            p = self.mdl.getPressure()
            self.lineEditP0.setText(str(p))

        v = self.mdl.getVelocity()
        self.lineEditV0.setText(str(v))

        init_length_choice = self.mdl.getLengthChoice()
        self.modelLength.setItem(str_model=init_length_choice)
        if init_length_choice == 'automatic':
            self.lineEditL0.setText(str())
            self.lineEditL0.hide()
            self.labelUnitL0.hide()
        else:
            self.lineEditL0.show()
            self.labelUnitL0.show()
            l = self.mdl.getLength()
            self.lineEditL0.setText(str(l))

        model = self.mdl.getParticularPhysical()
        if model == "atmo":
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(str(t))
        elif model != "off":
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(str(t))
            m = self.mdl.getMassemol()
            self.lineEditMassMolar.setText(str(m))
        else:
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(str(t))

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotPressure(self,  text):
        """
        Input PRESS.
        """
        if self.lineEditP0.validator().state == QValidator.Acceptable:
            p = from_qvariant(text, float)
            self.mdl.setPressure(p)


    @pyqtSlot(str)
    def slotVelocity(self,  text):
        """
        Input Velocity.
        """
        if self.lineEditV0.validator().state == QValidator.Acceptable:
            v = from_qvariant(text, float)
            self.mdl.setVelocity(v)


    @pyqtSlot(str)
    def slotLengthChoice(self,text):
        """
        Set value for parameterNTERUP
        """
        choice = self.modelLength.dicoV2M[str(text)]
        self.mdl.setLengthChoice(choice)
        if choice == 'automatic':
            self.lineEditL0.setText(str())
            self.lineEditL0.hide()
            self.labelUnitL0.hide()
        else:
            self.lineEditL0.show()
            self.labelUnitL0.show()
            value = self.mdl.getLength()
            self.lineEditL0.setText(str(value))
        log.debug("slotlengthchoice-> %s" % choice)


    @pyqtSlot(str)
    def slotLength(self,  text):
        """
        Input reference length.
        """
        if self.lineEditL0.validator().state == QValidator.Acceptable:
            l = from_qvariant(text, float)
            self.mdl.setLength(l)


    @pyqtSlot(str)
    def slotTemperature(self,  text):
        """
        Input TEMPERATURE.
        """
        if self.lineEditT0.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTemperature(t)


    @pyqtSlot(str)
    def slotTempOxydant(self,  text):
        """
        Input oxydant TEMPERATURE.
        """
        if self.lineEditOxydant.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTempOxydant(t)


    @pyqtSlot(str)
    def slotTempFuel(self,  text):
        """
        Input fuel TEMPERATURE.
        """
        if self.lineEditFuel.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTempFuel(t)


    @pyqtSlot(str)
    def slotMassemol(self,  text):
        """
        Input Mass molar.
        """
        if self.lineEditMassMolar.validator().state == QValidator.Acceptable:
            m = from_qvariant(text, float)
            self.mdl.setMassemol(m)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
