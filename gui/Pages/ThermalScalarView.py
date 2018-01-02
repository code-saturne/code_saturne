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
This module defines the thermal scalar management.

This module contains the following classes and function:
- ThermalScalarView
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
from code_saturne.Base.QtPage import ComboModel
from code_saturne.Pages.ThermalScalarForm import Ui_ThermalScalarForm
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.ElectricalModel import ElectricalModel
from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
from code_saturne.Pages.GasCombustionModel import GasCombustionModel
from code_saturne.Pages.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.Pages.CompressibleModel import CompressibleModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermalScalarView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ThermalScalarView(QWidget, Ui_ThermalScalarForm):
    """
    Class to open Thermal Scalar Transport Page.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ThermalScalarForm.__init__(self)
        self.setupUi(self)

        self.browser = tree
        self.case = case
        self.case.undoStopGlobal()
        self.thermal = ThermalScalarModel(self.case)

        # combo Model

        self.modelThermal = ComboModel(self.comboBoxThermal, 4, 1)

        self.modelThermal.addItem(self.tr("No thermal scalar"), 'off')
        self.modelThermal.addItem(self.tr("Temperature (Celsius)"), 'temperature_celsius')
        self.modelThermal.addItem(self.tr("Temperature (Kelvin)"), 'temperature_kelvin')
        self.modelThermal.addItem(self.tr("Enthalpy (J/kg)"), 'enthalpy')

        self.modelThermal.addItem(self.tr("Potential temperature"), 'potential_temperature')
        self.modelThermal.addItem(self.tr("Liquid potential temperature"), 'liquid_potential_temperature')
        self.modelThermal.addItem(self.tr("Total energy (J/kg)"), 'total_energy')

        self.comboBoxThermal.activated[str].connect(self.slotThermalScalar)

        # Update the thermal scalar list with the calculation features

        for sca in self.thermal.thermalModel:
            if sca not in self.thermal.thermalScalarModelsList():
                self.modelThermal.disableItem(str_model=sca)

        if ElectricalModel(self.case).getElectricalModel() != 'off':
            self.comboBoxThermal.setEnabled(False)

        if CoalCombustionModel(self.case).getCoalCombustionModel() != 'off':
            self.comboBoxThermal.setEnabled(False)

        if GasCombustionModel(self.case).getGasCombustionModel() != 'off':
            self.comboBoxThermal.setEnabled(False)

        if CompressibleModel(self.case).getCompressibleModel() != 'off':
            self.comboBoxThermal.setEnabled(False)
        else:
            self.modelThermal.delItem(6)

        if AtmosphericFlowsModel(self.case).getAtmosphericFlowsModel() != 'off':
            self.comboBoxThermal.setEnabled(False)
        else:
            self.modelThermal.delItem(5)
            self.modelThermal.delItem(4)

        # Select the thermal scalar model

        model = self.thermal.getThermalScalarModel()
        self.modelThermal.setItem(str_model=model)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotThermalScalar(self, text):
        """
        Update the thermal scalar markup.
        """
        th = self.modelThermal.dicoV2M[str(text)]
        self.thermal.setThermalModel(th)

        self.browser.configureTree(self.case)


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
