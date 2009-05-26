# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from ThermalScalarForm import Ui_ThermalScalarForm
import Base.QtPage as QtPage
from ThermalScalarModel import ThermalScalarModel

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
        self.thermal = ThermalScalarModel(self.case)

        # combo Model

        self.modelThermal = QtPage.ComboModel(self.comboBoxThermal, 4, 1)
        self.modelThermal.addItem(self.tr("No thermal scalar"), 'off')
        self.modelThermal.addItem(self.tr("Temperature (Celsius)"), 'temperature_celsius')
        self.modelThermal.addItem(self.tr("Temperature (Kelvin)"), 'temperature_kelvin')
        self.modelThermal.addItem(self.tr("Enthalpy (J/kg)"), 'enthalpy')

        self.connect(self.comboBoxThermal, SIGNAL("activated(const QString&)"), self.slotThermalScalar)

        # Update the thermal scalar list with the calculation features

        for sca in self.thermal.thermalModel:
            if sca not in self.thermal.thermalScalarModelsList():
                self.modelThermal.disableItem(str_model=sca)

        # Select the thermal scalar model

        model = self.thermal.getThermalScalarModel()
        self.modelThermal.setItem(str_model=model)


    @pyqtSignature("const QString &")
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