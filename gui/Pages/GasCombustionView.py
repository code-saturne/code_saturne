# -*- coding: iso-8859-1 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
This module contains the following classes and function:
- GasCombustionView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging, os

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------
import sys
if sys.version_info[0] == 2:
    import sip
    sip.setapi('QString', 2)

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from GasCombustionForm import Ui_GasCombustionForm
import Base.QtPage as QtPage
from Pages.GasCombustionModel import GasCombustionModel
from Base.QtPage import setGreenColor

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("GasCombustionView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class GasCombustionView(QWidget, Ui_GasCombustionForm):
    """
    Class to open the Gas Combustion option Page.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_GasCombustionForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = GasCombustionModel(self.case)

        # Set models and number of elements for combo boxes
        self.modelGasCombustionOption = QtPage.ComboModel(self.comboBoxGasCombustionOption,1,1)

        # Connections
        self.connect(self.comboBoxGasCombustionOption, SIGNAL("activated(const QString&)"), self.slotGasCombustionOption)
        self.connect(self.pushButtonThermochemistryData, SIGNAL("pressed()"), self.__slotSearchThermochemistryData)

        # Initialize Widgets
        model = self.mdl.getGasCombustionModel()

        if model == 'd3p':
            self.modelGasCombustionOption.addItem(self.tr("adiabatic model"), "adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("non adiabatic model"), "extended")
        elif model == 'ebu':
            self.modelGasCombustionOption.addItem(self.tr("reference Spalding model"), "spalding")
            self.modelGasCombustionOption.addItem(self.tr("extended model with enthalpy source term"), "enthalpy_st")
            self.modelGasCombustionOption.addItem(self.tr("extended model with mixture fraction transport"), "mixture_st")
            self.modelGasCombustionOption.addItem(self.tr("extended model with enthalpy and mixture fraction transport"), "enthalpy_mixture_st")
        elif model == 'lwp':
            self.modelGasCombustionOption.addItem(self.tr("reference two-peak model with adiabatic condition"), "2-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference two-peak model with enthalpy source term"), "2-peak_enthalpy")
            self.modelGasCombustionOption.addItem(self.tr("reference three-peak model with adiabatic condition"), "3-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference three-peak model with enthalpy source term"), "3-peak_enthalpy")
            self.modelGasCombustionOption.addItem(self.tr("reference four-peak model with adiabatic condition"), "4-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference four-peak model with enthalpy source term"), "4-peak_enthalpy")

        option = self.mdl.getGasCombustionOption()
        self.modelGasCombustionOption.setItem(str_model= option)

        name = self.mdl.getThermoChemistryDataFileName()
        if name != None:
            self.labelThermochemistryFile.setText(str(name))
            setGreenColor(self.pushButtonThermochemistryData, False)
        else:
            setGreenColor(self.pushButtonThermochemistryData, True)

        self.case.undoStartGlobal()


    @pyqtSignature("const QString&")
    def slotGasCombustionOption(self, text):
        """
        Private slot.
        Binding method for gas combustion models.
        """
        option = self.modelGasCombustionOption.dicoV2M[str(text)]
        self.mdl.setGasCombustionOption(option)


    @pyqtSignature("")
    def __slotSearchThermochemistryData(self):
        """
        Select a properties file of data for electric arc
        """
        data = self.case['data_path']
        title = self.tr("Thermochemistry file of data.")
        filetypes = self.tr("Thermochemistry (*dp_*);;All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, data, filetypes)
        file = str(file)
        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(data):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the DATA directory")
            QMessageBox.information(self, title, msg)
        else:
            self.labelThermochemistryFile.setText(str(file))
            self.mdl.setThermoChemistryDataFileName(file)
            setGreenColor(self.pushButtonThermochemistryData, False)


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
