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
- BoundaryConditionsEnergyView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BoundaryConditionsEnergy import Ui_BoundaryConditionsEnergy

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsEnergyView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsEnergyView(QWidget, Ui_BoundaryConditionsEnergy) :
    """
    Boundary condition for energy
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsEnergy.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxEnergy.activated[str].connect(self.__slotChoiceEnergy)

        self.__modelEnergy = ComboModel(self.comboBoxEnergy, 1, 1)

        self.lineEditEnergy.textChanged[str].connect(self.__slotEnergy)

        validatorEner = DoubleValidator(self.lineEditEnergy)

        self.lineEditEnergy.setValidator(validatorEner)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None
        self.__currentField = fieldId


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        for nb in range(len(self.__modelEnergy.getItems())):
            self.__modelEnergy.delItem(0)

        if self.__boundary.getNature() == "inlet"  or  self.__boundary.getNature() == "outlet":
            self.__modelEnergy.addItem(self.tr("Imposed enthalpy"), 'dirichlet')
            self.__modelEnergy.addItem(self.tr("Thermal flux"), 'flux')
            self.__modelEnergy.addItem(self.tr("Temperature"), 'timp_K')
            self.__modelEnergy.addItem(self.tr("Saturation enthalpy"), 'hsat_P')
        elif self.__boundary.getNature() == "wall" :
            self.__modelEnergy.addItem(self.tr("Flux"), 'flux')
            self.__modelEnergy.addItem(self.tr("Temperature"), 'temperature')

        if self.__boundary.getNature() == "inlet" or  self.__boundary.getNature() == "outlet":
            if MainFieldsModel(self.__case).getEnergyResolution(self.__currentField) == 'on' :
                energychoice = self.__boundary.getEnthalpyChoice(self.__currentField)
                self.__modelEnergy.setItem(str_model=energychoice)
                if energychoice == "hsat_P" :
                    self.lineEditEnergy.hide()
                    self.labelEnergy.hide()
                else :
                    self.lineEditEnergy.show()
                    self.labelEnergy.show()
                    val = self.__boundary.getEnthalpy(self.__currentField)
                    self.lineEditEnergy.setText(str(val))
                    if energychoice == 'dirichlet' :
                        self.labelEnergy.setText('J/kg')
                    elif energychoice == 'flux' :
                        self.labelEnergy.setText('W/m2')
                    else :
                        self.labelEnergy.setText('K')
                self.show()
                if ThermodynamicsModel(self.__case).getMaterials(self.__currentField) == 'user_material' :
                    self.__modelEnergy.disableItem(2)
                    self.__modelEnergy.disableItem(3)
                else :
                    self.__modelEnergy.enableItem(2)
                    self.__modelEnergy.enableItem(3)
            else :
                self.hideWidget()
        elif self.__boundary.getNature() == "wall" :
            if len(MainFieldsModel(self.__case).getEnthalpyResolvedField()) != 0 :
                energychoice = self.__boundary.getEnthalpyChoice("none")
                self.__modelEnergy.setItem(str_model=energychoice)
                val = self.__boundary.getEnthalpy("none")
                self.lineEditEnergy.setText(str(val))
                if energychoice == 'flux' :
                    self.labelEnergy.setText('W/m2')
                else :
                    self.labelEnergy.setText('K')
                self.lineEditEnergy.show()
                self.labelEnergy.show()
                self.show()
            else :
                self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotChoiceEnergy(self, text):
        """
        INPUT choice of method of calculation of the energy
        """
        energy_choice = self.__modelEnergy.dicoV2M[str(text)]
        self.__boundary.setEnthalpyChoice(self.__currentField, energy_choice)

        if energy_choice == "hsat_P" :
            self.lineEditEnergy.hide()
            self.labelEnergy.hide()
        else :
            self.lineEditEnergy.show()
            self.labelEnergy.show()
            val = self.__boundary.getEnthalpy(self.__currentField)
            self.lineEditEnergy.setText(str(val))
            if energy_choice == 'dirichlet':
                self.labelEnergy.setText('J/kg')
            elif energy_choice == 'flux' :
                self.labelEnergy.setText('W/m2')
            else :
                self.labelEnergy.setText('K')


    @pyqtSlot(str)
    def __slotEnergy(self, text):
        """
        INPUT energy value
        """
        if self.lineEditEnergy.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setEnthalpy(self.__currentField, value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
