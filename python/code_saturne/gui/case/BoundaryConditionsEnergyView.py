# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsEnergy import Ui_BoundaryConditionsEnergy

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
from code_saturne.model.ConjugateHeatTransferModel import ConjugateHeatTransferModel

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

        self.lineEditSyrthes.editingFinished.connect(self.__slotEnergySyrthes)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None
        self.__currentField = fieldId
        self.cht_model = ConjugateHeatTransferModel(self.case)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        # By default hide the syrthes instance line
        self.lineEditSyrthes.hide()

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
            self.__modelEnergy.addItem(self.tr("SYRTHES coupling"), "syrthes_coupling")

        if self.__boundary.getNature() == "inlet" or  self.__boundary.getNature() == "outlet":
            if MainFieldsModel(self.case).getEnergyResolution(self.__currentField) == 'on' :
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
                if ThermodynamicsModel(self.case).getMaterials(self.__currentField) == 'user_material' :
                    self.__modelEnergy.disableItem(2)
                    self.__modelEnergy.disableItem(3)
                else :
                    self.__modelEnergy.enableItem(2)
                    self.__modelEnergy.enableItem(3)
            else :
                self.hideWidget()
        elif self.__boundary.getNature() == "wall" :
            if len(MainFieldsModel(self.case).getEnthalpyResolvedField()) != 0 :
                energychoice = self.__boundary.getEnthalpyChoice("none")
                self.__modelEnergy.setItem(str_model=energychoice)
                val = self.__boundary.getEnthalpy("none")

                if energychoice == 'syrthes_coupling':
                    self.lineEditEnergy.hide()
                    self.lineEditSyrthes.show()
                    syrCompleter = \
                    QCompleter(self.cht_model.getSyrthesInstancesList())
                    self.lineEditSyrthes.setCompleter(syrCompleter)

                    self.lineEditSyrthes.setText(str(val))

                else:
                    self.lineEditSyrthes.hide()
                    self.lineEditEnergy.show()
                    self.lineEditEnergy.setText(str(val))

                if energychoice == 'flux' :
                    self.labelEnergy.setText('W/m2')
                elif energychoice == 'temperature':
                    self.labelEnergy.setText('K')
                else:
                    self.labelEnergy.setText('')

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
        old_choice = self.__boundary.getEnthalpyChoice(self.__currentField)
        # If we switch from coupling to BC condition, delete all previous
        # data
        if old_choice != energy_choice and old_choice == "syrthes_coupling":
            syrthes_name = self.__boundary.getEnthalpy(self.__currentField)
            bnd_label = self.__boundary.getLabel()
            self.cht_model.deleteSyrthesCoupling(syrthes_name, bnd_label)
            self.__boundary.setEnthalpy(self.__currentField, 0.)

        self.__boundary.setEnthalpyChoice(self.__currentField, energy_choice)


        if energy_choice == "hsat_P" :
            self.lineEditEnergy.hide()
            self.lineEditSyrthes.hide()
            self.labelEnergy.hide()
        else :
            self.labelEnergy.show()
            val = self.__boundary.getEnthalpy(self.__currentField)
            if energy_choice == "syrthes_coupling":
                self.lineEditSyrthes.show()
                self.lineEditEnergy.hide()
                if val not in self.cht_model.getSyrthesInstancesList():
                    val = ""
                    self.__boundary.setEnthalpy(self.__currentField, val)
                self.lineEditSyrthes.setText(str(val))
            else:
                self.lineEditSyrthes.hide()
                self.lineEditEnergy.show()
                self.lineEditEnergy.setText(str(val))

            if energy_choice == 'dirichlet':
                self.labelEnergy.setText('J/kg')
            elif energy_choice == 'flux' :
                self.labelEnergy.setText('W/m2')
            elif energy_choice == 'temperature' :
                self.labelEnergy.setText('K')
            else:
                self.labelEnergy.setText('')


    @pyqtSlot(str)
    def __slotEnergy(self, text):
        """
        INPUT energy value
        """
        if self.lineEditEnergy.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setEnthalpy(self.__currentField, value)


    @pyqtSlot()
    def __slotEnergySyrthes(self):
        """
        Input syrthes instance name
        """

        value = str(self.lineEditSyrthes.text())
        if value:
            value_p = val = self.__boundary.getEnthalpy("none")
            bnd_label = self.__boundary.getLabel()
            if value != value_p:
                self.cht_model.deleteSyrthesCoupling(value_p, bnd_label)


            if value not in self.cht_model.getSyrthesInstancesList():
                self.cht_model.addSyrthesCoupling(value)

            self.cht_model.addBoundaryLabel(value, bnd_label)

            self.__boundary.setEnthalpy(self.__currentField, value)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
