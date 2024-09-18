# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
- ImmersedBoundariesBoundaryConditionsWallView
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

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.gui.case.ImmersedBoundariesBoundaryConditionsWall import Ui_ImmersedBoundariesBoundaryConditionsWall
from code_saturne.model.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesBoundaryConditionsWall")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Internal dictionary for units
#-------------------------------------------------------------------------------
bc_energy_units = {'flux':'W/m2',
                   'flux_formula':'W/m2',
                   'temperature':'K',
                   'temperature_formula':'K'}

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesBoundaryConditionsWallView(QWidget, Ui_ImmersedBoundariesBoundaryConditionsWall) :
    """
    Wall boundary condition for the immersed object
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesBoundaryConditionsWall.__init__(self)
        self.setupUi(self)


    def setup(self, case, ibm, current_obj):
        """
        Setup the widget
        """
        self.case = case
        self.ibm = ibm
        self.current_obj = current_obj
        self.notebook = NotebookModel(self.case)
        self.mfm = MainFieldsModel(self.case)
        self.bc_thermal = False

        is_thermal = False
        for fId in self.mfm.getFieldIdList():
            if self.mfm.getFieldFromId(fId).enthalpy_model != 'off':
                is_thermal = True

        is_CHT = self.ibm.getObjectCHT(self.current_obj)
        is_FSI = self.ibm.getObjectFSI(self.current_obj)

        if ((is_CHT == 'off' and is_FSI == 'off') and is_thermal == True):
            self.bc_thermal = True

            self.__modelEnergy = ComboModel(self.comboBoxEnergy, 1, 1)

            self.__modelEnergy.addItem(self.tr("Thermal flux"), 'flux')
            self.__modelEnergy.addItem(self.tr("Thermal flux (user law)"), 'flux_formula')
            self.__modelEnergy.addItem(self.tr("Temperature"), 'temperature')
            self.__modelEnergy.addItem(self.tr("Temperature (user law)"), 'temperature_formula')
            #Thermal BC connections
            self.comboBoxEnergy.activated[str].connect(self.__slotChoiceEnergy)
            self.lineEditEnergy.textChanged[str].connect(self.__slotEnergy)
            validatorEner = DoubleValidator(self.lineEditEnergy)
            self.lineEditEnergy.setValidator(validatorEner)

        #Velocity BC Connections
        self.radioButtonSlip.clicked.connect(self.__slotChoiceVelocity)
        self.radioButtonNoSlip.clicked.connect(self.__slotChoiceVelocity)
        self.radioButtonWallLaw.clicked.connect(self.__slotChoiceVelocity)

        # MEG formula
        self.pushButtonEnergy.clicked.connect(self.slotThermalFormula)

        self.update()


    def update(self):

        self.groupBoxIBMThermal.hide()
        self.groupBoxVelocity.hide()

        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            return

        self.groupBoxVelocity.show()

        #Velocity boundary condition
        if self.ibm.getObjectBoundaryVelocityMode(self.current_obj) == "slip":
            self.radioButtonNoSlip.setChecked(False)
            self.radioButtonWallLaw.setChecked(False)
            self.radioButtonSlip.setChecked(True)
        elif self.ibm.getObjectBoundaryVelocityMode(self.current_obj) == "no_slip":
            self.radioButtonSlip.setChecked(False)
            self.radioButtonWallLaw.setChecked(False)
            self.radioButtonNoSlip.setChecked(True)
        elif self.ibm.getObjectBoundaryVelocityMode(self.current_obj) == "wall_law":
            self.radioButtonNoSlip.setChecked(False)
            self.radioButtonSlip.setChecked(False)
            self.radioButtonWallLaw.setChecked(True)

        if (self.bc_thermal == True):
            self.groupBoxIBMThermal.show()

            energy_mode = self.ibm.getObjectBoundaryEnergyMode(self.current_obj)

            if (energy_mode == "off"):
                energy_mode = "temperature" #by default
                self.__modelEnergy.setItem(str_model=energy_mode)
                self.ibm.setObjectBoundaryEnergyMode(self.current_obj, energy_mode)
            else:
                self.__modelEnergy.setItem(str_model=energy_mode)

            self.labelEnergy.setText(bc_energy_units.get(energy_mode, ''))

            if energy_mode in ['flux', 'temperature']:
                self.lineEditEnergy.setEnabled(True)
                self.pushButtonEnergy.setStyleSheet("background-color: grey")
                self.pushButtonEnergy.setEnabled(False) #no user law
                val = ''
                if energy_mode == 'temperature':
                    val = self.ibm.getObjectBoundaryTemperature(self.current_obj)
                    self.lineEditEnergy.setText(str(val))
                elif energy_mode == 'flux':
                    val = self.ibm.getObjectBoundaryThermalFlux(self.current_obj)
                    self.lineEditEnergy.setText(str(val))
                else:
                    self.lineEditEnergy.setEnabled(False)
                    self.pushButtonEnergy.setEnabled(True)

            if energy_mode in ['flux_formula', 'temperature_formula']:
                exp = ''
                if energy_mode == "flux_formula":
                    exp = self.ibm.getObjectBoundaryThermalFluxFormula(self.current_obj-1)
                elif energy_mode == "temperature_formula":
                    exp = self.ibm.getObjectBoundaryTemperatureFormula(self.current_obj-1)

                if exp:
                    self.pushButtonEnergy.setStyleSheet("background-color: green")
                    self.pushButtonEnergy.setToolTip(exp)
                else:
                    self.pushButtonEnergy.setStyleSheet("background-color: red")


    @pyqtSlot()
    def __slotChoiceVelocity(self):
        """
        Private slot.

        Select if the velocity boundary condition is slip, no-slip or wall law.
        """
        if self.radioButtonSlip.isChecked():
            self.ibm.setObjectBoundaryVelocityMode(self.current_obj,'slip')
        elif self.radioButtonNoSlip.isChecked():
            self.ibm.setObjectBoundaryVelocityMode(self.current_obj,'no_slip')
        elif self.radioButtonWallLaw.isChecked():
            self.ibm.setObjectBoundaryVelocityMode(self.current_obj,'wall_law')


    @pyqtSlot(str)
    def __slotChoiceEnergy(self, text):
        """
        INPUT choice of method of calculation of the energy
        """

        energy_mode = self.__modelEnergy.dicoV2M[str(text)]
        self.ibm.setObjectBoundaryEnergyMode(self.current_obj, energy_mode)

        if energy_mode[-7:] != "formula":
            val = ''
            if energy_mode == 'temperature':
                val = self.ibm.getObjectBoundaryTemperature(self.current_obj)
            elif energy_mode == 'flux':
                val = self.ibm.getObjectBoundaryThermalFlux(self.current_obj)

            self.pushButtonEnergy.setStyleSheet("background-color: grey")
            self.pushButtonEnergy.setEnabled(False)
            self.lineEditEnergy.setEnabled(True)
            self.lineEditEnergy.setText(str(val))
        else:
            self.lineEditEnergy.setText(str(''))
            self.lineEditEnergy.setEnabled(False)
            self.pushButtonEnergy.setEnabled(True)

            exp = ''
            if (energy_mode == 'temperature_formula'):
                exp = self.ibm.getObjectBoundaryTemperatureFormula(self.current_obj-1)
            elif (energy_mode == 'flux_formula'):
                exp = self.ibm.getObjectBoundaryThermalFluxFormula(self.current_obj-1)

            if exp:
                self.pushButtonEnergy.setToolTip(exp)
                self.pushButtonEnergy.setStyleSheet("background-color: green")
            else:
                self.pushButtonEnergy.setStyleSheet("background-color: red")

        self.labelEnergy.setText(bc_energy_units.get(energy_mode, ''))



    @pyqtSlot(str)
    def __slotEnergy(self, text):
        """
        INPUT energy value
        """
        if self.lineEditEnergy.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            energy_mode = self.ibm.getObjectBoundaryEnergyMode(self.current_obj)

            if (energy_mode == 'temperature'):
                self.ibm.setObjectBoundaryTemperature(self.current_obj, value)
            elif (energy_mode == 'flux'):
                self.ibm.setObjectBoundaryThermalFlux(self.current_obj, value)



    @pyqtSlot()
    def slotThermalFormula(self):
        objId = self.current_obj
        energy_mode = self.ibm.getObjectBoundaryEnergyMode(objId)

        exp, req, sym = self.ibm.getFormulaBoundaryEnergy(objId-1, energy_mode)

        var_name = ""
        exa = ""
        if energy_mode == "temperature_formula":
            var_name = "boundary_temperature"
            exa = """
temperature = 293.15;"""
        elif energy_mode == "flux_formula":
            var_name = "boundary_heat_flux"
            exa = """
flux = 1000.;"""

        obj_name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = "ibm_vol",
                                zone_name     = obj_name,
                                variable_name = var_name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotBoundaryThermalFormula -> %s" % str(result))
            self.pushButtonEnergy.setToolTip(result)
            self.pushButtonEnergy.setStyleSheet("background-color: green")

            if energy_mode == "temperature_formula":
                self.ibm.setObjectBoundaryTemperatureFormula(objId-1, result)
            elif energy_mode == "flux_formula":
                self.ibm.setObjectBoundaryThermalFluxFormula(objId-1, result)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
