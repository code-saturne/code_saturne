# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, IntValidator, DoubleValidator, from_qvariant
from code_saturne.gui.case.ThermalForm import Ui_ThermalForm
from code_saturne.model.ElectricalModel import ElectricalModel
from code_saturne.model.CoalCombustionModel import CoalCombustionModel
from code_saturne.model.GasCombustionModel import GasCombustionModel
from code_saturne.model.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.model.CompressibleModel import CompressibleModel
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.gui.case.ThermalRadiationAdvancedDialogForm import Ui_ThermalRadiationAdvancedDialogForm
from code_saturne.model.ThermalRadiationModel import ThermalRadiationModel
from code_saturne.model.ThermalParticlesRadiationModel import ThermalParticlesRadiationModel
from code_saturne.model.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermalView")
log.setLevel(GuiParam.DEBUG)

#--------------------------------------------------------------------------------
# Popup Class
#--------------------------------------------------------------------------------

class ThermalRadiationAdvancedDialogView(QDialog, Ui_ThermalRadiationAdvancedDialogForm):
    """
    Building of popup window for advanced options.
    """
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_ThermalRadiationAdvancedDialogForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result  = self.default.copy()

        # Combo models

        self.modelTSRay  = ComboModel(self.comboBoxTSRay, 4, 1)
        self.modelPrintT = ComboModel(self.comboBoxPrintT, 3, 1)
        self.modelPrintL = ComboModel(self.comboBoxPrintL, 3, 1)

        self.modelTSRay.addItem('-1', '-1')
        self.modelTSRay.addItem('0', '0')
        self.modelTSRay.addItem('1', '1')
        self.modelTSRay.addItem('2', '2')

        self.modelPrintT.addItem('0', '0')
        self.modelPrintT.addItem('1', '1')
        self.modelPrintT.addItem('2', '2')

        self.modelPrintL.addItem('0', '0')
        self.modelPrintL.addItem('1', '1')
        self.modelPrintL.addItem('2', '2')

        self.frequ     = self.default['frequency']
        self.tsr       = self.default['idiver']
        self.printTemp = self.default['tempP']
        self.printLum  = self.default['intensity']
        model          = self.default['model']

        # Initialization

        self.lineEditFreq.setText(str(self.frequ))
        self.modelTSRay.setItem(str_model=str(self.tsr))
        self.modelPrintT.setItem(str_model=str(self.printTemp))
        self.modelPrintL.setItem(str_model=str(self.printLum))

        if model == 'dom':
            self.labelPrintL.show()
            self.comboBoxPrintL.show()
        else:
            self.labelPrintL.hide()
            self.comboBoxPrintL.hide()

        # Validator

        validatorFreq = IntValidator(self.lineEditFreq, min=1)
        self.lineEditFreq.setValidator(validatorFreq)

        self.case.undoStartGlobal()


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """
        if self.lineEditFreq.validator().state == QValidator.Acceptable:
            self.result['frequency'] = from_qvariant(self.lineEditFreq.text(), int)
        self.result['idiver']    = from_qvariant(self.comboBoxTSRay.currentText(), int)
        self.result['tempP']     = from_qvariant(self.comboBoxPrintT.currentText(), int)
        self.result['intensity'] = from_qvariant(self.comboBoxPrintL.currentText(), int)

        QDialog.accept(self)


    def reject(self):
        """
        Method called when 'Cancel' button is clicked.
        """
        QDialog.reject(self)


    def get_result(self):
        """
        Method to get the result.
        """
        return self.result


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ThermalView(QWidget, Ui_ThermalForm):
    """
    Class to open Thermal Scalar Transport Page.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ThermalForm.__init__(self)
        self.setupUi(self)

        self.browser = tree
        self.case = case
        self.case.undoStopGlobal()

        # Current physics

        self.coal_or_gas = "off"

        if self.case.module_name() == 'code_saturne':
            from code_saturne.model.CoalCombustionModel import CoalCombustionModel
            self.coal_or_gas = CoalCombustionModel(self.case).getCoalCombustionModel("only")
            del CoalCombustionModel
            if self.coal_or_gas == "off":
                from code_saturne.model.GasCombustionModel import GasCombustionModel
                self.coal_or_gas = GasCombustionModel(self.case).getGasCombustionModel()
                del GasCombustionModel

        # Thermal scalar
        #---------------

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

        if self.case.module_name() != 'code_saturne':
            self.comboBoxThermal.setEnabled(False)
            if (MainFieldsModel(self.case).getSolidPhaseList() != []) :
                self.groupBoxParticlesRadiation.show()
            else:
                self.groupBoxParticlesRadiation.hide()
        else:
            self.groupBoxParticlesRadiation.hide()
            if ElectricalModel(self.case).getElectricalModel() != 'off':
                self.comboBoxThermal.setEnabled(False)

            elif self.coal_or_gas != 'off':
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

        # Fluid Radiation model
        #----------------
        self.__setFluidRadiation__(model)

        # inter-particles radiation model
        self.partRadiationModel = ThermalParticlesRadiationModel(self.case)
        self.__setParticlesRadiation__()

        # Soot model
        #----------------

        # Only activate soot for code_saturne (for which it is available).
        # Otherwise it may lead to the spurious apparition of pages for NCFD
        if self.case.module_name() == 'code_saturne':
            self.__setSoot__()

        # Undo/redo part

        self.case.undoStartGlobal()


    def __setParticlesRadiation__(self):
        """
        Setup particles radiation view
        """

        # Validators
        validator = DoubleValidator(self.lineEditEmissivity, min=0.0)
        validator.setExclusiveMin(False)
        self.lineEditEmissivity.setValidator(validator)

        self.groupBoxParticlesRadiation.toggled[bool].connect(self.slotActivateParticlesRadiation)
        self.lineEditEmissivity.textChanged[str].connect(self.slotSetEmissivity)

        # Initialization from XML model
        if self.partRadiationModel.isActivated == "on":
            self.groupBoxParticlesRadiation.setChecked(True)
        else:
            self.groupBoxParticlesRadiation.setChecked(False)
            self.slotActivateParticlesRadiation(False)
        self.lineEditEmissivity.setText(self.partRadiationModel.emissivity)


    def slotActivateParticlesRadiation(self, status):
        if status:
            self.partRadiationModel.isActivated = "on"
            self.lineEditEmissivity.show()
            self.labelEmissivity.show()
        else:
            self.partRadiationModel.isActivated = "off"
            self.lineEditEmissivity.hide()
            self.labelEmissivity.hide()


    def slotSetEmissivity(self, value):
        if self.lineEditEmissivity.validator().state == QValidator.Acceptable:
            self.partRadiationModel.emissivity = value

    def __setFluidRadiation__(self, model):
        """
        Update for radiation model
        """
        if self.case.module_name() == "code_saturne" and model == 'off':
            self.ThermalRadiationGroupBox.hide()
            return

        self.ThermalRadiationGroupBox.show()
        self.rmdl = ThermalRadiationModel(self.case)

        # Combo models

        self.modelRadModel   = ComboModel(self.comboBoxRadModel, 3, 1)
        self.modelDirection  = ComboModel(self.comboBoxQuadrature, 8, 1)
        self.modelAbsorption = ComboModel(self.comboBoxAbsorption, 3, 1)

        self.modelRadModel.addItem("No radiative transfers", 'off')
        self.modelRadModel.addItem("Discrete ordinates method", 'dom')
        self.modelRadModel.addItem("P-1 Model", 'p-1')

        self.modelDirection.addItem("24 directions (S4)",   "1")
        self.modelDirection.addItem("48 directions (S6)",   "2")
        self.modelDirection.addItem("80 directions (S8)",   "3")
        self.modelDirection.addItem("32 directions (T2)",   "4")
        self.modelDirection.addItem("128 directions (T4)",  "5")
        self.modelDirection.addItem("8n^2 directions (Tn)", "6")
        self.modelDirection.addItem("120 directions (LC11)", "7")
        self.modelDirection.addItem("48 directions (DCT020-2468)", "8")

        # Connections

        self.comboBoxRadModel.activated[str].connect(self.slotFluidRadiativeTransfer)
        self.checkBoxRadRestart.stateChanged.connect(self.slotRadRestart)
        self.comboBoxQuadrature.activated[str].connect(self.slotDirection)
        self.lineEditNdirec.textChanged[str].connect(self.slotNdirec)
        self.comboBoxAbsorption.activated[str].connect(self.slotTypeCoefficient)
        self.lineEditCoeff.textChanged[str].connect(self.slotAbsorptionCoefficient)
        self.toolButtonAdvanced.clicked.connect(self.slotAdvancedOptions)

        # Validator

        validatorCoeff = DoubleValidator(self.lineEditCoeff, min=0.0)
        self.lineEditCoeff.setValidator(validatorCoeff)

        validatorNdir = IntValidator(self.lineEditNdirec, min=2)
        self.lineEditNdirec.setValidator(validatorNdir)

        self.modelAbsorption.addItem('constant', 'constant')
        self.modelAbsorption.addItem('user function (cs_user_rad_transfer_absorption)', 'variable')
        self.modelAbsorption.addItem('H2O and CO2 mixing (Modak)', 'modak')

        if self.coal_or_gas != "off":
            self.modelAbsorption.disableItem(str_model='variable')
            self.modelAbsorption.enableItem(str_model='modak')
        else:
            self.modelAbsorption.disableItem(str_model='modak')
            self.modelAbsorption.enableItem(str_model='variable')

        # Initialization

        self.modelRadModel.setItem(str_model=self.rmdl.getRadiativeModel())

        # For multiphase flows, only droplet laden gas flows are accepted
        if self.case.module_name() != 'code_saturne':
            mfm = MainFieldsModel(self.case)
            predefined_flow = mfm.getPredefinedFlow()
            if predefined_flow not in ["droplet_flow", "multiregime_flow", "user_flow"]:
                self.modelRadModel.setItem(str_model='off')
                self.modelRadModel.disableItem(str_model='dom')
                self.modelRadModel.disableItem(str_model='p-1')

        self.slotFluidRadiativeTransfer()

        if self.rmdl.getRestart() == 'on':
            self.checkBoxRadRestart.setChecked(True)
        else:
            self.checkBoxRadRestart.setChecked(False)

        value = self.rmdl.getTypeCoeff()
        self.modelAbsorption.setItem(str_model=value)
        self.slotTypeCoefficient(self.modelAbsorption.dicoM2V[value])

        self.lineEditCoeff.setText(str(self.rmdl.getAbsorCoeff()))


    def __setSoot__(self):
        """
        Update for soot model
        """
        self.rmdl = ThermalRadiationModel(self.case)
        self.modelrad = self.rmdl.getRadiativeModel()

        self.gas = GasCombustionModel(self.case)
        self.model = self.gas.getGasCombustionModel()

        if self.model == "off" or self.modelrad =="off":
            self.groupBoxSoot.hide()
            return

        self.groupBoxSoot.show()

        # Combo models

        self.modelSoot   = ComboModel(self.comboBoxSoot, 3, 1)

        self.modelSoot.addItem("None", 'off')
        self.modelSoot.addItem("Constant soot yield", 'constant_soot_yield')
        self.modelSoot.addItem("2 equations model of Moss et al.", 'moss')

        # Connections

        self.comboBoxSoot.activated[str].connect(self.slotSoot)
        self.lineEditSootDensity.textChanged[str].connect(self.slotSootDensity)
        self.lineEditSootFraction.textChanged[str].connect(self.slotSootFraction)

        # Validator

        validatorSootDensity = DoubleValidator(self.lineEditSootDensity, min=0.0)
        self.lineEditSootDensity.setValidator(validatorSootDensity)

        validatorSootFraction = DoubleValidator(self.lineEditSootFraction, min=0.0, max=1.0)
        self.lineEditSootFraction.setValidator(validatorSootFraction)

        # Initialization

        self.modelSoot.setItem(str_model=self.gas.getSootModel())

        if self.rmdl.getRadiativeModel() == "off":
            self.modelSoot.setItem(str_model='off')

        self.slotSoot()


    @pyqtSlot(str)
    def slotThermalScalar(self, text):
        """
        Update the thermal scalar markup.
        """
        th = self.modelThermal.dicoV2M[str(text)]
        self.thermal.setThermalModel(th)

        self.__setFluidRadiation__(th)
        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotFluidRadiativeTransfer(self):
        """
        """
        self.gas = GasCombustionModel(self.case)
        model = self.modelRadModel.dicoV2M[str(self.comboBoxRadModel.currentText())]
        self.rmdl.setRadiativeModel(model)
        if model == 'off':
            self.frameOptions.hide()
            self.groupBoxDirection.hide()
            self.groupBoxAbsorptionCoeff.hide()
            self.gas.setSootModel("off")
            self.groupBoxSoot.hide()
        else:
            self.groupBoxSoot.show()
            self.groupBoxAbsorptionCoeff.show()
            self.frameOptions.show()

            if model == 'p-1':
                self.groupBoxDirection.hide()
            elif model == 'dom':
                self.groupBoxDirection.show()
                n = self.rmdl.getQuadrature()
                self.modelDirection.setItem(str_model=str(n))

                if str(n) == "6":
                    self.label_2.show()
                    self.lineEditNdirec.show()
                    self.lineEditNdirec.setText(str(self.rmdl.getNbDir()))
                else:
                    self.label_2.hide()
                    self.lineEditNdirec.hide()
        self.__setSoot__()


    @pyqtSlot(int)
    def slotRadRestart(self, val):
        """
        """
        if val == 0:
            self.rmdl.setRestart("off")
        else:
            self.rmdl.setRestart("on")


    @pyqtSlot(str)
    def slotDirection(self, text):
        """
        """
        n = int(self.modelDirection.dicoV2M[str(text)])
        self.rmdl.setQuadrature(n)

        if n == 6:
            self.label_2.show()
            self.lineEditNdirec.show()
            self.lineEditNdirec.setText(str(self.rmdl.getNbDir()))
        else:
            self.label_2.hide()
            self.lineEditNdirec.hide()


    @pyqtSlot(str)
    def slotNdirec(self, text):
        """
        """
        if self.lineEditNdirec.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, int)
            self.rmdl.setNbDir(n)


    @pyqtSlot(str)
    def slotTypeCoefficient(self, text):
        """
        """
        typeCoeff = self.modelAbsorption.dicoV2M[str(text)]
        self.rmdl.setTypeCoeff(typeCoeff)

        self.gas = GasCombustionModel(self.case)
        self.model = self.gas.getGasCombustionModel()

        if typeCoeff == 'constant':
            self.lineEditCoeff.show()
            self.lineEditCoeff.setEnabled(True)
            self.label.show()
            # we don't show the absorption coeff when a gas combustion model is used because the value is read from the Janaf File.
            if self.model != 'off':
                self.lineEditCoeff.hide()
                self.label.hide()
        elif typeCoeff == 'modak':
            self.lineEditCoeff.hide()
            self.lineEditCoeff.setDisabled(True)
            self.label.hide()
        else:
            self.lineEditCoeff.setDisabled(True)


    @pyqtSlot(str)
    def slotAbsorptionCoefficient(self, text):
        """
        """
        if self.lineEditCoeff.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.rmdl.setAbsorCoeff(c)


    @pyqtSlot()
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        default = {}
        default['frequency'] = self.rmdl.getFrequency()
        default['idiver']    = self.rmdl.getTrs()
        default['tempP']     = self.rmdl.getTemperatureListing()
        default['intensity'] = self.rmdl.getIntensityResolution()
        default['model']     = self.rmdl.getRadiativeModel()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = ThermalRadiationAdvancedDialogView(self, self.case, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.rmdl.setFrequency(result['frequency'])
            self.rmdl.setTrs(result['idiver'])
            self.rmdl.setTemperatureListing(result['tempP'])
            self.rmdl.setIntensityResolution(result['intensity'])


    @pyqtSlot(str)
    def slotSoot(self):
        """
        """
        model = self.modelSoot.dicoV2M[str(self.comboBoxSoot.currentText())]
        self.gas.setSootModel(model)
        if model == 'off':
            self.labelSootDensity.hide()
            self.lineEditSootDensity.hide()
            self.label_4.hide()

            self.labelSootFraction.hide()
            self.labelSootFractionbis.hide()
            self.lineEditSootFraction.hide()

        elif model == 'constant_soot_yield':

            self.labelSootDensity.show()
            self.lineEditSootDensity.show()
            self.lineEditSootDensity.setText(str(self.gas.getSootDensity()))
            self.label_4.show()

            self.labelSootFraction.show()
            self.labelSootFractionbis.show()
            self.lineEditSootFraction.show()
            self.lineEditSootFraction.setText(str(self.gas.getSootFraction()))

        elif model == 'moss':

            self.labelSootDensity.show()
            self.lineEditSootDensity.show()
            self.lineEditSootDensity.setText(str(self.gas.getSootDensity()))
            self.label_4.show()

            self.labelSootFraction.hide()
            self.labelSootFractionbis.hide()
            self.lineEditSootFraction.hide()

    @pyqtSlot(str)
    def slotSootDensity(self, text):
        """
        """
        if self.lineEditSootDensity.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.gas.setSootDensity(c)

    @pyqtSlot(str)
    def slotSootFraction(self, text):
        """
        """
        if self.lineEditSootFraction.validator().state == QValidator.Acceptable:
            c  = from_qvariant(text, float)
            self.gas.setSootFraction(c)

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
