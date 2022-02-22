# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2015 EDF S.A.
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
This module defines the 'Global numerical parameters' page.

This module contains the following classes:
- GlobalNumericalParametersAdvancedOptionsDialogView
- GlobalNumericalParametersView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.case.GlobalNumericalParameters import Ui_GlobalNumericalParameters
from code_saturne.gui.case.GlobalNumericalParametersAdvancedOptionsDialog import Ui_GlobalNumericalParametersAdvancedOptionsDialog
from code_saturne.model.GlobalNumericalParametersModel import GlobalNumericalParametersModel
from code_saturne.model.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("GlobalNumericalParametersView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# GlobalNumericalParametersAdvancedOptionsDialogView class
#-------------------------------------------------------------------------------

class GlobalNumericalParametersAdvancedOptionsDialogView(QDialog, Ui_GlobalNumericalParametersAdvancedOptionsDialog):
    """
    Advanced global numerical parameters layout.
    """
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_GlobalNumericalParametersAdvancedOptionsDialog.__init__(self)
        self.setupUi(self)
        title = self.tr("Advanced options for Global Numerical Parameters")
        self.setWindowTitle(title)
        self.default = default
        self.result  = self.default.copy()
        self.case    = case

        self.case.undoStopGlobal()

        # Combo model
        self.modelVelocityUpdate = ComboModel(self.comboBoxVelocityUpdate, 4, 1)
        self.modelVelocityUpdate.addItem(self.tr("By the pressure gradient increment"), "pressure_gradient_increment")
        self.modelVelocityUpdate.addItem(self.tr("RT0 by flow rate (internal+boundary)"), "flow_rate")
        self.modelVelocityUpdate.addItem(self.tr("RT0 by flumas(b) increment"), "flumas_increment")
        self.modelVelocityUpdate.addItem(self.tr("By means of Conv+Diff+ST equation"), "conv_diff_equation")

        self.modelGradientPressure = ComboModel(self.comboBoxGradientPressure, 5, 1)
        self.modelGradientPressure.addItem(self.tr("Mass ponderation"), "mass_ponderation")
        self.modelGradientPressure.addItem(self.tr("Standard"), "standard")
        self.modelGradientPressure.addItem(self.tr("Controlled mass ponderation"), "controlled_mass_pound")
        self.modelGradientPressure.addItem(self.tr("Momentum ponderation"), "momentum_pound")
        self.modelGradientPressure.addItem(self.tr("Gravity+HL+ST momentum"), "gravity_momentum")

        # Validator
        validatorSumAlpha = DoubleValidator(self.lineEditMaxSumAlpha, min = 0., max = 1.)
        validatorAlphaP   = IntValidator(self.lineEditNumberAlphaPCycle, min = 1)

        validatorSumAlpha.setExclusiveMin(False)
        validatorAlphaP.setExclusiveMin(False)

        self.lineEditMaxSumAlpha.setValidator(validatorSumAlpha)
        self.lineEditNumberAlphaPCycle.setValidator(validatorAlphaP)

        # Initialization
        if self.result['pressure_symetrisation'] == 'on' :
            self.checkBoxSymetPressure.setChecked(True)
        else :
            self.checkBoxSymetPressure.setChecked(False)

        self.modelVelocityUpdate.setItem(str_model=self.result['velocity_update'])
        self.modelGradientPressure.setItem(str_model=self.result['pressure_gradient'])
        self.lineEditMaxSumAlpha.setText(str(self.result['max_sum_alpha']))
        self.lineEditNumberAlphaPCycle.setText(str(self.result['alpha_p_cycle']))

        self.case.undoStartGlobal()


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        if self.checkBoxSymetPressure.isChecked():
            self.result['pressure_symetrisation'] = 'on'
        else :
            self.result['pressure_symetrisation'] = 'off'

        self.result['velocity_update'] = \
                 self.modelVelocityUpdate.dicoV2M[str(self.comboBoxVelocityUpdate.currentText())]

        self.result['pressure_gradient'] = \
                 self.modelGradientPressure.dicoV2M[str(self.comboBoxGradientPressure.currentText())]

        self.result['max_sum_alpha']     = float(str(self.lineEditMaxSumAlpha.text()))
        self.result['alpha_p_cycle']     = int(str(self.lineEditNumberAlphaPCycle.text()))

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


#-------------------------------------------------------------------------------
# GlobalNumericalParametersView class
#-------------------------------------------------------------------------------

class GlobalNumericalParametersView(QWidget, Ui_GlobalNumericalParameters):
    """
    Global numerical parameters layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_GlobalNumericalParameters.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = GlobalNumericalParametersModel(self.case)

        # Combo model
        self.modelVelocityAlgorithm = ComboModel(self.comboBoxVelocityAlgorithm, 3, 1)
        self.modelVelocityAlgorithm.addItem(self.tr("Standard"), "standard_difvit")
        self.modelVelocityAlgorithm.addItem(self.tr("Coupled"), "coupled_difvitc")
        self.modelVelocityAlgorithm.addItem(self.tr("Mean velocity - relative velocity"), "mean_velocity_relative_velocity")

        mfm = MainFieldsModel(self.case)
        if len(mfm.getFieldIdList()) < 2 :
            self.modelVelocityAlgorithm.disableItem(2)
        else :
            self.modelVelocityAlgorithm.enableItem(2)

        # Validator
        validatorMaxRestart = IntValidator(self.lineEditMaxRestart, min = 0)
        validatorSplitting  = DoubleValidator(self.lineEditTimeSplitting, min = 0.)
        validatorPRelax     = DoubleValidator(self.lineEditPressureRelaxation, min = 0.)
        validatorMinP       = DoubleValidator(self.lineEditMinimumPressure)
        validatorMaxP       = DoubleValidator(self.lineEditMaximumPressure)

        validatorMaxRestart.setExclusiveMin(False)
        validatorSplitting.setExclusiveMin(True)
        validatorPRelax.setExclusiveMin(True)

        self.lineEditMaxRestart.setValidator(validatorMaxRestart)
        self.lineEditTimeSplitting.setValidator(validatorSplitting)
        self.lineEditPressureRelaxation.setValidator(validatorPRelax)
        self.lineEditMinimumPressure.setValidator(validatorMinP)
        self.lineEditMaximumPressure.setValidator(validatorMaxP)

        # Connections
        self.checkBoxRestart.clicked.connect(self.slotRestart)
        self.checkBoxPotentialState.clicked.connect(self.slotPotentialState)
        self.checkBoxFacesReconstruction.clicked.connect(self.slotFacesReconstruction)
        self.checkBoxMultigrid.clicked.connect(self.slotMultigrid)
        self.lineEditMinimumPressure.textChanged[str].connect(self.slotMinimumPressure)
        self.lineEditMaximumPressure.textChanged[str].connect(self.slotMaximumPressure)
        self.pushButtonAdvanced.clicked.connect(self.slotAdvancedOptions)
        self.lineEditMaxRestart.textChanged[str].connect(self.slotMaxRestart)
        self.lineEditTimeSplitting.textChanged[str].connect(self.slotTimeSplitting)
        self.lineEditPressureRelaxation.textChanged[str].connect(self.slotPressureRelaxation)
        self.checkBoxUpwindAlphaEnergy.clicked.connect(self.slotUpwindAlphaEnergy)
        self.checkBoxStopRestart.clicked.connect(self.slotStopRestart)
        self.comboBoxVelocityAlgorithm.activated[str].connect(self.slotVelocityAlgorithm)

        self.checkBoxRegulBadCells.clicked.connect(self.slotRegulateBadCells)

        # Initialize widget
        status = self.mdl.getRestartTimeStep()
        if status == 'on':
            self.checkBoxRestart.setChecked(1)
            self.groupBoxRestartOption.show()

            value = self.mdl.getMaxNumberOfRestart()
            self.lineEditMaxRestart.setText(str(value))
            value = self.mdl.getTimeSplit()
            self.lineEditTimeSplitting.setText(str(value))
            value = self.mdl.getPressureRelaxation()
            self.lineEditPressureRelaxation.setText(str(value))
            status = self.mdl.getUpwindScheme()
            if status == 'on':
                self.checkBoxUpwindAlphaEnergy.setChecked(True)
            else :
                self.checkBoxUpwindAlphaEnergy.setChecked(False)
            status = self.mdl.getStopNoConvergence()
            if status == 'on':
                self.checkBoxStopRestart.setChecked(True)
            else :
                self.checkBoxStopRestart.setChecked(False)
        else :
            self.checkBoxRestart.setChecked(0)
            self.groupBoxRestartOption.hide()

        is_compressible = False
        for fid in mfm.getFieldIdList():
            if mfm.getCompressibleStatus(int(fid)) == 'on':
                is_compressible = True
                break

        if is_compressible:
            self.mdl.setPotentielState('off')
            self.checkBoxPotentialState.setChecked(0)
            self.checkBoxPotentialState.setEnabled(False)
        else:
            status = self.mdl.getPotentielState()
            if status == 'on':
                self.checkBoxPotentialState.setChecked(1)
            else :
                self.checkBoxPotentialState.setChecked(0)

        status = self.mdl.getFacesReconstruction()
        if status == 'on':
            self.checkBoxFacesReconstruction.setChecked(1)
        else :
            self.checkBoxFacesReconstruction.setChecked(0)

        status = self.mdl.getRegulateBadCElls()
        self.checkBoxRegulBadCells.setChecked(status == 'on')

        status = self.mdl.getMultigridStatus()
        if status == 'on':
            self.checkBoxMultigrid.setChecked(1)
        else :
            self.checkBoxMultigrid.setChecked(0)

        value = self.mdl.getMinPressure()
        self.lineEditMinimumPressure.setText(str(value))
        value = self.mdl.getMaxPressure()
        self.lineEditMaximumPressure.setText(str(value))

        model = self.mdl.getVelocityPredictorAlgo()
        self.modelVelocityAlgorithm.setItem(str_model=model)

        self.case.undoStartGlobal()

        predefined_flow = MainFieldsModel(self.case).getPredefinedFlow()
        if predefined_flow in ["free_surface", "droplet_flow", "multiregime"]:
            self.modelVelocityAlgorithm.setItem(str_model="coupled_difvitc")
            self.comboBoxVelocityAlgorithm.setEnabled(False)
        elif predefined_flow == "boiling_flow":
            self.modelVelocityAlgorithm.setItem(str_model="mean_velocity_relative_velocity")
            self.comboBoxVelocityAlgorithm.setEnabled(False)

    @pyqtSlot()
    def slotRestart(self):
        """
        Input if restart time step if not converged
        """
        if self.checkBoxRestart.isChecked():
            self.mdl.setRestartTimeStep('on')
            self.groupBoxRestartOption.show()

            value = self.mdl.getMaxNumberOfRestart()
            self.lineEditMaxRestart.setText(str(value))
            value = self.mdl.getTimeSplit()
            self.lineEditTimeSplitting.setText(str(value))
            value = self.mdl.getPressureRelaxation()
            self.lineEditPressureRelaxation.setText(str(value))
            status = self.mdl.getUpwindScheme()
            if status == 'on':
                self.checkBoxUpwindAlphaEnergy.setChecked(True)
            else :
                self.checkBoxUpwindAlphaEnergy.setChecked(False)
            status = self.mdl.getStopNoConvergence()
            if status == 'on':
                self.checkBoxStopRestart.setChecked(True)
            else :
                self.checkBoxStopRestart.setChecked(False)
        else:
            self.mdl.setRestartTimeStep('off')
            self.groupBoxRestartOption.hide()


    @pyqtSlot()
    def slotPotentialState(self):
        """
        Input if restart time step if not converged
        """
        if self.checkBoxPotentialState.isChecked():
            self.mdl.setPotentielState('on')
        else:
            self.mdl.setPotentielState('off')


    @pyqtSlot()
    def slotFacesReconstruction(self):
        """
        Input if faces reconstruction
        """
        if self.checkBoxFacesReconstruction.isChecked():
            self.mdl.setFacesReconstruction('on')
            self.checkBoxFacesReconstruction.setChecked(1)
        else:
            self.mdl.setFacesReconstruction('off')


    @pyqtSlot()
    def slotMultigrid(self):
        """
        Input if multigrid for pressure
        """
        if self.checkBoxMultigrid.isChecked():
            self.mdl.setMultigridStatus('on')
        else:
            self.mdl.setMultigridStatus('off')


    @pyqtSlot(str)
    def slotMinimumPressure(self, text):
        """
        Input value of minimum pressure
        """
        if self.lineEditMinimumPressure.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMinPressure(value)


    @pyqtSlot(str)
    def slotMaximumPressure(self, text):
        """
        Input value of maximum pressure
        """
        if self.lineEditMaximumPressure.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setMaxPressure(value)


    @pyqtSlot()
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        default = {}
        default['velocity_update'] = self.mdl.getVelocityUpdate()
        default['pressure_symetrisation'] = self.mdl.getPressureSymetrisation()
        default['pressure_gradient'] = self.mdl.getPressureGradient()
        default['max_sum_alpha'] = self.mdl.getSumAlpha()
        default['alpha_p_cycle'] = self.mdl.getAlphaPressureCycles()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = GlobalNumericalParametersAdvancedOptionsDialogView(self, self.case, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.mdl.setVelocityUpdate(result['velocity_update'])
            self.mdl.setPressureSymetrisation(result['pressure_symetrisation'])
            self.mdl.setPressureGradient(result['pressure_gradient'])
            self.mdl.setSumAlpha(result['max_sum_alpha'])
            self.mdl.setAlphaPressureCycles(result['alpha_p_cycle'])


    @pyqtSlot(str)
    def slotMaxRestart(self, text):
        """
        Input value of Maximum number of restart
        """
        if self.lineEditMaxRestart.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.mdl.setMaxNumberOfRestart(value)


    @pyqtSlot(str)
    def slotTimeSplitting(self, text):
        """
        Input value of time-step splitting
        """
        if self.lineEditTimeSplitting.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setTimeSplit(value)


    @pyqtSlot(str)
    def slotPressureRelaxation(self, text):
        """
        Input value of pressure increment relaxation
        """
        if self.lineEditPressureRelaxation.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setPressureRelaxation(value)


    @pyqtSlot()
    def slotUpwindAlphaEnergy(self):
        """
        Input if upwind scheme for mass and energy
        """
        if self.checkBoxUpwindAlphaEnergy.isChecked():
            self.mdl.setUpwindScheme('on')
        else:
            self.mdl.setUpwindScheme('off')


    @pyqtSlot()
    def slotStopRestart(self):
        """
        Input if stop if no convergence
        """
        if self.checkBoxStopRestart.isChecked():
            self.mdl.setStopNoConvergence('on')
        else:
            self.mdl.setStopNoConvergence('off')


    @pyqtSlot(str)
    def slotVelocityAlgorithm(self, text):
        """
        Input velocity algorithm model
        """
        model = self.modelVelocityAlgorithm.dicoV2M[str(text)]
        self.mdl.setVelocityPredictorAlgo(model)


    @pyqtSlot()
    def slotRegulateBadCells(self):
        """
        Activate bad cells regulations.
        """
        if self.checkBoxRegulBadCells.isChecked():
            self.mdl.setRegulateBadCells('on')
        else:
            self.mdl.setRegulateBadCells('off')

