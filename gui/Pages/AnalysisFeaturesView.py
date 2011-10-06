# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
This module defines the Page in which the user defines the physical options
of the treated case.

This module contains the following class :
- AnalysisFeaturesView
"""

#-------------------------------------------------------------------------------
# Standard modules
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
from AnalysisFeaturesForm import Ui_AnalysisFeaturesForm

import Base.QtPage as QtPage

from Pages.TurbulenceModel import TurbulenceModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.LagrangianModel import LagrangianModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.CoalCombustionModel import CoalCombustionModel
from Pages.ElectricalModelsModel import ElectricalModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.ThermalRadiationModel import ThermalRadiationModel
from Pages.SteadyManagementModel import SteadyManagementModel
from Pages.AtmosphericFlowsModel import AtmosphericFlowsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("AnalysisFeaturesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Analysis Features View class
#-------------------------------------------------------------------------------

class AnalysisFeaturesView(QWidget, Ui_AnalysisFeaturesForm):
    """
    Class to open Calculation Features Page.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)

        Ui_AnalysisFeaturesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.browser = tree

        self.lagr  = LagrangianModel(self.case)
        self.turb  = TurbulenceModel(self.case)
        self.therm = ThermalScalarModel(self.case)
        self.gas   = GasCombustionModel(self.case)
        self.pcoal = CoalCombustionModel(self.case)
        self.elect = ElectricalModel(self.case)
        self.scal  = DefineUserScalarsModel(self.case)
        self.std   = SteadyManagementModel(self.case)
        self.atmo  = AtmosphericFlowsModel(self.case)

        # Set models and number of elements for combo boxes

        self.modelSteadyFlow         = QtPage.ComboModel(self.comboBoxSteadyFlow,2,1)
        self.modelLagrangian         = QtPage.ComboModel(self.comboBoxLagrangian,2,1)
        self.modelAtmospheric        = QtPage.ComboModel(self.comboBoxAtmospheric,4,1)
        self.modelGasCombustionModel = QtPage.ComboModel(self.comboBoxGasCombustionModel,3,1)
        self.modelPulverizedCoal     = QtPage.ComboModel(self.comboBoxPulverizedCoal,3,1)
        self.modelJouleEffect        = QtPage.ComboModel(self.comboBoxJouleEffect,3,1)

        self.modelSteadyFlow.addItem(self.tr("steady flow"), "on")
        self.modelSteadyFlow.addItem(self.tr("unsteady flow"), "off")

        self.modelLagrangian.addItem(self.tr("single Phase Flow"), "single_phase")
        self.modelLagrangian.addItem(self.tr("particles and droplets tracking"), "lagrangian")

        self.modelAtmospheric.addItem(self.tr("off"             ), "off")
        self.modelAtmospheric.addItem(self.tr("constant density"), "constant")
        self.modelAtmospheric.addItem(self.tr("dry atmosphere"  ), "dry")
        self.modelAtmospheric.addItem(self.tr("humid atmosphere"), "humid")

        self.modelGasCombustionModel.addItem(self.tr("off"), "off")
        self.modelGasCombustionModel.addItem(self.tr("premixed flame (Eddy Break-Up)"), "ebu")
        self.modelGasCombustionModel.addItem(self.tr("equilibrium chemistry diffusion flame (Presumed PDF)"), "d3p")

        self.modelPulverizedCoal.addItem(self.tr("off"), "off")
        self.modelPulverizedCoal.addItem(self.tr("homogeneous approach"), "coal_homo")
        self.modelPulverizedCoal.addItem(self.tr("homogeneous approach (wet)"), "coal_homo2")
        # WARNING: the 'coal_lagr' model is deprecated
        #self.modelPulverizedCoal.addItem(self.tr("Gaseous phase coupling with Lagrangian coal transport"), "coal_lagr")

        self.modelJouleEffect.addItem(self.tr("off"), "off")
        self.modelJouleEffect.addItem(self.tr("Joule Effect"), "joule")
        self.modelJouleEffect.addItem(self.tr("Joule Effect and Lorentz Forces"), "arc")

        # Connect signals to slots

        self.connect(self.comboBoxSteadyFlow, SIGNAL("activated(const QString&)"), self.slotSteadyFlow)
        self.connect(self.comboBoxLagrangian, SIGNAL("activated(const QString&)"), self.slotLagrangian)
        self.connect(self.comboBoxAtmospheric, SIGNAL("activated(const QString&)"), self.slotAtmospheric)
        self.connect(self.comboBoxGasCombustionModel, SIGNAL("activated(const QString&)"), self.slotGasCombustionModel)
        self.connect(self.comboBoxPulverizedCoal, SIGNAL("activated(const QString&)"), self.slotPulverizedCoal)
        self.connect(self.comboBoxJouleEffect, SIGNAL("activated(const QString&)"), self.slotJouleEffect)

        # Initialize Widgets

        val = self.std.getSteadyFlowManagement()
        if val == 'on':
            self.modelSteadyFlow.setItem(str_model='on')
            self.modelLagrangian.disableItem(str_model='lagrangian')
        else:
            self.modelSteadyFlow.setItem(str_model='off')
            self.modelLagrangian.enableItem(str_model='lagrangian')

        val = self.lagr.getLagrangianStatus()
        if val == 'off':
            self.modelLagrangian.setItem(str_model='single_phase')
            self.modelSteadyFlow.enableItem(str_model='on')
        else:
            self.modelLagrangian.setItem(str_model='lagrangian')
            self.modelSteadyFlow.disableItem(str_model='on')

        #self.modelLagrangian.disableItem(str_model='lagrangian')

        val = self.atmo.getAtmosphericFlowsModel()
        self.modelAtmospheric.setItem(str_model=val)
        self.modelAtmospheric.disableItem(str_model='humid')

        self.modelGasCombustionModel.setItem(str_model='off') # to delete
        self.modelGasCombustionModel.disableItem(str_model='ebu') # to delete
        self.modelGasCombustionModel.disableItem(str_model='d3p') # to delete
        model = self.gas.getGasCombustionModel()
        self.modelGasCombustionModel.setItem(str_model=model)

        self.modelJouleEffect.setItem(str_model='off') # to delete
        self.modelJouleEffect.disableItem(str_model='joule') # to delete
        self.modelJouleEffect.disableItem(str_model='arc') # to delete
        elec = self.elect.getElectricalModel()
        self.modelJouleEffect.setItem(str_model=elec)

        #self.modelPulverizedCoal.disableItem(str_model='coal_homo') # to delete
        #self.modelPulverizedCoal.disableItem(str_model='coal_homo2') # to delete
        coal = self.pcoal.getCoalCombustionModel()
        self.modelPulverizedCoal.setItem(str_model=coal)

        # Multi-phase flow and coal combustion
        # WARNING: the 'coal_lagr' model is deprecated

        #if self.lagr.getLagrangianStatus() == 'off':
        #    self.modelLagrangian.setItem(str_model='single_phase')
        #    self.modelPulverizedCoal.disableItem(str_model='coal_lagr')
        #else:
        #    self.modelLagrangian.setItem(str_model='lagrangian')
        #    self.modelPulverizedCoal.enableItem(str_model='coal_lagr')

        #if self.pcoal.getCoalCombustionModel() == 'coal_lagr':
        #    self.modelLagrangian.setItem(str_model='lagrangian')
        #    self.modelLagrangian.disableItem(str_model='single_phase')

        # Compatibility between turbulence model and multi-phases flow model

        if self.turb.getTurbulenceModel() not in \
                ('off', 'k-epsilon', 'k-epsilon-PL',
                 'Rij-epsilon', 'Rij-SSG', 'v2f-phi', 'k-omega-SST'):
            self.modelLagrangian.setItem(str_model='single_phase')
            self.comboBoxLagrangian.setEnabled(False)

        # Compatibility between turbulence model and reactive flow models

        if self.turb.getTurbulenceModel() not in ('k-epsilon',
                                                  'k-epsilon-PL',
                                                  'Rij-epsilon',
                                                  'Rij-SSG',
                                                  'v2f-phi',
                                                  'k-omega-SST'):

            self.modelGasCombustionModel.setItem(str_model='off')
            self.modelPulverizedCoal.setItem(str_model='off')
            self.modelJouleEffect.setItem(str_model='off')

            self.modelGasCombustionModel.disableItem(str_model='ebu')
            self.modelGasCombustionModel.disableItem(str_model='d3p')

            self.modelJouleEffect.disableItem(str_model='joule')
            self.modelJouleEffect.disableItem(str_model='arc')

            self.modelPulverizedCoal.disableItem(str_model='coal_homo')
            self.modelPulverizedCoal.disableItem(str_model='coal_homo2')
            # WARNING: the 'coal_lagr' model is deprecated
            #self.modelPulverizedCoal.disableItem(str_model='coal_lagr')

            self.comboBoxGasCombustionModel.setEnabled(False)
            self.comboBoxPulverizedCoal.setEnabled(False)
            self.comboBoxJouleEffect.setEnabled(False)

        # Update the QComboBox

        flame = self.gas.getGasCombustionModel()
        coal  = self.pcoal.getCoalCombustionModel()
        joule = self.elect.getElectricalModel()
        atmospheric = self.atmo.getAtmosphericFlowsModel()

        self.modelGasCombustionModel.setItem(str_model=flame)
        self.modelPulverizedCoal.setItem(str_model=coal)
        self.modelJouleEffect.setItem(str_model=joule)
        self.modelAtmospheric.setItem(str_model=atmospheric)

        # If one model is turned on, the others are turned off

        if (flame, coal, joule, atmospheric) != ('off', 'off', 'off', 'off'):

            if flame == 'off':
                self.comboBoxGasCombustionModel.setEnabled(False)

            if coal == 'off':
                self.comboBoxPulverizedCoal.setEnabled(False)

            if joule == 'off':
                self.comboBoxJouleEffect.setEnabled(False)

            if atmospheric == 'off':
                self.comboBoxAtmospheric.setEnabled(False)

        # Update the Tree files and folders

        self.browser.configureTree(self.case)


    def __activateComboBox(self):
        """
        Private Method.
        Change to NORMAL the state of the reactive flow OptionMenu buttons.
        """
        self.comboBoxSteadyFlow.setEnabled(True)
        #self.comboBoxLagrangian.setEnabled(True)
        self.comboBoxGasCombustionModel.setEnabled(True)
        self.comboBoxPulverizedCoal.setEnabled(True)
        self.comboBoxJouleEffect.setEnabled(True)
        self.comboBoxAtmospheric.setEnabled(True)


    def __disableComboBox(self):
        """
        Private Method.
        Change to DISABLED the state of the reactive flow OptionMenu buttons.
        """
        #self.comboBoxSteadyFlow.setEnabled(False)
        #self.comboBoxLagrangian.setEnabled(False)
        self.comboBoxGasCombustionModel.setEnabled(False)
        self.comboBoxPulverizedCoal.setEnabled(False)
        self.comboBoxJouleEffect.setEnabled(False)
        self.comboBoxAtmospheric.setEnabled(False)


    def __stringModelFromCombo(self, name):
        """
        Private Method.
        Method to get the current item from a QComboBox and returns
        the correct string for the model
        """
        if not name in ['SteadyFlow',
                        'Lagrangian',
                        'Atmospheric',
                        'GasCombustionModel',
                        'PulverizedCoal',
                        'JouleEffect']:
            log.debug("__stringModelFromCombo() Incorrect name for QComboBox name")
            string = ""
        else:
            combo   = eval('self.comboBox' + name)
            dico    = eval('self.model' + name + '.dicoV2M')
            string  = dico[str(combo.currentText())]

        return string


    @pyqtSignature("const QString&")
    def slotSteadyFlow(self, text):
        """
        Private slot.
        Configure tree and update xmlfile beyond the steady or unsteady flow type.
        """
        log.debug("slotSteadyFlow")
        steady = self.__stringModelFromCombo('SteadyFlow')

        if steady == 'on':
            self.modelLagrangian.disableItem(str_model='lagrangian')
        else:
            self.modelLagrangian.enableItem(str_model='lagrangian')

        self.std.setSteadyFlowManagement(steady)
        self.browser.configureTree(self.case)


    @pyqtSignature("const QString&")
    def slotLagrangian(self, text):
        """
        Private slot.
        Put value beyond the multi-phase flow treatment is choosen or not.
        """
        if self.__stringModelFromCombo('Lagrangian') == 'single_phase':
            val = 'off'
        else:
            val = 'on'

        if val == 'off':
            self.modelSteadyFlow.enableItem(str_model='on')
            # WARNING: the 'coal_lagr' model is deprecated
            #self.modelPulverizedCoal.disableItem(str_model='coal_lagr')
        else:
            self.modelSteadyFlow.disableItem(str_model='on')
            # WARNING: the 'coal_lagr' model is deprecated
            #self.modelPulverizedCoal.enableItem(str_model='coal_lagr')

        self.lagr.setLagrangianStatus(val)
        self.browser.configureTree(self.case)


    @pyqtSignature("const QString&")
    def slotAtmospheric(self, text):
        """
        Called when the comboBoxAtmospheric changed
        """
        self.__activateComboBox()

        model = self.__stringModelFromCombo('Atmospheric')

        if model != 'off':
            # we inform that thermal scalar will be removed if it exists
            th_label = self.scal.getThermalScalarLabel()
            if th_label != '':
                title = self.tr("Warning")
                msg   = self.tr("This selection implies the destruction of the thermal scalar")
                QMessageBox.warning(self, title, msg)

            self.therm.setThermalModel('off')
            self.__disableComboBox()
            self.comboBoxAtmospheric.setEnabled(True)

        self.atmo.setAtmosphericFlowsModel(model)
        self.browser.configureTree(self.case)


    @pyqtSignature("const QString&")
    def slotGasCombustionModel(self, text):
        """
        Private slot.
        Binding method for gas combustion models.
        """
        self.__activateComboBox()

        model = self.__stringModelFromCombo('GasCombustionModel')

        if model != 'off':
            # we inform that thermal scalar will be removed if it exists
            th_label = self.scal.getThermalScalarLabel()
            if th_label != '':
                title = self.tr("Warning")
                msg   = self.tr("This selection implies the destruction of the thermal scalar")
                QMessageBox.warning(self, title, msg)

            self.therm.setThermalModel('off')
            self.__disableComboBox()
            self.comboBoxGasCombustionModel.setEnabled(True)

        self.gas.setGasCombustionModel(model)
        self.browser.configureTree(self.case)


    @pyqtSignature("const QString&")
    def slotPulverizedCoal(self, text):
        """
        Private slot.
        Binding method for pulverized coal combustion models
        """
        self.__activateComboBox()

        model = self.__stringModelFromCombo('PulverizedCoal')

        # WARNING: the 'coal_lagr' model is deprecated
#        if model == 'coal_lagr':
#            self.modelLagrangian.disableItem(str_model='single_phase')
#        else:
#            self.modelLagrangian.enableItem(str_model='single_phase')

        if model != 'off':
            # we inform that thermal scalar will be removed if it exists
            th_label = self.scal.getThermalScalarLabel()
            if th_label != '':
                title = self.tr("Warning")
                msg   = self.tr("This selection implies the destruction of the thermal scalar")
                QMessageBox.warning(self, title, msg)

            self.therm.setThermalModel('off')
            self.__disableComboBox()
            self.comboBoxPulverizedCoal.setEnabled(True)

        self.pcoal.setCoalCombustionModel(model)
        self.browser.configureTree(self.case)


    @pyqtSignature("const QString&")
    def slotJouleEffect(self, text):
        """
        Private slot.
        Binding method for electrical models
        """
        self.__activateComboBox()

        model = self.__stringModelFromCombo('JouleEffect')

        if model != 'off':
            # we inform that thermal scalar will be removed if it exists
            th_label = self.scal.getThermalScalarLabel()
            if th_label != '':
                title = self.tr("Warning")
                msg   = self.tr("This selection implies the destruction of the thermal scalar")
                QMessageBox.warning(self, title, msg)

            self.therm.setThermalModel('off')
            self.__disableComboBox()
            self.comboBoxJouleEffect.setEnabled(True)

        self.elect.setElectricalModel(model)
        self.browser.configureTree(self.case)


    def tr(self, text):
        """
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
