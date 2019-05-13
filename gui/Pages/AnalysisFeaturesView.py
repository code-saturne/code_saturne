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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam

import code_saturne.Base.QtPage as QtPage

from code_saturne.Pages.AnalysisFeaturesForm import Ui_AnalysisFeaturesForm
from code_saturne.model.TurbulenceModel import TurbulenceModel
from code_saturne.model.GasCombustionModel import GasCombustionModel
from code_saturne.model.CompressibleModel import CompressibleModel
from code_saturne.model.CoalCombustionModel import CoalCombustionModel
from code_saturne.model.ElectricalModel import ElectricalModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.ThermalRadiationModel import ThermalRadiationModel
from code_saturne.model.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.model.MainFieldsModel import MainFieldsModel

from code_saturne.model.LagrangianModel import LagrangianModel
from code_saturne.model.TurboMachineryModel import TurboMachineryModel
from code_saturne.model.MobileMeshModel import MobileMeshModel
from code_saturne.model.FansModel import FansStatus

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

        self.parent = parent

        import os
        self.icondir = os.path.dirname(os.path.abspath(__file__)) + '/../Base/'

        QWidget.__init__(self, parent)

        Ui_AnalysisFeaturesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.browser = tree

        self.case.undoStopGlobal()

        self.checkPrev = None

        # Set models and number of elements for combo boxes

        self.modelSinglePhase    = QtPage.ComboModel(self.comboBoxSinglePhase,3,1)
        self.modelAtmospheric    = QtPage.ComboModel(self.comboBoxAtmospheric,3,1)
        self.modelReactiveFlows  = QtPage.ComboModel(self.comboBoxReactiveFlows,2,1)
        self.modelGasCombustion  = QtPage.ComboModel(self.comboBoxGasCombustion,3,1)
        self.modelCoalCombustion = QtPage.ComboModel(self.comboBoxCoalCombustion,2,1)
        self.modelJouleEffect    = QtPage.ComboModel(self.comboBoxJouleEffect,2,1)
        self.modelGroundwater    = QtPage.ComboModel(self.comboBoxGroundwater,1,1)

        self.modelNeptuneCFD     = QtPage.ComboModel(self.comboBoxNeptuneCFD,2,1)

        self.modelSinglePhase.addItem(self.tr("Incompressible"), 'off')
        self.modelSinglePhase.addItem(self.tr("Compressible, perfect gas with constant gamma"),
                                      'constant_gamma')
        self.modelSinglePhase.addItem(self.tr("Compressible, perfect gas with variable gamma"),
                                      'variable_gamma')
        self.modelSinglePhase.disableItem(str_model='variable_gamma')

        self.modelAtmospheric.addItem(self.tr("constant density"), "constant")
        self.modelAtmospheric.addItem(self.tr("dry atmosphere"  ), "dry")
        self.modelAtmospheric.addItem(self.tr("humid atmosphere"), "humid")

        self.modelReactiveFlows.addItem(self.tr("Gas combustion"),
                                        "gas_combustion")
        self.modelReactiveFlows.addItem(self.tr("Pulverized Coal"),
                                        "pulverized_coal")

        self.modelGasCombustion.addItem(self.tr("perfect premixed flame (Eddy Break-Up)"),
                                        "ebu")
        self.modelGasCombustion.addItem(self.tr("infinitely fast chemistry diffusion flame"),
                                        "d3p")
        self.modelGasCombustion.addItem(self.tr("partial premixed flame (Libby_Williams)"),
                                        "lwp")

        self.modelCoalCombustion.addItem(self.tr("homogeneous approach"),
                                         "homogeneous_fuel")
        self.modelCoalCombustion.addItem(self.tr("homogeneous approach with moisture"),
                                         "homogeneous_fuel_moisture")

        self.modelJouleEffect.addItem(self.tr("Joule Effect"), "joule")
        self.modelJouleEffect.addItem(self.tr("Joule Effect and Laplace Forces"), "arc")

        self.modelGroundwater.addItem(self.tr("Groundwater flows"), 'groundwater')

        self.modelNeptuneCFD.addItem(self.tr("User-defined"), "None")
        self.modelNeptuneCFD.addItem(self.tr("Thermal free surface flow"), "free_surface")
        self.modelNeptuneCFD.addItem(self.tr("Boiling flow"), "boiling_flow")
        self.modelNeptuneCFD.addItem(self.tr("Droplet-laden flow"), "droplet_flow")
        self.modelNeptuneCFD.addItem(self.tr("Particle-laden flow"), "particles_flow")

        # Lagrangian combobox

        self.modelLagrangian = QtPage.ComboModel(self.comboBoxLagrangian,4,1)
        self.modelLagrangian.addItem(self.tr("off"),                 "off")
        self.modelLagrangian.addItem(self.tr("One-way coupling"),    "one_way")
        self.modelLagrangian.addItem(self.tr("Two-way coupling"),    "two_way")
        self.modelLagrangian.addItem(self.tr("Frozen carrier flow"), "frozen")

        # Turbomachinery combobox

        self.modelTurboMachinery = QtPage.ComboModel(self.comboBoxTurboMachinery, 5, 1)
        self.modelTurboMachinery.addItem(self.tr("None"), "off")
        self.modelTurboMachinery.addItem(self.tr("Full transient simulation"), "transient")
        self.modelTurboMachinery.addItem(self.tr("Transient with explicit coupling"), "transient_coupled")
        self.modelTurboMachinery.addItem(self.tr("Frozen rotor model"), "frozen")
        self.modelTurboMachinery.addItem(self.tr("Frozen rotor with explicit coupling"), "frozen_coupled")

        self.__uncheckRadioButtons()

        # Connect signals to slots

        self.comboBoxAtmospheric.activated[str].connect(self.slotAtmospheric)
        self.comboBoxReactiveFlows.activated[str].connect(self.slotReactiveFlows)
        self.comboBoxGasCombustion.activated[str].connect(self.slotGasCombustion)
        self.comboBoxCoalCombustion.activated[str].connect(self.slotCoalCombustion)
        self.comboBoxJouleEffect.activated[str].connect(self.slotJouleEffect)
        self.comboBoxSinglePhase.activated[str].connect(self.slotSinglePhase)
        self.comboBoxGroundwater.activated[str].connect(self.slotGroundwater)
        self.comboBoxNeptuneCFD.activated[str].connect(self.slotNeptuneCFD)
        self.checkBoxALE.stateChanged.connect(self.slotALE)
        self.checkBoxFans.stateChanged.connect(self.slotFans)

        self.comboBoxLagrangian.activated[str].connect(self.slotLagrangian)
        self.comboBoxTurboMachinery.activated[str].connect(self.slotTurboModel)

        for ind in ['SinglePhase', 'Atmospheric', 'JouleEffect',
                    'Groundwater', 'ReactiveFlows',
                    'NeptuneCFD']:
            eval('self.radioButton'+ind+'.toggled.connect(self.slotRadioButton)')

        # Initializations based on code

        if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
            self.init_saturne()

            # Enable NEPTUNE_CFD based on presence and/or environment
            # variable (environment variable has precedence
            # for ease of testing)
            enable_neptune_cfd = False
            try:
                from neptune_cfd.nc_package import package as nc_package
                enable_neptune_cfd = True
            except Exception:
                pass
            ev_enable = os.getenv('CS_GUI_ENABLE_NEPTUNE_CFD')
            if ev_enable:
                try:
                    if int(ev_enable) != 0:
                        enable_neptune_cfd = True
                    else:
                        enable_neptune_cfd = False
                except Exception:
                    pass
            self.radioButtonNeptuneCFD.setEnabled(enable_neptune_cfd)

        elif self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            self.init_neptune()

        self.tbm = TurboMachineryModel(self.case)
        self.ale = MobileMeshModel(self.case)
        self.fans = FansStatus(self.case)

        self.init_common()

        # Update the Tree files and folders

        self.browser.configureTree(self.case)

        self.case.undoStartGlobal()


    def __uncheckRadioButtons(self):

        self.__hideComboBox()

        for ind in ['SinglePhase', 'Atmospheric', 'JouleEffect',
                    'Groundwater', 'ReactiveFlows',
                    'NeptuneCFD']:
            eval('self.radioButton'+ind+'.setChecked(False)')


    def __hideComboBox(self):
        """
        Private Method.
        Change to DISABLED the state of the reactive flow OptionMenu buttons.
        """
        self.comboBoxSinglePhase.hide()
        self.comboBoxAtmospheric.hide()
        self.comboBoxReactiveFlows.hide()
        self.comboBoxGasCombustion.hide()
        self.comboBoxCoalCombustion.hide()
        self.comboBoxJouleEffect.hide()
        self.comboBoxGroundwater.hide()

        self.comboBoxNeptuneCFD.hide()


    def __stringModelFromCombo(self, name):
        """
        Private Method.
        Method to get the current item from a QComboBox and returns
        the correct string for the model
        """

        if not name in ['SinglePhase',
                        'Atmospheric',
                        'GasCombustion',
                        'CoalCombustion',
                        'JouleEffect',
                        'Groundwater',
                        'NeptuneCFD',
                        'Lagrangian',
                        'TurboMachinery']:
            log.debug("__stringModelFromCombo() Incorrect name for QComboBox name")
            string = ""
        else:
            combo   = eval('self.comboBox' + name)
            dico    = eval('self.model' + name + '.dicoV2M')
            string  = dico[str(combo.currentText())]

        return string


    def switch_case_to_saturne(self):

        if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
            return

        from code_saturne.cs_package import package as cs_package
        self.case['package'] = cs_package()

        from code_saturne.model.XMLinitialize import XMLinit

        self.case.root().xmlRemoveChild("additional_scalars")
        self.case.root().xmlRemoveChild("closure_modeling")
        self.case.root().xmlRemoveChild("thermophysical_models")
        self.case.root().xmlRemoveChild("numerical_parameters")

        self.case.xmlRootNode().tagName = "Code_Saturne_GUI"
        self.case.root()["solver_version"] = ""

        XMLinit(self.case).initialize()


    def init_saturne(self):

        self.__hideComboBox()

        self.turb  = TurbulenceModel(self.case)
        self.gas   = GasCombustionModel(self.case)
        self.pcoal = CoalCombustionModel(self.case)
        self.elect = ElectricalModel(self.case)
        self.scal  = DefineUserScalarsModel(self.case)
        self.atmo  = AtmosphericFlowsModel(self.case)
        self.comp  = CompressibleModel(self.case)
        self.darc  = GroundwaterModel(self.case)

        from code_saturne.model.TimeStepModel import TimeStepModel

        joule = self.elect.getElectricalModel()
        atmospheric = self.atmo.getAtmosphericFlowsModel()
        darcy = self.darc.getGroundwaterModel()
        coal = self.pcoal.getCoalCombustionModel()
        gas = self.gas.getGasCombustionModel()
        compressible = self.comp.getCompressibleModel()

        # Compatibility between turbulence model and reactive flow models
        if self.turb.getTurbulenceModel() not in self.turb.RANSmodels():
            if coal != 'off':
                self.pcoal.setCoalCombustionModel('off')
                coal = 'off'
            elif gas != 'off':
                self.gas.setGasCombustionModel('off')
                gas = 'off'

            self.radioButtonReactiveFlows.setEnabled(False)
        else:
            self.radioButtonReactiveFlows.setEnabled(True)

        # Set combobox values

        if atmospheric != 'off':
            self.modelAtmospheric.setItem(str_model=atmospheric)

        elif gas != 'off':
            self.modelReactiveFlows.setItem(str_model='gas_combustion')
            self.modelGasCombustion.setItem(str_model=gas)
            self.comboBoxGasCombustion.show()

        elif coal != 'off':
            self.modelReactiveFlows.setItem(str_model='pulverized_coal')
            self.modelCoalCombustion.setItem(str_model=coal)
            self.comboBoxCoalCombustion.show()

        elif joule != 'off':
            self.modelJouleEffect.setItem(str_model=joule)

        elif darcy != 'off':
            self.modelGroundwater.setItem(str_model=darcy)

            self.labelLagrangian.hide()
            self.comboBoxLagrangian.hide()
            self.labelTurboMachinery.hide()
            self.comboBoxTurboMachinery.hide()
            self.checkBoxALE.hide()
            self.checkBoxFans.hide()

        else:
            self.modelSinglePhase.setItem(str_model=compressible)

        # Set radiobutton, which also leads to combobox updates
        # through the matching signals and slot.

        Atmospheric = self.atmo.getAtmosphericFlowsModel() != 'off'
        JouleEffect = self.elect.getElectricalModel() != 'off'
        Groundwater = self.darc.getGroundwaterModel() != 'off'

        ReactiveFlows = self.gas.getGasCombustionModel()  != 'off' \
                        or self.pcoal.getCoalCombustionModel() != 'off'

        self.checkPrev = 'SinglePhase'
        combo = None

        for ind in ['Atmospheric',
                    'JouleEffect',
                    'Groundwater',
                    'ReactiveFlows']:

            radioButton = eval('self.radioButton'+ind)
            model_on = eval(ind)

            if model_on:
                base_model = False
                self.checkPrev = ind
                radioButton.setChecked(True)
            else:
                radioButton.setChecked(False)

        if self.checkPrev == 'SinglePhase':
            self.radioButtonSinglePhase.setChecked(True)

        # Lagrangian model features

        self.lagr  = LagrangianModel(self.case)

        lagr = self.lagr.getLagrangianModel()
        self.modelLagrangian.setItem(str_model=lagr)

        idtvar = TimeStepModel(self.case).getTimePassing()
        idtvar_p = idtvar

        if lagr != 'off':
            if idtvar not in [0, 1]:
                idtvar = 0

        if coal == 'homogeneous_fuel_moisture':
            self.modelLagrangian.disableItem(str_model='two_way')
        self.modelLagrangian.enableItem(str_model='frozen')

        if self.turb.getTurbulenceModel() not in \
                ('off', 'k-epsilon', 'k-epsilon-PL',
                 'Rij-epsilon', 'Rij-SSG', 'Rij-EBRSM', 'v2f-BL-v2/k',
                 'k-omega-SST', 'Spalart-Allmaras'):
            self.modelLagrangian.setItem(str_model='off')
            self.comboBoxLagrangian.setEnabled(False)
        else:
            self.comboBoxLagrangian.setEnabled(True)


    def switch_case_to_neptune(self):

        if self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            return

        from code_saturne.model.XMLinitializeNeptune import XMLinitNeptune

        self.case.root().xmlRemoveChild("thermophysical_models")
        self.case.root().xmlRemoveChild("numerical_parameters")
        self.case.root().xmlRemoveChild("time_parameters")
        self.case.root().xmlRemoveChild("analysis_control")
        self.case.root().xmlRemoveChild("physical_properties")
        self.case.root().xmlInitNode("physical_properties")

        self.case.xmlRootNode().tagName = "NEPTUNE_CFD_GUI"
        self.case.root()["solver_version"] = ""

        XMLinitNeptune(self.case).initialize()


    def init_neptune(self):

        self.__hideComboBox()

        self.comboBoxNeptuneCFD.show()

        self.nept = MainFieldsModel(self.case)

        predefined_flow = self.nept.getPredefinedFlow()
        self.modelNeptuneCFD.setItem(str_model=predefined_flow)

        self.radioButtonNeptuneCFD.setChecked(True)

        self.comboBoxLagrangian.setEnabled(True)

        self.checkPrev = 'NeptuneCFD'

        # Lagrangian model features

        self.lagr  = LagrangianModel(self.case)

        lagr = self.lagr.getLagrangianModel()
        self.modelLagrangian.setItem(str_model=lagr)

        self.modelLagrangian.disableItem(str_model='two_way')
        self.modelLagrangian.disableItem(str_model='frozen')

        self.comboBoxLagrangian.setEnabled(True)

        # Other features

        self.checkBoxALE.hide()
        self.checkBoxFans.hide()


    def init_common(self):

        from code_saturne.model.TimeStepModel import TimeStepModel

        # Turbomachinery features

        mdl = self.tbm.getTurboMachineryModel()
        self.modelTurboMachinery.setItem(str_model = mdl)

        # ALE

        ale_status = self.ale.getMethod()
        if ale_status == 'on':
            self.checkBoxALE.setChecked(True)
        else:
            self.checkBoxALE.setChecked(False)
            if self.ale.isMobileMeshCompatible():
                self.checkBoxALE.show()
            else:
                self.checkBoxALE.hide()

        # Fans

        count = self.fans.getFanCount()
        if count >= 0:
            self.checkBoxFans.setChecked(True)
            if count > 0:
                self.checkBoxFans.setEnabled(False)
        else:
            self.checkBoxFans.setEnabled(True)


        # Configure tree

        self.browser.configureTree(self.case)


    @pyqtSlot()
    def slotRadioButton(self):

        checkCur = None
        model_expr = None

        for ind in ['SinglePhase',
                    'Atmospheric',
                    'JouleEffect',
                    'Groundwater',
                    'ReactiveFlows',
                    'NeptuneCFD']:

            radioButton = eval('self.radioButton'+ind)
            combo = eval('self.comboBox'+ind)
            model = eval('self.model'+ind)

            if radioButton.isChecked():
                combo.show()
                model_expr = 'self.slot'+ind+'(model)'
                checkCur = ind
            else:
                combo.hide()
                model.setItem(0)

        # 2 signals may be sent on button toggles, so update
        # what is necessary.
        if self.checkPrev == checkCur:
            return

        # Restore visibility of features which may be masked with
        # some models.

        self.labelLagrangian.show()
        self.comboBoxLagrangian.show()
        self.labelTurboMachinery.show()
        self.comboBoxTurboMachinery.show()
        self.checkBoxALE.show()
        self.checkBoxFans.show()

        # Update

        if checkCur == 'NeptuneCFD':

            self.switch_case_to_neptune()

            self.init_neptune()

            try:
                from neptune_cfd.nc_package import package as nc_package
                self.case['package'] = nc_package()
            except Exception:
                from code_saturne.cs_package import package as cs_package
                self.case['package'] = cs_package()
                self.case['package'].code_name = 'NEPTUNE_CFD'

            if hasattr(self.parent, 'updateTitleBar'):
                self.parent.updateTitleBar()

        elif self.checkPrev == 'NeptuneCFD':

            self.switch_case_to_saturne()

            self.init_saturne()

            if hasattr(self.parent, 'updateTitleBar'):
                self.parent.updateTitleBar()

            eval('self.radioButton'+checkCur+'.setChecked(True)')

        else: # Both selected and previous models are Code_Saturne models.

            # setting of the previous model saturne to 'off'
            if self.checkPrev == 'SinglePhase':
                self.comp.setCompressibleModel('off')
            if self.checkPrev == 'Atmospheric':
                self.atmo.setAtmosphericFlowsModel('off')
            if self.checkPrev == 'JouleEffect':
                self.elect.setElectricalModel('off')
            if self.checkPrev == 'Groundwater':
                self.darc.setGroundwaterModel('off')
                self.labelLagrangian.show()
                self.comboBoxLagrangian.show()
                self.labelTurboMachinery.show()
                self.comboBoxTurboMachinery.show()
                self.checkBoxALE.show()
                self.checkBoxFans.show()
            if self.checkPrev == 'ReactiveFlows':
                self.gas.setGasCombustionModel('off')
                self.pcoal.setCoalCombustionModel('off')

                self.comboBoxGasCombustion.hide()
                self.comboBoxCoalCombustion.hide()

                self.modelLagrangian.enableItem(str_model='two_way')

            if checkCur == 'ReactiveFlows':
                if self.pcoal.getCoalCombustionModel() != 'off':
                    combo = self.comboBoxCoalCombustion
                    model = self.modelCoalCombustion
                else:
                    combo = self.comboBoxGasCombustion
                    model = self.modelGasCombustion
                combo.show()
                model_expr = 'self.slotReactiveFlows(model)'

            eval(model_expr)

        self.browser.configureTree(self.case)
        self.checkPrev = checkCur


    @pyqtSlot(str)
    def slotLagrangian(self, text):
        """
        Private slot.
        Put value beyond the multi-phase flow treatment is choosen or not.
        """
        model = self.__stringModelFromCombo('Lagrangian')

        self.lagr.setLagrangianModel(model)
        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotAtmospheric(self, text):
        """
        Called when the comboBoxAtmospheric changed
        """
        model = self.__stringModelFromCombo('Atmospheric')
        self.atmo.setAtmosphericFlowsModel(model)

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotJouleEffect(self, text):
        """
        Private slot.
        Binding method for electrical models
        """

        model = self.__stringModelFromCombo('JouleEffect')
        self.elect.setElectricalModel(model)

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotSinglePhase(self, text):
        """
        Private slot.
        Binding method for gas SinglePhase models.
        """

        model = self.__stringModelFromCombo('SinglePhase')
        self.comp.setCompressibleModel(model)

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotGroundwater(self, text):
        """
        Called when the comboBoxGroundwater changed
        """

        model = self.__stringModelFromCombo('Groundwater')
        self.darc.setGroundwaterModel(model)

        self.labelLagrangian.hide()
        self.comboBoxLagrangian.hide()
        self.labelTurboMachinery.hide()
        self.comboBoxTurboMachinery.hide()
        self.checkBoxALE.hide()
        self.checkBoxFans.hide()

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotNeptuneCFD(self, text):
        """
        Called when the comboBoxNeptuneCFD changed
        """
        model = self.__stringModelFromCombo('NeptuneCFD')
        self.nept.setPredefinedFlow(model)

        self.checkBoxALE.hide()
        self.checkBoxFans.hide()

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotReactiveFlows(self, text):
        """
        Called when the ReactiveFlows changed
        """
        model = str(self.comboBoxReactiveFlows.currentText())

        if model == "Gas combustion":
            self.pcoal.setCoalCombustionModel("off")
            self.comboBoxCoalCombustion.hide()
            self.comboBoxGasCombustion.show()
            model = self.__stringModelFromCombo('GasCombustion')
            self.slotGasCombustion(model)

        elif model == "Pulverized Coal":
            self.gas.setGasCombustionModel("off")
            self.comboBoxGasCombustion.hide()
            self.comboBoxCoalCombustion.show()
            model = self.__stringModelFromCombo('CoalCombustion')
            self.slotCoalCombustion(model)


    @pyqtSlot(str)
    def slotGasCombustion(self, text):
        """
        Private slot.
        Binding method for gas combustion models.
        """

        model = self.__stringModelFromCombo('GasCombustion')

        self.gas.setGasCombustionModel(model)

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotCoalCombustion(self, text):
        """
        Private slot.
        Binding method for coal combustion models.
        """

        model = self.__stringModelFromCombo('CoalCombustion')

        self.pcoal.setCoalCombustionModel(model)

        if model == 'homogeneous_fuel_moisture':
            self.modelLagrangian.disableItem(str_model='two_way')
        else:
            self.modelLagrangian.enableItem(str_model='two_way')

        self.browser.configureTree(self.case)


    @pyqtSlot(str)
    def slotTurboModel(self, text):
        """
        Input turbomachinery model.
        """
        mdl = self.modelTurboMachinery.dicoV2M[str(text)]
        self.tbm.setTurboMachineryModel(mdl)

        self.browser.configureTree(self.case)


    @pyqtSlot(int)
    def slotALE(self, val):

        if val == 0:
            self.ale.setMethod ("off")
        else:
            self.ale.setMethod ("on")

        self.browser.configureTree(self.case)


    @pyqtSlot(int)
    def slotFans(self, val):

        if val == 0:
            self.fans.cleanFans()
        else:
            self.fans.initFans()

        self.browser.configureTree(self.case)


    def tr(self, text):
        """
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
