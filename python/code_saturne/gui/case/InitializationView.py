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
This module contains the following class:
- InitializationView
"""

#-------------------------------------------------------------------------------
# Standard modules
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
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, ComboModel

from code_saturne.gui.case.InitializationForm import Ui_InitializationForm
from code_saturne.model.TurbulenceModel import TurbulenceModel
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from code_saturne.model.InitializationModel import InitializationModel
from code_saturne.model.CompressibleModel import CompressibleModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.model.HgnModel import HgnModel
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.GasCombustionModel import GasCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InitializationView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class InitializationView(QWidget, Ui_InitializationForm):
    """
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InitializationForm.__init__(self)
        self.setupUi(self)

        self.case = None
        self.zone = None
        self.parent = parent
        self.is_solid  = False

        # create group to control hide/show options
        self.hydraulic_head_group = [self.labelHydraulicHead, self.checkBoxHydraulicHead,
                                     self.pushButtonHydraulicHead]
        self.velocity_group = [self.labelVelocity, self.checkBoxVelocity,
                               self.pushButtonVelocity]
        self.void_fraction_group = [self.labelVoidFraction, self.checkBoxVoidFraction,
                               self.pushButtonVoidFraction]
        self.turb_group = [self.labelTurbulence, self.pushButtonTurbulence,
                           self.comboBoxTurbulence]
        self.thermal_group = [self.labelThermal, self.comboBoxThermal,
                              self.pushButtonThermal]
        self.species_group = [self.labelSpecies, self.comboBoxSpecies,
                              self.checkBoxSpecies, self.pushButtonSpecies,]
        self.meteo_group = [self.labelMeteo, self.comboBoxMeteo,
                            self.checkBoxMeteo, self.pushButtonMeteo]
        self.combustion_group = [self.labelCombustion, self.comboBoxCombustion,
                                 self.checkBoxCombustion, self.pushButtonCombustion]

        self.thermodynamic_list = ['Pressure', 'Density', 'Temperature', 'Energy']

        self.modelThermal = ComboModel(self.comboBoxThermal, 2, 1)
        self.modelThermal.addItem(self.tr("Automatic initialization"), 'automatic')
        self.modelThermal.addItem(self.tr("Initialization by formula"), 'formula')

        self.modelTurbulence = ComboModel(self.comboBoxTurbulence, 2, 1)
        self.modelTurbulence.addItem(self.tr("Initialization by formula"), 'formula')
        self.modelTurbulence.addItem(self.tr("Initialization by reference value(s)"),
                                     'reference_value')


    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()
        self.comp = CompressibleModel(self.case)
        self.th_sca = DefineUserScalarsModel(self.case)
        self.init = InitializationModel(self.case)
        self.turb = TurbulenceModel(self.case)
        self.therm = ThermalScalarModel(self.case)
        self.hgn = HgnModel(self.case)
        self.notebook = NotebookModel(self.case)

        # Identify selected zone_id
        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone

        if self.zone.isNatureActivated("initialization"):
            if self.zone.isNatureActivated("solid"):
                self.is_solid = True
            self.setViewFromCase()
        else:
            # Should not arrive here, as the whole containing tab is normally
            # hidden in this case.
            self.setEnabled(False)
        self.case.undoStartGlobal()


    def setViewFromCase(self):

        zone_id = str(self.zone.getCodeNumber())
        choice = self.init.getInitialTurbulenceChoice(zone_id)

        # Thermodynamic for compressible flows
        if self.comp.getCompressibleModel() != 'off':
            self.groupBoxThermodynamic.show()
        else:
            self.groupBoxThermodynamic.hide()

        # Velocity or groundwater flow
        if GroundwaterModel(self.case).getGroundwaterModel() == "groundwater":
            for item in self.velocity_group:
                item.hide()
            self.updateViewHydraulicHead(zone_id)
        else:
            for item in self.hydraulic_head_group:
                item.hide()
            self.updateViewVelocity(zone_id)

        # Void fraction
        self.initViewVoidFraction(zone_id)

        # Thermal
        self.initViewThermal(zone_id)

        # Turbulence
        self.initViewTurbulence(zone_id)

        # Species treatment
        self.modelSpecies = ComboModel(self.comboBoxSpecies, 1, 1)
        self.scalar = ""
        scalar_list = self.th_sca.getUserScalarNameList()
        for s in self.th_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)
        if scalar_list != []:
            self.scalar = scalar_list[0]
            for item in self.species_group:
                item.show()
            for scalar in scalar_list:
                self.modelSpecies.addItem(self.tr(scalar), scalar)
            self.modelSpecies.setItem(str_model=self.scalar)
            self.updateViewSpecies(zone_id)
        else:
            for item in self.species_group:
                item.hide()

        # Meteo scalars
        self.modelMeteo = ComboModel(self.comboBoxMeteo, 1, 1)
        self.scalar_meteo = ""
        scalar_meteo_list = DefineUserScalarsModel(self.case).getMeteoScalarsNameList()
        if scalar_meteo_list != None and scalar_meteo_list != []:
            self.scalar_meteo = scalar_meteo_list[0]
            for item in self.meteo_group:
                item.show()
            for scalar in scalar_meteo_list:
                self.modelMeteo.addItem(self.tr(scalar), scalar)
            self.modelMeteo.setItem(str_model=self.scalar_meteo)
            self.updateViewMeteo(zone_id)
        else:
            for item in self.meteo_group:
                item.hide()

        # Combustion
        self.modelCombustion = ComboModel(self.comboBoxCombustion, 1, 1)
        self.scalar_combustion = ""
        scalar_combustion_list = DefineUserScalarsModel( self.case).getGasCombScalarsNameList()
        if GasCombustionModel(self.case).getGasCombustionModel() == "d3p":
            # For the d3p model (option extended),
            # we allow only the Automatic Initialization for the enthalpy
            if not self.is_solid:
                option = GasCombustionModel(self.case).getGasCombustionOption()
                if option == 'extended':
                    self.modelThermal.disableItem(str_model = 'formula')
            self.scalar_combustion = scalar_combustion_list[0]
            for item in self.combustion_group:
                item.show()
            for scalar in scalar_combustion_list:
                self.modelCombustion.addItem(self.tr(scalar), scalar)
            self.modelCombustion.setItem(str_model = self.scalar_combustion)
            self.updateViewCombustion(zone_id)
        else:
            for item in self.combustion_group:
                item.hide()

        # Compressible
        self.initializeCompressibleVariables()

        # Initialize widget
        self.defineConnections()


    def defineConnections(self):
        self.comboBoxThermal.activated[str].connect(self.slotThermalChoice)
        self.comboBoxTurbulence.activated[str].connect(self.slotTurbulenceChoice)
        self.comboBoxSpecies.activated[str].connect(self.slotSpeciesChoice)
        self.comboBoxMeteo.activated[str].connect(self.slotMeteoChoice)
        self.comboBoxCombustion.activated[str].connect(self.slotCombustionChoice)
        self.checkBoxPressure.clicked.connect(self.slotPressure)
        self.checkBoxDensity.clicked.connect(self.slotDensity)
        self.checkBoxTemperature.clicked.connect(self.slotTemperature)
        self.checkBoxEnergy.clicked.connect(self.slotEnergy)
        self.checkBoxSpecies.clicked.connect(self.slotSpeciesCheck)
        self.checkBoxMeteo.clicked.connect(self.slotMeteoCheck)
        self.checkBoxCombustion.clicked.connect(self.slotCombustionCheck)
        self.checkBoxHydraulicHead.clicked.connect(self.slotHydraulicHeadCheck)
        self.pushButtonVelocity.clicked.connect(self.slotVelocityFormula)
        self.checkBoxVelocity.clicked.connect(self.slotVelocityChoice)
        self.pushButtonVoidFraction.clicked.connect(self.slotVoidFractionFormula)
        self.checkBoxVoidFraction.clicked.connect(self.slotVoidFractionChoice)
        self.pushButtonThermal.clicked.connect(self.slotThermalFormula)
        self.pushButtonTurbulence.clicked.connect(self.slotTurbulenceFormula)
        self.pushButtonSpecies.clicked.connect(self.slotSpeciesFormula)
        self.pushButtonMeteo.clicked.connect(self.slotMeteoFormula)
        self.pushButtonCombustion.clicked.connect(self.slotCombustionFormula)
        self.pushButtonPressure.clicked.connect(self.slotPressureFormula)
        self.pushButtonDensity.clicked.connect(self.slotDensityFormula)
        self.pushButtonTemperature.clicked.connect(self.slotTemperatureFormula)
        self.pushButtonEnergy.clicked.connect(self.slotEnergyFormula)
        self.pushButtonHydraulicHead.clicked.connect(self.slotHydraulicHeadFormula)


    def updateViewHydraulicHead(self, zone_id):
        """
        Upate species view
        """
        exp = self.init.getHydraulicHeadFormula(zone_id)
        if exp:
            self.checkBoxHydraulicHead.setChecked(True)
            self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
            self.pushButtonHydraulicHead.show()
            self.pushButtonHydraulicHead.setToolTip(exp)
        else:
            self.checkBoxHydraulicHead.setChecked(False)
            self.pushButtonHydraulicHead.hide()


    def updateViewVelocity(self, zone_id):
        """
        Upate species view
        """
        exp = self.init.getVelocityFormula(zone_id)
        if exp:
            self.checkBoxVelocity.setChecked(True)
            self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.pushButtonVelocity.show()
            self.pushButtonVelocity.setToolTip(exp)
        else:
            self.checkBoxVelocity.setChecked(False)
            self.pushButtonVelocity.hide()


    def initViewVoidFraction(self, zone_id):
        """
        Update void fraction view
        """
        model = self.hgn.getHgnModel()

        if model == "off":
            for item in self.void_fraction_group:
                item.hide()
            return

        exp = self.init.getVoidFractionFormula(zone_id)
        if exp:
            self.checkBoxVoidFraction.setChecked(True)
            self.pushButtonVoidFraction.setStyleSheet("background-color: green")
            self.pushButtonVoidFraction.show()
            self.pushButtonVoidFraction.setToolTip(exp)
        else:
            self.checkBoxVoidFraction.setChecked(False)
            self.pushButtonVoidFraction.hide()


    def initViewTurbulence(self, zone_id):
        """
        Update turbulence view
        """
        turb_model = self.turb.getTurbulenceModel()

        if turb_model not in ('k-epsilon',
                              'k-epsilon-PL',
                              'Rij-epsilon',
                              'Rij-SSG',
                              'Rij-EBRSM',
                              'v2f-BL-v2/k',
                              'k-omega-SST',
                              'Spalart-Allmaras'):
            for item in self.turb_group:
                item.hide()
            return

        for item in self.turb_group:
            item.show()

        turb_init = self.init.getInitialTurbulenceChoice(zone_id)
        self.modelTurbulence.setItem(str_model=turb_init)

        if turb_init == 'formula':
            self.pushButtonTurbulence.show()
            exp = self.init.getTurbFormula(zone_id, turb_model)
            if not exp:
                exp = self.init.getDefaultTurbFormula(turb_model)
                self.pushButtonTurbulence.setEnabled(True)
                self.pushButtonTurbulence.setStyleSheet("background-color: red")
            else:
                self.pushButtonTurbulence.setEnabled(True)
                self.pushButtonTurbulence.setStyleSheet("background-color: green")
            self.init.setTurbFormula(zone_id, exp)
            self.pushButtonTurbulence.setToolTip(exp)
        else:
            self.pushButtonTurbulence.hide()


    def initViewThermal(self, zone_id):
        """
        INPUT choice of method of initialization
        """
        model = self.therm.getThermalScalarModel()

        if model == "off" or self.comp.getCompressibleModel() != 'off':
            for item in self.thermal_group:
                item.hide()
            return

        th_formula = self.init.getThermalFormula(zone_id)
        # If no formula, then automatic initialization
        if not th_formula:
            self.pushButtonThermal.hide()
            self.modelThermal.setItem(str_model="automatic")
        else:
            self.modelThermal.setItem(str_model="formula")
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(th_formula)


    def updateViewSpecies(self, zone_id):
        """
        Upate species view
        """
        exp = self.init.getSpeciesFormula(zone_id, self.scalar)
        if exp:
            self.checkBoxSpecies.setChecked(True)
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(exp)
            self.pushButtonSpecies.show()
        else:
            self.checkBoxSpecies.setChecked(False)
            self.pushButtonSpecies.hide()


    def updateViewMeteo(self, zone_id):
        """
        Upate meteo scalars view
        """
        exp = self.init.getMeteoFormula(zone_id, self.scalar_meteo)
        if exp:
            self.checkBoxMeteo.setChecked(True)
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(exp)
            self.pushButtonMeteo.show()
        else:
            self.checkBoxMeteo.setChecked(False)
            self.pushButtonMeteo.hide()


    def updateViewCombustion(self, zone_id):
        """
        Upate meteo scalars view
        """
        exp = self.init.getCombustionFormula(zone_id, self.scalar_combustion)
        if exp:
            self.checkBoxCombustion.setChecked(True)
            self.pushButtonCombustion.setStyleSheet("background-color: green")
            self.pushButtonCombustion.setToolTip(exp)
            self.pushButtonCombustion.show()
        else:
            self.checkBoxCombustion.setChecked(False)
            self.pushButtonCombustion.hide()


    @pyqtSlot(str)
    def slotTurbulenceChoice(self, text):
        """
        INPUT choice of method of initialization
        """
        choice = self.modelTurbulence.dicoV2M[str(text)]
        log.debug("slotTurbulenceChoice choice =  %s " % str(choice))
        zone_id = str(self.zone.getCodeNumber())
        self.init.setInitialTurbulenceChoice(zone_id, choice)
        turb_model = self.turb.getTurbulenceModel()
        if choice == 'formula':
            self.pushButtonTurbulence.show()
            turb_formula = self.init.getTurbFormula(zone_id, turb_model)
            if not turb_formula:
                turb_formula = self.init.getDefaultTurbFormula(turb_model)
                self.pushButtonTurbulence.setStyleSheet("background-color: red")
            else:
                self.pushButtonTurbulence.setStyleSheet("background-color: green")
                self.init.setTurbFormula(zone_id, turb_formula)
                self.pushButtonTurbulence.setToolTip(turb_formula)
        else:
            self.pushButtonTurbulence.hide()
        if choice == 'reference_value':
            self.pushButtonTurbulence.setDisabled(True)
            self.pushButtonTurbulence.setStyleSheet("background-color: None")
            self.init.setTurbFormula(zone_id, None)


    @pyqtSlot(str)
    def slotThermalChoice(self, text):
        """
        INPUT choice of method of initialization
        """
        zone_id = str(self.zone.getCodeNumber())
        choice = self.modelThermal.dicoV2M[str(text)]
        log.debug("slotThermalChoice choice =  %s " % str(choice))
        if choice == 'formula':
            th_formula = self.init.getThermalFormula(zone_id)
            if not th_formula:
                th_formula = self.init.getDefaultThermalFormula()
                self.pushButtonThermal.setStyleSheet("background-color: red")
            else:
                self.pushButtonThermal.setStyleSheet("background-color: green")
            self.init.setThermalFormula(zone_id, th_formula)
            self.pushButtonThermal.show()
        else:
            self.init.setThermalFormula(zone_id, None)
            self.pushButtonThermal.hide()


    @pyqtSlot(str)
    def slotSpeciesChoice(self, text):
        """
        INPUT label for choice of zone_id
        """
        self.scalar = self.modelSpecies.dicoV2M[str(text)]
        zone_id = str(self.zone.getCodeNumber())
        self.updateViewSpecies(zone_id)


    @pyqtSlot(str)
    def slotMeteoChoice(self, text):
        """
        INPUT label for choice of zone_id
        """
        self.scalar_meteo = self.modelMeteo.dicoV2M[str(text)]
        zone_id = str(self.zone.getCodeNumber())
        self.updateViewMeteo(zone_id)


    @pyqtSlot(str)
    def slotCombustionChoice(self, text):
        """
        INPUT label for choice of zone_id
        """
        self.scalar_combustion = self.modelCombustion.dicoV2M[str(text)]
        zone_id = str(self.zone.getCodeNumber())
        self.updateViewCombustion(zone_id)


    @pyqtSlot()
    def slotVelocityFormula(self):
        """
        INPUT choice of method of initialization
        """
        exa = """#example: \n""" + self.init.getDefaultVelocityFormula()

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getVelocityFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="velocity",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.init.setVelocityFormula(zone_id, str(result))
            self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.pushButtonVelocity.setToolTip(result)


    @pyqtSlot()
    def slotVelocityChoice(self):
        zone_id = str(self.zone.getCodeNumber())
        exp = self.init.getVelocityFormula(zone_id)
        if self.checkBoxVelocity.isChecked():
            self.pushButtonVelocity.show()
            self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.pushButtonVelocity.setToolTip(exp)
            self.init.setVelocityFormula(zone_id, exp)
        else:
            self.pushButtonVelocity.hide()
            self.init.setVelocityFormula(zone_id, None)


    @pyqtSlot()
    def slotVoidFractionFormula(self):
        """
        INPUT choice of method of initialization
        """
        exa = """#example: \n""" + self.init.getDefaultVoidFractionFormula()

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getVoidFractionFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="void_fraction",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVoidFraction -> %s" % str(result))
            self.init.setVoidFractionFormula(zone_id, str(result))
            self.pushButtonVoidFraction.setStyleSheet("background-color: green")
            self.pushButtonVoidFraction.setToolTip(result)


    @pyqtSlot()
    def slotVoidFractionChoice(self):
        zone_id = str(self.zone.getCodeNumber())
        exp = self.init.getVoidFractionFormula(zone_id)
        if self.checkBoxVoidFraction.isChecked():
            self.pushButtonVoidFraction.show()
            self.pushButtonVoidFraction.setStyleSheet("background-color: green")
            self.pushButtonVoidFraction.setToolTip(exp)
            self.init.setVoidFractionFormula(zone_id, exp)
        else:
            self.pushButtonVoidFraction.hide()
            self.init.setVoidFractionFormula(zone_id, None)


    @pyqtSlot()
    def slotHydraulicHeadCheck(self):
        zone_id = str(self.zone.getCodeNumber())

        if self.checkBoxHydraulicHead.isChecked():
            exp = self.init.getHydraulicHeadFormula(zone_id)
            self.pushButtonHydraulicHead.show()
            self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
            self.pushButtonHydraulicHead.setToolTip(exp)
            self.init.setHydraulicHeadFormula(zone_id, exp)
        else:
            self.pushButtonHydraulicHead.hide()
            self.init.setHydraulicHeadFormula(zone_id, None)


    @pyqtSlot()
    def slotSpeciesCheck(self):
        zone_id = str(self.zone.getCodeNumber())
        self.scalar = self.modelSpecies.dicoV2M[self.comboBoxSpecies.currentText()]
        if self.checkBoxSpecies.isChecked():
            exp = self.init.getSpeciesFormula(zone_id, self.scalar)
            self.pushButtonSpecies.show()
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(exp)
            self.init.setSpeciesFormula(zone_id, self.scalar, exp)
        else:
            self.pushButtonSpecies.hide()
            self.init.setSpeciesFormula(zone_id, self.scalar, None)


    @pyqtSlot()
    def slotMeteoCheck(self):
        zone_id = str(self.zone.getCodeNumber())
        self.scalar = self.modelMeteo.dicoV2M[self.comboBoxMeteo.currentText()]
        if self.checkBoxMeteo.isChecked():
            exp = self.init.getMeteoFormula(zone_id, self.scalar_meteo)
            self.pushButtonMeteo.show()
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(exp)
            self.init.setMeteoFormula(zone_id, self.scalar_meteo, exp)
        else:
            self.pushButtonMeteo.hide()
            self.init.setMeteoFormula(zone_id, self.scalar, None)


    @pyqtSlot()
    def slotCombustionCheck(self):
        zone_id = str(self.zone.getCodeNumber())
        self.scalar = self.modelCombustion.dicoV2M[self.comboBoxCombustion.currentText()]
        if self.checkBoxCombustion.isChecked():
            species_formula = self.init.getCombustionFormula(zone_id, self.scalar)
            self.pushButtonCombustion.show()
            self.pushButtonCombustion.setStyleSheet("background-color: green")
            self.init.setCombustionFormula(zone_id, self.scalar, species_formula)
        else:
            self.pushButtonCombustion.hide()
            self.pushButtonCombustion.setStyleSheet("background-color: None")
            self.init.setCombustionFormula(zone_id, self.scalar, None)


    @pyqtSlot()
    def slotTurbulenceFormula(self):
        """
        INPUT user formula
        """
        turb_model = self.turb.getTurbulenceModel()
        exa = """#example \n""" + self.init.getDefaultTurbFormula(turb_model)

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getTurbFormulaComponents(zone_id, turb_model)

        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            turb_vname = 'turbulence_ke'
        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            turb_vname = 'turbulence_rije'
        elif turb_model == 'Rij-EBRSM':
            turb_vname = 'turbulence_rij_ebrsm'
        elif turb_model == 'v2f-BL-v2/k':
            turb_vname = 'turbulence_v2f'
        elif turb_model == 'k-omega-SST':
            turb_vname = 'turbulence_kw'
        elif turb_model == 'Spalart-Allmaras':
            turb_vname = 'turbulence_spalart'

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name=turb_vname,
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTurb -> %s" % str(result))
            self.init.setTurbFormula(zone_id, str(result))
            self.pushButtonTurbulence.setStyleSheet("background-color: green")
            self.pushButtonTurbulence.setToolTip(result)


    @pyqtSlot()
    def slotThermalFormula(self):
        """
        Input the initial formula of thermal scalar
        """
        exa = """#example \n""" + self.init.getDefaultThermalFormula()

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym, knf = self.init.getThermalFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="thermal",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                known_fields=knf,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermal -> %s" % str(result))
            self.init.setThermalFormula(zone_id, str(result))
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(result)


    @pyqtSlot()
    def slotSpeciesFormula(self):
        """
        Input the initial formula of species
        """
        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym, knf = self.init.getSpeciesFormulaComponents(zone_id, self.scalar)

        name = DefineUserScalarsModel(self.case).getScalarName(self.scalar)
        exa = """#example: \n""" + str(name) + """ = 0;\n"""

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name=name,
                                expression=exp,
                                required=req,
                                symbols=sym,
                                known_fields=knf,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaSpecies -> %s" % str(result))
            self.init.setSpeciesFormula(zone_id, self.scalar, str(result))
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(result)


    @pyqtSlot()
    def slotMeteoFormula(self):
        """
        """
        name = self.scalar_meteo
        exa = """#example: \n""" + str(name) + """ = 0;\n"""

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getMeteoFormulaComponents(zone_id, self.scalar_meteo)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name=name,
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMeteo -> %s" % str(result))
            self.init.setMeteoFormula(zone_id, self.scalar_meteo, str(result))
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(result)


    @pyqtSlot()
    def slotCombustionFormula(self):
        """
        """
        name = self.scalar_combustion
        exa = """#example: \n""" + str(name) + """ = 0;\n"""

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getCombustionFormulaComponents(zone_id,
                                                                 self.scalar_combustion)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name=name,
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaCombustion -> %s" % str(result))
            self.init.setCombustionFormula(zone_id, self.scalar_combustion, str(result))
            self.pushButtonCombustion.setStyleSheet("background-color: green")
            self.pushButtonCombustion.setToolTip(result)


    @pyqtSlot()
    def slotPressure(self):
        """
        Pressure selected or not for the initialisation.
        """
        zone_id = str(self.zone.getCodeNumber())

        if self.checkBoxPressure.isChecked():
            self.init.setPressureStatus(zone_id, "on")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonPressure.setEnabled(True)
            exp = self.init.getPressureFormula(zone_id)
            if exp:
                self.pushButtonPressure.setStyleSheet("background-color: green")
                self.pushButtonPressure.setToolTip(exp)
            else:
                self.pushButtonPressure.setStyleSheet("background-color: red")
            if len(box_list) == 2:
                for name in self.thermodynamic_list:
                    if name not in box_list:
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(False)
        else:
            self.init.setPressureStatus(zone_id, "off")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonPressure.setEnabled(False)
            self.pushButtonPressure.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Pressure':
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(True)
                if box_list[0] == 'Energy':
                    self.checkBoxTemperature.setEnabled(False)
                if box_list[0] == 'Temperature':
                    self.checkBoxEnergy.setEnabled(False)

    @pyqtSlot()
    def slotDensity(self):
        """
        Density selected or not for the initialisation.
        """
        zone_id = str(self.zone.getCodeNumber())

        if self.checkBoxDensity.isChecked():
            self.init.setDensityStatus(zone_id, "on")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonDensity.setEnabled(True)
            exp = self.init.getDensityFormula(zone_id)
            if exp:
                self.pushButtonDensity.setStyleSheet("background-color: green")
                self.pushButtonDensity.setToolTip(exp)
            else:
                self.pushButtonDensity.setStyleSheet("background-color: red")
            if len(box_list) == 2:
                for name in self.thermodynamic_list:
                    if name not in box_list:
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(False)
        else:
            self.init.setDensityStatus(zone_id, "off")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonDensity.setEnabled(False)
            self.pushButtonDensity.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Density':
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(True)
                if box_list[0] == 'Energy':
                    self.checkBoxTemperature.setEnabled(False)
                if box_list[0] == 'Temperature':
                    self.checkBoxEnergy.setEnabled(False)

    @pyqtSlot()
    def slotTemperature(self):
        """
        Temperature selected or not for the initialisation.
        """
        zone_id = str(self.zone.getCodeNumber())

        if self.checkBoxTemperature.isChecked():
            self.init.setTemperatureStatus(zone_id, "on")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonTemperature.setEnabled(True)
            exp = self.init.getTemperatureFormula(zone_id)
            if exp:
                self.pushButtonTemperature.setStyleSheet("background-color: green")
                self.pushButtonTemperature.setToolTip(exp)
            else:
                self.pushButtonTemperature.setStyleSheet("background-color: red")
            if len(box_list) == 2:
                for name in self.thermodynamic_list:
                    if name not in box_list:
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(False)
            self.checkBoxEnergy.setEnabled(False)
        else:
            self.init.setTemperatureStatus(zone_id, "off")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonTemperature.setEnabled(False)
            self.pushButtonTemperature.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Temperature':
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(True)
            self.checkBoxEnergy.setEnabled(True)

    @pyqtSlot()
    def slotEnergy(self):
        """
        Energy selected or not for the initialisation.
        """
        zone_id = str(self.zone.getCodeNumber())

        if self.checkBoxEnergy.isChecked():
            self.init.setEnergyStatus(zone_id, "on")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonEnergy.setEnabled(True)
            exp = self.init.getEnergyFormula(zone_id)
            if exp:
                self.pushButtonEnergy.setStyleSheet("background-color: green")
                self.pushButtonEnergy.setToolTip(exp)
            else:
                self.pushButtonEnergy.setStyleSheet("background-color: red")
            if len(box_list) == 2:
                for name in self.thermodynamic_list:
                    if name not in box_list:
                        __checkBox = getattr(self, "checkBox" + name)
                        __Button = getattr(self, "pushButton" + name)
                        __checkBox.setEnabled(False)
                        __Button.setEnabled(False)
                        __Button.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                self.checkBoxTemperature.setEnabled(False)
        else:
            self.init.setEnergyStatus(zone_id, "off")
            box_list = self.init.getCheckedBoxList(zone_id)
            self.pushButtonEnergy.setEnabled(False)
            self.pushButtonEnergy.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Energy':
                        __checkBox = getattr(self, "checkBox" + name)
                        __Button = getattr(self, "pushButton" + name)
                        __checkBox.setEnabled(True)
                        __Button.setEnabled(False)
                        __Button.setStyleSheet("background-color: None")
            self.checkBoxTemperature.setEnabled(True)


    @pyqtSlot()
    def slotPressureFormula(self):
        """
        Input the initial Pressure formula
        """
        exa = """#example: """

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getPressureFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="pressure",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotPressureFormula -> %s" % str(result))
            self.init.setPressureFormula(zone_id, str(result))
            self.pushButtonPressure.setStyleSheet("background-color: green")
            self.pushButtonPressure.setToolTip(result)


    @pyqtSlot()
    def slotHydraulicHeadFormula(self):
        """
        Input the initial Hydraulic Head formula
        """
        exa = """#example: """

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getHydraulicHeadFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="hydraulic_head",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotHydraulicHeadFormula -> %s" % str(result))
            self.init.setHydraulicHeadFormula(zone_id, str(result))
            self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
            self.pushButtonHydraulicHead.setToolTip(result)


    @pyqtSlot()
    def slotDensityFormula(self):
        """
        Input the initial Density formula
        """
        exa = """#example: """

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getDensityFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="density",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotDensityFormula -> %s" % str(result))
            self.init.setDensityFormula(zone_id, str(result))
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(result)


    @pyqtSlot()
    def slotTemperatureFormula(self):
        """
        Input the initial Temperature formula
        """
        exa = """#example: """

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getTemperatureFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="temperature",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotTemperatureFormula -> %s" % str(result))
            self.init.setTemperatureFormula(zone_id, str(result))
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(result)


    @pyqtSlot()
    def slotEnergyFormula(self):
        """
        Input the initial Energy formula
        """
        exa = """#example: """

        zone_id = str(self.zone.getCodeNumber())
        exp, req, sym = self.init.getEnergyFormulaComponents(zone_id)
        exp, req, sym = self.init.getEnergyFormulaComponents(zone_id)

        dialog = QMegEditorView(parent=self,
                                function_type="ini",
                                zone_name=self.zone.getLabel(),
                                variable_name="energy",
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotEnergyFormula -> %s" % str(result))
            self.init.setEnergyFormula(zone_id, str(result))
            self.pushButtonEnergy.setStyleSheet("background-color: green")
            self.pushButtonEnergy.setToolTip(result)


    def initializeCompressibleVariables(self):
        """
        Initialisation of the thermodynamics values for the compressible model
        """
        if self.comp.getCompressibleModel() != 'off':
            zone_id = str(self.zone.getCodeNumber())
            nb_box = 0
            box_list = self.init.getCheckedBoxList(zone_id)
            if box_list == []:
                for name in self.thermodynamic_list:
                    __checkBox = getattr(self, "checkBox" + name)
                    __Button = getattr(self, "pushButton" + name)
                    __checkBox.setChecked(False)
                    __Button.setEnabled(False)
                    __Button.setStyleSheet("background-color: None")
            elif len(box_list) == 1:
                box = box_list[0]
                for name in self.thermodynamic_list:
                    if name != box:
                        __checkBox = getattr(self, "checkBox" + name)
                        __Button = getattr(self, "pushButton" + name)
                        __checkBox.setChecked(False)
                        __Button.setEnabled(False)
                        __Button.setStyleSheet("background-color: None")
                if box == 'Temperature':
                    self.checkBoxEnergy.setEnabled(False)
                elif box == 'Energy':
                    self.checkBoxTemperature.setEnabled(False)
                __checkBox = getattr(self, "checkBox" + box)
                __checkBox.setChecked(True)
                __Button = getattr(self, "pushButton" + box)
                __Button.setEnabled(True)
                __Button.setStyleSheet("background-color: red")
            elif len(box_list) == 2:
                box1 = box_list[0]
                box2 = box_list[1]
                for name in self.thermodynamic_list:
                    if name not in box_list:
                        __checkBox = getattr(self, "checkBox" + name)
                        __Button = getattr(self, "pushButton" + name)
                        __checkBox.setChecked(False)
                        __checkBox.setEnabled(False)
                        __Button.setEnabled(False)
                for name in box_list:
                    __checkBox = getattr(self, "checkBox" + name)
                    __Button = getattr(self, "pushButton" + name)
                    __checkBox.setChecked(True)
                    __Button.setEnabled(True)
                    __Button.setStyleSheet("background-color: red")


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
