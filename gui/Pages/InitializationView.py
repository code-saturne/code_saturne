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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, ComboModel

from code_saturne.Pages.InitializationForm import Ui_InitializationForm
from code_saturne.Pages.TurbulenceModel import TurbulenceModel
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from code_saturne.Pages.InitializationModel import InitializationModel
from code_saturne.Pages.CompressibleModel import CompressibleModel
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.GroundwaterModel import GroundwaterModel
from code_saturne.Pages.NotebookModel import NotebookModel

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
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InitializationForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.case.undoStopGlobal()

        self.init    = InitializationModel(self.case)
        self.turb    = TurbulenceModel(self.case)
        self.therm   = ThermalScalarModel(self.case)
        self.th_sca  = DefineUserScalarsModel(self.case)
        self.comp    = CompressibleModel(self.case)
        self.volzone = LocalizationModel('VolumicZone', self.case)
        self.notebook = NotebookModel(self.case)

        # create group to control hide/show options
        self.turb_group = [self.labelTurbulence, self.pushButtonTurbulence,
                           self.comboBoxTurbulence]
        self.thermal_group = [self.labelThermal, self.pushButtonThermal]
        self.velocity_group = [self.labelVelocity, self.pushButtonVelocity]
        self.species_group = [self.labelSpecies, self.comboBoxSpecies, self.pushButtonSpecies]
        self.meteo_group =   [self.labelMeteo, self.comboBoxMeteo, self.pushButtonMeteo]
        self.thermodynamic_list = ['Pressure', 'Density', 'Temperature', 'Energy']

        # 1/ Combo box models

        self.modelZone = ComboModel(self.comboBoxZone, 1, 1)
        if self.comp.getCompressibleModel() != 'off':
            self.groupBoxThermodynamic.show()
        else:
            self.groupBoxThermodynamic.hide()

        self.zone = ""
        zones = self.volzone.getZones()
        for zone in zones:
            if zone.getNature()['initialization'] == "on":
                label = zone.getLabel()
                name = str(zone.getCodeNumber())
                self.modelZone.addItem(self.tr(label), name)
                if label == "all_cells":
                    self.zone = name
                if not self.zone:
                    self.zone = name

        self.modelZone.setItem(str_model = self.zone)

        self.modelTurbulence = ComboModel(self.comboBoxTurbulence, 2, 1)
        self.modelTurbulence.addItem(self.tr("Initialization by formula"), 'formula')
        self.modelTurbulence.addItem(self.tr("Initialization by reference value(s)"), 'reference_value')

        # 2/ Connections

        self.comboBoxZone.activated[str].connect(self.slotZone)
        self.comboBoxTurbulence.activated[str].connect(self.slotChoice)
        self.comboBoxSpecies.activated[str].connect(self.slotSpeciesChoice)
        self.comboBoxMeteo.activated[str].connect(self.slotMeteoChoice)
        self.checkBoxPressure.clicked.connect(self.slotPressure)
        self.checkBoxDensity.clicked.connect(self.slotDensity)
        self.checkBoxTemperature.clicked.connect(self.slotTemperature)
        self.checkBoxEnergy.clicked.connect(self.slotEnergy)
        self.pushButtonVelocity.clicked.connect(self.slotVelocityFormula)
        self.pushButtonThermal.clicked.connect(self.slotThermalFormula)
        self.pushButtonTurbulence.clicked.connect(self.slotTurbulenceFormula)
        self.pushButtonSpecies.clicked.connect(self.slotSpeciesFormula)
        self.pushButtonMeteo.clicked.connect(self.slotMeteoFormula)
        self.pushButtonPressure.clicked.connect(self.slotPressureFormula)
        self.pushButtonDensity.clicked.connect(self.slotDensityFormula)
        self.pushButtonTemperature.clicked.connect(self.slotTemperatureFormula)
        self.pushButtonEnergy.clicked.connect(self.slotEnergyFormula)
        self.pushButtonHydraulicHead.clicked.connect(self.slotHydraulicHeadFormula)

        choice = self.init.getInitialTurbulenceChoice(self.zone)
        self.modelTurbulence.setItem(str_model = choice)

        # species treatment
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
            self.modelSpecies.setItem(str_model = self.scalar)
            exp = self.init.getSpeciesFormula(self.zone, self.scalar)
            if exp:
                self.pushButtonSpecies.setStyleSheet("background-color: green")
                self.pushButtonSpecies.setToolTip(exp)
            else:
                self.pushButtonSpecies.setStyleSheet("background-color: red")
        else:
            for item in self.species_group:
                item.hide()

        # meteo
        self.modelMeteo = ComboModel(self.comboBoxMeteo, 1, 1)
        self.scalar_meteo = ""
        scalar_meteo_list = DefineUserScalarsModel( self.case).getMeteoScalarsNameList()
        if scalar_meteo_list != None and scalar_meteo_list != []:
            self.scalar_meteo = scalar_meteo_list[0]
            for item in self.meteo_group:
                item.show()
            for scalar in scalar_meteo_list:
                self.modelMeteo.addItem(self.tr(scalar), scalar)
            self.modelMeteo.setItem(str_model = self.scalar_meteo)
            exp = self.init.getMeteoFormula(self.zone, self.scalar_meteo)
            if exp:
                self.pushButtonMeteo.setStyleSheet("background-color: green")
                self.pushButtonMeteo.setToolTip(exp)
            else:
                self.pushButtonMeteo.setStyleSheet("background-color: red")
        else:
            for item in self.meteo_group:
                item.hide()

        if GroundwaterModel(self.case).getGroundwaterModel() == "off":
            self.labelHydraulicHead.hide()
            self.pushButtonHydraulicHead.hide()
        else:
            exp = self.init.getHydraulicHeadFormula(self.zone)
            if exp:
                self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
                self.pushButtonHydraulicHead.setToolTip(exp)
            else:
                self.pushButtonHydraulicHead.setStyleSheet("background-color: red")

        # Initialize widget
        self.initializeVariables(self.zone)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotZone(self, text):
        """
        INPUT label for choice of zone
        """
        self.zone = self.modelZone.dicoV2M[str(text)]
        self.initializeVariables(self.zone)


    @pyqtSlot(str)
    def slotChoice(self, text):
        """
        INPUT choice of method of initialization
        """
        choice = self.modelTurbulence.dicoV2M[str(text)]
        log.debug("slotChoice choice =  %s "%str(choice))
        self.init.setInitialTurbulenceChoice(self.zone, choice)
        turb_model = self.turb.getTurbulenceModel()

        self.initializeVariables(self.zone)


    @pyqtSlot(str)
    def slotMeteoChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.scalar_meteo= self.modelMeteo.dicoV2M[str(text)]
        self.initializeVariables(self.zone)
        exp = self.init.getMeteoFormula(self.zone, self.scalar_meteo)
        if exp:
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(exp)
        else:
            self.pushButtonMeteo.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def slotSpeciesChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.scalar= self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables(self.zone)
        exp = self.init.getSpeciesFormula(self.zone, self.scalar)
        if exp:
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(exp)
        else:
            self.pushButtonSpecies.setStyleSheet("background-color: red")


    @pyqtSlot()
    def slotVelocityFormula(self):
        """
        """
        exp = self.init.getVelocityFormula(self.zone)
        if not exp:
            exp = self.init.getDefaultVelocityFormula()
        exa = """#example: \n""" + self.init.getDefaultVelocityFormula()
        req = [('velocity[0]', "velocity"),
               ('velocity[1]', "velocity"),
               ('velocity[2]', "velocity")]
        sym = [('uref', 'reference velocity'),
               ('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.init.setVelocityFormula(self.zone, str(result))
            self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.pushButtonVelocity.setToolTip(result)


    @pyqtSlot()
    def slotTurbulenceFormula(self):
        """
        INPUT user formula
        """
        turb_model = self.turb.getTurbulenceModel()
        exa = """#example \n""" + self.init.getDefaultTurbFormula(turb_model)
        exp = self.init.getTurbFormula(self.zone, turb_model)
        sym = [('rho0', 'density (reference value)'),
               ('mu0', 'viscosity (reference value)'),
               ('cp0', 'specific heat (reference value)'),
               ('lambda0', 'thermal conductivity (reference value)'),
               ('x','cell center coordinate'),
               ('y','cell center coordinate'),
               ('z','cell center coordinate'),
               ('uref','reference velocity'),
               ('almax','reference length')]
        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            req = [('k', "turbulent energy"),
                   ('epsilon', "turbulent dissipation")]
        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            req = [('r11', "Reynolds stress R11"),
                   ('r22', "Reynolds stress R22"),
                   ('r33', "Reynolds stress R33"),
                   ('r12', "Reynolds stress R12"),
                   ('r23', "Reynolds stress R23"),
                   ('r13', "Reynolds stress R13"),
                   ('epsilon', "turbulent dissipation")]
        elif turb_model == 'Rij-EBRSM':
            req = [('r11', "Reynolds stress R11"),
                   ('r22', "Reynolds stress R22"),
                   ('r33', "Reynolds stress R33"),
                   ('r12', "Reynolds stress R12"),
                   ('r23', "Reynolds stress R23"),
                   ('r13', "Reynolds stress R13"),
                   ('epsilon', "turbulent dissipation"),
                   ('alpha', "alpha")]
        elif turb_model == 'v2f-BL-v2/k':
            req = [('k', "turbulent energy"),
                   ('epsilon', "turbulent dissipation"),
                   ('phi', "variable phi in v2f model"),
                   ('alpha', "variable alpha in v2f model")]
        elif turb_model == 'k-omega-SST':
            req = [('k', "turbulent energy"),
                   ('omega', "specific dissipation rate")]
        elif turb_model == 'Spalart-Allmaras':
            req = [('nu_tilda', "nusa")]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTurb -> %s" % str(result))
            self.init.setTurbFormula(self.zone, str(result))
            self.pushButtonTurbulence.setStyleSheet("background-color: green")
            self.pushButtonTurbulence.setToolTip(result)


    @pyqtSlot()
    def slotThermalFormula(self):
        """
        Input the initial formula of thermal scalar
        """
        exp = self.init.getThermalFormula(self.zone)
        if not exp:
            exp = self.init.getDefaultThermalFormula()
        exa = """#example \n""" + self.init.getDefaultThermalFormula()
        if self.therm.getThermalScalarModel() == "enthalpy":
            req = [('enthalpy', 'enthalpy')]
        elif self.therm.getThermalScalarModel() == "total_energy":
            req = [('total_energy', 'total energy')]
        else:
            req = [('temperature', 'temperature')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermal -> %s" % str(result))
            self.init.setThermalFormula(self.zone, str(result))
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(result)


    @pyqtSlot()
    def slotSpeciesFormula(self):
        """
        Input the initial formula of species
        """
        exp = self.init.getSpeciesFormula(self.zone, self.scalar)
        name = self.th_sca.getScalarName(self.scalar)
        if not exp:
            exp = str(name)+""" = 0;\n"""
        exa = """#example: \n""" + str(name)+""" = 0;\n"""
        req = [(str(name), str(name))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaSpecies -> %s" % str(result))
            self.init.setSpeciesFormula(self.zone, self.scalar, str(result))
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(result)


    @pyqtSlot()
    def slotMeteoFormula(self):
        """
        """
        exp = self.init.getMeteoFormula(self.zone, self.scalar_meteo)
        name = self.scalar_meteo
        if not exp:
            exp = str(name)+""" = 0;\n"""
        exa = """#example: \n""" + str(name)+""" = 0;\n"""
        req = [(str(name), str(name))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMeteo -> %s" % str(result))
            self.init.setMeteoFormula(self.zone, self.scalar_meteo, str(result))
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(result)


    @pyqtSlot()
    def slotPressure(self):
        """
        Pressure selected or not for the initialisation.
        """
        if self.checkBoxPressure.isChecked():
            self.init.setPressureStatus(self.zone,"on")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonPressure.setEnabled(True)
            exp = self.init.getPressureFormula(self.zone)
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
            self.init.setPressureStatus(self.zone,"off")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonPressure.setEnabled(False)
            self.pushButtonPressure.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Pressure':
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(True)
                if box_list[0] =='Energy':
                    self.checkBoxTemperature.setEnabled(False)
                if box_list[0] =='Temperature':
                    self.checkBoxEnergy.setEnabled(False)


    @pyqtSlot()
    def slotDensity(self):
        """
        Density selected or not for the initialisation.
        """
        if self.checkBoxDensity.isChecked():
            self.init.setDensityStatus(self.zone,"on")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonDensity.setEnabled(True)
            exp = self.init.getDensityFormula(self.zone)
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
            self.init.setDensityStatus(self.zone,"off")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonDensity.setEnabled(False)
            self.pushButtonDensity.setStyleSheet("background-color: None")
            if len(box_list) == 1:
                for name in self.thermodynamic_list:
                    if name != 'Density':
                        __checkBox = getattr(self, "checkBox" + name)
                        __checkBox.setEnabled(True)
                if box_list[0] =='Energy':
                    self.checkBoxTemperature.setEnabled(False)
                if box_list[0] =='Temperature':
                    self.checkBoxEnergy.setEnabled(False)


    @pyqtSlot()
    def slotTemperature(self):
        """
        Temperature selected or not for the initialisation.
        """
        if self.checkBoxTemperature.isChecked():
            self.init.setTemperatureStatus(self.zone,"on")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonTemperature.setEnabled(True)
            exp = self.init.getTemperatureFormula(self.zone)
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
            self.init.setTemperatureStatus(self.zone,"off")
            box_list = self.init.getCheckedBoxList(self.zone)
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
        if self.checkBoxEnergy.isChecked():
            self.init.setEnergyStatus(self.zone,"on")
            box_list = self.init.getCheckedBoxList(self.zone)
            self.pushButtonEnergy.setEnabled(True)
            exp = self.init.getEnergyFormula(self.zone)
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
            self.init.setEnergyStatus(self.zone,"off")
            box_list = self.init.getCheckedBoxList(self.zone)
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
        exp = self.init.getPressureFormula(self.zone)
        if not exp:
            exp = """p0 = 0.;
g = 9.81;
ro = 1.17862;
pressure = p0 + g * ro * z;\n"""
        exa = """#example: """
        req = [('pressure', 'pressure')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotPressureFormula -> %s" % str(result))
            self.init.setPressureFormula(self.zone, str(result))
            self.pushButtonPressure.setStyleSheet("background-color: green")
            self.pushButtonPressure.setToolTip(result)



    @pyqtSlot()
    def slotHydraulicHeadFormula(self):
        """
        Input the initial Hydraulic Head formula
        """
        exp = self.init.getHydraulicHeadFormula(self.zone)
        if not exp:
            exp = """H = z;\n"""
        exa = """#example: """
        req = [('H', 'hydraulic head')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotHydraulicHeadFormula -> %s" % str(result))
            self.init.setHydraulicHeadFormula(self.zone, str(result))
            self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
            self.pushButtonHydraulicHead.setToolTip(result)



    @pyqtSlot()
    def slotDensityFormula(self):
        """
        Input the initial Density formula
        """
        exp = self.init.getDensityFormula(self.zone)
        if not exp:
            exp = """density = 0;\n"""
        exa = """#example: """
        req = [('density', 'density')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotDensityFormula -> %s" % str(result))
            self.init.setDensityFormula(self.zone, str(result))
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(result)



    @pyqtSlot()
    def slotTemperatureFormula(self):
        """
        Input the initial Temperature formula
        """
        exp = self.init.getTemperatureFormula(self.zone)
        if not exp:
            exp = """temperature = 0;\n"""
        exa = """#example: """
        req = [('temperature', 'temperature')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotTemperatureFormula -> %s" % str(result))
            self.init.setTemperatureFormula(self.zone, str(result))
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(result)



    @pyqtSlot()
    def slotEnergyFormula(self):
        """
        Input the initial Energy formula
        """
        exp = self.init.getEnergyFormula(self.zone)
        if not exp:
            exp = """total_energy = 0;\n"""
        exa = """#example: """
        req = [('total_energy', 'Energy')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotEnergyFormula -> %s" % str(result))
            self.init.setEnergyFormula(self.zone, str(result))
            self.pushButtonEnergy.setStyleSheet("background-color: green")
            self.pushButtonEnergy.setToolTip(result)


    def initializeVariables(self, zone):
        """
        Initialize variables when a new volumic zone is choosen
        """
        # Initialisation of Turbulence

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
        else:
            for item in self.turb_group:
                item.show()

            turb_init = self.init.getInitialTurbulenceChoice(self.zone)
            self.modelTurbulence.setItem(str_model = turb_init)

            if turb_init == 'formula':
                self.pushButtonTurbulence.setEnabled(True)
                turb_formula = self.init.getTurbFormula(zone, turb_model)
                if not turb_formula:
                    turb_formula = self.init.getDefaultTurbFormula(turb_model)
                    self.pushButtonTurbulence.setStyleSheet("background-color: red")
                else:
                    self.pushButtonTurbulence.setStyleSheet("background-color: green")
                self.init.setTurbFormula(zone, turb_formula)
                self.pushButtonTurbulence.setToolTip(turb_formula)
            else:
                self.pushButtonTurbulence.setEnabled(False)
                self.pushButtonTurbulence.setStyleSheet("background-color: None")

        #velocity
        if GroundwaterModel(self.case).getGroundwaterModel() == "groundwater":
            for item in self.velocity_group:
                item.hide()
        else:
            velocity_formula = self.init.getVelocityFormula(zone)
            if not velocity_formula:
                velocity_formula = self.init.getDefaultVelocityFormula()
                self.pushButtonVelocity.setStyleSheet("background-color: red")
            else:
                self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.init.setVelocityFormula(zone, velocity_formula)
            self.pushButtonVelocity.setToolTip(velocity_formula)

        # Initialisation of Model Variables if thermal model is selectionned
        for item in self.thermal_group:
            item.hide()

        model = self.therm.getThermalScalarModel()

        if model != "off" and self.comp.getCompressibleModel() == 'off':
            for item in self.thermal_group:
                item.show()
            th_formula = self.init.getThermalFormula(zone)
            if not th_formula:
                th_formula = self.init.getDefaultThermalFormula()
                self.pushButtonThermal.setStyleSheet("background-color: red")
            else:
                self.pushButtonThermal.setStyleSheet("background-color: green")
            self.init.setThermalFormula(zone, th_formula)
            self.pushButtonThermal.setToolTip(th_formula)

        # Initialisation of the termodynamics values for the compressible model
        if self.comp.getCompressibleModel() != 'off':
            nb_box = 0
            box_list = self.init.getCheckedBoxList(self.zone)
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
