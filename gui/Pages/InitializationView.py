# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.InitializationForm import Ui_InitializationForm

from Base.Toolbox import GuiParam
from Base.QtPage import IntValidator, DoubleValidator, ComboModel, setGreenColor
from Pages.TurbulenceModel import TurbulenceModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from Pages.InitializationModel import InitializationModel
from Pages.QMeiEditorView import QMeiEditorView

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

        self.init    = InitializationModel(self.case)
        self.turb    = TurbulenceModel(self.case)
        self.therm   = ThermalScalarModel(self.case)
        self.th_sca  = DefineUserScalarsModel(self.case)
        self.volzone = LocalizationModel('VolumicZone', self.case)

        # create group to control hide/show options
        self.turb_group = [self.labelTurbulence, self.pushButtonTurbulence,
                           self.comboBoxTurbulence]
        self.thermal_group = [self.labelThermal, self.pushButtonThermal]
        self.species_group = [self.labelSpecies, self.comboBoxSpecies, self.pushButtonSpecies]

        # 1/ Combo box models

        self.modelZone = ComboModel(self.comboBoxZone, 1, 1)

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

        self.connect(self.comboBoxZone,         SIGNAL("activated(const QString&)"),   self.slotZone)
        self.connect(self.comboBoxTurbulence,   SIGNAL("activated(const QString&)"),   self.slotChoice)
        self.connect(self.comboBoxSpecies,      SIGNAL("activated(const QString&)"),   self.slotSpeciesChoice)
        self.connect(self.pushButtonVelocity,   SIGNAL("clicked()"),                   self.slotVelocityFormula)
        self.connect(self.pushButtonThermal,    SIGNAL("clicked()"),                   self.slotThermalFormula)
        self.connect(self.pushButtonTurbulence, SIGNAL("clicked()"),                   self.slotTurbulenceFormula)
        self.connect(self.pushButtonSpecies,    SIGNAL("clicked()"),                   self.slotSpeciesFormula)

        # Define thermal variable if needed
        th_sca_label = ''
        model = self.therm.getThermalScalarModel()
        if model != 'off':
            th_sca_label = self.therm.getThermalScalarLabel()

        self.th_sca_label = th_sca_label

        choice = self.init.getInitialTurbulenceChoice(self.zone)
        self.modelTurbulence.setItem(str_model = choice)

        # species treatment
        self.modelSpecies = ComboModel(self.comboBoxSpecies, 1, 1)
        self.scalar = ""
        scalar_list = self.th_sca.getUserScalarLabelsList()
        for s in self.th_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)

        if scalar_list != []:
            self.scalar = scalar_list[0]
            for item in self.species_group:
                item.show()
            for scalar in scalar_list:
                self.modelSpecies.addItem(self.tr(scalar),self.scalar)
            self.modelSpecies.setItem(str_model = self.scalar)
            setGreenColor(self.pushButtonSpecies, True)
        else:
            for item in self.species_group:
                item.hide()

        # Initialize widget
        self.initializeVariables(self.zone)


    @pyqtSignature("const QString&")
    def slotZone(self, text):
        """
        INPUT label for choice of zone
        """
        self.zone = self.modelZone.dicoV2M[str(text)]
        self.initializeVariables(self.zone)


    @pyqtSignature("const QString&")
    def slotChoice(self, text):
        """
        INPUT choice of method of initialization
        """
        choice = self.modelTurbulence.dicoV2M[str(text)]
        log.debug("slotChoice choice =  %s "%str(choice))
        self.init.setInitialTurbulenceChoice(self.zone, choice)
        turb_model = self.turb.getTurbulenceModel()

        self.initializeVariables(self.zone)


    @pyqtSignature("const QString&")
    def slotSpeciesChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.scalar= self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables(self.zone)
        setGreenColor(self.pushButtonSpecies, True)


    @pyqtSignature("const QString&")
    def slotVelocityFormula(self):
        """
        """
        exp = self.init.getVelocityFormula(self.zone)
        if not exp:
            exp = """u = 0;\nv = 0;\nw = 0;\n"""
        exa = """#example: """
        req = [('u', "x velocity"),
        ('v', "y velocity"),
        ('w', "z velocity")]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')] #quel symbol
        dialog = QMeiEditorView(self,expression = exp,
                                 required   = req,
                                 symbols    = sym,
                                 examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.init.setVelocityFormula(self.zone, result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotTurbulenceFormula(self):
        """
        INPUT user formula
        """
        turb_model = self.turb.getTurbulenceModel()
        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#example

cmu = 0.42;
k = 1.5*(0.02*uref)^2;
eps = k^1.5*cmu/almax;

"""
            req = [('k', "turbulent energy"),
            ('eps', "turbulent dissipation")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)

        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#exemple :
trii   = (0.02*uref)^2;

cmu = 0.42;

r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.d;
r23 = 0.;
k = 0.5*(r11+r22+r33);
eps = k^1.5*cmu/almax;"""
            req = [('r11', "Reunolds stress R11"),
            ('r22', "Reynolds stress R22"),
            ('r33', "Reynolds stress R33"),
            ('r12', "Reynolds stress R12"),
            ('r23', "Reynolds stress R23"),
            ('r13', "Reynolds stress R13"),
            ('eps', "turbulent dissipation")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)

        elif turb_model == 'Rij-EBRSM':
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#exemple :
trii   = (0.02*uref)^2;

cmu = 0.42;

r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.d;
r23 = 0.;
k = 0.5*(r11+r22+r33);
eps = k^1.5*cmu/almax;
alpha = 1.;"""
            req = [('r11', "Reunolds stress R11"),
            ('r22', "Reynolds stress R22"),
            ('r33', "Reynolds stress R33"),
            ('r12', "Reynolds stress R12"),
            ('r23', "Reynolds stress R23"),
            ('r13', "Reynolds stress R13"),
            ('eps', "turbulent dissipation"),
            ('alpha', "alpha")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)

        elif turb_model == 'v2f-phi':
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#example

cmu = 0.22;
k = 1.5*(0.02*uref)^2;
eps = k^1.5*cmu/almax;

phi = 2./3.;
fb = 0.;
"""
            req = [('k', "turbulent energy"),
            ('eps', "turbulent dissipation"),
            ('phi', "variable phi in v2f model"),
            ('fb', "variable f in v2f model")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)

        elif turb_model == 'k-omega-SST':
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#exemple :
k = 1.5*(0.02*uref)^2;
omega = k^0.5/almax;"""
            req = [('k', "turbulent energy"),
            ('omega', "specific dissipation rate")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)

        elif turb_model == 'Spalart-Allmaras':
            exp = self.init.getTurbFormula(self.zone, turb_model)
            exa = """#exemple :
nusa = (cmu * k)/eps;;"""
            req = [('nusa', "nusa")]
            sym = [('rho0', 'density (reference value)'),
                   ('mu0', 'viscosity (reference value)'),
                   ('cp0', 'specific heat (reference value)'),
                   ('lambda0', 'thermal conductivity (reference value)'),
                   ('x','cell center coordinate'),
                   ('y','cell center coordinate'),
                   ('z','cell center coordinate'),
                   ('uref','reference velocity'),
                   ('almax','reference length')]
            dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.init.setTurbFormula(self.zone, result)
                setGreenColor(self.sender(), False)


    @pyqtSignature("const QString&")
    def slotThermalFormula(self):
        """
        Input the initial formula of thermal scalar
        """
        exp = self.init.getThermalFormula(self.zone, self.th_sca_label)
        if not exp:
            exp = self.th_sca_label+""" = 0;\n"""
        exa = """#example: """
        req = [(self.th_sca_label, str(self.th_sca_label))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]
        dialog = QMeiEditorView(self,expression = exp,
                                 required   = req,
                                 symbols    = sym,
                                 examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermal -> %s" % str(result))
            self.init.setThermalFormula(self.zone, self.th_sca_label, result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("const QString&")
    def slotSpeciesFormula(self):
        """
        Input the initial formula of species
        """
        exp = self.init.getSpeciesFormula(self.zone, self.scalar)
        if not exp:
            exp = str(self.scalar)+""" = 0;\n"""
        exa = """#example: """
        req = [(str(self.scalar), str(self.scalar))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]
        dialog = QMeiEditorView(self,expression = exp,
                                 required   = req,
                                 symbols    = sym,
                                 examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaSpecies -> %s" % str(result))
            self.init.setSpeciesFormula(self.zone, self.scalar, result)
            setGreenColor(self.sender(), False)


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
                              'v2f-phi',
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
                self.init.setTurbFormula(zone, turb_formula)
                setGreenColor(self.pushButtonTurbulence, True)
            else:
                self.pushButtonTurbulence.setEnabled(False)
                setGreenColor(self.pushButtonTurbulence, False)

        #velocity
        velocity_formula = self.init.getVelocityFormula(zone)
        if not velocity_formula:
            velocity_formula = """u = 0;\nv = 0;\nw = 0;\n"""
        self.init.setVelocityFormula(zone, velocity_formula)
        setGreenColor(self.pushButtonVelocity, True)

        # Initialisation of Model Variables if thermal model is selectionned

        model = self.therm.getThermalScalarModel()

        if model == "off":
            for item in self.thermal_group:
                item.hide()
        else:
            for item in self.thermal_group:
                item.show()
            th_formula = self.init.getThermalFormula(zone, self.th_sca_label)
            if not th_formula:
                th_formula = self.th_sca_label+""" = 0;\n"""
            self.init.setThermalFormula(zone, self.th_sca_label, th_formula)
            setGreenColor(self.pushButtonThermal, True)


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
