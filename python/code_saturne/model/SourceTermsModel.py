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
This module initialize model dynamics variables and model scalars

This module contents the following classes:
- SourceTermsModel
- SourceTermsTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Variables and Scalar model initialization modelling class
#-------------------------------------------------------------------------------

class SourceTermsModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.models = self.case.xmlGetNode('thermophysical_models')
        self.node_userscalar = self.case.xmlGetNode('additional_scalars')
        self.node_veloce   = self.models.xmlGetNode('velocity_pressure')
        self.node_turb     = self.models.xmlGetNode('turbulence', 'model')
        self.node_therm    = self.models.xmlGetNode('thermal_scalar', 'model')
        self.node_sterm   = self.models.xmlInitNode('source_terms')

        self.therm   = ThermalScalarModel(self.case)
        self.th_sca  = DefineUserScalarsModel(self.case)
        self.notebook = NotebookModel(self.case)


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())


    def getMomentumFormulaComponents(self, zone):

        exp = self.getMomentumFormula(zone)
        if not exp:
            exp = """Su = 0;\nSv = 0;\nSw = 0;\n
dSudu = 0;\ndSudv = 0;\ndSudw = 0;\n
dSvdu = 0;\ndSvdv = 0;\ndSvdw = 0;\n
dSwdu = 0;\ndSwdv = 0;\ndSwdw = 0;\n"""

        req = [('Su', "x component of the momentum source term"),
               ('Sv', "y component of the momentum source term"),
               ('Sw', "z component of the momentum source term"),
               ('dSudu', "x component x velocity derivative"),
               ('dSudv', "x component y velocity derivative"),
               ('dSudw', "x component z velocity derivative"),
               ('dSvdu', "y component x velocity derivative"),
               ('dSvdv', "y component y velocity derivative"),
               ('dSvdw', "y component z velocity derivative"),
               ('dSwdu', "z component x velocity derivative"),
               ('dSwdv', "z component y velocity derivative"),
               ('dSwdw', "z component z velocity derivative")]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume'),
               ('fluid_volume', 'Source terms zone fluid volume')]

        sym.append( ("u", 'x velocity component'))
        sym.append( ("v", 'y velocity component'))
        sym.append( ("w", 'z velocity component'))
        sym.append( ("rho", 'local density (kg/m^3)'))

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setMomentumFormula(self, zone, formula):
        """
        Public method.
        Set the formula for the velocity.
        """
        self.__verifyZone(zone)
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        n = node.xmlInitChildNode('momentum_formula', zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getMomentumFormula(self, zone):
        """
        Public method.
        Return the formula for the velocity.
        """
        self.__verifyZone(zone)
        node = self.node_sterm

        formula = node.xmlGetString('momentum_formula', zone_id=zone)

        return formula


    def getSpeciesFormulaComponents(self, zone, species):

        exp = self.getSpeciesFormula(zone, species)
        if not exp:
            exp = """S = 0;\ndS = 0;\n"""
        req = [('S', 'Explicit species source term ([species]*kg/m^3/s)'),
                ('dS', 'Species source term derivative (kg/m^3/s)')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume'),
               ('fluid_volume', 'Source terms zone fluid volume')]

        sym.append( ("rho", 'local density (kg/m^3)'))

        name = self.th_sca.getScalarName(species)
        sym.append((name, 'current species'))

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        # Known fields
        knf = [(str(name), str(name)), ('rho', 'density')]

        return exp, req, sym, knf


    @Variables.undoGlobal
    def setSpeciesFormula(self, zone, species, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        name_species = DefineUserScalarsModel(self.case).getScalarName(species)
        n = node.xmlInitChildNode('scalar_formula', name = name_species, label = species, zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getSpeciesFormula(self, zone, species):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        name_species = DefineUserScalarsModel(self.case).getScalarName(species)
        formula = node.xmlGetString('scalar_formula', name = name_species, label = species, zone_id=zone)

        return formula


    def getGroundWaterSpeciesFormulaComponents(self, zone, species):

        exp = self.getGroundWaterSpeciesFormula(zone, species)
        if not exp:
            exp = """Q = 0;\nlambda = 0;"""

        req = [('Q',      'species source term'),
               ('lambda', 'radioactive decay')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume'),
               ('fluid_volume', 'Source terms zone fluid volume')]

        name = self.th_sca.getScalarName(species)
        sym.append((name, 'current species'))

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        # Known fields
        knf = [(str(name), str(name))]

        return exp, req, sym, knf

    @Variables.undoGlobal
    def setGroundWaterSpeciesFormula(self, zone, species, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        name_species = DefineUserScalarsModel(self.case).getScalarName(species)
        n = node.xmlInitChildNode('scalar_formula', name = name_species, label = species, zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getGroundWaterSpeciesFormula(self, zone, species):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        name_species = DefineUserScalarsModel(self.case).getScalarName(species)
        formula = node.xmlGetString('scalar_formula', name = name_species, label = species, zone_id=zone)

        return formula


    def getRichardsFormulaComponents(self, zone):
        exp = self.getRichardsFormula(zone)
        if not exp:
            exp = """Qs = 0;\n"""

        req = [('Qs', 'volumetric source term')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume'),
               ('fluid_volume', 'Source terms zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoGlobal
    def setRichardsFormula(self, zone, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        n = node.xmlInitChildNode('volumetric_source_term', zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getRichardsFormula(self, zone):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        formula = node.xmlGetString('volumetric_source_term', zone_id=zone)

        return formula


    def getThermalFormulaComponents(self, zone, scalar):

        exp = self.getThermalFormula(zone, scalar)
        if not exp:
            exp = self.getDefaultThermalFormula(scalar)

        req = [('S', 'Explicit thermal source term (W/m^3)'),
               ('dS', 'Thermal source term derivative (W/m^3/[thermal scalar])')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume'),
               ('fluid_volume', 'Source terms zone fluid volume'),
               ('rho', 'density (kg/m^3)')]

        name = 'temperature'
        if self.case.module_name() == 'code_saturne':
            if self.therm.getThermalScalarModel() == 'enthalpy':
                name = 'enthalpy'
                sym.append(('enthalpy', 'thermal scalar'))
            if self.therm.getThermalScalarModel() == 'total_energy':
                name = 'total_energy'
                sym.append(('total_energy', 'thermal scalar'))
            else:
                sym.append(('temperature', 'thermal scalar'))

        elif self.case.module_name() == 'neptune_cfd':
            name = 'enthalpy'
            sym.append(('enthalpy', 'Enthalpy'))

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        # Known fields
        knf = [(str(name), str(name)), ('rho', 'density')]

        return exp, req, sym, knf


    @Variables.undoGlobal
    def setThermalFormula(self, zone, scalar, formula):
        """
        Public method.
        Set the formula for tharmal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, ['enthalpy', 'total_energy', 'temperature'])
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        n = node.xmlInitChildNode('thermal_formula', name = scalar, zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getThermalFormula(self, zone, scalar):
        """
        Public method.
        Return the formula for thermal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, ['enthalpy', 'total_energy', 'temperature'])
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)

        formula = node.xmlGetString('thermal_formula', name = scalar, zone_id=zone)

        return formula


    def getDefaultThermalFormula(self, scalar):
        """
        Public method.
        Return the default formula for thermal scalars.
        """
        formula = """S = 0;\ndS = 0;"""

        return formula


#-------------------------------------------------------------------------------
# InitializationModel test case
#-------------------------------------------------------------------------------


class SourceTermsTestCase(ModelTest):
    """
    Unittest.
    """
    def checkSourceTermsInstantiation(self):
        """Check whether the SourceTerms class could be instantiated."""
        model = None
        model = SourceTermsModel(self.case)
        assert model != None, 'Could not instantiate SourceTerms'


def suite():
    testSuite = unittest.makeSuite(SourceTermsTestCase, "check")
    return testSuite


def runTest():
    print("SourceTermsTestCase - OK !!!!")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
