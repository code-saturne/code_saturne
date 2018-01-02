# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

import sys, string, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.LocalizationModel import LocalizationModel

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


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())


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
