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
This module defines the electrical thermal flow modelling management.

This module contains the following classes and function:
- ElectricalModel
- ElectricalTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Variables, Model
from Pages.ThermalRadiationModel import ThermalRadiationModel
from Pages.FluidCharacteristicsModel import FluidCharacteristicsModel

#-------------------------------------------------------------------------------
# Coal combustion model class
#-------------------------------------------------------------------------------

class ElectricalModel(Variables, Model):
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        nModels         = self.case.xmlGetNode('thermophysical_models')
        self.node_turb  = nModels.xmlGetNode('turbulence',      'model')
        self.node_gas   = nModels.xmlInitNode('gas_combustion',  'model')
        self.node_coal  = nModels.xmlInitNode('pulverized_coal', 'model')
        self.node_joule = nModels.xmlInitNode('joule_effect',    'model')
        self.node_therm = nModels.xmlGetNode('thermal_scalar',  'model')

        self.electricalModel = ('off', 'joule', 'arc')


    def defaultElectricalValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['model'] = "off"

        return default


    def getAllElectricalModels(self):
        """
        Return all defined electrical models in a tuple.
        """
        return self.electricalModel


    def electricalModelsList(self):
        """
        Create a tuple with the electrical models allowed
        by the calculation features.
        """
        electricalList = self.electricalModel

        n, m = FluidCharacteristicsModel(self.case).getThermalModel()
        if m != "off" and m not in electricalList:
            electricalList = ('off',)

        if self.node_turb['model'] not in ('off','k-epsilon',
                                           'k-epsilon-PL',
                                           'Rij-epsilon',
                                           'Rij-SSG',
                                           'v2f-phi',
                                           'k-omega-SST'):
            electricalList = ('off',)

        return electricalList


    def setElectricalModel(self, model):
        """
        Update the electrical model markup from the XML document.
        """
        self.isInList(model, self.electricalModelsList())

        if model == 'off':
            self.node_joule['model']   = 'off'
            ThermalRadiationModel(self.case).setRadiativeModel('off')
        else:
            self.node_gas['model']   = 'off'
            self.node_coal['model']  = 'off'
            self.node_joule['model'] = model
            self.node_therm['model'] = 'off'


    def getElectricalModel(self):
        """
        Return the current electrical model.
        """
        model = self.node_joule['model']
        if model not in self.electricalModelsList():
            model = self.defaultElectricalValues()['model']
            self.setElectricalModel(model)

        return model


#-------------------------------------------------------------------------------
# Electrical model test case
#-------------------------------------------------------------------------------


class ElectricalTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from Base.XMLengine import Case, XMLDocument
        from Base.XMLinitialize import XMLinit
        Tool.GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case).initialize()
        self.doc = XMLDocument()

    def tearDown(self):
        """This method is executed after all "check" methods."""
        del self.case
        del self.doc

    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()

    def checkElectricalInstantiation(self):
        """
        Check whether the ElectricalModel class could be instantiated
        """
        model = None
        model = ElectricalModel(self.case)
        assert model != None, 'Could not instantiate ElectricalModel'


def suite():
    testSuite = unittest.makeSuite(ElectricalTestCase, "check")
    return testSuite


def runTest():
    print("ElectricalTestCase - TODO**************")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
