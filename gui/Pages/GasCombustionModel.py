# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module defines the gas combustion thermal flow modelling management.

This module contains the following classes and function:
- GasCombustionModel
- GasCombustionTestCase
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
# Gas combustion model class
#-------------------------------------------------------------------------------

class GasCombustionModel(Variables, Model):
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        nModels         = self.case.xmlGetNode('thermophysical_models')
        self.node_turb  = nModels.xmlGetNode('turbulence',      'model')
        self.node_gas   = nModels.xmlGetNode('gas_combustion',  'model')
        self.node_coal  = nModels.xmlGetNode('pulverized_coal', 'model')
        self.node_joule = nModels.xmlGetNode('joule_effect',    'model')
        self.node_therm = nModels.xmlGetNode('thermal_scalar',  'model')

        self.gasCombustionModel = ('off', 'ebu', 'd3p')


    def defaultGasCombustionValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['model'] = "off"

        return default


    def getAllGasCombustionModels(self):
        """
        Return all defined gas combustion models in a tuple.
        """
        return self.gasCombustionModel


    def gasCombustionModelsList(self):
        """
        Create a tuple with the gas combustion models allowed
        by the calculation features.
        """
        gasCombustionList = self.gasCombustionModel

        n, m = FluidCharacteristicsModel(self.case).getThermalModel()
        if m != "off" and m not in gasCombustionList:
            gasCombustionList = ('off',)

        if self.node_turb['model'] not in ('k-epsilon',
                                           'k-epsilon-PL',
                                           'Rij-epsilon',
                                           'Rij-SSG',
                                           'v2f-phi',
                                           'k-omega-SST'):
            gasCombustionList = ('off',)

        return gasCombustionList


    def setGasCombustionModel(self, model):
        """
        Update the gas combustion model markup from the XML document.
        """
        self.isInList(model, self.gasCombustionModelsList())

        if model == 'off':
            self.node_gas['model'] = model
            ThermalRadiationModel(self.case).setRadiativeModel('off')
        else:
            self.node_gas['model'] = model
            self.node_coal['model']  = 'off'
            self.node_joule['model'] = 'off'
            self.node_therm['model'] = 'off'


    def getGasCombustionModel(self):
        """
        Return the current gas combustion model.
        """
        model = self.node_gas['model']
        if model not in self.gasCombustionModelsList():
            model = self.defaultGasCombustionValues()['model']
            self.setGasCombustionModel(model)

        return model


#-------------------------------------------------------------------------------
# Gas combustion test case
#-------------------------------------------------------------------------------


class GasCombustionTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from Base.XMLengine import Case, XMLDocument
        from Base.XMLinitialize import XMLinit
        Tool.GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case)
        self.doc = XMLDocument()

    def tearDown(self):
        """This method is executed after all "check" methods."""
        del self.case
        del self.doc

    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()

    def checkGasCombustionInstantiation(self):
        """
        Check whether the gasCombustionModel class could be instantiated
        """
        model = None
        model = GasCombustionModel(self.case)
        assert model != None, 'Could not instantiate GasCombustionModel'


def suite():
    testSuite = unittest.makeSuite(GasCombustionTestCase, "check")
    return testSuite


def runTest():
    print("GasCombustionTestCase - TODO**************")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
