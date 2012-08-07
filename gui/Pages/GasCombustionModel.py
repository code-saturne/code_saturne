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
from Pages.NumericalParamEquationModel import NumericalParamEquatModel
from Pages.LocalizationModel import LocalizationModel
from Pages.Boundary import Boundary

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

        nModels          = self.case.xmlGetNode('thermophysical_models')
        self.node_turb   = nModels.xmlInitNode('turbulence',        'model')
        self.node_gas    = nModels.xmlInitNode('gas_combustion',    'model')
        self.node_coal   = nModels.xmlInitNode('solid_fuels',       'model')
        self.node_joule  = nModels.xmlInitNode('joule_effect',      'model')
        self.node_therm  = nModels.xmlInitNode('thermal_scalar',    'model')
        self.node_atmo   = nModels.xmlInitNode('atmospheric_flows', 'model')
        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_reference = self.node_models.xmlInitNode('reference_values')

        self.gasCombustionModel = ('off', 'ebu', 'd3p','lwp')
        self.d3p_list = ("adiabatic", "extended")
        self.ebu_list = ("spalding", "enthalpy_st", "mixture_st", "enthalpy_mixture_st")
        self.lwp_list = ("2-peak_adiabatic", "2-peak_enthalpy",
                         "3-peak_adiabatic", "3-peak_enthalpy",
                         "4-peak_adiabatic", "4-peak_enthalpy")


    def defaultGasCombustionValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['model'] = "off"

        model = self.getGasCombustionModel()
        if model == 'd3p':
            default['option'] = "adiabatic"
        elif model == 'ebu':
            default['option'] = "spalding"
        elif model == 'lwp':
            default['option'] = "2-peak_adiabatic"
        elif model == 'off':
            default['option'] = "off"

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
                                           'Rij-EBRSM',
                                           'v2f-phi',
                                           'k-omega-SST',
                                           'Spalart-Allmaras'):
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
            for tag in ('scalar',
                        'property',
                        'reference_mass_molar',
                        'reference_temperature'):
                for node in self.node_gas.xmlGetNodeList(tag):
                    node.xmlRemoveNode()
            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "inlet":
                    Boundary("gas_comb_inlet", zone.getLabel(), self.case).deleteGas()

        else:
            self.node_gas['model'] = model
            self.node_coal['model']  = 'off'
            self.node_joule['model'] = 'off'
            self.node_therm['model'] = 'off'

        if model != 'd3p':
            self.node_reference.xmlRemoveChild('oxydant_temperature')
            self.node_reference.xmlRemoveChild('fuel_temperature')

        self.createModel()


    def getGasCombustionModel(self):
        """
        Return the current gas combustion model.
        """
        model = self.node_gas['model']
        if model not in self.gasCombustionModelsList():
            model = self.defaultGasCombustionValues()['model']
            self.setGasCombustionModel(model)

        return model


    def getGasCombustionOption(self):
        """
        Return the current gas combustion option.
        """
        option = self.node_gas['option']
        if option == None:
            option = self.defaultGasCombustionValues()['option']
            self.setGasCombustionOption(option)

        model = self.getGasCombustionModel()
        if model == 'd3p':
            if option not in self.d3p_list:
                option = self.defaultGasCombustionValues()['option']
                self.setGasCombustionOption(option)
        elif model == 'ebu':
            if option not in self.ebu_list:
                option = self.defaultGasCombustionValues()['option']
                self.setGasCombustionOption(option)
        elif model == 'lwp':
            if option not in self.lwp_list:
                option = self.defaultGasCombustionValues()['option']
                self.setGasCombustionOption(option)
        elif model == 'off':
            option = 'off'
        return option


    def setGasCombustionOption(self, option):
        """
        Return the current gas combustion option.
        """
        model = self.getGasCombustionModel()
        if model == 'd3p':
            self.isInList(option, self.d3p_list)
        elif model == 'ebu':
            self.isInList(option, self.ebu_list)
        elif model == 'lwp':
            self.isInList(option, self.lwp_list)
        elif model == 'off':
            self.isInList(option, ('off'))
        self.node_gas['option'] = option
        option = self.node_gas['option']
        self.createModel()


    def __createModelScalarsList(self , model):
        """
        Private method
        Create model scalar list
        """
        option = self.getGasCombustionOption()
        list_options = ["3-peak_adiabatic", "3-peak_enthalpy",
                        "4-peak_adiabatic", "4-peak_enthalpy"]
        acceptable_options = ["2-peak_enthalpy", "3-peak_enthalpy",
                              "4-peak_enthalpy"]
        lst = []

        if model == 'd3p':
            lst.append("Fra_MEL")
            lst.append("Var_FMe")
            if option == 'extended':
                lst.append("Enthalpy")
        elif model == 'ebu':
            lst.append("Fra_GF")
            if option == "mixture_st" or option =="enthalpy_misture_st":
                lst.append("Fra_MEL")
            elif option == "enthalpy_st" or option =="enthalpy_mixture_st":
                lst.append("Enthalpy")
        elif model == 'lwp':
            lst.append("Fra_MEL")
            lst.append("Var_FMe")
            lst.append("Fra_Mas")
            lst.append("COYF_PP4")
            if option in list_options:
                lst.append("Var_FMa")
            if option in acceptable_options:
                lst.append("Enthalpy")
        return lst


    def __createModelPropertiesList(self, model):
        """
        Private method
        Create model properties
        """
        lst = []
        lst.append("Temperature")
        lst.append("YM_Fuel")
        lst.append("YM_Oxyd")
        lst.append("YM_Prod")
        if model == 'lwp':
            lst.append("T.SOURCE")
            lst.append("Mas_Mol")
            ndirac = self.getNdirac()
            for idirac in range(ndirac):
                lst.append("RHOL0" + str(idirac + 1))
                lst.append("TEML0" + str(idirac + 1))
                lst.append("FMEL0" + str(idirac + 1))
                lst.append("FMAL0" + str(idirac + 1))
                lst.append("AMPL0" + str(idirac + 1))
                lst.append("TSCL0" + str(idirac + 1))
                lst.append("MAML0" + str(idirac + 1))
        return lst


    def __createModelScalars(self , model):
        """
        Private method
        Create model scalar
        """
        previous_list = []
        nodes = self.node_gas.xmlGetChildNodeList('scalar')
        for node in nodes:
            previous_list.append(node['name'])

        if model == "off":
            for node in nodes:
                node.xmlRemoveNode()
        else:
            new_list = self.__createModelScalarsList(model)
            for name in previous_list:
                if name not in new_list:
                    self.node_gas.xmlRemoveChild('scalar',  name = name)

            for name in new_list:
                if name not in previous_list:
                    self.setNewModelScalar(self.node_gas, name)

            NPE = NumericalParamEquatModel(self.case)
            for node in self.node_gas.xmlGetChildNodeList('scalar'):
                NPE.setBlendingFactor(node['label'], 0.)
                NPE.setScheme(node['label'], 'upwind')
                NPE.setFluxReconstruction(node['label'], 'off')


    def getNdirac(self):
        """
        """
        option = self.getGasCombustionOption()
        self.isInList(option, self.lwp_list)
        if option == '2-peak_adiabatic' or option == '2-peak_enthalpy':
            ndirac = 2
        if option == '3-peak_adiabatic' or option == '3-peak_enthalpy':
            ndirac = 3
        if option == '4-peak_adiabatic' or option == '4-peak_enthalpy':
            ndirac = 4
        return ndirac


    def __createModelProperties(self, model):
        """
        Private method
        Create model properties
        """
        previous_list = []
        nodes = self.node_gas.xmlGetChildNodeList('property')
        if model == "off":
            for node in nodes:
                node.xmlRemoveNode()
        else:
            for node in nodes:
                previous_list.append(node['name'])

            new_list = self.__createModelPropertiesList(model)
            for name in previous_list:
                if name not in new_list:
                    self.node_gas.xmlRemoveChild('property',  name = name)

            for name in new_list:
                if name not in previous_list:
                    self.setNewProperty(self.node_gas, name)


    def createModel (self) :
        """
        Private method
        Create scalars and properties when gas combustion is selected
        """
        model = self.getGasCombustionModel()
        self.__createModelScalars(model)
        self.__createModelProperties(model)


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
        XMLinit(self.case).initialize()
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
