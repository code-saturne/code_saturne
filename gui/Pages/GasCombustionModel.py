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

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.Pages.NumericalParamEquationModel import NumericalParamEquationModel
from code_saturne.Pages.LocalizationModel import LocalizationModel
from code_saturne.Pages.Boundary import Boundary

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


    @Variables.noUndo
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

        if self.node_turb['model'] not in ('k-epsilon',
                                           'k-epsilon-PL',
                                           'Rij-epsilon',
                                           'Rij-SSG',
                                           'Rij-EBRSM',
                                           'v2f-BL-v2/k',
                                           'k-omega-SST',
                                           'Spalart-Allmaras'):
            gasCombustionList = ('off',)

        return gasCombustionList


    @Variables.undoGlobal
    def setGasCombustionModel(self, model):
        """
        Update the gas combustion model markup from the XML document.
        """
        self.isInList(model, self.gasCombustionModelsList())
        node_prop   = self.case.xmlGetNode('physical_properties')
        node_fluid  = node_prop.xmlInitNode('fluid_properties')

        old_model = self.node_gas['model']
        ThermalScalarModel(self.case).setThermalModel('off')

        if model == 'off':
            self.node_gas['model'] = model
            self.node_gas['option'] = "off"
            ThermalRadiationModel(self.case).setRadiativeModel('off')
            for tag in ('variable',
                        'property',
                        'reference_mass_molar',
                        'reference_temperature'):
                for node in self.node_gas.xmlGetNodeList(tag):
                    node.xmlRemoveNode()
            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "inlet":
                    Boundary("inlet", zone.getLabel(), self.case).deleteGas()

            node_fluid.xmlRemoveChild('property', name='dynamic_diffusion')

        else:
            self.node_gas['model'] = model
            self.node_coal['model']  = 'off'
            self.node_joule['model'] = 'off'
            self.setNewFluidProperty(node_fluid, 'dynamic_diffusion')

            if old_model != model:
                for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                    if zone.getNature() == "inlet":
                        Boundary("inlet", zone.getLabel(), self.case).deleteGas()

        if model != 'd3p':
            self.node_reference.xmlRemoveChild('oxydant_temperature')
            self.node_reference.xmlRemoveChild('fuel_temperature')

        self.createModel()


    @Variables.noUndo
    def getGasCombustionModel(self):
        """
        Return the current gas combustion model.
        """
        model = self.node_gas['model']
        if model not in self.gasCombustionModelsList():
            model = 'off'
            self.setGasCombustionModel(model)

        return model


    @Variables.noUndo
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


    @Variables.undoGlobal
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

        ThermalScalarModel(self.case).setThermalModel('off')

        if model == 'd3p':
            lst.append("mixture_fraction")
            lst.append("mixture_fraction_variance")
            if option == 'extended':
                ThermalScalarModel(self.case).setThermalModel('enthalpy')
        elif model == 'ebu':
            lst.append("fresh_gas_fraction")
            if option == "mixture_st" or option =="enthalpy_misture_st":
                lst.append("mixture_fraction")
            elif option == "enthalpy_st" or option =="enthalpy_mixture_st":
                ThermalScalarModel(self.case).setThermalModel('enthalpy')
        elif model == 'lwp':
            lst.append("mixture_fraction")
            lst.append("mixture_fraction_variance")
            lst.append("mass_fraction")
            lst.append("mass_fraction_covariance")
            if option in list_options:
                lst.append("mass_fraction_variance")
            if option in acceptable_options:
                ThermalScalarModel(self.case).setThermalModel('enthalpy')
        return lst


    def __createModelPropertiesList(self, model):
        """
        Private method
        Create model properties
        """
        lst = []
        lst.append("temperature")
        lst.append("ym_fuel")
        lst.append("ym_oxyd")
        lst.append("ym_prod")
        if model == 'lwp':
            lst.append("source_term")
            lst.append("molar_mass")
            ndirac = self.getNdirac()
            for idirac in range(ndirac):
                lst.append("rho_local_" + str(idirac + 1))
                lst.append("temperature_local_" + str(idirac + 1))
                lst.append("ym_local_" + str(idirac + 1))
                lst.append("w_local_" + str(idirac + 1))
                lst.append("amplitude_local_" + str(idirac + 1))
                lst.append("chemical_st_local_" + str(idirac + 1))
                lst.append("molar_mass_local_" + str(idirac + 1))
        return lst


    def __createModelScalars(self , model):
        """
        Private method
        Create model scalar
        """
        previous_list = []
        nodes = self.node_gas.xmlGetChildNodeList('variable')
        for node in nodes:
            previous_list.append(node['name'])

        if model == "off":
            for node in nodes:
                node.xmlRemoveNode()
        else:
            new_list = self.__createModelScalarsList(model)
            for name in previous_list:
                if name not in new_list:
                    self.node_gas.xmlRemoveChild('variable',  name = name)

            for name in new_list:
                if name not in previous_list:
                    self.setNewVariable(self.node_gas, name, tpe="model", label=name)

            NPE = NumericalParamEquationModel(self.case)
            for node in self.node_gas.xmlGetChildNodeList('variable'):
                NPE.setBlendingFactor(node['label'], 0.)
                NPE.setScheme(node['label'], 'upwind')
                NPE.setFluxReconstruction(node['label'], 'off')


    @Variables.noUndo
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


    @Variables.noUndo
    def getThermoChemistryDataFileName(self):
        """
        Get name for properties data (return None if not defined)i
        """
        f = self.node_gas.xmlGetString('data_file')
        return f


    @Variables.undoLocal
    def setThermoChemistryDataFileName(self, name):
        """
        Put name for properties data and load file for number gaz and radiative model
        """
        self.node_gas.xmlSetData('data_file', name)


#-------------------------------------------------------------------------------
# Gas combustion test case
#-------------------------------------------------------------------------------


class GasCombustionTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from code_saturne.Base.XMLengine import Case, XMLDocument
        from code_saturne.Base.XMLinitialize import XMLinit
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
