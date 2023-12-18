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
This module defines the gas combustion thermal flow modelling management.

This module contains the following classes and function:
- GasCombustionModel
- GasCombustionTestCase
- ThermochemistryData
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.TurbulenceModel import TurbulenceModel
from code_saturne.model.ThermalRadiationModel import ThermalRadiationModel
from code_saturne.model.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.model.NumericalParamEquationModel import NumericalParamEquationModel
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.Boundary import Boundary

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
        self.node_prop   = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = self.node_prop.xmlInitNode('fluid_properties')

        self.gasCombustionModel = ('off', 'ebu', 'd3p','lwp')
        self.d3p_list = ("adiabatic", "extended")
        self.ebu_list = ("spalding", "enthalpy_st", "mixture_st", "enthalpy_mixture_st")
        self.lwp_list = ("2-peak_adiabatic", "2-peak_enthalpy",
                         "3-peak_adiabatic", "3-peak_enthalpy",
                         "4-peak_adiabatic", "4-peak_enthalpy")

        self.sootModels = ('off', 'constant_soot_yield', 'moss')


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
        Create a tuple with the available gas combustion models.
        """
        return self.gasCombustionModel


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
                        'reference_temperature',
                        'soot_model'):
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
            self.node_fluid.xmlRemoveChild('reference_oxydant_temperature')
            self.node_fluid.xmlRemoveChild('reference_fuel_temperature')

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
        if option is None:
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
            if option == "mixture_st" or option =="enthalpy_mixture_st":
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
                    self.setNewVariable(self.node_gas, name, tpe='var_model', label=name)

            NPE = NumericalParamEquationModel(self.case)
            for node in self.node_gas.xmlGetChildNodeList('variable'):
                name = node['name']
                NPE.setBlendingFactor(name, 0.)
                NPE.setScheme(name, 'centered')
                NPE.setFluxReconstruction(name, 'off')

                if self.getGasCombustionModel() == "d3p":
                    if name == "mixture_fraction":
                        NPE.setMinValue(name, 0.)
                        NPE.setMaxValue(name, 1.)
                    elif name == "mixture_fraction_variance":
                        NPE.setMinValue(name, 0.)
                        NPE.setMaxValue(name, 1.e+12)
                        node.xmlSetData('variance', "mixture_fraction")


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

    def _defaultValues(self):
        """
        default values
        """
        self.default = {}
        self.default['thermodynamical_pressure'] = 'off'
        self.default['soot_model']               = 'off'
        self.default['soot_density']             = 0.0
        self.default['soot_fraction']            = 0.0
        return self.default

    @Variables.noUndo
    def getUniformVariableThermodynamicalPressure(self):
        """
        Return status of uniform variable thermodynamical pressure
        """
        node = self.node_gas.xmlInitNode('thermodynamical_pressure', 'status')
        status = node['status']
        if not status:
            status = self._defaultValues()['thermodynamical_pressure']
            self.setUniformVariableThermodynamicalPressure(status)
        return status

    @Variables.undoLocal
    def setUniformVariableThermodynamicalPressure(self, status):
        """
        Put status of uniform variable thermodynamical pressure
        """
        self.isOnOff(status)
        node = self.node_gas.xmlInitNode('thermodynamical_pressure', 'status')
        node['status'] = status

    @Variables.noUndo
    def getSootModel(self):
        """
        Return value of attribute model
        """
        node = self.node_gas.xmlInitChildNode('soot_model', 'model')
        model = node['model']
        if model not in self.sootModels:
            model = self._defaultValues()['soot_model']
            self.setSootModel(model)
        return model

    @Variables.undoGlobal
    def setSootModel(self, model):
        """
        Put value of attribute model to soot model
        """
        self.isInList(model, self.sootModels)
        node  = self.node_gas.xmlInitChildNode('soot_model', 'model')
        node['model'] = model
        if model == 'moss':
            self.node_gas.xmlRemoveChild('soot_fraction')
        if model == 'off':
            self.node_gas.xmlRemoveChild('soot_density')
            self.node_gas.xmlRemoveChild('soot_fraction')

    @Variables.noUndo
    def getSootDensity(self):
        """
        Return value of soot density
        """
        val = self.node_gas.xmlGetDouble('soot_density')
        if val is None:
            val = self._defaultValues()['soot_density']
            self.setSootDensity(val)
        return val

    @Variables.undoGlobal

    @Variables.undoGlobal
    def setSootDensity(self, val):
        """
        Put value of soot density
        """
        self.isPositiveFloat(val)
        self.node_soot = self.node_gas.xmlGetNode('soot_model')
        self.node_soot.xmlSetData('soot_density', val)

    @Variables.noUndo
    def getSootFraction(self):
        """
        Return value of soot fraction
        """
        val = self.node_gas.xmlGetDouble('soot_fraction')
        if val is None:
            val = self._defaultValues()['soot_fraction']
            self.setSootFraction(val)
        return val

    @Variables.undoGlobal
    def setSootFraction(self, val):
        """
        Put value of soot fraction
        """
        self.isPositiveFloat(val)
        self.node_soot = self.node_gas.xmlGetNode('soot_model')
        self.node_soot.xmlSetData('soot_fraction', val)

#-------------------------------------------------------------------------------
# Gas combustion test case
#-------------------------------------------------------------------------------


class GasCombustionTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from code_saturne.model.XMLengine import Case, XMLDocument
        from code_saturne.model.XMLinitialize import XMLinit
        GuiParam.lang = 'en'
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
# Thermochemistry data class
#-------------------------------------------------------------------------------

class ThermochemistryData(Model):
    """
    Useful methods to create a Janaf File with the GUI
    """
    def __init__(self, case):

        """
        Constructor.
        """
        self.case = case

        nModels              = self.case.xmlGetNode('thermophysical_models')
        self.node_thermodata = nModels.xmlInitNode('Thermochemistry_data')

        # Parameters to create the Thermochemistry Data File
        # Available Chemical Elements (WARNING must be in upper case)
        self.ChemicalElem = ['C', 'H', 'O', 'N', 'CL']
        # Molar mass for the Chemical Elements
        # (WARNING must be in upper case)
        self.MolarMass = {}
        self.MolarMass['C']  = 0.012
        self.MolarMass['H']  = 0.001
        self.MolarMass['O']  = 0.016
        self.MolarMass['N']  = 0.014
        self.MolarMass['CL'] = 0.0354
        # Number of known global species, (ngazg in colecd.f90)
        # = 2 -> Automatic equilibrium of the reaction
        # = 3 -> Equilibrium defined by the user
        # NumberOfKnownGlobalSpecies is updated in WriteThermochemistryDataFile
        self.NumberOfKnownGlobalSpecies = 3
        # Always 1 reaction (see colecd.f90)
        self.NbReaction = 1
        self.Error_GUI = False


    def defaultParamforTabu(self):
        """
        Return in a dictionnary which contains default values
        for parameters necessarry to create the Janaf File
        """
        default = {}
        default['CreateThermoDataFile'] = 'off'
        default['user_mode_chem'] = 'auto'
        default['NbPointsTabu'] = 10
        default['MaximumTemp'] = 3000.0
        default['MinimumTemp'] = 273.0

        return default


    @Variables.noUndo
    def getCreateThermoDataFile(self):
        """
        Return status of CreateThermoDataFile
        """
        node = self.node_thermodata.xmlInitNode('CreateThermoDataFile', 'status')
        status = node['status']
        if not status:
            status = self.defaultParamforTabu()['CreateThermoDataFile']
            self.setCreateThermoDataFile(status)
        return status


    @Variables.undoLocal
    def setCreateThermoDataFile(self, status):
        """
        Put status of CreateThermoDataFile
        """
        self.isOnOff(status)
        node = self.node_thermodata.xmlInitNode('CreateThermoDataFile', 'status')
        node['status'] = status
        if status != 'on':
            self.node_thermodata.xmlRemoveChild('NbPointsTabu')
            self.node_thermodata.xmlRemoveChild('MaximumTemp')
            self.node_thermodata.xmlRemoveChild('MinimumTemp')
            self.node_thermodata.xmlRemoveChild('user_mode_chem')


    @Variables.noUndo
    def getUserModeForChemicalReaction(self):
        """
        Return the mode to build the Janaf File
        """
        node = self.node_thermodata.xmlInitChildNode('user_mode_chem', 'choice')
        model = node['choice']
        if not model:
            model = self.defaultParamforTabu()['user_mode_chem']
            self.setUserModeForChemicalReaction(model)
        return model


    @Variables.undoGlobal
    def setUserModeForChemicalReaction(self, model):
        """
        Put the mode to build the Janaf File
        """
        self.node_mode_chem = self.node_thermodata.xmlInitChildNode('user_mode_chem', 'choice')
        self.node_mode_chem['choice'] = model

        if model == 'user':
            self.node_mode_chem.xmlRemoveChild('Fuel')
            self.node_mode_chem.xmlRemoveChild('VolPropO2')
            self.node_mode_chem.xmlRemoveChild('VolPropN2')
            self.node_mode_chem.xmlRemoveChild('COyield')
            self.node_mode_chem.xmlRemoveChild('CSyield')
        if model == 'auto':
            nodes = self.node_mode_chem.xmlGetChildNodeList('variable')
            for node in nodes:
                node.xmlRemoveNode()


    @Variables.noUndo
    def getNbPointsTabu(self):
        """
        Return value of NbPointsTabu
        """
        NbPointsTabu = self.node_thermodata.xmlGetInt('NbPointsTabu')
        if NbPointsTabu is None:
            NbPointsTabu = self.defaultParamforTabu()['NbPointsTabu']
            self.setNbPointsTabu(NbPointsTabu)

        return NbPointsTabu


    @Variables.undoGlobal
    def setNbPointsTabu(self, value):
        """
        Put value of NbPointsTabu
        """
        self.isInt(value)
        self.node_thermodata.xmlSetData('NbPointsTabu', value)


    @Variables.noUndo
    def getMaximumTemp(self):
        """
        Return value of MaximumTemp
        """
        MaximumTemp = self.node_thermodata.xmlGetDouble('MaximumTemp')
        if MaximumTemp is None:
            MaximumTemp = self.defaultParamforTabu()['MaximumTemp']
            self.setMaximumTemp(MaximumTemp)

        return MaximumTemp


    @Variables.undoGlobal
    def setMaximumTemp(self, value):
        """
        Put value of MaximumTemp
        """
        self.isFloat(value)
        self.node_thermodata.xmlSetData('MaximumTemp', value)


    @Variables.noUndo
    def getMinimumTemp(self):
        """
        Return value of MinimumTemp
        """
        MinimumTemp = self.node_thermodata.xmlGetDouble('MinimumTemp')
        if MinimumTemp is None:
            MinimumTemp = self.defaultParamforTabu()['MinimumTemp']
            self.setMinimumTemp(MinimumTemp)

        return MinimumTemp


    @Variables.undoGlobal
    def setMinimumTemp(self, value):
        """
        Put value of MinimumTemp
        """
        self.isFloat(value)
        self.node_thermodata.xmlSetData('MinimumTemp', value)


    def defaultSpeciesProperties(self):
        """
        Return the default properties for a species (mode "auto" and "user")
        """
        default = {}
        default['species_label']       = "species"
        default['chemical_formula']    = "CHON"
        default['fuel_composition']    = 0.0
        default['oxi_composition']     = 0.0
        default['prod_composition']    = 0.0
        default['coeff_absorption']    = 0.35
        default['fuel']                = "C4H10"
        default['volPropO2']           = 0.21
        default['volPropN2']           = 0.79
        default['COyield']             = 0.0
        default['CSyield']             = 0.0

        return default


    @Variables.noUndo
    def getChemicalFormulaFuel(self):
        """
        Return the chemical formula for the Fuel
        """
        ChemicalFormula = self.node_thermodata.xmlGetString('Fuel')
        if ChemicalFormula == "":
            ChemicalFormula = self.defaultSpeciesProperties()['fuel']
            self.setChemicalFormulaFuel(ChemicalFormula)

        return ChemicalFormula


    @Variables.undoGlobal
    def setChemicalFormulaFuel(self, chemical_formula):
        """
        Put the chemical formula for the Fuel
        """
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        self.node_mode_chem.xmlSetData('Fuel', chemical_formula)


    @Variables.noUndo
    def getVolPropO2(self):
        """
        Return the value of VolPropO2
        """
        VolPropO2 = self.node_thermodata.xmlGetDouble('VolPropO2')
        if VolPropO2 is None:
            VolPropO2 = self.defaultSpeciesProperties()['volPropO2']
            self.setVolPropO2(VolPropO2)

        return VolPropO2


    @Variables.undoGlobal
    def setVolPropO2(self, value):
        """
        Put the value of VolPropO2
        """
        self.isFloat(value)
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        self.node_mode_chem.xmlSetData('VolPropO2', value)


    @Variables.noUndo
    def getVolPropN2(self):
        """
        Return the value of VolPropN2
        """
        VolPropN2 = self.node_thermodata.xmlGetDouble('VolPropN2')
        if VolPropN2 is None:
            VolPropN2 = self.defaultSpeciesProperties()['volPropN2']
            self.setVolPropN2(VolPropN2)

        return VolPropN2


    @Variables.undoGlobal
    def setVolPropN2(self, value):
        """
        Put the value of VolPropN2
        """
        self.isFloat(value)
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        self.node_mode_chem.xmlSetData('VolPropN2', value)


    @Variables.noUndo
    def getCOyield(self):
        """
        Return the value of COyield
        """
        COyield = self.node_thermodata.xmlGetDouble('COyield')
        if COyield is None:
            COyield = self.defaultSpeciesProperties()['COyield']
            self.setCOyield(COyield)

        return COyield


    @Variables.undoGlobal
    def setCOyield(self, value):
        """
        Put the value of COyield
        """
        self.isFloat(value)
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        self.node_mode_chem.xmlSetData('COyield', value)


    @Variables.noUndo
    def getCSyield(self):
        """
        Return the value of CSyield
        """
        CSyield = self.node_thermodata.xmlGetDouble('CSyield')
        if CSyield is None:
            CSyield = self.defaultSpeciesProperties()['CSyield']
            self.setCSyield(CSyield)

        return CSyield


    @Variables.undoGlobal
    def setCSyield(self, value):
        """
        Put the value of CSyield
        """
        self.isFloat(value)
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        self.node_mode_chem.xmlSetData('CSyield', value)


    def __defaultSpeciesLabel(self, species_name=None):
        """
        Private method.
        Return a default label for a new species.
        """
        __coef = {}
        for l in self.getSpeciesNamesList():
            __coef[l] = l
        length = len(__coef)
        Lspe = self.defaultSpeciesProperties()['species_label']

        # new species: default value

        if not species_name:
            if length != 0:
                i = 1
                while (Lspe + str(i)) in list(__coef.values()):
                    i = i + 1
                num = str(i)
            else:
                num = str(1)
            species_name = Lspe + num
        return species_name


    @Variables.undoGlobal
    def addSpecies(self, name=None):
        """
        Public method.
        Input a new species I{name}
        """

        c = self.__defaultSpeciesLabel(name)

        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')

        if c not in self.getSpeciesNamesList():
            self.node_mode_chem.xmlInitNode('variable', label=c, name=c)
            self.setSpeciesChemicalFormula(c, self.defaultSpeciesProperties()['chemical_formula'])
            self.setCompFuel(c, self.defaultSpeciesProperties()['fuel_composition'])
            self.setCompOxi(c, self.defaultSpeciesProperties()['oxi_composition'])
            self.setCompProd(c, self.defaultSpeciesProperties()['prod_composition'])
            self.setCoeffAbsorp(c, self.defaultSpeciesProperties()['coeff_absorption'])

        return c


    @Variables.undoGlobal
    def updateSpecies(self, data):
        """
        Public method.
        Update the XML from the tableView
        """
        c = data[0]

        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')

        self.node_mode_chem.xmlInitNode('variable', label=c, name=c)
        self.setSpeciesChemicalFormula(c, data[1])
        self.setCompFuel(c, data[2])
        self.setCompOxi(c, data[3])
        self.setCompProd(c, data[4])
        self.setCoeffAbsorp(c, data[5])


    @Variables.undoGlobal
    def deleteSpecies(self, species_label):
        """
        Public method.
        Delete species I{name}
        """
        self.isInList(species_label, self.getSpeciesNamesList())
        node = self.node_thermodata.xmlGetNode('variable', label=species_label)
        node.xmlRemoveNode()


    @Variables.noUndo
    def getSpeciesNamesList(self):
        """
        Public method.
        Return the Species name list
        """
        self.node_mode_chem = self.node_thermodata.xmlGetNode('user_mode_chem')
        lst = []
        for node in self.node_mode_chem.xmlGetChildNodeList('variable'):
                lst.append(node['label'])
        return lst


    @Variables.noUndo
    def getSpeciesChemicalFormula(self, l):
        """
        Public method.
        Return the Chemical Formula
        """
        n = self.case.xmlGetNode('variable', label=l)
        name = n.xmlGetString('chemical_formula')
        return name


    @Variables.undoGlobal
    def setSpeciesChemicalFormula(self, species_label, species_name):
        """
        Put the Chemical Formula
        """
        n = self.case.xmlGetNode('variable', label=species_label)
        n.xmlSetData('chemical_formula', species_name)


    @Variables.undoGlobal
    def getCompFuel(self, l):
        """
        Get Fuel Composition
        """
        n = self.case.xmlGetNode('variable', label=l)
        val = n.xmlGetDouble('fuel_composition')
        if val is None:
            val = self.defaultSpeciesProperties()['fuel_composition']
            self.setCompFuel(l, val)

        return val


    @Variables.undoGlobal
    def setCompFuel(self, species_label, CompFuel):
        """
        Put Fuel Composition
        """
        n = self.case.xmlGetNode('variable', label=species_label)
        n.xmlSetData('fuel_composition', CompFuel)


    @Variables.undoGlobal
    def getCompOxi(self, l):
        """
        Get Oxi Composition
        """
        n = self.case.xmlGetNode('variable', label=l)
        val = n.xmlGetDouble('oxi_composition')
        if val is None:
            val = self.defaultSpeciesProperties()['oxi_composition']
            self.setCompOxi(l, val)

        return val


    @Variables.undoGlobal
    def setCompOxi(self, species_label, CompOxi):
        """
        Put Oxi Composition
        """
        n = self.case.xmlGetNode('variable', label=species_label)
        n.xmlSetData('oxi_composition', CompOxi)


    @Variables.undoGlobal
    def getCompProd(self, l):
        """
        Get Product Composition
        """
        n = self.case.xmlGetNode('variable', label=l)
        val = n.xmlGetDouble('prod_composition')
        if val is None:
            val = self.defaultSpeciesProperties()['prod_composition']
            self.setCompProd(l, val)

        return val


    @Variables.undoGlobal
    def setCompProd(self, species_label, CompProd):
        """
        Put Product Composition
        """
        n = self.case.xmlGetNode('variable', label=species_label)
        n.xmlSetData('prod_composition', CompProd)


    @Variables.undoGlobal
    def getCoeffAbsorp(self, l):
        """
        Get absorption coefficient
        """
        n = self.case.xmlGetNode('variable', label=l)
        val = n.xmlGetDouble('coeff_absorption')
        if val is None:
            val = self.defaultSpeciesProperties()['coeff_absorption']
            self.setCoeffAbsorp(l, val)

        return val


    @Variables.undoGlobal
    def setCoeffAbsorp(self, species_label, CoeffAbsorp):
        """
        Put absorption coefficient
        """
        n = self.case.xmlGetNode('variable', label=species_label)
        n.xmlSetData('coeff_absorption', CoeffAbsorp)


    @Variables.undoGlobal
    def getNumberOfChemicalElem(self, ChemicalFormula):
        """
        Get the number of Chemical Element
        """

        ChemicalFormula = str(ChemicalFormula.upper())

        NumberOfChemElem = {}
        for Elem in self.ChemicalElem:
            NumberOfChemElem[Elem] = "0"

        ChemicalElem_sorted = sorted(self.ChemicalElem, key=len, reverse=True)

        for Elem in ChemicalElem_sorted :
            ReadFirstElem = True
            if Elem in ChemicalFormula :
                NumberOfChemElem[Elem] = "1"
                ReadFirstElem = True
                position = ChemicalFormula.find(Elem)+len(Elem)
                ChemicalFormula = ChemicalFormula.replace(str(Elem), ' ')
                while ChemicalFormula[position:position+1].isdigit():
                    if ReadFirstElem:
                        NumberOfChemElem[Elem] = ChemicalFormula[position:position+1]
                    else :
                        NumberOfChemElem[Elem] = NumberOfChemElem[Elem]+ChemicalFormula[position:position+1]

                    temp = list(ChemicalFormula)
                    temp[position] = ' '
                    ChemicalFormula = "".join(temp)

                    ReadFirstElem = False
                    position=position+1

        return NumberOfChemElem

    @Variables.undoGlobal
    def WriteThermochemistryDataFile(self, file_path):
        """
        Write the thermochemistry Data File
        """
        # Assumption : only one reaction (in colecd : ir = 1) (parameter NbReaction)
        # in colecd : ngazg can equal to 2 or 3 (parameter NumberOfKnownGlobalSpecies)
        # if ngazg = 3 we add the parameters:
        # in colecd : igfuel and igoxy (parameters position_fuel and position_oxi) and in colecd : stoeg (parameter global_stoeg)
        # see colecd.f90 for more information

        #Comment to add at the end of each lines FR
        InfoLine = {}
        InfoLine['NumberOfSpecies'] = "Nb especes courantes"
        InfoLine['NbPointsTabu'] = "Nb de points de tabulation ENTH-TEMP"
        InfoLine['MinimumTemp'] = "TMIN"
        InfoLine['MaximumTemp'] = "TMAX"
        InfoLine['LineInfo-GaseousSpecies'] = "Especes Gazeuses"
        InfoLine['CoeffAbsorp'] = "Coeff absorption (ray)"
        InfoLine['GlobalElemCompo'] = "Nb especes elementaires"
        InfoLine['LineInfo-ChemElem'] = "Composition CHON"
        InfoLine['NumberOfKnownGlobalSpecies'] = "Nb d'especes globales connues (ici : / Fuel / Oxydant / Produits)"
        InfoLine['FuelComposition'] = "Numero espece reactive / Composition Fuel     en especes elementaires"
        InfoLine['OxiComposition'] = "Numero espece reactive / Composition Oxydant  en especes elementaires"
        InfoLine['ProdComposition'] = "Numero espece reactive / Composition Produits en especes elementaires"
        InfoLine['NbReaction'] = "Nb de reactions mises en jeu pour les especes globales"
        InfoLine['Stoeg'] = "Stoechiometrie en nb de mole especes globales"
        #----------------------------------------
        ErrorLine = {}
        ErrorLine['Error'] = "Error in the GUI\n"
        ErrorLine['No_species'] = "There are no species defined\n"
        ErrorLine['FLAG_No_Oxi'] = "There is no Oxidiser specified\n"
        ErrorLine['FLAG_No_Fuel'] = "There is no Fuel specified\n"
        ErrorLine['FLAG_No_Prod'] = "There is no Product specified\n"
        ErrorLine['FLAG_No_compo'] = " has no Nb of moles specified\n"

        #Prepare the data depending from the user mode (automatic or defined by user)
        ListOfSpecies       = []
        CompFuelDict        = {}
        CompOxiDict         = {}
        CompProdDict        = {}
        CoeffAbsorpDict     = {}
        ChemicalFormulaDict = {}
        Option_UserMode = self.getUserModeForChemicalReaction()
        if Option_UserMode == 'user':
            self.NumberOfKnownGlobalSpecies = 3
            ListOfSpecies       = self.getSpeciesNamesList()
            for label in ListOfSpecies:
                CompFuelDict[label]        = self.getCompFuel(label)
                CompOxiDict[label]         = self.getCompOxi(label)
                CompProdDict[label]        = self.getCompProd(label)
                CoeffAbsorpDict[label]     = self.getCoeffAbsorp(label)
                ChemicalFormulaDict[label] = self.getSpeciesChemicalFormula(label)
        elif Option_UserMode == 'auto':
            self.NumberOfKnownGlobalSpecies = 2
            #Fuel
            ListOfSpecies.append('fuel')
            CompFuelDict['fuel']        = 1.0
            CompOxiDict['fuel']         = 0.0
            CompProdDict['fuel']        = 0.0
            CoeffAbsorpDict['fuel']     = self.defaultSpeciesProperties()['coeff_absorption']
            ChemicalFormulaDict['fuel'] = self.getChemicalFormulaFuel()
            #O2
            ListOfSpecies.append('O2')
            CompFuelDict['O2']          = 0.0
            CompOxiDict['O2']           = 1.0
            CompProdDict['O2']          = 0.0
            CoeffAbsorpDict['O2']       = self.defaultSpeciesProperties()['coeff_absorption']
            ChemicalFormulaDict['O2']   = 'O2'
            #N2
            if self.getVolPropN2() != 0.0:
                ListOfSpecies.append('N2')
                CompFuelDict['N2']          = 0.0
                CompOxiDict['N2']           = round(self.getVolPropN2()/self.getVolPropO2(), 3)
                CompProdDict['N2']          = 0.0
                CoeffAbsorpDict['N2']       = self.defaultSpeciesProperties()['coeff_absorption']
                ChemicalFormulaDict['N2']   = 'N2'
            #CO yield
            if self.getCOyield() != 0.0:
                ListOfSpecies.append('CO')
                CompFuelDict['CO']          = 0.0
                CompOxiDict['CO']           = 0.0
                CompProdDict['CO']          = self.getCOyield()
                CoeffAbsorpDict['CO']       = self.defaultSpeciesProperties()['coeff_absorption']
                ChemicalFormulaDict['CO']   = 'CO'
            #CS yield
            if self.getCSyield() != 0.0:
                ListOfSpecies.append('CS')
                CompFuelDict['CS']          = 0.0
                CompOxiDict['CS']           = 0.0
                CompProdDict['CS']          = self.getCSyield()
                CoeffAbsorpDict['CS']       = self.defaultSpeciesProperties()['coeff_absorption']
                ChemicalFormulaDict['CS']   = 'C(S)'
            #H2O
            ListOfSpecies.append('H2O')
            CompFuelDict['H2O']          = 0.0
            CompOxiDict['H2O']           = 0.0
            CompProdDict['H2O']          = 0.0
            CoeffAbsorpDict['H2O']       = self.defaultSpeciesProperties()['coeff_absorption']
            ChemicalFormulaDict['H2O']   = 'H2O'
            #CO2
            if 'C' in str(self.getChemicalFormulaFuel()).upper():
                ListOfSpecies.append('CO2')
                CompFuelDict['CO2']          = 0.0
                CompOxiDict['CO2']           = 0.0
                CompProdDict['CO2']          = 0.0
                CoeffAbsorpDict['CO2']       = self.defaultSpeciesProperties()['coeff_absorption']
                ChemicalFormulaDict['CO2']   = 'CO2'
            #HCl
            if 'CL' in str(self.getChemicalFormulaFuel()).upper():
                ListOfSpecies.append('HCl')
                CompFuelDict['HCl']          = 0.0
                CompOxiDict['HCl']           = 0.0
                CompProdDict['HCl']          = 0.0
                CoeffAbsorpDict['HCl']       = self.defaultSpeciesProperties()['coeff_absorption']
                ChemicalFormulaDict['HCl']   = 'HCl'


        #Open the Thermochemistry Data File
        f = open(file_path, "w")

        #Check if all the species have a given composition with the "user" mode
        self.Error_GUI = False
        if Option_UserMode == 'user':
            CheckComp = {}
            FLAG_No_Oxi = False
            FLAG_No_Fuel = False
            FLAG_No_Prod = False
            FLAG_No_compo = False

            if not ListOfSpecies :
                f.write(ErrorLine['Error'])
                f.write(ErrorLine['No_species'])
                f.close()
                self.Error_GUI = True
                return

            for label in ListOfSpecies:
                CheckComp[label] = CompFuelDict[label]+CompOxiDict[label]+CompProdDict[label]

            if max(CompOxiDict.values()) == 0 : FLAG_No_Oxi = True
            if max(CompFuelDict.values()) == 0 : FLAG_No_Fuel = True
            if max(CompProdDict.values()) == 0 : FLAG_No_Prod = True
            if 0 in CheckComp.values(): FLAG_No_compo = True

            if FLAG_No_Oxi | FLAG_No_Fuel | FLAG_No_Prod | FLAG_No_compo :
                f.write(ErrorLine['Error'])
                if FLAG_No_Oxi : f.write(ErrorLine['FLAG_No_Oxi'])
                if FLAG_No_Fuel : f.write(ErrorLine['FLAG_No_Fuel'])
                if FLAG_No_Prod : f.write(ErrorLine['FLAG_No_Prod'])
                if FLAG_No_compo :
                    for label in ListOfSpecies:
                        if CheckComp.get(label) == 0 :
                            f.write("The " + label + " " + ChemicalFormulaDict[label]+ ErrorLine['FLAG_No_compo'])
                f.close()
                self.Error_GUI = True
                return

        # Write the number of species, and parameters for the tabulation ENTH TEMP
        NumberOfSpecies = len(ListOfSpecies)
        NbPointsTabu = self.getNbPointsTabu()
        MinimumTemp  = self.getMinimumTemp()
        MaximumTemp  = self.getMaximumTemp()
        f.write('{:<10}'.format(str(NumberOfSpecies))+InfoLine['NumberOfSpecies']+"\n")
        f.write('{:<10}'.format(str(NbPointsTabu))+InfoLine['NbPointsTabu']+"\n")
        f.write('{:<10}'.format(str(MinimumTemp))+InfoLine['MinimumTemp']+"\n")
        f.write('{:<10}'.format(str(MaximumTemp))+InfoLine['MaximumTemp']+"\n")
        f.write(InfoLine['LineInfo-GaseousSpecies']+"\n")
        f.write(" "*10)

        # Reorder the species to have the Fuel in first, then Oxidiser and
        # finish by the Product and sorted the compositions (see colecd.f90)
        maxCompFuel = 0.0
        maxCompOxi  = 0.0
        maxCompProd = 0.0
        for label in ListOfSpecies:
            maxCompFuel = max(maxCompFuel, CompFuelDict[label])
            maxCompOxi = max(maxCompOxi, CompOxiDict[label])
            maxCompProd = max(maxCompProd, CompOxiDict[label])

        Order = []
        for i in range(NumberOfSpecies):
            Order.append(3)

        i = 0
        for label in ListOfSpecies:
            if CompFuelDict[label] > 0.0:
                Order[i] = 0 + CompFuelDict[label]/maxCompFuel
            if CompOxiDict[label] > 0.0:
                Order[i] = 1 + CompOxiDict[label]/maxCompOxi
            if CompProdDict[label] > 0.0:
                Order[i] = 3 + CompProdDict[label]/maxCompProd
            i = i + 1

        getSpeciesNamesList_Sorted = ListOfSpecies
        Order, getSpeciesNamesList_Sorted = zip(*sorted(zip(Order, getSpeciesNamesList_Sorted)))

        # Write the chemical formula of all the species
        for label in getSpeciesNamesList_Sorted:
            f.write('{:>15}'.format(str(ChemicalFormulaDict[label])))
        f.write("\n")
        f.write(" "*10)

        # Write the absoption coefficient for all the species
        for label in getSpeciesNamesList_Sorted:
            f.write('{:>15}'.format(str(CoeffAbsorpDict[label])))
        f.write(" "*5+InfoLine['CoeffAbsorp']+"\n")

        # List of all chemical formula
        LstChemicalFormula = []
        for label in getSpeciesNamesList_Sorted:
            LstChemicalFormula.append(str(ChemicalFormulaDict[label]))

        # Write the number of chemical element (CHON)
        GlobalElemCompo = {}
        ChemicalElem_sorted = sorted(self.ChemicalElem, key=len, reverse=True)
        for Elem in ChemicalElem_sorted:
            GlobalElemCompo[Elem] = False
        for Elem in ChemicalElem_sorted:
            if Elem in str(LstChemicalFormula).upper():
                GlobalElemCompo[Elem] = True
        f.write('{:<15}'.format(str(str(GlobalElemCompo).count("True")))+" "*15*NumberOfSpecies+InfoLine['GlobalElemCompo']+"\n")

        # Create a dictionnary with the composition in chemical elements for each species
        Composition = {}
        for label in getSpeciesNamesList_Sorted:
            Composition[label] = self.getNumberOfChemicalElem(ChemicalFormulaDict[label])

        # Write the number of chemical element (CHON) for each species
        for Elem in self.ChemicalElem:
            if GlobalElemCompo[Elem] == True:
                f.write('{:<10}'.format(str(self.MolarMass.get(Elem))))
                for label in getSpeciesNamesList_Sorted:
                    f.write('{:>15}'.format(str(Composition[label][Elem])))
                if Elem == self.ChemicalElem[0]:
                    f.write(" "*5+InfoLine['LineInfo-ChemElem'])
                f.write('\n')

        # Write the number of known global species (always 2 (Fuel and Oxy) for the moment)
        f.write('{:<15}'.format(str(self.NumberOfKnownGlobalSpecies))+" "*15*NumberOfSpecies+InfoLine['NumberOfKnownGlobalSpecies']+"\n")
        f.write(" "*10)

        #Write the Fuel composition
        for label in getSpeciesNamesList_Sorted:
            f.write('{:>15}'.format(str(CompFuelDict[label])))
        f.write(" "*5+InfoLine['FuelComposition']+"\n")
        f.write(" "*10)

        #Write the Oxidiser composition
        for label in getSpeciesNamesList_Sorted:
            f.write('{:>15}'.format(str(CompOxiDict[label])))
        f.write(" "*5+InfoLine['OxiComposition']+"\n")
        f.write(" "*10)

        #Write the Product composition
        for label in getSpeciesNamesList_Sorted:
            f.write('{:>15}'.format(str(CompProdDict[label])))
        f.write(" "*5+InfoLine['ProdComposition'])

        #If all the species have a given composition : ngazg == 3
        if self.NumberOfKnownGlobalSpecies == 3 :
            position_fuel = 1
            position_oxi = 2
            f.write("\n")
            f.write('{:<15}'.format(str(self.NbReaction))+" "*15*NumberOfSpecies+InfoLine['NbReaction']+"\n")
            f.write('{:<15}'.format(str(position_fuel)))
            f.write('{:<15}'.format(str(position_oxi)))

            Prodsum = {}
            Prodsum['Fuel'] = 0.0
            Prodsum['Oxi'] = 0.0
            Prodsum['Prod'] = 0.0
            global_stoeg = {}
            global_stoeg['Fuel'] = 0.0
            global_stoeg['Oxi'] = 0.0
            global_stoeg['Prod'] = 0.0

            for label in getSpeciesNamesList_Sorted:
                Prodsum['Fuel'] = Prodsum['Fuel'] + float(Composition[label]['O'])*CompFuelDict[label]
                Prodsum['Oxi'] = Prodsum['Oxi'] + float(Composition[label]['O'])*CompOxiDict[label]
                Prodsum['Prod'] = Prodsum['Prod'] + float(Composition[label]['O'])*CompProdDict[label]
                global_stoeg['Fuel'] = global_stoeg['Fuel'] + CompFuelDict[label]
                global_stoeg['Oxi'] = global_stoeg['Oxi'] + CompOxiDict[label]
                global_stoeg['Prod'] = global_stoeg['Prod'] + CompProdDict[label]

            nreact_oxi = (Prodsum['Prod']-Prodsum['Fuel'])/Prodsum['Oxi']
            global_stoeg['Oxi'] = -global_stoeg['Oxi']*nreact_oxi
            f.write('{:<15}'.format(str(-global_stoeg['Fuel'])))
            f.write('{:<15}'.format(str(global_stoeg['Oxi'])))
            f.write('{:<15}'.format(str(global_stoeg['Prod'])))
            f.write(" "*15*(max(1,NumberOfSpecies-4))+InfoLine['Stoeg']+"\n")

        #Close the Thermochemistry Data File
        f.close()

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
