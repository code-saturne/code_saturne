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
This module defines the coal combustion thermal flow modelling management.

This module contains the following classes and function:
- CoalCombustionModel
- CoalCombustionTestCase
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
import code_saturne.Pages.IdentityAndPathesModel as IdentityAndPathesModel
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
from code_saturne.Pages.LocalizationModel import LocalizationModel
from code_saturne.Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# Solid fuel combustion model class
#-------------------------------------------------------------------------------

class CoalCombustionModel(Variables, Model):
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_lagr   = self.case.xmlGetNode('lagrangian', 'model')
        self.node_turb   = self.node_models.xmlGetNode('turbulence',   'model')
        self.node_fuel   = self.node_models.xmlInitNode('solid_fuels', 'model')

        self.coalCombustionModel = ('off', 'homogeneous_fuel', 'homogeneous_fuel_moisture',
                                    'homogeneous_fuel_moisture_lagr')


    def defaultValues(self):
        """
        Private method
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['model']                                 = "off"
        default['diameter']                              = 0.000122
        default['mass_percent']                          = 0.1
        default['C_compo_dry']                           = 70.9
        default['H_compo_dry']                           = 4.6
        default['O_compo_dry']                           = 10.8
        default['N_compo_dry']                           = 0
        default['S_compo_dry']                           = 0
        default['C_coke_compo_dry']                      = 100
        default['H_coke_compo_dry']                      = 0.
        default['O_coke_compo_dry']                      = 0.
        default['N_coke_compo_dry']                      = 0.
        default['S_coke_compo_dry']                      = 0.
        default['O2_oxi']                                = 1
        default['N2_oxi']                                = 3.76
        default['H2O_oxi']                               = 0
        default['CO2_oxi']                               = 0
        default['PCIChoice']                             = 'LHV'
        default['PCIType']                               = 'dry_basis'
        default['diameter_type']                         = 'automatic'
        default['stoichiometric_coefficient']            = 'user_define'
        default['PCI']                                   = 0
        default['density']                               = 1200
        default['volatile_matter']                       = 0.
        default['ashes_enthalpy']                        = 0
        default['ashes_thermal_capacity']                = 1800
        default['rate_of_ashes_on_mass']                 = 11.5
        default['specific_heat_average']                 = 1800
        default['thermal_conductivity']                  = 1.e-5
        default['moisture']                              = 0.
        default['Y1']                                    = 0.416
        default['Y2']                                    = 0.582
        default['A1_pre-exponential_factor']             = 370000
        default['A2_pre-exponential_factor']             = 1.3e13
        default['E1_energy_of_activation']               = 74000
        default['E2_energy_of_activation']               = 250000
        default['HCN_NH3_partitionning_reaction_1']      = 0.5
        default['HCN_NH3_partitionning_reaction_2']      = 0.5
        default['pre-exp_constant']                      = 38.
        default['E_activation']                          = 15.96
        default['order_reaction']                        = "1"
        default['nitrogen_fraction']                     = 1
        default['nitrogen_concentration']                = 0.015
        default['nitrogen_in_char_at_low_temperatures']  = 0.8
        default['nitrogen_in_char_at_high_temperatures'] = 0.2
        default['percentage_HCN_char_combustion']        = 1.
        default['mean_devolatilisation_rate']            = 0
        default['refusal_value']                         = 0.001
        default['oxidant_type']                          = 'molar'
        default['fuel_type']                             = 'coal'
        default['NOx_formation']                         = 'on'
        default['CO2_Kinetics']                          = 'off'
        default['H2O_Kinetics']                          = 'off'
        default['improved_NOx_model']                    = 'off'
        default['reburning']                             = 'unused'

        return default


    def __coalCombustionModelsList(self):
        """
        Private method
        Create a tuple with the coal combustion models allowed
        by the calculation features.
        """
        coalCombustionList = self.coalCombustionModel

        if self.node_lagr and self.node_lagr['model'] != 'off':
            coalCombustionList = ('off', 'homogeneous_fuel_moisture_lagr')

        if self.node_turb != None:
            if self.node_turb['model'] not in ('k-epsilon',
                                               'k-epsilon-PL',
                                               'Rij-epsilon',
                                               'Rij-SSG',
                                               'Rij-EBRSM',
                                               'v2f-BL-v2/k',
                                               'k-omega-SST',
                                               'Spalart-Allmaras'):
                coalCombustionList = ('off',)
        else:
            coalCombustionList = ('off',)

        return coalCombustionList


    def __getVariableList(self):
        """
        Private method
        Create list of variables for a class
        """
        modelVariables =  ["n_p_", "x_p_coal_", "x_p_char_", "x_p_h_"]
        if self.getCoalCombustionModel() == 'homogeneous_fuel_moisture' or self.getCoalCombustionModel() == 'homogeneous_fuel_moisture_lagr':
            modelVariables.append("x_p_wt_")

        return modelVariables


    def __getPropertiesList(self):
        """
        Private method
        Create list of properties for a class
        """
        modelProperties = ["t_p_", "x_p_", "rho_p_", "diam_p_", "dissapear_rate_p_",
                           "m_transfer_v1_p_",  "m_transfer_v2_p_", "het_ts_o2_p_"]
        if self.getCoalCombustionModel() == 'homogeneous_fuel_moisture' or self.getCoalCombustionModel() == 'homogeneous_fuel_moisture_lagr':
            modelProperties.append("dry_ts_p")

        if self.getCO2KineticsStatus() == 'on':
            modelProperties.append("het_ts_co2_p")

        if self.getH2OKineticsStatus() == 'on':
            modelProperties.append("het_ts_h2o_p")

        return modelProperties


    def __createModelVariableList(self):
        """
        Private method
        Create list of scalar for all fuels and classes
        """
        coalsNumber = self.getCoalNumber()      # total number of solid fuel
        classesNumber = self.getClassesNumber() # total number of class (diameter)

        lst = []

        # list of coal variables
        baseNames = ["fr_mv1_", "fr_mv2_"]
        for baseName in baseNames:
            for coal in range(0, coalsNumber):
                name = '%s%2.2i' % (baseName, coal+1)
                lst.append(name)

        # list of class variables
        baseNames = self.__getVariableList()
        for baseName in baseNames:
            for classe in range(0, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                lst.append(name)

        lst.append("fr_het_o2")
        lst.append("x_c_h")

        if self.getCO2KineticsStatus() == "on":
            lst.append("fr_het_co2")

        if self.getH2OKineticsStatus() == "on":
            lst.append("fr_het_h2o")

        if self.getNOxFormationStatus() == "on":
            lst.append("x_c_hcn")
            lst.append("x_c_no")
            lst.append("x_c_h_ox")
            lst.append("x_c_nh3")

        if self.getCoalCombustionModel() == 'homogeneous_fuel_moisture' or self.getCoalCombustionModel() == 'homogeneous_fuel_moisture_lagr':
            lst.append("fr_h2o")

        if self.getOxidantNumber() >= 2:
            lst.append("fr_oxyd2")

        if self.getOxidantNumber() == 3:
            lst.append("fr_oxyd3")

        # ieqco2 fix to true
        lst.append("x_c_co2")

        lst.append("f1f2_variance")

        return lst


    def __createModelVariableMinMaxList(self):
        """
        Private method
        Create list of scalar for all fuels and classes for special min/max initialization
        """
        coalsNumber = self.getCoalNumber()      # total number of solid fuel
        classesNumber = self.getClassesNumber() # total number of class (diameter)

        lst = []

        # list of class variables
        baseNames =  ["x_p_h_"]
        for baseName in baseNames:
            for classe in range(0, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                lst.append(name)

        if self.getNOxFormationStatus() == "on":
            lst.append("x_c_h_ox")

        return lst


    def __createModelScalars(self):
        """
        Private method
        Create and update model scalar
        """
        previous_list = []
        nodes = self.node_fuel.xmlGetChildNodeList('variable')
        for node in nodes:
            previous_list.append(node['name'])

        new_list = self.__createModelVariableList()
        for name in previous_list:
            if name not in new_list:
                self.node_fuel.xmlRemoveChild('variable',  name = name)

        for name in new_list:
            if name not in previous_list:
                self.setNewVariable(self.node_fuel, name, tpe="model", label=name)


    def __createModelPropertiesList(self):
        """
        Private method
        Create list of properties for all fuels and classes
        """
        classesNumber = self.getClassesNumber()

        lst = ["t_gas", "rho_gas", "ym_chx1m", "ym_chx2m",
               "ym_co", "ym_o2", "ym_co2", "ym_h2o", "ym_n2",
               "ym_h2s", "ym_h2", "ym_hcn", "ym_nh3", "ym_so2",
               "xm"]

        baseNames = self.__getPropertiesList()

        for baseName in baseNames:
            for classe in range(0, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                lst.append(name)

        lst.append("intensity")
        lst.append("x_carbone")
        lst.append("x_oxygen")
        lst.append("x_hydrogen")

        if self.getNOxFormationStatus() == "on":
            lst.append("exp1")
            lst.append("exp2")
            lst.append("exp3")
            lst.append("exp4")
            lst.append("exp5")
            lst.append("f_hcn_dev")
            lst.append("f_hcn_het")
            lst.append("f_nh3_dev")
            lst.append("f_nh3_het")
            lst.append("f_no_hcn")
            lst.append("f_no_nh3")
            lst.append("f_no_het")
            lst.append("f_no_the")
            lst.append("c_no_hcn")
            lst.append("c_no_nh3")
            lst.append("f_hcn_rb")
            lst.append("c_no_rb")
            lst.append("exp_rb")

        return lst


    def __createModelProperties(self):
        """
        Private method
        Create and update model properties
        """
        previous_list = []
        nodes = self.node_fuel.xmlGetChildNodeList('property')
        for node in nodes:
            previous_list.append(node['name'])

        new_list = self.__createModelPropertiesList()
        for name in previous_list:
            if name not in new_list:
                self.node_fuel.xmlRemoveChild('property',  name = name)

        for name in new_list:
            if name not in previous_list:
                self.setNewProperty(self.node_fuel, name)


    def createModel(self) :
        """
        Private method
        Create scalars and properties when coal combustion is selected
        """
        if not self.node_fuel.xmlGetNode('CO2_kinetics'):
            self.node_fuel.xmlInitNode('CO2_kinetics',
                                       status = self.defaultValues()['CO2_Kinetics'])
        if not self.node_fuel.xmlGetNode('H2O_kinetics'):
            self.node_fuel.xmlInitNode('H2O_kinetics',
                                       status = self.defaultValues()['H2O_Kinetics'])
        if not self.node_fuel.xmlGetNode('NOx_formation'):
            self.node_fuel.xmlInitNode('NOx_formation',
                                       status = self.defaultValues()['NOx_formation'])
        self.__createModelScalars()
        self.__createModelProperties()


    def __updateWetScalarsAndProperty(self, model):
        """
        Private method
        Delete scalars x_p_wt_ and fr_h2o and property dry_ts_coal
        if model isn't 'homogeneous_fuel_moisture'
        """
        # TODO a supprimer doit etre appele si on change le modele uniquement
        if model != 'homogeneous_fuel_moisture' and self.getCoalCombustionModel() != 'homogeneous_fuel_moisture_lagr':
            nod = self.node_fuel.xmlGetNode('variable', type="model", name="fr_h2o")
            if nod:
                nod.xmlRemoveNode()
            for node in self.node_fuel.xmlGetNodeList('variable', type="model"):
                if node['name'][:6] == "x_p_wt_":
                    node.xmlRemoveNode()
            for node in self.node_fuel.xmlGetNodeList('property'):
                if node['name'][:6] == "dry_ts_coal":
                    node.xmlRemoveNode()


    def __updateSolidFuelCombustionDensity(self, model):
        """
        Private method
        Update the coal combustion model markup from the XML document.
        """
        self.isInList(model, self.__coalCombustionModelsList())

        mdl = FluidCharacteristicsModel(self.case)

        if model != 'off':
            mdl.setPropertyMode('density', 'variable')
            if mdl.getPropertyMode('density') == 'constant':
                mdl.setPropertyMode('density', 'variable')


    def __createCoalModelScalars(self, coalsNumber, coalClassesNumber, classesNumber):
        """
        Private method
        Create new scalars for one coal
        """
        # add new scalars
        baseNames = self.__getVariableList()

        for baseName in baseNames:
            for classe in range(classesNumber - coalClassesNumber, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                self.setNewVariable(self.node_fuel, name, tpe="model", label=name)

        baseNames = ["fr_mv1_", "fr_mv2_"]
        for baseName in baseNames:
            name = '%s%2.2i' % (baseName, coalsNumber)
            self.setNewVariable(self.node_fuel, name, tpe="model", label=name)


    def __createCoalModelProperties(self, coalsNumber, coalClassesNumber, classesNumber):
        """
        Private method
        Create new properties for one coal
        """
        # create new properties
        baseNames = self.__getPropertiesList()

        for baseName in baseNames:
            for classe in range(classesNumber - coalClassesNumber, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                self.setNewProperty(self.node_fuel, name)


    def __createClassModelProperties(self, classNum, classesNumber):
        """
        Private method
        Create class of model properties
        """
        baseNames = self.__getPropertiesList()

        # Rename other classes
        nodeList = self.node_fuel.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum + 1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name
        #
        # create new properties
        for i in range(len(baseNames)):
            name = '%s%2.2i' % (baseNames[i], classNum)
            self.setNewProperty(self.node_fuel, name)


    def __createClassModelScalars(self, classNum, classesNumber):
        """
        Private method
        Create a new coal and associated scalars
        """
        baseNames = self.__getVariableList()

        # Rename other classes
        nodeList = self.node_fuel.xmlGetNodeList('variable')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum+1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name

        # create new scalars
        for i in range(len(baseNames)):
            name = '%s%2.2i' % (baseNames[i], classNum)
            self.setNewVariable(self.node_fuel, name, tpe="model", label=name)


    @Variables.undoGlobal
    def setCoalCombustionModel(self, model):
        """
        Update the coal combustion model markup from the XML document.
        """
        self.isInList(model, self.__coalCombustionModelsList())

        if self.node_fuel['model'] == model:
            return

        # Avoid changing model settings for regular and Lagrangian
        # variants of homogeneous fuel moisture combustion model
        if model in ('homogeneous_fuel_moisture',
                     'homogeneous_fuel_moisture_lagr'):
            if self.node_fuel['model'] in ('homogeneous_fuel_moisture',
                                           'homogeneous_fuel_moisture_lagr'):
                self.node_fuel['model']  = model
                return

        self.node_fuel['model']  = model

        self.__updateScalarAndProperty(model)

        if model == 'off':
            ThermalScalarModel(self.case).setThermalModel('off')
        else:
            ThermalScalarModel(self.case).setThermalModel('enthalpy')


    def __updateScalarAndProperty(self, model):
        """
        Update scalars and properties depending on model
        """
        if model == 'off':
            for tag in ('variable',
                        'property',
                        'reference_mass_molar',
                        'reference_temperature'):
                for node in self.node_fuel.xmlGetNodeList(tag):
                    node.xmlRemoveNode()

            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "inlet":
                    Boundary("coal_inlet", zone.getLabel(), self.case).deleteCoals()

            self.node_fuel.xmlRemoveChild('oxidants')
            self.node_fuel.xmlRemoveChild('solid_fuel')
            self.node_fuel.xmlRemoveChild('CO2_kinetics')
            self.node_fuel.xmlRemoveChild('H2O_kinetics')
            self.node_fuel.xmlRemoveChild('NOx_formation')
            self.node_ray = self.node_models.xmlInitNode('radiative_transfer')
            self.node_ray.xmlRemoveChild('absorption_coefficient')

        else:
            self.createModel()
            self.createCoal()
            self.createOxidant()
            self.createCoalModelScalarsAndProperties()

            if model != 'homogeneous_fuel_moisture' and self.getCoalCombustionModel() != 'homogeneous_fuel_moisture_lagr':
                self.__updateWetScalarsAndProperty(model)


    @Variables.noUndo
    def getCoalCombustionModel(self, status="full"):
        """
        Return the current coal combustion model.
        """
        model = self.node_fuel['model']

        if model not in self.__coalCombustionModelsList():
            model = self.defaultValues()['model']
            if status == "only":
                self.node_fuel.xmlRemoveNode()
            else:
                self.setCoalCombustionModel(model)
        else:
            self.__updateSolidFuelCombustionDensity(model)

        return model


    def createCoalModelScalarsAndProperties(self):
        """
        Create new scalars and new properties for one coal
        """
        coalsNumber = self.getCoalNumber()
        coalClassesNumber = self.getClassNumber(coalsNumber)
        classesNumber = self.getClassesNumber()

        # add new scalars and properties
        self.__createCoalModelScalars(coalsNumber, coalClassesNumber, classesNumber)
        self.__createCoalModelProperties(coalsNumber, coalClassesNumber, classesNumber)


    def createClassModelScalarsAndProperties(self, coalNumber):
        """
        Create class of model scalars and properties for one given coal
        """
        self.isInt(coalNumber)

        classNum = 0
        for coal in range(0, coalNumber):
            node= self.node_fuel.xmlGetNode('solid_fuel', fuel_id = str(coal + 1))
            diameter_type = self.getDiameterType(coal + 1)

            if diameter_type == 'automatic':
                classNum += len(node.xmlGetNodeList('diameter', 'class_id'))
            else:
                classNum += len(node.xmlGetNodeList('mass_percent', 'class_id'))

        classesNumber = self.getClassesNumber()

        # create new scalars
        self.__createClassModelScalars(classNum, classesNumber)
        self.__createClassModelProperties(classNum, classesNumber)


    def deleteCoalModelScalarsAndProperties(self, coalNumber):
        """
        Delete scalars and properties for one coal
        """
        self.isInt(coalNumber)

        for classId in range(self.getClassNumber(coalNumber), 0, -1):
            self.deleteClass(coalNumber, classId)

        # Remove fuel scalars
        baseNames = [ "fr_mv1_", "fr_mv2_"]
        nodeList = self.node_fuel.xmlGetNodeList('variable')
        if nodeList != None:
            for node in nodeList :
                nameNode = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, coalNumber)
                    if (nameNode == name):
                        node.xmlRemoveNode()

        # Rename other fuels
        nodeList = self.node_fuel.xmlGetNodeList('variable')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(coalNumber + 1, self.getCoalNumber() + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum - 1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name


    def deleteClassModelScalars(self, coalNumber, classeNumber):
        """
        delete scalar model for a define coalNumber and class number
        """
        self.isInt(coalNumber)
        self.isInt(classeNumber)

        classNum = 0
        if (coalNumber >= 1) :
            for coal in range(0, coalNumber - 1):
                classNum += self.getClassNumber(coal + 1)

        # global class number
        classNum += classeNumber

        # total class number
        classesNumber = self.getClassesNumber()

        # list of variables for a class
        baseNames = self.__getVariableList()

        # Remove coal classes
        nodeList = self.node_fuel.xmlGetNodeList('variable')
        if nodeList != None:
            for node in nodeList :
                nodeName = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, classNum)
                    if (nodeName == name):
                        node.xmlRemoveNode()

        # Rename other classes
        nodeList = self.node_fuel.xmlGetNodeList('variable')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum + 1, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum - 1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name


    def deleteClassModelProperties(self, coalNumber, classeNumber):
        """
        delete properties model for a define coalNumber and class number
        """
        self.isInt(coalNumber)
        self.isInt(classeNumber)

        classNum = 0
        if (coalNumber >= 1) :
            for coal in range(0, coalNumber - 1):
                classNum += self.getClassNumber(coal + 1)
        # global class number
        classNum += classeNumber

        # total class number
        classesNumber = self.getClassesNumber()

        # list of properties for a class
        baseNames = self.__getPropertiesList()

        # Remove coal classes
        nodeList = self.node_fuel.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList :
                nodeName = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, classNum)
                    if (nodeName == name):
                        node.xmlRemoveNode()

        # Rename other classes
        nodeList = self.node_fuel.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum + 1, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum - 1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name


    @Variables.noUndo
    def getFuelNameList(self):
        """
        Return the fuel name list
        """
        fuel = []
        for node in self.node_fuel.xmlGetNodeList('solid_fuel'):
            fuel.append(node['name'])
        return fuel


    @Variables.noUndo
    def getFuelIdList(self):
        """
        return list of fuel Id's
        """
        fuel = []
        for node in self.node_fuel.xmlGetNodeList('solid_fuel'):
            fuel.append(node['fuel_id'])
        return fuel


    @Variables.noUndo
    def getLabelIdList(self):
        """
        return list of fuel label
        """
        fuel = []
        for node in self.node_fuel.xmlGetNodeList('solid_fuel'):
            fuel.append(node['name'])
        return fuel


    @Variables.noUndo
    def getCoalNumber(self):
        """
        return number of solid fuel
        """
        nb = len(self.node_fuel.xmlGetNodeList('solid_fuel'))
        return nb


    @Variables.noUndo
    def getClassIdList(self, fuelId):
        """
        return list of class_id for define fuel id
        """
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_classe = solid_fuel.xmlInitNode('class')
        diameter_type = self.getDiameterType(fuelId)
        class_list= []
        if diameter_type == 'automatic':
            for node in node_classe.xmlGetNodeList('diameter'):
                class_list.append(node['class_id'])
        elif diameter_type == 'rosin-rammler_law':
            for node in node_classe.xmlGetNodeList('mass_percent'):
                class_list.append(node['class_id'])
        return class_list


    @Variables.noUndo
    def getClassNumber(self, fuelId):
        """
        return number of class for define fuel id
        """
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlInitNode('class')
        diameter_type = self.getDiameterType(fuelId)
        nb = 0
        if diameter_type == 'automatic':
            nb = len(node_class.xmlGetNodeList('diameter'))
        elif diameter_type == 'rosin-rammler_law':
            nb = len(node_class.xmlGetNodeList('mass_percent'))
        return nb


    @Variables.noUndo
    def getClassesNumber(self):
        """
        return global number of class for fuel(s)
        """
        return len(self.case.xmlGetNodeList('diameter', 'class_id'))


    @Variables.noUndo
    def getRefusalIdList(self, fuelId):
        """
        return number of refusal for define fuel id
        """
        v = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        class_list= []
        for node in v.xmlGetNodeList('refusal'):
            class_list.append(node['id'])
        return class_list


    @Variables.noUndo
    def getRefusalNumber(self, id):
        """
        return global number of refusal
        """
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = id)
        nb = len(solid_fuel.xmlGetNodeList('refusal'))
        return nb


    @Variables.noUndo
    def getOxidantNumber(self):
        """
        return global number oxidant
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        nb = len(node_oxi.xmlGetNodeList('oxidant'))
        return nb


    @Variables.noUndo
    def getOxidantIdList(self):
        """
        return list of oxidant Id's
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        oxidant = []
        for node in node_oxi.xmlGetNodeList('oxidant'):
            oxidant.append(node['ox_id'])
        return oxidant


    @Variables.noUndo
    def createCoal(self):
        """
        create a new solid fuel
        """
        number = self.getCoalNumber()
        new = number + 1
        solid_fuel = self.node_fuel.xmlInitNode('solid_fuel', fuel_id = new)
        self.isLowerOrEqual(number, 3)
        self.createClass(new)
        char_comb = solid_fuel.xmlInitNode('char_combustion')
        char_comb.xmlInitNode('specie', nature= "O2")
        char_comb.xmlInitNode('specie', nature= "CO2")
        char_comb.xmlInitNode('specie', nature= "H2O")
        solid_fuel.xmlInitNode('nox_formation')


    def deleteSolidFuel(self, fuelId):
        """
        delete a solid fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.deleteCoalModelScalarsAndProperties(fuelId)
        node = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node.xmlRemoveNode()

        # update boundary conditions (1/2)
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            if nature == "inlet":
                bc = Boundary("coal_inlet", label, self.case)
                bc.deleteCoalFlow(fuelId, self.getCoalNumber())

        self.__updateFuelId()


    def __updateFuelId(self):
        """
        """
        n = 0
        for node in self.node_fuel.xmlGetNodeList('solid_fuel'):
            if int(node['fuel_id']) > 0 :
                n = n + 1
                node['fuel_id'] = str(n)


    def createOxidant(self):
        """
        create a new oxidant
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        number = self.getOxidantNumber()
        new = number + 1
        node_oxi.xmlInitNode('oxidant', ox_id = new)
        self.isLowerOrEqual(number, 3)
        # update list of scalar
        self.__createModelScalars()
        self.__createModelProperties()


    def deleteOxidant(self, number):
        """
        delete an oxidant
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        self.isInList(str(number), self.getOxidantIdList())
        node = node_oxi.xmlGetNode('oxidant', ox_id = number)
        node.xmlRemoveNode()

        # Update boundary conditions
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            if nature == "inlet":
                bc = Boundary("coal_inlet", label, self.case)
                oxi_max = bc.getOxidantNumber()
                if oxi_max >= number:
                    bc.setOxidantNumber(oxi_max-1)

        self.__updateOxidantId()


    def __updateOxidantId(self):
        """
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        n = 0
        for node in node_oxi.xmlGetNodeList('oxidant'):
            if int(node['ox_id']) > 0 :
                n = n + 1
                node['ox_id'] = str(n)


    @Variables.noUndo
    def getElementComposition(self, oxId, element):
        """
        return contribution of an element for a define oxidant
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        self.isInList(str(oxId), self.getOxidantIdList())
        oxidant = node_oxi.xmlGetNode('oxidant', ox_id = oxId)
        name = element + "_composition"
        value = oxidant.xmlGetDouble(name)
        if value == None:
            defName = element + "_oxi"
            value = self.defaultValues()[defName]
            self.setElementComposition(oxId, element, value)
        return value


    @Variables.undoLocal
    def setElementComposition(self, oxId, element, value):
        """
        set contribution of an element for a define oxidant
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        self.isInList(str(oxId), self.getOxidantIdList())
        self.isPositiveFloat(value)
        oxidant = node_oxi.xmlGetNode('oxidant', ox_id = oxId)
        name = element + "_composition"
        oxidant.xmlSetData(name, value)


    @Variables.noUndo
    def getOxidant(self, oxId):
        """
        return an oxidant
        """
        self.isInList(str(oxId), self.getOxidantIdList())
        O2  = self.getElementComposition(oxId, "O2")
        N2  = self.getElementComposition(oxId, "N2")
        H2O = self.getElementComposition(oxId, "H2O")
        CO2 = self.getElementComposition(oxId, "CO2")
        oxi = [id, O2, N2, H2O, CO2]
        return oxi


    @Variables.noUndo
    def getDiameterType(self, fuelId):
        """
        return diameter model for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)

        value = solid_fuel.xmlGetString('diameter_type')
        if value == '':
            value = self.defaultValues()['diameter_type']
            self.setDiameterType(fuelId, value)
        return value


    @Variables.undoGlobal
    def setDiameterType(self, fuelId, choice):
        """
        put diameter model for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(choice, ('automatic', 'rosin-rammler_law'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        old_type = solid_fuel.xmlGetString('diameter_type')
        solid_fuel.xmlSetData('diameter_type', choice)
        if old_type != '':
            if old_type != choice:
                if old_type == 'rosin-rammler_law':
                    for node in solid_fuel.xmlGetNodeList('refusal'):
                        node.xmlRemoveNode()
                solid_fuel.xmlRemoveChild('class')
                # one class needed
                self.createClass(fuelId)

                if choice == 'rosin-rammler_law':
                    self.createRefusal(fuelId)


    def createClass(self, fuelId):
        """
        create a new diameter class for define fuel Id
        """
        number = self.getClassNumber(fuelId)
        new = number + 1
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlInitNode('class')
        diameter_type = self.getDiameterType(fuelId)
        if diameter_type == 'automatic':
            node_class.xmlInitNode('diameter', class_id = new)
            self.getDiameter(fuelId, new)
        elif diameter_type == 'rosin-rammler_law':
            node_class.xmlInitNode('mass_percent', class_id = new)
            self.getMassPercent(fuelId, new)

        # Update boundary conditions
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            if zone.getNature() == "inlet":
                b = Boundary("coal_inlet", zone.getLabel(), self.case)
                b.updateCoalRatios(fuelId)


    def deleteClass(self, fuelId,  number):
        """
        delete a diameter class for a define fuel Id
        """
        self.isInList(str(number), self.getClassIdList(fuelId))

        self.deleteClassModelProperties(fuelId, number)
        self.deleteClassModelScalars(fuelId, number)

        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlInitNode('class')
        diameter_type = self.getDiameterType(fuelId)
        if diameter_type == 'automatic':
            node = node_class.xmlGetNode('diameter', class_id = number)
            node.xmlRemoveNode()
        elif diameter_type == 'rosin-rammler_law':
            node = node_class.xmlGetNode('mass_percent', class_id = number)
            node.xmlRemoveNode()

            # delete refusal if needed
            refusal_number = self.getRefusalNumber(fuelId)
            class_number = self.getClassNumber(fuelId)
            if refusal_number > class_number:
                for num in range (class_number + 1, refusal_number + 1):
                    self.deleteRefusal(fuelId, num)

        self.__updateClassId(fuelId)


    def __updateClassId(self, fuelId):
        """
        """
        n = 0
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlInitNode('class')
        diameter_type = self.getDiameterType(fuelId)
        if diameter_type == 'automatic':
            for node in node_class.xmlGetNodeList('diameter'):
                if int(node['class_id']) > 0 :
                    n = n + 1
                    node['class_id'] = str(n)
        elif diameter_type == 'rosin-rammler_law':
            for node in node_class.xmlGetNodeList('mass_percent'):
                if int(node['class_id']) > 0 :
                    n = n + 1
                    node['class_id'] = str(n)


    def createRefusal(self, fuelId):
        """
        create de new refusal for define fuel Id
        """
        number = 0
        number = self.getRefusalNumber(fuelId)
        new = number + 1
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel.xmlInitNode('refusal', id = new)
        self.getRefusalDiameter(fuelId, new)
        self.getRefusalValue(fuelId, new)


    def deleteRefusal(self, fuelId,  number):
        """
        delete a refusal for define fuel Id
        """
        self.isInList(str(number), self.getRefusalIdList(fuelId))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlGetNode('refusal', id = number)
        node.xmlRemoveNode()
        self.__updateRefusalId(fuelId)


    def __updateRefusalId(self, fuelId):
        """
        """
        n = 0
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        for node in solid_fuel.xmlGetNodeList('refusal'):
            if int(node['id']) > 0 :
                n = n + 1
                node['id'] = str(n)


    @Variables.noUndo
    def getRefusalDiameter(self, fuelId, refusal_number):
        """
        Return the refusal diameter for a define fuel Id  and refusal Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(str(refusal_number), self.getRefusalIdList(fuelId))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_refusal = solid_fuel.xmlGetNode('refusal', id = refusal_number)
        value = node_refusal.xmlGetDouble('diameter')
        if value == None:
            value = self.defaultValues()['diameter']
            self.setRefusalDiameter(fuelId, refusal_number, value)
        return value


    @Variables.undoLocal
    def setRefusalDiameter(self, fuelId, refusal_number, value):
        """
        Put the refusal diameter for a define fuel Id  and refusal Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(str(refusal_number), self.getRefusalIdList(fuelId))
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_refusal = solid_fuel.xmlGetNode('refusal', id = refusal_number)
        node_refusal.xmlSetData('diameter', value)


    @Variables.noUndo
    def getRefusalValue(self, fuelId, refusal_number):
        """
        Return the refusal value for a define fuel Id  and refusal Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(str(refusal_number), self.getRefusalIdList(fuelId))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_refusal = solid_fuel.xmlGetNode('refusal', id = refusal_number)
        value = node_refusal.xmlGetDouble('value')
        if value == None:
            value = self.defaultValues()['refusal_value']
            self.setRefusalValue(fuelId, refusal_number, value)
        return value


    @Variables.noUndo
    def getRefusal(self, fuelId, refusal_number):
        """
        Return all characteristics of a refusal
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(str(refusal_number), self.getRefusalIdList(fuelId))
        diameter = self.getRefusalDiameter(fuelId, refusal_number)
        value = self.getRefusalValue(fuelId, refusal_number)
        return [refusal_number, diameter, value]


    @Variables.undoLocal
    def setRefusalValue(self, fuelId, refusal_number, value):
        """
        Put the refusal value for a define fuel Id  and refusal Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(str(refusal_number), self.getRefusalIdList(fuelId))
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_refusal = solid_fuel.xmlGetNode('refusal', id = refusal_number)
        node_refusal.xmlSetData('value', value)


    @Variables.noUndo
    def getDiameter(self, fuelId, class_number):
        """
        Return diameter for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlGetNode('class')
        value = node_class.xmlGetDouble('diameter', class_id = class_number)
        if value == None:
            value = self.defaultValues()['diameter']
            self.setDiameter(fuelId, class_number, value)
        return value


    @Variables.undoLocal
    def setDiameter(self, fuelId, class_number, value):
        """
        Put diameter for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlGetNode('class')
        node_class.xmlSetData('diameter', value, class_id = class_number)


    @Variables.noUndo
    def getMassPercent(self, fuelId, class_number):
        """
        Return mass percent for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_classe = solid_fuel.xmlGetNode('class')
        value = node_classe.xmlGetDouble('mass_percent', class_id = class_number)
        if value == None:
            value = self.defaultValues()['mass_percent']
            self.setMassPercent(fuelId, class_number, value)
        return value


    @Variables.undoLocal
    def setMassPercent(self, fuelId, class_number, value):
        """
        Put mass percent for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node_class = solid_fuel.xmlGetNode('class')
        node_class.xmlSetData('mass_percent', value, class_id = class_number)


    @Variables.noUndo
    def getFuelLabel(self, fuelId):
        """
        Return label for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        label = solid_fuel['name']
        if label == None:
            label = self.__getDefaultLabel()
            self.setFuelLabel(fuelId, label)
        return label


    def __getDefaultLabel(self):
        """
        Determine a default label for a fuel
        """
        name = "SolidFuel_" + str(self.getCoalNumber())
        if name in self.getFuelNameList():
            labelNumber = 1
            name = "SolidFuel_" + str(labelNumber)
            while name in self.getFuelNameList():
                labelNumber += 1
                name = "SolidFuel_" + str(labelNumber)
        return name


    @Variables.undoLocal
    def setFuelLabel(self, fuelId, label):
        """
        Set a fuel label
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel['name']= label


    @Variables.noUndo
    def getFuelType(self, fuelId):
        """
        Return the type for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        fuel_type = solid_fuel['type']
        if fuel_type == None:
            fuel_type = self.defaultValues()['fuel_type']
            self.setFuelType(fuelId, fuel_type)
        return fuel_type


    @Variables.undoLocal
    def setFuelType(self, fuelId, fuel_type):
        """
        Set the type for a define fuel Id
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(fuel_type, ('biomass', 'coal' ))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel['type']= fuel_type


    @Variables.noUndo
    def getComposition(self, fuelId, element):
        """
        Return composition for a define fuel Id and an element
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        name = element + "_composition_on_dry"
        composition = solid_fuel.xmlGetDouble(name)
        if composition == None:
            defName = element + "_compo_dry"
            composition = self.defaultValues()[defName]
            self.setComposition(fuelId, element, composition)
        return composition


    @Variables.undoLocal
    def setComposition(self, fuelId, element, composition):
        """
        Set composition for a define fuel Id and an element
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(composition)
        name = element + "_composition_on_dry"
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel.xmlSetData(name, composition)


    @Variables.noUndo
    def getCokeComposition(self, fuelId, element):
        """
        Return composition for a define fuel Id and an element
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        name = element + "_coke_composition_on_dry"
        composition = solid_fuel.xmlGetDouble(name)
        if composition == None:
            defName = element + "_coke_compo_dry"
            composition = self.defaultValues()[defName]
            self.setCokeComposition(fuelId, element, composition)
        return composition


    @Variables.undoLocal
    def setCokeComposition(self, fuelId, element, composition):
        """
        Set composition for a define fuel Id and an element
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(composition)
        name = element + "_coke_composition_on_dry"
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel.xmlSetData(name, composition)


    @Variables.noUndo
    def getPCIValue(self, fuelId):
        """
        Return PCI Value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        value = node.xmlGetDouble('value')
        if value == None:
            value = self.defaultValues()['PCI']
            self.setPCIValue(fuelId, value)
        return value


    @Variables.undoLocal
    def setPCIValue(self, fuelId, value):
        """
        Set PCI Value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        node.xmlSetData('value', value)


    @Variables.noUndo
    def getPCIChoice(self, fuelId):
        """
        Return PCI choice for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        PCIChoice = node['choice']
        if PCIChoice == None:
            PCIChoice = self.defaultValues()['PCIChoice']
            self.setPCIChoice(fuelId, PCIChoice)
        return PCIChoice


    @Variables.undoLocal
    def setPCIChoice(self, fuelId, choice):
        """
        Set PCI choice for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(choice, ('LHV', 'HHV', 'IGT_correlation'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        node['choice'] = choice
        if choice == 'IGT_correlation':
            node.xmlRemoveChild('value')
            node.xmlRemoveChild('type')


    @Variables.noUndo
    def getPCIType(self, fuelId):
        """
        Return PCI type for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        value = node.xmlGetString('type')
        if value == None:
            value = self.defaultValues()['PCIType']
            self.setPCIType(fuelId, value)
        return value


    @Variables.undoLocal
    def setPCIType(self, fuelId, choice):
        """
        Set PCI type for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(choice, ('dry_basis', 'dry_ash_free', 'as_received'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('Heating_model')
        node.xmlSetData('type', choice)


    @Variables.noUndo
    def getProperty(self, fuelId, name):
        """
        Return value for a define fuel Id and property
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        value = solid_fuel.xmlGetDouble(name)
        if value == None:
            value = self.defaultValues()[name]
            self.setProperty(fuelId, name, value)
        return value


    @Variables.undoLocal
    def setProperty(self, fuelId, name, value):
        """
        Set value for a define fuel Id and property
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        solid_fuel.xmlSetData(name, value)


    @Variables.noUndo
    def getY1Y2(self, fuelId):
        """
        Return Y1Y2 value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        node = devolatilisation.xmlInitNode('stoichiometric_coefficient')
        choice = node['type']
        if choice == None:
            choice = self.defaultValues()['stoichiometric_coefficient']
            self.setY1Y2(fuelId, choice)
        return choice


    @Variables.undoLocal
    def setY1Y2(self, fuelId, choice):
        """
        Set Y1Y2 value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(choice, ('user_define', 'automatic_CHONS', 'automatic_formula'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('stoichiometric_coefficient')
        node['type'] = choice
        node.xmlRemoveChild('Y1')
        node.xmlRemoveChild('Y2')


    @Variables.noUndo
    def getY1StoichiometricCoefficient(self, fuelId):
        """
        Return Y1 stoichiometric coefficient for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(self.getY1Y2(fuelId),('user_define', 'automatic_formula'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        node = devolatilisation.xmlInitNode('stoichiometric_coefficient')
        value = devolatilisation.xmlGetDouble('Y1')
        if value == None:
            choice = self.getY1Y2(fuelId)
            if choice == 'automatic_formula':
                value = self.__Y1AutomaticFormula(fuelId)
            else:
                value = self.defaultValues()['Y1']
            self.setY1StoichiometricCoefficient(fuelId, value)
        return value


    @Variables.undoLocal
    def setY1StoichiometricCoefficient(self, fuelId, value):
        """
        Set Y1 stoichiometric coefficient for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(self.getY1Y2(fuelId),('user_define', 'automatic_formula'))
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        node = solid_fuel.xmlInitNode('stoichiometric_coefficient')
        node.xmlSetData('Y1', value)


    @Variables.noUndo
    def getY2StoichiometricCoefficient(self, fuelId):
        """
        Return Y2 stoichiometric coefficient for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(self.getY1Y2(fuelId),('user_define', 'automatic_formula'))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        node = devolatilisation.xmlInitNode('stoichiometric_coefficient')
        value = devolatilisation.xmlGetDouble('Y2')
        if value == None:
            choice = self.getY1Y2(fuelId)
            if choice == 'automatic_formula':
                value = self.__Y2AutomaticFormula(fuelId)
            else:
                value = self.defaultValues()['Y2']
            self.setY2StoichiometricCoefficient(fuelId, value)
        return value


    @Variables.undoLocal
    def setY2StoichiometricCoefficient(self, fuelId, value):
        """
        Set Y2 stoichiometric coefficient for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(self.getY1Y2(fuelId),('user_define', 'automatic_formula'))
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        node = solid_fuel.xmlInitNode('stoichiometric_coefficient')
        node.xmlSetData('Y2', value)

    def __Y1AutomaticFormula(self, fuelId):
        """
        Return Y1 automatic value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        volatile_matter = self.getProperty(fuelId, "volatile_matter")
        ash_ratio = self.getProperty(fuelId, "rate_of_ashes_on_mass")
        moisture = self.getProperty(fuelId, "moisture")
        Y1 = (volatile_matter / 100) / (1 - ash_ratio / 100 - moisture / 100)
        return Y1

    def __Y2AutomaticFormula(self, fuelId):
        """
        Return Y2 automatic value for a fuel
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        Y1 = self.__Y1AutomaticFormula(fuelId)
        Y2 = 1.4 * Y1
        return Y2


    @Variables.noUndo
    def getHCNParameter(self, fuelId, param):
        """
        Return value for a define fuel Id and nitrogen partition parameter
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        value = devolatilisation.xmlGetDouble(param)
        if value == None:
            value = self.defaultValues()[param]
            self.setHCNParameter(fuelId, param, value)
        return value


    @Variables.undoLocal
    def setHCNParameter(self, fuelId, param, value):
        """
        Set value for a define fuel Id and nitrogen partition parameter
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        devolatilisation.xmlSetData(param, value)


    @Variables.noUndo
    def getDevolatilisationParameter(self, fuelId, param):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        value = devolatilisation.xmlGetDouble(param)
        if value == None:
            value = self.defaultValues()[param]
            self.setDevolatilisationParameter(fuelId, param, value)
        return value


    @Variables.undoLocal
    def setDevolatilisationParameter(self, fuelId, param, value):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        devolatilisation = solid_fuel.xmlInitNode('devolatilisation_parameters')
        devolatilisation.xmlSetData(param, value)


    @Variables.noUndo
    def getPreExponentialConstant(self, fuelId, specie):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature = specie)
        value = specie_node.xmlGetDouble('pre-exponential_constant')
        if value == None:
            value = self.defaultValues()['pre-exp_constant']
            self.setPreExponentialConstant(fuelId, specie, value)
        return value


    @Variables.undoLocal
    def setPreExponentialConstant(self, fuelId, specie, value):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature = specie)
        specie_node.xmlSetData('pre-exponential_constant', value)


    @Variables.noUndo
    def getEnergyOfActivation(self, fuelId, specie):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature = specie)
        value = specie_node.xmlGetDouble('energy_of_activation')
        if value == None:
            value = self.defaultValues()['E_activation']
            self.setEnergyOfActivation(fuelId, specie, value)
        return value


    @Variables.undoLocal
    def setEnergyOfActivation(self, fuelId, specie, value):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature= specie)
        specie_node.xmlSetData('energy_of_activation', value)


    @Variables.noUndo
    def getOrderOfReaction(self, fuelId, specie):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature = specie)
        node = specie_node.xmlInitNode('order_of_reaction')
        choice = node['choice']
        if choice == None:
            choice = self.defaultValues()['order_reaction']
            self.setOrderOfReaction(fuelId, specie, choice)
        return choice


    @Variables.undoLocal
    def setOrderOfReaction(self, fuelId, specie, choice):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isInList(choice, ("1", "0.5"))
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        char_comb = solid_fuel.xmlGetNode('char_combustion')
        specie_node = char_comb.xmlGetNode('specie', nature = specie)
        node = specie_node.xmlInitNode('order_of_reaction')
        node['choice'] = choice


    @Variables.noUndo
    def getNOxFormationParameter(self, fuelId, parameter):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlInitNode('nox_formation')
        value = node.xmlGetDouble(parameter)
        if value == None:
            value = self.defaultValues()[parameter]
            self.setNOxFormationParameter(fuelId, parameter, value)
        return value


    @Variables.undoLocal
    def setNOxFormationParameter(self, fuelId, parameter, value):
        """
        """
        self.isInList(str(fuelId), self.getFuelIdList())
        self.isPositiveFloat(value)
        solid_fuel = self.node_fuel.xmlGetNode('solid_fuel', fuel_id = fuelId)
        node = solid_fuel.xmlGetNode('nox_formation')
        node.xmlSetData(parameter, value)


    @Variables.noUndo
    def getOxidantType(self):
        """
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        value = node_oxi.xmlGetString('oxidant_type')
        if value == '':
            value = self.defaultValues()['oxidant_type']
            self.setOxidantType(value)
        return value


    @Variables.undoLocal
    def setOxidantType(self, choice):
        """
        """
        node_oxi = self.node_fuel.xmlInitNode('oxidants')
        self.isInList(choice, ('volumic_percent', 'molar' ))
        node_oxi.xmlSetData('oxidant_type', choice)


    @Variables.noUndo
    def getReburning(self, fuelId):
        """
        """
        fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
        node = fuelNode.xmlGetNode("nox_formation")
        value = node.xmlGetString('reburning_model')
        if value == '':
            value = self.defaultValues()['reburning']
            self.setReburning(fuelId, value)
        return value


    @Variables.undoLocal
    def setReburning(self, fuelId, choice):
        """
        """
        self.isInList(choice, ('unused', 'chen', 'dimitriou'))
        fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
        node = fuelNode.xmlGetNode("nox_formation")
        node.xmlSetData('reburning_model', choice)


    @Variables.undoGlobal
    def setNOxFormationStatus(self, status):
        """
        put NOx formation status
        """
        self.isOnOff(status)

        node = self.node_fuel.xmlGetNode('NOx_formation')
        node['status'] = status
        if status == "off":
            for fuelId in range (1, self.getCoalNumber() + 1):
                fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
                devol_node = fuelNode.xmlGetNode("devolatilisation_parameters")
                nox_node = fuelNode.xmlGetNode("nox_formation")
                nox_node.xmlRemoveChild('nitrogen_fraction')
                nox_node.xmlRemoveChild('nitrogen_concentration')
                devol_node.xmlRemoveChild('HCN_NH3_partitionning_reaction_1')
                devol_node.xmlRemoveChild('HCN_NH3_partitionning_reaction_2')
        else:
            for fuelId in range (1, self.getCoalNumber() + 1):
                fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
                self.getNOxFormationParameter(fuelId, 'nitrogen_fraction')
                self.getNOxFormationParameter(fuelId, 'nitrogen_concentration')
                self.getNOxFormationParameter(fuelId, 'nitrogen_in_char_at_low_temperatures')
                self.getNOxFormationParameter(fuelId, 'nitrogen_in_char_at_high_temperatures')
                self.getNOxFormationParameter(fuelId, 'percentage_HCN_char_combustion')
                self.getHCNParameter(fuelId, 'HCN_NH3_partitionning_reaction_1')
                self.getHCNParameter(fuelId, 'HCN_NH3_partitionning_reaction_2')
        # update list of scalar
        self.__createModelScalars()
        self.__createModelProperties()


    @Variables.noUndo
    def getNOxFormationStatus(self):
        """
        get NOx formation status
        """
        node = self.node_fuel.xmlGetNode('NOx_formation')
        status = node['status']

        return status


    @Variables.undoGlobal
    def setNOxFormationFeature(self, fuelId, status):
        """
        put NOx formation status
        """
        self.isOnOff(status)

        fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
        node_nox = fuelNode.xmlGetNode("nox_formation")
        node = node_nox.xmlInitNode('improved_NOx_model')
        node['status'] = status


    @Variables.noUndo
    def getNOxFormationFeature(self, fuelId):
        """
        get NOx formation feature status
        """
        fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
        node_nox = fuelNode.xmlGetNode("nox_formation")
        node = node_nox.xmlGetNode('improved_NOx_model')
        if node == None:
            self.setNOxFormationFeature(fuelId, self.defaultValues()['improved_NOx_model'])
            node = node_nox.xmlGetNode('improved_NOx_model')
        status = node['status']

        return status


    @Variables.undoGlobal
    def setCO2KineticsStatus(self, status):
        """
        put CO2 Kinetics status
        """
        self.isOnOff(status)

        node = self.node_fuel.xmlGetNode('CO2_kinetics')
        node['status'] = status
        if status == "off":
            for fuelId in range (1, self.getCoalNumber() + 1):
                fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
                combu_node = fuelNode.xmlGetNode("char_combustion")
                co2_node = combu_node.xmlGetNode("specie", nature = "CO2")
                co2_node.xmlRemoveChild('pre-exponential_constant')
                co2_node.xmlRemoveChild('energy_of_activation')
                co2_node.xmlRemoveChild('order_of_reaction')
        else:
            for fuelId in range (1, self.getCoalNumber() + 1):
                self.getPreExponentialConstant(fuelId, "CO2")
                self.getEnergyOfActivation(fuelId, "CO2")
                self.getOrderOfReaction(fuelId, "CO2")
        # update list of scalar
        self.__createModelScalars()
        self.__createModelProperties()


    @Variables.noUndo
    def getCO2KineticsStatus(self):
        """
        get CO2 Kinetics status
        """
        node = self.node_fuel.xmlGetNode('CO2_kinetics')
        status = node['status']

        return status


    @Variables.undoGlobal
    def setH2OKineticsStatus(self, status):
        """
        put H2O Kinetics status
        """
        self.isOnOff(status)

        node = self.node_fuel.xmlGetNode('H2O_kinetics')
        node['status'] = status
        if status == "off":
            for fuelId in range (1, self.getCoalNumber() + 1):
                fuelNode = self.node_fuel.xmlGetNode("solid_fuel", fuel_id = fuelId)
                combu_node = fuelNode.xmlGetNode("char_combustion")
                h2o_node = combu_node.xmlGetNode("specie", nature = "H2O")
                h2o_node.xmlRemoveChild('pre-exponential_constant')
                h2o_node.xmlRemoveChild('energy_of_activation')
                h2o_node.xmlRemoveChild('order_of_reaction')
        else:
            for fuelId in range (1, self.getCoalNumber() + 1):
                self.getPreExponentialConstant(fuelId, "H2O")
                self.getEnergyOfActivation(fuelId, "H2O")
                self.getOrderOfReaction(fuelId, "H2O")
        # update list of scalar
        self.__createModelScalars()
        self.__createModelProperties()


    @Variables.noUndo
    def getH2OKineticsStatus(self):
        """
        get H2O Kinetics status
        """
        node = self.node_fuel.xmlGetNode('H2O_kinetics')
        status = node['status']

        return status


#-------------------------------------------------------------------------------
# SolidFuel combustion test case
#-------------------------------------------------------------------------------

class CoalCombustionModelTestCase(ModelTest):
    """
    """

    def checkCoalCombustionModelInstantiation(self):
        """Check whether the CoalCombustionModel class could be instantiated"""
        model = None
        model = CoalCombustionModel(self.case)
        assert model != None, 'Could not instantiate coalCombustionModel'

    def checkSetAndGetCoalCombustionModel(self):
        """Check whether the CoalCombustionModel class could be set and get combustion model"""
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        doc = '''<solid_fuels model="coal_homo">
                    <variable label="NP_CP01" name="n_p_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP01" name="x_p_coal_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP01" name="x_p_char_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV101" name="fr_mv1_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV201" name="fr_mv2_01" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_HET" name="het_fraction" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET_O201" name="het_ts_o2_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
            </solid_fuels>'''

        assert model.node_fuel == self.xmlNodeFromString(doc),\
           'Could not set a model of combustion'
        assert model.getCoalCombustionModel() == 'coal_homo',\
           'Could not get a model of combustion'


    def checkCreateSolidFuelModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        created new variables and properties for one new coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        model.createCoalModelScalarsAndProperties(coalThermoChModel)

        doc = '''<solid_fuels model="coal_homo">
                    <variable label="NP_CP01" name="n_p_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCH_CP01" name="x_p_coal_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCK_CP01" name="x_p_char_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV101" name="fr_mv1_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV201" name="fr_mv2_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_HET" name="het_fraction" type="model"><flux_reconstruction status="off"/></variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET_O201" name="het_ts_o2_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <variable label="NP_CP02" name="n_p_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP02" name="x_p_coal_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP02" name="x_p_char_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent02" name="x_p_h_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV102" name="fr_mv1_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV202" name="fr_mv2_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP02" name="t_coal02"/>
                    <property label="Frm_CP02" name="w_solid_coal02"/>
                    <property label="Rho_CP02" name="rho_coal02"/>
                    <property label="Dia_CK02" name="diameter_coal02"/>
                    <property label="Ga_DCH02" name="dissapear_rate_coal02"/>
                    <property label="Ga_DV102" name="m_transfer_v1_coal02"/>
                    <property label="Ga_DV202" name="m_transfer_v2_coal02"/>
                    <property label="Ga_HET_O202" name="het_ts_o2_coal02"/>
            </solid_fuels>'''

        assert model.node_fuel == self.xmlNodeFromString(doc),\
            'Could not create newvariables and properties for new coal'


    def checkCreateClassModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        created a new class of variables and properties for one coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)

        doc = '''<solid_fuels model="coal_homo">
                    <variable label="NP_CP01" name="n_p_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCH_CP01" name="x_p_coal_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCK_CP01" name="x_p_char_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV101" name="fr_mv1_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV201" name="fr_mv2_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_HET" name="het_fraction" type="model"><flux_reconstruction status="off"/></variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET_O201" name="het_ts_o2_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <variable label="NP_CP03" name="n_p_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP03" name="x_p_coal_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP03" name="x_p_char_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent03" name="x_p_h_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV102" name="fr_mv1_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV202" name="fr_mv2_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP02" name="t_coal02"/>
                    <property label="Frm_CP02" name="w_solid_coal02"/>
                    <property label="Rho_CP02" name="rho_coal02"/>
                    <property label="Dia_CK02" name="diameter_coal02"/>
                    <property label="Ga_DCH02" name="dissapear_rate_coal02"/>
                    <property label="Ga_DV102" name="m_transfer_v1_coal02"/>
                    <property label="Ga_DV202" name="m_transfer_v2_coal02"/>
                    <property label="Ga_HET_O202" name="het_ts_o2_coal02"/>
                    <variable label="NP_CP04" name="n_p_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP04" name="x_p_coal_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP04" name="x_p_char_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent04" name="x_p_h_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV103" name="fr_mv1_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV203" name="fr_mv2_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP03" name="t_coal03"/>
                    <property label="Frm_CP03" name="w_solid_coal03"/>
                    <property label="Rho_CP03" name="rho_coal03"/>
                    <property label="Dia_CK03" name="diameter_coal03"/>
                    <property label="Ga_DCH03" name="dissapear_rate_coal03"/>
                    <property label="Ga_DV103" name="m_transfer_v1_coal03"/>
                    <property label="Ga_DV203" name="m_transfer_v2_coal03"/>
                    <property label="Ga_HET_O203" name="het_ts_o2_coal03"/>
                    <variable label="NP_CP02" name="n_p_02" type="model"/>
                    <variable label="XCH_CP02" name="x_p_coal_02" type="model"/>
                    <variable label="XCK_CP02" name="x_p_char_02" type="model"/>
                    <variable label="Xp_Ent02" name="x_p_h_02" type="model"/>
                    <property label="Temp_CP04" name="t_coal04"/>
                    <property label="Frm_CP04" name="w_solid_coal04"/>
                    <property label="Rho_CP04" name="rho_coal04"/>
                    <property label="Dia_CK04" name="diameter_coal04"/>
                    <property label="Ga_DCH04" name="dissapear_rate_coal04"/>
                    <property label="Ga_DV104" name="m_transfer_v1_coal04"/>
                    <property label="Ga_DV204" name="m_transfer_v2_coal04"/>
                    <property label="Ga_HET_O204" name="het_ts_o2_coal04"/>
                </solid_fuels>'''

        assert model.node_fuel == self.xmlNodeFromString(doc),\
           'Could not create a new class of variables and properties for one coal'


    def checkDeleteSolidFuelModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        deleted class of variables and properties for one given coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)
        model.deleteCoalModelScalarsAndProperties(coalThermoChModel, 1)

        doc = '''<solid_fuels model="coal_homo">
                    <variable label="NP_CP01" name="n_p_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCH_CP01" name="x_p_coal_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCK_CP01" name="x_p_char_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV101" name="fr_mv1_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV201" name="fr_mv2_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_HET" name="het_fraction" type="model"><flux_reconstruction status="off"/></variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET_O201" name="het_ts_o2_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <variable label="NP_CP02" name="n_p_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP02" name="x_p_coal_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP02" name="x_p_char_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent02" name="x_p_h_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP02" name="t_coal02"/>
                    <property label="Frm_CP02" name="w_solid_coal02"/>
                    <property label="Rho_CP02" name="rho_coal02"/>
                    <property label="Dia_CK02" name="diameter_coal02"/>
                    <property label="Ga_DCH02" name="dissapear_rate_coal02"/>
                    <property label="Ga_DV102" name="m_transfer_v1_coal02"/>
                    <property label="Ga_DV202" name="m_transfer_v2_coal02"/>
                    <property label="Ga_HET_O202" name="het_ts_o2_coal02"/>
                    <variable label="NP_CP04" name="n_p_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP04" name="x_p_coal_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP04" name="x_p_char_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent04" name="x_p_h_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV102" name="fr_mv1_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV202" name="fr_mv2_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP04" name="t_coal04"/>
                    <property label="Frm_CP04" name="w_solid_coal04"/>
                    <property label="Rho_CP04" name="rho_coal04"/>
                    <property label="Dia_CK04" name="diameter_coal04"/>
                    <property label="Ga_DCH04" name="dissapear_rate_coal04"/>
                    <property label="Ga_DV104" name="m_transfer_v1_coal04"/>
                    <property label="Ga_DV204" name="m_transfer_v2_coal04"/>
                    <property label="Ga_HET_O204" name="het_ts_o2_coal04"/>
                </solid_fuels>'''

        assert model.node_fuel == self.xmlNodeFromString(doc),\
           "Could not delete one coal's variables and properties"


    def checkDeleteClassModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        deleted class of variables and properties for one given coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)
        model.deleteClassModelScalars(coalThermoChModel, 2, 1)

        doc = '''<solid_fuels model="coal_homo">
                    <variable label="NP_CP01" name="n_p_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCH_CP01" name="x_p_coal_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="XCK_CP01" name="x_p_char_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV101" name="fr_mv1_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_MV201" name="fr_mv2_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_HET" name="het_fraction" type="model"><flux_reconstruction status="off"/></variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET_O201" name="het_ts_o2_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <variable label="Fr_MV102" name="fr_mv1_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV202" name="fr_mv2_02" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP03" name="t_coal03"/>
                    <property label="Frm_CP03" name="w_solid_coal03"/>
                    <property label="Rho_CP03" name="rho_coal03"/>
                    <property label="Dia_CK03" name="diameter_coal03"/>
                    <property label="Ga_DCH03" name="dissapear_rate_coal03"/>
                    <property label="Ga_DV103" name="m_transfer_v1_coal03"/>
                    <property label="Ga_DV203" name="m_transfer_v2_coal03"/>
                    <property label="Ga_HET_O203" name="het_ts_o2_coal03"/>
                    <variable label="NP_CP04" name="n_p_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCH_CP04" name="x_p_coal_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="XCK_CP04" name="x_p_char_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Xp_Ent04" name="x_p_h_04" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV103" name="fr_mv1_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <variable label="Fr_MV203" name="fr_mv2_03" type="model">
                            <flux_reconstruction status="off"/>
                    </variable>
                    <property label="Temp_CP04" name="t_coal04"/>
                    <property label="Frm_CP04" name="w_solid_coal04"/>
                    <property label="Rho_CP04" name="rho_coal04"/>
                    <property label="Dia_CK04" name="diameter_coal04"/>
                    <property label="Ga_DCH04" name="dissapear_rate_coal04"/>
                    <property label="Ga_DV104" name="m_transfer_v1_coal04"/>
                    <property label="Ga_DV204" name="m_transfer_v2_coal04"/>
                    <property label="Ga_HET_O204" name="het_ts_o2_coal04"/>
                    <variable label="NP_CP02" name="n_p_02" type="model"/>
                    <variable label="XCH_CP02" name="x_p_coal_02" type="model"/>
                    <variable label="XCK_CP02" name="x_p_char_02" type="model"/>
                    <variable label="Xp_Ent02" name="x_p_h_02" type="model"/>
                    <property label="Temp_CP02" name="t_coal02"/>
                    <property label="Frm_CP02" name="w_solid_coal02"/>
                    <property label="Rho_CP02" name="rho_coal02"/>
                    <property label="Dia_CK02" name="diameter_coal02"/>
                    <property label="Ga_DCH02" name="dissapear_rate_coal02"/>
                    <property label="Ga_DV102" name="m_transfer_v1_coal02"/>
                    <property label="Ga_DV202" name="m_transfer_v2_coal02"/>
                    <property label="Ga_HET_O202" name="het_ts_o2_coal02"/>
            </solid_fuels>'''

        assert model.node_fuel == self.xmlNodeFromString(doc),\
           "Could not delete one class of variables and properties for one given coal"


def suite():
    testSuite = unittest.makeSuite(CoalCombustionModelTestCase, "check")
    return testSuite


def runTest():
    print("CoalCombustionModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
