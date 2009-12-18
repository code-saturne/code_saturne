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

from Base.Common import *
import Base.Toolbox as Tool
import Pages.CoalThermoChemistry as CoalThermoChemistry
import Pages.IdentityAndPathesModel as IdentityAndPathesModel
from Base.XMLvariables import Variables, Model
from Base.XMLmodel import ModelTest
from Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from Pages.NumericalParamEquationModel import NumericalParamEquatModel
from ThermalRadiationModel import ThermalRadiationModel
from ConjugateHeatTransferModel import ConjugateHeatTransferModel
from LocalizationModel import LocalizationModel
from Boundary import Boundary

#-------------------------------------------------------------------------------
# Coal combustion model class
#-------------------------------------------------------------------------------

class CoalCombustionModel(Variables, Model):
    """
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        nModels         = self.case.xmlGetNode('thermophysical_models')
        self.node_lagr  = self.case.xmlGetNode('lagrangian', 'model')
        self.node_turb  = nModels.xmlGetNode('turbulence',       'model')
        self.node_gas   = nModels.xmlInitNode('gas_combustion',  'model')
        self.node_coal  = nModels.xmlInitNode('pulverized_coal', 'model')
        self.node_joule = nModels.xmlInitNode('joule_effect',    'model')
        self.node_therm = nModels.xmlInitNode('thermal_scalar',  'model')

        self.coalCombustionModel = ('off', 'coal_homo', 'coal_homo2')


    def defaultValues(self):
        """
        Private method
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['model'] = "off"
        default['diameter'] = 0.0001
        default['ihtco2'] = 0
        default['ieqco2'] = 0

        return default


    def __coalCombustionModelsList(self):
        """
        Private method
        Create a tuple with the coal combustion models allowed
        by the calculation features.
        """
        coalCombustionList = self.coalCombustionModel

        if self.node_lagr and self.node_lagr['model'] != 'off':
            coalCombustionList = ('off', 'coal_homo', 'coal_homo2')

        n, m = FluidCharacteristicsModel(self.case).getThermalModel()
        if m != "off" and m not in coalCombustionList:
            coalCombustionList = ('off',)

        if self.node_turb['model'] not in ('k-epsilon',
                                           'k-epsilon-PL',
                                           'Rij-epsilon',
                                           'Rij-SSG',
                                           'v2f-phi',
                                           'k-omega-SST'):
            coalCombustionList = ('off',)

        return coalCombustionList


    def __createModelScalarsList(self , thermoChemistryModel):
        """
        Private method
        Create model scalar list
        """
        coalsNumber = thermoChemistryModel.getCoals().getCoalNumber()
        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())

        list = []
        # add new scalars
        list.append("Enthalpy")

        baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "XWT_CP", "ENT_CP"]

        for baseName in baseNames:
            for classe in range(0,classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                list.append(name)

        baseNames = [ "Fr_MV1", "Fr_MV2"]
        for baseName in baseNames:
            for coal in range(0,coalsNumber):
                name = '%s%2.2i' % (baseName, coal+1)
                list.append(name)

        self.setNewModelScalar(self.node_coal, "Fr_HET_O2")
        list.append("Fr_HET_O2")

        if self.defaultValues()['ihtco2'] == 1:
            list.append("Fr_HET_CO2")

        list.append("Var_AIR")

        if self.getCoalCombustionModel() == 'coal_homo2':
            list.append("FR_H20")

        if thermoChemistryModel.getOxydants().getNumber() >= 2:
            list.append("FR_OXYD2")

        if thermoChemistryModel.getOxydants().getNumber() == 3:
            list.append("FR_OXYD3")

        if self.defaultValues()['ieqco2'] == 1:
            list.append("FR_CO2")

        return list


    def __createModelScalars(self , thermoChemistryModel):
        """
        Private method
        Create model scalar
        """
        previous_list = []
        nodes = self.node_coal.xmlGetChildNodeList('scalar')
        for node in nodes:
            previous_list.append(node['name'])

        new_list = self.__createModelScalarsList(thermoChemistryModel)
        for name in previous_list:
            if name not in new_list:
                self.node_coal.xmlRemoveChild('scalar',  name = name)

        for name in new_list:
            if name not in previous_list:
                self.setNewModelScalar(self.node_coal, name)

        NPE = NumericalParamEquatModel(self.case)
        for node in self.node_coal.xmlGetChildNodeList('scalar'):
            NPE.setBlendingFactor(node['label'], 0.)
            NPE.setScheme(node['label'], 'upwind')
            NPE.setFluxReconstruction(node['label'], 'off')


    def __createModelPropertiesList(self, thermoChemistryModel):
        """
        Private method
        Create model properties
        """
        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())

        list = []
        list.append("Temp_GAZ")
        list.append("ROM_GAZ")
        list.append("YM_CHx1m")
        list.append("YM_CHx2m")
        list.append("YM_CO")
        list.append("YM_O2")
        list.append("YM_CO2")
        list.append("YM_H2O")
        list.append("YM_N2")
        list.append("XM")

        baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                     "Ga_DV1", "Ga_DV2", "Ga_HET_O2"]

        if self.defaultValues()['ihtco2'] == 1:
            baseNames.append("Ga_HET_CO2")

        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames.append("Ga_SEC")

        for baseName in baseNames:
            for classe in range(0,classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                list.append(name)

        list.append("IntLuminance_4PI")

        return list


    def __createModelProperties(self, thermoChemistryModel):
        """
        Private method
        Create model properties
        """
        previous_list = []
        nodes = self.node_coal.xmlGetChildNodeList('property')
        for node in nodes:
            previous_list.append(node['name'])

        new_list = self.__createModelPropertiesList(thermoChemistryModel)
        for name in previous_list:
            if name not in new_list:
                self.node_coal.xmlRemoveChild('property',  name = name)

        for name in new_list:
            if name not in previous_list:
                self.setNewProperty(self.node_coal, name)


    def createModel (self) :
        """
        Private method
        Create scalars and properties when coal combustion is selected
        """
        model = CoalThermoChemistry.CoalThermoChemistryModel("dp_FCP", self.case)
        self.__createModelScalars(model)
        self.__createModelProperties(model)


    def __deleteWetScalarsAndProperty(self):
        """
        Private method
        Delete scalars XWT_CP and Fr_H20 and property Ga_SEC
        if model is'nt 'coal_homo2'
        """
        if self.getCoalCombustionModel() != 'coal_homo2':
            nod = self.node_coal.xmlGetNode('scalar', type="model", name="FR_H20")
            if nod:
                nod.xmlRemoveNode()
            for node in self.node_coal.xmlGetNodeList('scalar', type="model"):
                if node['name'][:6] == "XWT_CP":
                    node.xmlRemoveNode()
            for node in self.node_coal.xmlGetNodeList('property'):
                if node['name'][:6] == "Ga_SEC":
                    node.xmlRemoveNode()


    def __deleteCoalModelProperties(self, classMin, classMax, classesNumber):
        """
        Private method
        Delete properties for one coal
        """
        baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                     "Ga_DV1", "Ga_DV2", "Ga_HET"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                         "Ga_DV1", "Ga_DV2", "Ga_HET", "Ga_SEC"]
        #
        # Remove coal classes
        nodeList = self.node_coal.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList :
                nameNode =node['name']
                for baseName in baseNames:
                    for classe in range(classMin, classMax):
                        name = '%s%2.2i' % (baseName, classe+1)
                        if ( nameNode == name):
                            node.xmlRemoveNode()
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('property')
        delta = classMax - classMin
        if nodeList != None:
            for node in nodeList:
                try:
                    oldName = node['name']
                    if oldName[:-2] in baseNames :
                        oldNum = int(oldName[-2:])
                        if oldNum in range(classMax+1, classesNumber+1):
                            name = '%s%2.2i' % (oldName[:-2], oldNum-delta)
                            node['name'] = name
                            if node['label'] == oldName:
                                node['label'] = name

                except:
                    pass


    def __deleteCoalModelScalars(self, classMin, classMax, classesNumber, coalNumber, coalsNumber):
        """
        Private method
        Delete scalars for one coal
        """
        baseNames = ["NP_CP",  "XCH_CP", "XCK_CP", "ENT_CP"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP", "XWT_CP"]
        #
        # Remove coal classes
        nodeList = self.node_coal.xmlGetNodeList('scalar', 'name')
        if nodeList != None:
            for node in nodeList :
                nameNode = node['name']
                for baseName in baseNames:
                    for classe in range(classMin, classMax):
                        name = '%s%2.2i' % (baseName, classe+1)
                        if (nameNode == name):
                            node.xmlRemoveNode()
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('scalar')
        delta = classMax - classMin
        if nodeList != None:
            for node in nodeList:
                try:
                    oldName = node['name']
                    if oldName[:-2] in baseNames :
                        oldNum = int(oldName[-2:])
                        if oldNum in range(classMax+1, classesNumber+1):
                            name = '%s%2.2i' % (oldName[:-2], oldNum-delta)
                            node['name'] = name
                            if node['label'] == oldName:
                                node['label'] = name
                except:
                    pass
        #
        # Remove coal scalars
        baseNames = [ "Fr_MV1", "Fr_MV2"]
        nodeList = self.node_coal.xmlGetNodeList('scalar')
        if nodeList != None:
            for node in nodeList :
                nameNode = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, coalNumber+1)
                    if (nameNode == name):
                        node.xmlRemoveNode()
        #
        # Rename other coals
        nodeList = self.node_coal.xmlGetNodeList('scalar')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(coalNumber+1, coalsNumber+1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum-1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name


    def __updateCoalCombustionDensity(self, model):
        """
        Private method
        Update the coal combustion model markup from the XML document.
        """
        self.isInList(model, self.__coalCombustionModelsList())

        mdl = FluidCharacteristicsModel(self.case)

        if model == 'off':
            w = mdl.node_density.xmlGetDouble('initial_value')

            if w == None:
                v = mdl.defaultFluidCharacteristicsValues()['density']
                mdl.node_density.xmlSetData('initial_value',v)
                mdl.setPropertyMode('density', 'constant')
                mdl.node_density.xmlInitNode('listing_printing', status='off')
                mdl.node_density.xmlInitNode('postprocessing_recording', status='off')

                v = mdl.defaultFluidCharacteristicsValues()['thermal_conductivity']
                mdl.node_cond.xmlSetData('initial_value',v)
                mdl.setPropertyMode('thermal_conductivity', 'constant')
                mdl.node_cond.xmlInitNode('listing_printing', status='off')
                mdl.node_cond.xmlInitNode('postprocessing_recording', status='off')

        else:
            mdl.setPropertyMode('density', 'variable')
            mdl.node_density.xmlRemoveChild('initial_value')
            mdl.node_density.xmlRemoveChild('listing_printing')
            mdl.node_density.xmlRemoveChild('postprocessing_recording')
            mdl.node_cond.xmlRemoveNode()


    def __createCoalModelScalars(self , coalsNumber, coalClassesNumber, classesNumber):
        """
        Private method
        Create new scalars for one coal
        """
        # add new scalars
        baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP", "XWT_CP"]
        else:
            self.__deleteWetScalarsAndProperty()

        for baseName in baseNames:
            for classe in range(classesNumber - coalClassesNumber, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                self.setNewModelScalar(self.node_coal, name)

        baseNames = [ "Fr_MV1", "Fr_MV2"]
        for baseName in baseNames:
            name = '%s%2.2i' % (baseName, coalsNumber)
            self.setNewModelScalar(self.node_coal, name)

        NPE = NumericalParamEquatModel(self.case)
        for node in self.node_coal.xmlGetChildNodeList('scalar'):
            NPE.setBlendingFactor(node['label'], 0.)
            NPE.setScheme(node['label'], 'upwind')
            NPE.setFluxReconstruction(node['label'], 'off')


    def __createCoalModelProperties(self, coalsNumber, coalClassesNumber, classesNumber):
        """
        Private method
        Create new properties for one coal
        """
        # create new properties
        baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                     "Ga_DV1", "Ga_DV2", "Ga_HET"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                         "Ga_DV1", "Ga_DV2", "Ga_HET", "Ga_SEC"]
        for baseName in baseNames:
            for classe in range(classesNumber - coalClassesNumber, classesNumber):
                name = '%s%2.2i' % (baseName, classe+1)
                self.setNewProperty(self.node_coal, name)


    def __createClassModelProperties(self, classNum, classesNumber):
        """
        Private method
        Create class of model properties
        """
        baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                     "Ga_DV1", "Ga_DV2", "Ga_HET"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                         "Ga_DV1", "Ga_DV2", "Ga_HET", "Ga_SEC"]
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('property')
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
            self.setNewProperty(self.node_coal, name)


    def __createClassModelScalars(self, classNum, classesNumber):
        """
        Private method
        Create a new coal and associated scalars
        """
        baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP", "XWT_CP"]
        else:
            self.__deleteWetScalarsAndProperty()
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('scalar')
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
        #
        # create new scalars
        for i in range(len(baseNames)):
            name = '%s%2.2i' % (baseNames[i], classNum)
            self.setNewModelScalar(self.node_coal, name)


    def setCoalCombustionModel(self, model):
        """
        Update the coal combustion model markup from the XML document.
        """
        self.isInList(model, self.__coalCombustionModelsList())

        if model == 'off':
            for tag in ('scalar',
                        'property',
                        'reference_mass_molar',
                        'reference_temperature'):
                for node in self.node_coal.xmlGetNodeList(tag):
                    node.xmlRemoveNode()

            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "inlet":
                    Boundary("coal_inlet", zone.getLabel(), self.case).deleteCoals()

            ThermalRadiationModel(self.case).setRadiativeModel('off')
            ConjugateHeatTransferModel(self.case).setConjugateHeatTransferStatus('off')
            self.node_coal['model'] = 'off'

        else:
            self.node_gas['model']   = 'off'
            self.node_coal['model']  = model
            self.node_joule['model'] = 'off'
            self.node_therm['model'] = 'off'
            self.createModel()
            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "inlet":
                    b = Boundary("coal_inlet", zone.getLabel(), self.case)
                    b.getOxydantTemperature()
                    b.getOxydantNumber()

        self.__updateCoalCombustionDensity(model)


    def getCoalCombustionModel(self):
        """
        Return the current coal combustion model.
        """
        model = self.node_coal['model']

        if model not in self.__coalCombustionModelsList():
            model = self.defaultValues()['model']
            self.setCoalCombustionModel(model)
        else:
            self.__updateCoalCombustionDensity(model)

        return model


    def createCoalModelScalarsAndProperties(self, thermoChemistryModel ):
        """
        Create new scalars and new properties for one coal
        """
        coalsNumber = thermoChemistryModel.getCoals().getCoalNumber()
        coalClassesNumber = thermoChemistryModel.getCoals().getClassesNumberList()[coalsNumber - 1]
        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())

        # add new scalars and properties
        self.__createCoalModelScalars(coalsNumber, coalClassesNumber, classesNumber)
        self.__createCoalModelProperties(coalsNumber, coalClassesNumber, classesNumber)


    def createClassModelScalarsAndProperties(self, thermoChemistryModel, coalNumber):
        """
        Create class of model scalars and properties for one given coal
        """
        self.isInt(coalNumber)

        classNum = 0
        for coal in range(0, coalNumber):
            classNum += thermoChemistryModel.getCoals().getClassesNumberList()[coal]

        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())
        #
        # create new scalars
        self.__createClassModelScalars(classNum, classesNumber)
        self.__createClassModelProperties(classNum, classesNumber)


    def deleteCoalModelScalarsAndProperties(self, thermoChemistryModel, coalNumber):
        """
        Delete scalars and properties for one coal
        """
        self.isInt(coalNumber)

        classMin = 0
        for coal in range(0, coalNumber):
            classMin += thermoChemistryModel.getCoals().getClassesNumberList()[coal]
        classMax = classMin + thermoChemistryModel.getCoals().getClassesNumberList()[coalNumber]
        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())
        coalsNumber = thermoChemistryModel.getCoals().getCoalNumber()

        self.__deleteCoalModelScalars(classMin, classMax, classesNumber, coalNumber, coalsNumber)
        self.__deleteCoalModelProperties(classMin, classMax, classesNumber)


    def deleteClassModelScalars(self, thermoChemistryModel, coalNumber, classeNumber):
        """
        delete class of model scalars
        """
        self.isInt(coalNumber)
        self.isInt(classeNumber)

        classNum = 0
        if (coalNumber >= 1) :
            for coal in range(0, coalNumber - 1):
                classNum += thermoChemistryModel.getCoals().getClassesNumberList()[coal]
        classNum += classeNumber


        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())
        baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["NP_CP", "XCH_CP", "XCK_CP", "ENT_CP", "XWT_CP"]

        #
        # Remove coal classes
        nodeList = self.node_coal.xmlGetNodeList('scalar')
        if nodeList != None:
            for node in nodeList :
                nodeName = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, classNum+1)
                    if (nodeName == name):
                        node.xmlRemoveNode()
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('scalar')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum + 1, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum-1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name


    def deleteClassModelProperties(self, thermoChemistryModel, coalNumber, classeNumber):
        """
        delete class of model properties
        """
        self.isInt(coalNumber)
        self.isInt(classeNumber)

        classNum = 0
        for coal in range(0, coalNumber):
            classNum += thermoChemistryModel.getCoals().getClassesNumberList()[coal]
        classNum += classeNumber

        classesNumber = sum(thermoChemistryModel.getCoals().getClassesNumberList())
        baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                     "Ga_DV1", "Ga_DV2", "Ga_HET"]
        if self.getCoalCombustionModel() == 'coal_homo2':
            baseNames = ["Temp_CP", "Frm_CP", "Rho_CP", "Dia_CK", "Ga_DCH",
                         "Ga_DV1", "Ga_DV2", "Ga_HET", "Ga_SEC"]
        #
        # Remove coal classes
        nodeList = self.node_coal.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList :
                nodeName = node['name']
                for baseName in baseNames:
                    name = '%s%2.2i' % (baseName, classNum+1)
                    if (nodeName == name):
                        node.xmlRemoveNode()
        #
        # Rename other classes
        nodeList = self.node_coal.xmlGetNodeList('property')
        if nodeList != None:
            for node in nodeList:
                oldName = node['name']
                if oldName[:-2] in baseNames :
                    oldNum = int(oldName[-2:])
                    if oldNum in range(classNum + 1, classesNumber + 1):
                        name = '%s%2.2i' % (oldName[:-2], oldNum-1)
                        node['name'] = name
                        if node['label'] == oldName:
                            node['label'] = name

#-------------------------------------------------------------------------------
# Coal combustion test case
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

        doc = '''<pulverized_coal model="coal_homo">
                    <scalar label="Enthalpy" name="Enthalpy" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="NP_CP01" name="NP_CP01" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP01" name="XCH_CP01" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP01" name="XCK_CP01" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP01" name="ENT_CP01" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV101" name="Fr_MV101" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV201" name="Fr_MV201" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_HET" name="Fr_HET" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Var_AIR" name="Var_AIR" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_GAZ" name="Temp_GAZ"/>
                    <property label="ROM_GAZ" name="ROM_GAZ"/>
                    <property label="YM_CHx1m" name="YM_CHx1m"/>
                    <property label="YM_CHx2m" name="YM_CHx2m"/>
                    <property label="YM_CO" name="YM_CO"/>
                    <property label="YM_O2" name="YM_O2"/>
                    <property label="YM_CO2" name="YM_CO2"/>
                    <property label="YM_H2O" name="YM_H2O"/>
                    <property label="YM_N2" name="YM_N2"/>
                    <property label="XM" name="XM"/>
                    <property label="Temp_CP01" name="Temp_CP01"/>
                    <property label="Frm_CP01" name="Frm_CP01"/>
                    <property label="Rho_CP01" name="Rho_CP01"/>
                    <property label="Dia_CK01" name="Dia_CK01"/>
                    <property label="Ga_DCH01" name="Ga_DCH01"/>
                    <property label="Ga_DV101" name="Ga_DV101"/>
                    <property label="Ga_DV201" name="Ga_DV201"/>
                    <property label="Ga_HET01" name="Ga_HET01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
            </pulverized_coal>'''

        assert model.node_coal == self.xmlNodeFromString(doc),\
           'Could not set a model of combustion'
        assert model.getCoalCombustionModel() == 'coal_homo',\
           'Could not get a model of combustion'


    def checkCreateCoalModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        created new scalars and properties for one new coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        from CoalThermoChemistry import CoalThermoChemistryModel, Coal
        coalThermoChModel = CoalThermoChemistryModel("dp_FCP",self.case)
        coalThermoChModel.getCoals().addCoal(Coal())
        del CoalThermoChemistryModel

        model.createCoalModelScalarsAndProperties(coalThermoChModel)

        doc = '''<pulverized_coal model="coal_homo">
                    <scalar label="Enthalpy" name="Enthalpy" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="NP_CP01" name="NP_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCH_CP01" name="XCH_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCK_CP01" name="XCK_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="ENT_CP01" name="ENT_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV101" name="Fr_MV101" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV201" name="Fr_MV201" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_HET" name="Fr_HET" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Var_AIR" name="Var_AIR" type="model"><flux_reconstruction status="off"/></scalar>
                    <property label="Temp_GAZ" name="Temp_GAZ"/>
                    <property label="ROM_GAZ" name="ROM_GAZ"/>
                    <property label="YM_CHx1m" name="YM_CHx1m"/>
                    <property label="YM_CHx2m" name="YM_CHx2m"/>
                    <property label="YM_CO" name="YM_CO"/>
                    <property label="YM_O2" name="YM_O2"/>
                    <property label="YM_CO2" name="YM_CO2"/>
                    <property label="YM_H2O" name="YM_H2O"/>
                    <property label="YM_N2" name="YM_N2"/>
                    <property label="XM" name="XM"/>
                    <property label="Temp_CP01" name="Temp_CP01"/>
                    <property label="Frm_CP01" name="Frm_CP01"/>
                    <property label="Rho_CP01" name="Rho_CP01"/>
                    <property label="Dia_CK01" name="Dia_CK01"/>
                    <property label="Ga_DCH01" name="Ga_DCH01"/>
                    <property label="Ga_DV101" name="Ga_DV101"/>
                    <property label="Ga_DV201" name="Ga_DV201"/>
                    <property label="Ga_HET01" name="Ga_HET01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <scalar label="NP_CP02" name="NP_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP02" name="XCH_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP02" name="XCK_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP02" name="ENT_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV102" name="Fr_MV102" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV202" name="Fr_MV202" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP02" name="Temp_CP02"/>
                    <property label="Frm_CP02" name="Frm_CP02"/>
                    <property label="Rho_CP02" name="Rho_CP02"/>
                    <property label="Dia_CK02" name="Dia_CK02"/>
                    <property label="Ga_DCH02" name="Ga_DCH02"/>
                    <property label="Ga_DV102" name="Ga_DV102"/>
                    <property label="Ga_DV202" name="Ga_DV202"/>
                    <property label="Ga_HET02" name="Ga_HET02"/>
            </pulverized_coal>'''

        assert model.node_coal == self.xmlNodeFromString(doc),\
            'Could not create newscalars and properties for new coal'


    def checkCreateClassModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        created a new class of scalars and properties for one coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        from CoalThermoChemistry import CoalThermoChemistryModel, Coal
        coalThermoChModel = CoalThermoChemistryModel("dp_FCP",self.case)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        del CoalThermoChemistryModel

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)

        doc = '''<pulverized_coal model="coal_homo">
                    <scalar label="Enthalpy" name="Enthalpy" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="NP_CP01" name="NP_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCH_CP01" name="XCH_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCK_CP01" name="XCK_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="ENT_CP01" name="ENT_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV101" name="Fr_MV101" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV201" name="Fr_MV201" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_HET" name="Fr_HET" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Var_AIR" name="Var_AIR" type="model"><flux_reconstruction status="off"/></scalar>
                    <property label="Temp_GAZ" name="Temp_GAZ"/>
                    <property label="ROM_GAZ" name="ROM_GAZ"/>
                    <property label="YM_CHx1m" name="YM_CHx1m"/>
                    <property label="YM_CHx2m" name="YM_CHx2m"/>
                    <property label="YM_CO" name="YM_CO"/>
                    <property label="YM_O2" name="YM_O2"/>
                    <property label="YM_CO2" name="YM_CO2"/>
                    <property label="YM_H2O" name="YM_H2O"/>
                    <property label="YM_N2" name="YM_N2"/>
                    <property label="XM" name="XM"/>
                    <property label="Temp_CP01" name="Temp_CP01"/>
                    <property label="Frm_CP01" name="Frm_CP01"/>
                    <property label="Rho_CP01" name="Rho_CP01"/>
                    <property label="Dia_CK01" name="Dia_CK01"/>
                    <property label="Ga_DCH01" name="Ga_DCH01"/>
                    <property label="Ga_DV101" name="Ga_DV101"/>
                    <property label="Ga_DV201" name="Ga_DV201"/>
                    <property label="Ga_HET01" name="Ga_HET01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <scalar label="NP_CP03" name="NP_CP03" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP03" name="XCH_CP03" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP03" name="XCK_CP03" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP03" name="ENT_CP03" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV102" name="Fr_MV102" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV202" name="Fr_MV202" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP02" name="Temp_CP02"/>
                    <property label="Frm_CP02" name="Frm_CP02"/>
                    <property label="Rho_CP02" name="Rho_CP02"/>
                    <property label="Dia_CK02" name="Dia_CK02"/>
                    <property label="Ga_DCH02" name="Ga_DCH02"/>
                    <property label="Ga_DV102" name="Ga_DV102"/>
                    <property label="Ga_DV202" name="Ga_DV202"/>
                    <property label="Ga_HET02" name="Ga_HET02"/>
                    <scalar label="NP_CP04" name="NP_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP04" name="XCH_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP04" name="XCK_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP04" name="ENT_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV103" name="Fr_MV103" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV203" name="Fr_MV203" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP03" name="Temp_CP03"/>
                    <property label="Frm_CP03" name="Frm_CP03"/>
                    <property label="Rho_CP03" name="Rho_CP03"/>
                    <property label="Dia_CK03" name="Dia_CK03"/>
                    <property label="Ga_DCH03" name="Ga_DCH03"/>
                    <property label="Ga_DV103" name="Ga_DV103"/>
                    <property label="Ga_DV203" name="Ga_DV203"/>
                    <property label="Ga_HET03" name="Ga_HET03"/>
                    <scalar label="NP_CP02" name="NP_CP02" type="model"/>
                    <scalar label="XCH_CP02" name="XCH_CP02" type="model"/>
                    <scalar label="XCK_CP02" name="XCK_CP02" type="model"/>
                    <scalar label="ENT_CP02" name="ENT_CP02" type="model"/>
                    <property label="Temp_CP04" name="Temp_CP04"/>
                    <property label="Frm_CP04" name="Frm_CP04"/>
                    <property label="Rho_CP04" name="Rho_CP04"/>
                    <property label="Dia_CK04" name="Dia_CK04"/>
                    <property label="Ga_DCH04" name="Ga_DCH04"/>
                    <property label="Ga_DV104" name="Ga_DV104"/>
                    <property label="Ga_DV204" name="Ga_DV204"/>
                    <property label="Ga_HET04" name="Ga_HET04"/>
                </pulverized_coal>'''

        assert model.node_coal == self.xmlNodeFromString(doc),\
           'Could not create a new class of scalars and properties for one coal'


    def checkDeleteCoalModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        deleted class of scalars and properties for one given coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        from CoalThermoChemistry import CoalThermoChemistryModel, Coal
        coalThermoChModel = CoalThermoChemistryModel("dp_FCP",self.case)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        del CoalThermoChemistryModel

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)
        model.deleteCoalModelScalarsAndProperties(coalThermoChModel, 1)

        doc = '''<pulverized_coal model="coal_homo">
                    <scalar label="Enthalpy" name="Enthalpy" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="NP_CP01" name="NP_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCH_CP01" name="XCH_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCK_CP01" name="XCK_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="ENT_CP01" name="ENT_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV101" name="Fr_MV101" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV201" name="Fr_MV201" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_HET" name="Fr_HET" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Var_AIR" name="Var_AIR" type="model"><flux_reconstruction status="off"/></scalar>
                    <property label="Temp_GAZ" name="Temp_GAZ"/>
                    <property label="ROM_GAZ" name="ROM_GAZ"/>
                    <property label="YM_CHx1m" name="YM_CHx1m"/>
                    <property label="YM_CHx2m" name="YM_CHx2m"/>
                    <property label="YM_CO" name="YM_CO"/>
                    <property label="YM_O2" name="YM_O2"/>
                    <property label="YM_CO2" name="YM_CO2"/>
                    <property label="YM_H2O" name="YM_H2O"/>
                    <property label="YM_N2" name="YM_N2"/>
                    <property label="XM" name="XM"/>
                    <property label="Temp_CP01" name="Temp_CP01"/>
                    <property label="Frm_CP01" name="Frm_CP01"/>
                    <property label="Rho_CP01" name="Rho_CP01"/>
                    <property label="Dia_CK01" name="Dia_CK01"/>
                    <property label="Ga_DCH01" name="Ga_DCH01"/>
                    <property label="Ga_DV101" name="Ga_DV101"/>
                    <property label="Ga_DV201" name="Ga_DV201"/>
                    <property label="Ga_HET01" name="Ga_HET01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <scalar label="NP_CP02" name="NP_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP02" name="XCH_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP02" name="XCK_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP02" name="ENT_CP02" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP02" name="Temp_CP02"/>
                    <property label="Frm_CP02" name="Frm_CP02"/>
                    <property label="Rho_CP02" name="Rho_CP02"/>
                    <property label="Dia_CK02" name="Dia_CK02"/>
                    <property label="Ga_DCH02" name="Ga_DCH02"/>
                    <property label="Ga_DV102" name="Ga_DV102"/>
                    <property label="Ga_DV202" name="Ga_DV202"/>
                    <property label="Ga_HET02" name="Ga_HET02"/>
                    <scalar label="NP_CP04" name="NP_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP04" name="XCH_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP04" name="XCK_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP04" name="ENT_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV102" name="Fr_MV102" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV202" name="Fr_MV202" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP04" name="Temp_CP04"/>
                    <property label="Frm_CP04" name="Frm_CP04"/>
                    <property label="Rho_CP04" name="Rho_CP04"/>
                    <property label="Dia_CK04" name="Dia_CK04"/>
                    <property label="Ga_DCH04" name="Ga_DCH04"/>
                    <property label="Ga_DV104" name="Ga_DV104"/>
                    <property label="Ga_DV204" name="Ga_DV204"/>
                    <property label="Ga_HET04" name="Ga_HET04"/>
                </pulverized_coal>'''

        assert model.node_coal == self.xmlNodeFromString(doc),\
           "Could not delete one coal's scalars and properties"


    def checkDeleteClassModelScalarsAndProperties(self):
        """
        Check whether the CoalCombustionModel class could be
        deleted class of scalars and properties for one given coal
        """
        model = CoalCombustionModel(self.case)
        model.setCoalCombustionModel('coal_homo')

        from CoalThermoChemistry import CoalThermoChemistryModel, Coal
        coalThermoChModel = CoalThermoChemistryModel("dp_FCP",self.case)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        coalThermoChModel.getCoals().addCoal(Coal())
        model.createCoalModelScalarsAndProperties(coalThermoChModel)
        del CoalThermoChemistryModel

        model.createClassModelScalarsAndProperties(coalThermoChModel, 2)
        model.deleteClassModelScalars(coalThermoChModel, 2, 1)

        doc = '''<pulverized_coal model="coal_homo">
                    <scalar label="Enthalpy" name="Enthalpy" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="NP_CP01" name="NP_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCH_CP01" name="XCH_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="XCK_CP01" name="XCK_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="ENT_CP01" name="ENT_CP01" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV101" name="Fr_MV101" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_MV201" name="Fr_MV201" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Fr_HET" name="Fr_HET" type="model"><flux_reconstruction status="off"/></scalar>
                    <scalar label="Var_AIR" name="Var_AIR" type="model"><flux_reconstruction status="off"/></scalar>
                    <property label="Temp_GAZ" name="Temp_GAZ"/>
                    <property label="ROM_GAZ" name="ROM_GAZ"/>
                    <property label="YM_CHx1m" name="YM_CHx1m"/>
                    <property label="YM_CHx2m" name="YM_CHx2m"/>
                    <property label="YM_CO" name="YM_CO"/>
                    <property label="YM_O2" name="YM_O2"/>
                    <property label="YM_CO2" name="YM_CO2"/>
                    <property label="YM_H2O" name="YM_H2O"/>
                    <property label="YM_N2" name="YM_N2"/>
                    <property label="XM" name="XM"/>
                    <property label="Temp_CP01" name="Temp_CP01"/>
                    <property label="Frm_CP01" name="Frm_CP01"/>
                    <property label="Rho_CP01" name="Rho_CP01"/>
                    <property label="Dia_CK01" name="Dia_CK01"/>
                    <property label="Ga_DCH01" name="Ga_DCH01"/>
                    <property label="Ga_DV101" name="Ga_DV101"/>
                    <property label="Ga_DV201" name="Ga_DV201"/>
                    <property label="Ga_HET01" name="Ga_HET01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <scalar label="Fr_MV102" name="Fr_MV102" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV202" name="Fr_MV202" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP03" name="Temp_CP03"/>
                    <property label="Frm_CP03" name="Frm_CP03"/>
                    <property label="Rho_CP03" name="Rho_CP03"/>
                    <property label="Dia_CK03" name="Dia_CK03"/>
                    <property label="Ga_DCH03" name="Ga_DCH03"/>
                    <property label="Ga_DV103" name="Ga_DV103"/>
                    <property label="Ga_DV203" name="Ga_DV203"/>
                    <property label="Ga_HET03" name="Ga_HET03"/>
                    <scalar label="NP_CP04" name="NP_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCH_CP04" name="XCH_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="XCK_CP04" name="XCK_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="ENT_CP04" name="ENT_CP04" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV103" name="Fr_MV103" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <scalar label="Fr_MV203" name="Fr_MV203" type="model">
                            <flux_reconstruction status="off"/>
                    </scalar>
                    <property label="Temp_CP04" name="Temp_CP04"/>
                    <property label="Frm_CP04" name="Frm_CP04"/>
                    <property label="Rho_CP04" name="Rho_CP04"/>
                    <property label="Dia_CK04" name="Dia_CK04"/>
                    <property label="Ga_DCH04" name="Ga_DCH04"/>
                    <property label="Ga_DV104" name="Ga_DV104"/>
                    <property label="Ga_DV204" name="Ga_DV204"/>
                    <property label="Ga_HET04" name="Ga_HET04"/>
                    <scalar label="NP_CP02" name="NP_CP02" type="model"/>
                    <scalar label="XCH_CP02" name="XCH_CP02" type="model"/>
                    <scalar label="XCK_CP02" name="XCK_CP02" type="model"/>
                    <scalar label="ENT_CP02" name="ENT_CP02" type="model"/>
                    <property label="Temp_CP02" name="Temp_CP02"/>
                    <property label="Frm_CP02" name="Frm_CP02"/>
                    <property label="Rho_CP02" name="Rho_CP02"/>
                    <property label="Dia_CK02" name="Dia_CK02"/>
                    <property label="Ga_DCH02" name="Ga_DCH02"/>
                    <property label="Ga_DV102" name="Ga_DV102"/>
                    <property label="Ga_DV202" name="Ga_DV202"/>
                    <property label="Ga_HET02" name="Ga_HET02"/>
            </pulverized_coal>'''

        assert model.node_coal == self.xmlNodeFromString(doc),\
           "Could not delete one class of scalars and properties for one given coal"


def suite():
    testSuite = unittest.makeSuite(CoalCombustionModelTestCase, "check")
    return testSuite


def runTest():
    print "CoalCombustionModelTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
