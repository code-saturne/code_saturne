# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the atmospheric flows modelling management.

This module contains the following classes and function:
- AtmosphericFlowsModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLmodel     import  ModelTest
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.model.NumericalParamGlobalModel import NumericalParamGlobalModel

#-------------------------------------------------------------------------------
# Atmospheric flows model class
#-------------------------------------------------------------------------------

class AtmosphericFlowsModel(Model):
    """
    Model for atmospheric flows
    """
    off             = 'off'
    constant        = 'constant'
    dry             = 'dry'
    humid           = 'humid'
    read_meteo_data = 'read_meteo_data'
    model           = 'model'
    status          = 'status'


    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.__fluidProp = FluidCharacteristicsModel(self.case)

        models = case.xmlGetNode('thermophysical_models')
        self.__node_atmos  = models.xmlInitChildNode('atmospheric_flows')
        self.__atmosphericModel = (AtmosphericFlowsModel.off,
                                   AtmosphericFlowsModel.constant,
                                   AtmosphericFlowsModel.dry,
                                   AtmosphericFlowsModel.humid)
        self.__default = {}
        self.__default[self.model] = AtmosphericFlowsModel.off
        self.__default[self.read_meteo_data] = AtmosphericFlowsModel.off
        self.__default['meteo_data'] = "meteo"


    @Variables.undoLocal
    def setAtmosphericFlowsModel(self, model):
        """
        Update the atmospheric flows model markup from the XML document.
        """
        self.isInList(model, self.__atmosphericModel)
        self.__node_atmos[self.model] = model
        self.__updateScalarAndProperty()
        if (model == "humid" or model == "dry"):
            NumericalParamGlobalModel(self.case).setHydrostaticPressure("on")
            if (model == "dry"):
                ThermalScalarModel(self.case).setThermalModel('potential_temperature')
            else:
                ThermalScalarModel(self.case).setThermalModel('liquid_potential_temperature')
        else:
            ThermalScalarModel(self.case).setThermalModel('off')
            NumericalParamGlobalModel(self.case).setHydrostaticPressure("off")


    @Variables.noUndo
    def getAtmosphericFlowsModel(self):
        """
        Return the current atmospherics flows model.
        """
        model = self.__node_atmos[self.model]
        if model not in self.__atmosphericModel:
            model = self.__default[self.model]
            self.setAtmosphericFlowsModel(model)
        return model


    @Variables.noUndo
    def getMeteoDataStatus(self):
        """
        Return if reading meteo data status is 'on' or 'off'.
        """
        node = self.__node_atmos.xmlInitChildNode(self.read_meteo_data)
        if not node[self.status]:
            status = self.__default[self.read_meteo_data]
            self.setMeteoDataStatus(status)
        return node[self.status]


    @Variables.undoLocal
    def setMeteoDataStatus(self, status):
        """
        Set meteo data status to 'on' / 'off'.
        """
        self.isOnOff(status)
        self.__node_atmos.xmlInitChildNode(self.read_meteo_data)[self.status] = status

        if status == 'off':
            for tag in ['read_meteo_data', 'meteo_automatic']:
                for node in self.case.xmlGetNodeList(tag):
                    node['status'] = "off"


    @Variables.noUndo
    def getMeteoDataFileName(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetString('meteo_data')
        if f == None:
            f = self.__default['meteo_data']
            self.setMeteoDataFile(f)
        return f


    @Variables.undoLocal
    def setMeteoDataFileName(self, tag):
        """
        Set the name of the meteo data file.
        """
        self.__node_atmos.xmlSetData('meteo_data', tag)


    def __updateScalarAndProperty(self):
        """
        Update scalars and properties depending on model
        """
        model = self.getAtmosphericFlowsModel()
        node = self.__node_atmos

        # Update only if getMeteoDataStatus is not off
        if model != AtmosphericFlowsModel.off:

            if model == AtmosphericFlowsModel.dry:
                self.__removeScalar(node, 'ym_water')
                self.__removeScalar(node, 'number_of_droplets')
                self.__removeProperty(node, 'liquid_water')
                self.__setProperty(node, 'RealTemp', 'real_temperature')
                if self.__fluidProp.getPropertyMode('density') == 'constant':
                    self.__fluidProp.setPropertyMode('density', 'predefined_law')

            elif model == AtmosphericFlowsModel.humid:
                self.__setScalar(node, 'TotWater', 'ym_water', 'model')
                self.__setScalar(node, 'TotDrop', 'number_of_droplets', 'model')
                self.__setProperty(node, 'RealTemp', 'real_temperature')
                self.__setProperty(node, 'LiqWater', 'liquid_water')
                if self.__fluidProp.getPropertyMode('density') == 'constant':
                    self.__fluidProp.setPropertyMode('density', 'predefined_law')

            elif model == AtmosphericFlowsModel.constant:
                self.__removeScalar(node, 'ym_water')
                self.__removeScalar(node, 'number_of_droplets')
                self.__removeProperty(node, 'liquid_water')
                FluidCharacteristicsModel(self.case).setPropertyMode('density', 'constant')

        else:
            self.__removeScalar(node, 'ym_water')
            self.__removeScalar(node, 'number_of_droplets')
            self.__removeProperty(node, 'liquid_water')


    def atmosphericFlowsNode(self):
        """
        Get the atmospheric flows node
        """
        return self.__node_atmos


    def __setScalar(self, parentNode, labelStr, nameStr, typeStr ):
        """
        Create xml scalar
        """
        scalar = parentNode.xmlInitChildNode('variable', name = nameStr)
        scalar['label'] = labelStr
        scalar['type']  = typeStr


    def __setProperty(self, parentNode, labelStr, nameStr):
        """
        Create xml property
        """
        prop = parentNode.xmlInitChildNode('property', name = nameStr)
        prop['label']  = labelStr


    def __removeScalar(self, parentNode, nameStr):
        """
        Delete scalar
        """
        parentNode.xmlRemoveChild('variable', name = nameStr)


    def __removeProperty(self, parentNode, nameStr):
        """
        Delete property
        """
        parentNode.xmlRemoveChild('property', name = nameStr)

#-------------------------------------------------------------------------------
# AtmosphericFlowsModel test case
#-------------------------------------------------------------------------------

class AtmosphericFlowsTestCase(ModelTest):
    """
    Test case for AtmosphericFlows
    """
    def checkAtmosphericFlowsInstantiation(self):
        """
        Check whether the AtmosphericFlowsModel class could be instantiated
        """
        model = None
        model = AtmosphericFlowsModel(self.case)
        assert model != None, 'Could not instantiate '


    def checkGetandSetAtmosphericFlowsModel(self):
        """Check whether the AtmosphericFlowsModel class could be set and get the model"""
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.dry)

        doc = """<atmospheric_flows model="dry">
                    <read_meteo_data status="off"/>
                 </atmospheric_flows>"""
        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set atmospheric flows model'
        assert mdl.getAtmosphericFlowsModel() == AtmosphericFlowsModel.dry, \
            'Could not get atmospheric flows model'


    def checkGetandSetMeteoDataStatus(self):
        """Check whether the AtmosphericFlowsModel class could be set and get the meteo data status"""
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.constant)
        mdl.setMeteoDataStatus('on')

        doc = """<atmospheric_flows model="constant">
                    <read_meteo_data status="on"/>
                 </atmospheric_flows>"""

        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set meteo data status'
        assert mdl.getMeteoDataStatus() == 'on', \
            'Could not get meteo data status'


    def checkDryModel(self):
        """
        Check whether the AtmosphericFlowsModel class could set the correct
        properties and scalar for dry model
        """
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.dry)
        mdl.setMeteoDataStatus('on')

        doc = """<atmospheric_flows model="dry">
                    <read_meteo_data status="on">
                        <property label="Real temp" name="real_temperature"/>
                    </read_meteo_data>
                </atmospheric_flows>
                """

        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set scalars and properties for dry model'


    def checkHumidModel(self):
        """
        Check whether the AtmosphericFlowsModel class could set the correct
        properties and scalar for humid model
        """
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.humid)
        mdl.setMeteoDataStatus('on')

        doc = """<atmospheric_flows model="humid">
                    <read_meteo_data status="on">
                        <variable label="total water" name="ym_water" type="model"/>
                        <variable label="number of droplets" name="number_of_droplets" type="model"/>
                        <property label="Real temp" name="real_temperature"/>
                        <property label="Liquid water" name="liquid_water"/>
                    </read_meteo_data>
                </atmospheric_flows>"""

        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set scalars and properties for humid model'


def suite():
    """
    Test Suite for AtmosphericFlows
    """
    testSuite = unittest.makeSuite(AtmosphericFlowsTestCase, "check")
    return testSuite


def runTest():
    """
    run test
    """
    print("AtmosphericFlowsTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
