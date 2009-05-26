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

from Base.XMLvariables import Model
from Base.XMLmodel     import  ModelTest

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
        self.__case = case

        models = case.xmlGetNode('thermophysical_models')
        self.__node_atmos  = models.xmlInitChildNode('atmospheric_flows')
        self.__atmosphericModel = (AtmosphericFlowsModel.off,
                                   AtmosphericFlowsModel.constant,
                                   AtmosphericFlowsModel.dry,
                                   AtmosphericFlowsModel.humid)
        self.__default = {}
        self.__default[self.model] = AtmosphericFlowsModel.off
        self.__default[self.read_meteo_data] = AtmosphericFlowsModel.off


    def setAtmosphericFlowsModel(self, model):
        """
        Update the atmospheric flows model markup from the XML document.
        """
        self.isInList(model, self.__atmosphericModel)
        self.__node_atmos[self.model]  = model
        self.__updateScalarAndProperty()


    def getAtmosphericFlowsModel(self):
        """
        Return the current atmospherics flows model.
        """
        model = self.__node_atmos[self.model]
        if model not in self.__atmosphericModel:
            model = self.__default[self.model]
            self.setAtmosphericFlowsModel(model)
        return model


    def getMeteoDataStatus(self):
        """
        Return if meteoData status is 'on' or 'off'
        """
        node = self.__node_atmos.xmlInitChildNode(self.read_meteo_data)
        if not node[self.status]:
            status = self.__default[self.read_meteo_data]
            self.setMeteoDataStatus(status)
        return node[self.status]


    def setMeteoDataStatus(self, status):
        """
        Set meteo data status to 'on' / 'off'
        """
        self.isOnOff(status)
        node = self.__node_atmos.xmlInitChildNode(self.read_meteo_data)
        if status != node[self.status]:
            node[self.status] = status
            self.__updateScalarAndProperty()


    def __updateScalarAndProperty(self):
        """
        Update scalars and properties depending on model
        """
        # Update only if getMeteoDataStatus is not off

        if self.getMeteoDataStatus() != AtmosphericFlowsModel.off:

            model = self.getAtmosphericFlowsModel()
            node = self.__node_atmos

            if model == AtmosphericFlowsModel.dry:
                self.__removeScalar(node, 'liquid_potential_temperature')
                self.__removeScalar(node, 'total_water')
                self.__removeScalar(node, 'number_of_droplets')
                self.__removeProperty(node, 'liquid_water')

                self.__setScalar(node, 'Potential temp', 
                                 'potential_temperature', 'model')
                self.__setProperty(node, 'Real temp', 'real_temperature' )

            elif model == AtmosphericFlowsModel.humid:
                self.__removeScalar(node, 'potential_temperature')

                self.__setScalar(node, 'Liq potential temp', 
                                 'liquid_potential_temperature', 'model')
                self.__setScalar(node, 'total water', 'total_water', 'model')
                self.__setScalar(node, 'number of droplets', 
                                'number_of_droplets', 'model')
                self.__setProperty(node, 'Real temp', 'real_temperature')
                self.__setProperty(node, 'Liquid water', 'liquid_water')


    def atmosphericFlowsNode(self):
        """
        Get the atmospheric flows node
        """
        return self.__node_atmos


    def __setScalar(self, parentNode, labelStr, nameStr, typeStr ):
        """
        Create xml scalar
        """
        scalar = parentNode.xmlInitChildNode('scalar', label = labelStr)
        scalar['name']  = nameStr
        scalar['type']  = typeStr


    def __setProperty(self, parentNode, labelStr, nameStr):
        """
        Create xml property
        """
        prop = parentNode.xmlInitChildNode('property', label = labelStr)
        prop['name']  = nameStr


    def __removeScalar(self, parentNode, nameStr):
        """
        Delete scalar
        """
        scalar = parentNode.xmlRemoveChild('scalar', name = nameStr)


    def __removeProperty(self, parentNode, nameStr):
        """
        Delete property
        """
        prop = parentNode.xmlRemoveChild('property', name = nameStr)

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
                        <scalar label="Potential temp" name="potential_temperature" type="model"/>
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
                        <scalar label="Liq potential temp" name="liquid_potential_temperature" type="model"/>
                        <scalar label="total water" name="total_water" type="model"/>
                        <scalar label="number of droplets" name="number_of_droplets" type="model"/>
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
    print "AtmosphericFlowsTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
