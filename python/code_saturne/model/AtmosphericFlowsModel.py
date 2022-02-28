# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
from datetime import datetime

#-------------------------------------------------------------------------------
# Atmospheric flows model class
#-------------------------------------------------------------------------------

class AtmosphericFlowsModel(Model):
    """
    Model for atmospheric flows
    """
    off                 = 'off'
    constant            = 'constant'
    dry                 = 'dry'
    humid               = 'humid'
    read_meteo_data     = 'read_meteo_data'
    large_scale_meteo   = 'large_scale_meteo'
    act_chemistry       = 'activate_chemistry'
    model               = 'model'
    status              = 'status'

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
        self.__default[self.large_scale_meteo] = AtmosphericFlowsModel.off
        self.__default[self.act_chemistry] = AtmosphericFlowsModel.off
        self.__default['meteo_data'] = "meteo"
        self.__default['longitude'] = 4.39
        self.__default['latitude'] = 45.44
        self.__default['domain_orientation'] = 0
        self.__default['wind_direction'] = 0
        self.__default['meteo_z0'] = 0.1
        self.__default['meteo_psea'] = 101325.0
        self.__default['meteo_t0'] = 284.15
        self.__default['meteo_angle'] = 0.
        self.__default['meteo_dlmo'] = 0.
        self.__default['meteo_zref'] = 10.
        self.__default['meteo_uref'] = 5.
        self.__default['meteo_ustar'] = -1.

    @Variables.undoLocal
    def setAtmosphericFlowsModel(self, model):
        """
        Update the atmospheric flows  markup from the XML document.
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

    #--------------------------------------------------------------------------
    #-----------------Large scale meteo block----------------------------------
    @Variables.noUndo
    def getLargeScaleMeteoStatus(self):
        """
        Return if reading meteo data status is 'on' or 'off'.
        """
        node = self.__node_atmos.xmlInitChildNode(self.large_scale_meteo)
        if not node[self.status]:
            status = self.__default[self.large_scale_meteo]
            self.setLargeScaleMeteoStatus(status)
        return node[self.status]
    @Variables.undoLocal
    def setLargeScaleMeteoStatus(self, status):
        """
        Set meteo data status to 'on' / 'off'.
        """
        self.isOnOff(status)
        self.__node_atmos.xmlInitChildNode(self.large_scale_meteo)[self.status] = status

        if status == 'off':
            for tag in ['large_scale_meteo']:
                for node in self.case.xmlGetNodeList(tag):
                    node['status'] = "off"
        elif status == 'on':
            self.setMeteoDataStatus('off');

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getLongitude(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetDouble('longitude')
        if f is None:
            f = self.__default['longitude']
            self.setLongitude(f)
        return f
    @Variables.undoLocal
    def setLongitude(self, tag):
        self.__node_atmos.xmlSetData('longitude', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getLatitude(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetDouble('latitude')
        if f is None:
            f = self.__default['latitude']
            self.setLatitude(f)
        return f
    @Variables.undoLocal
    def setLatitude(self, tag):
        self.__node_atmos.xmlSetData('latitude', tag)


    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getDomainOrientation(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetInt('domain_orientation')
        if f is None:
            f = self.__default['domain_orientation']
            self.setDomainOrientation(f)
        return f
    @Variables.undoLocal
    def setDomainOrientation(self, tag):
        """
        Set the name of the meteo data file.
        """
        self.__node_atmos.xmlSetData('domain_orientation', tag)

   #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoZ0(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetDouble('meteo_z0')
        if f is None:
            f = self.__default['meteo_z0']
            self.setMeteoZ0(f)
        return f
    @Variables.undoLocal
    def setMeteoZ0(self, tag):
        self.__node_atmos.xmlSetData('meteo_z0', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoPsea(self):
        f = self.__node_atmos.xmlGetDouble('meteo_psea')
        if f is None:
            f = self.__default['meteo_psea']
            self.setMeteoPsea(f)
        return f
    @Variables.undoLocal
    def setMeteoPsea(self, tag):
        self.__node_atmos.xmlSetData('meteo_psea', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoT0(self):
        f = self.__node_atmos.xmlGetDouble('meteo_t0')
        if f is None:
            f = self.__default['meteo_t0']
            self.setMeteoT0(f)
        return f
    @Variables.undoLocal
    def setMeteoT0(self, tag):
        self.__node_atmos.xmlSetData('meteo_t0', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getWindDir(self):
        f = self.__node_atmos.xmlGetInt('wind_direction')
        if f is None:
            f = self.__default['wind_direction']
            self.setWindDir(f)
        return f
    @Variables.undoLocal
    def setWindDir(self, tag):
        self.__node_atmos.xmlSetData('wind_direction', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoUref(self):
        f = self.__node_atmos.xmlGetDouble('meteo_uref')
        if f is None:
            f = self.__default['meteo_uref']
            self.setMeteoUref(f)
        return f
    @Variables.undoLocal
    def setMeteoUref(self, tag):
        self.__node_atmos.xmlSetData('meteo_uref', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoUstar(self):
        f = self.__node_atmos.xmlGetDouble('meteo_ustar')
        if f is None:
            f = self.__default['meteo_ustar']
            self.setMeteoUstar(f)
        return f
    @Variables.undoLocal
    def setMeteoUstar(self, tag):
        self.__node_atmos.xmlSetData('meteo_ustar', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoDlmo(self):
        f = self.__node_atmos.xmlGetDouble('meteo_dlmo')
        if f is None:
            f = self.__default['meteo_dlmo']
            self.setMeteoDlmo(f)
        return f
    @Variables.undoLocal
    def setMeteoDlmo(self, tag):
        self.__node_atmos.xmlSetData('meteo_dlmo', tag)

    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getMeteoZref(self):
        f = self.__node_atmos.xmlGetDouble('meteo_zref')
        if f is None:
            f = self.__default['meteo_zref']
            self.setMeteoZref(f)
        return f
    @Variables.undoLocal
    def setMeteoZref(self, tag):
        self.__node_atmos.xmlSetData('meteo_zref', tag)


    #-------------------------------------------------------------------------
    @Variables.noUndo
    def getStartTime(self):
        startYear = self.__node_atmos.xmlGetInt('start_year')
        startDay = self.__node_atmos.xmlGetInt('start_day')
        startHour = self.__node_atmos.xmlGetInt('start_hour')
        startMin = self.__node_atmos.xmlGetInt('start_min')
        startSec = self.__node_atmos.xmlGetInt('start_sec')
        dateTime = datetime.now()
        if (startYear is None) or (startDay is None) or (startHour is None) \
            or (startMin is None) or (startSec is None):
            self.setStartTime(dateTime)
        else:
            #FIXME after startDay.rjust(3,'0')
            dateTimeStr = str(startYear)+'-'+str(startDay)+\
                                    ' '+str(startHour)+':' +str(startMin)+':'+str(startSec);
            formatDateTime = "%Y-%j %H:%M:%S";
            dateTime = datetime.strptime(dateTimeStr, formatDateTime)
        #convert back to the string in order to read it by QDateTime
        dateTimeStr = dateTime.strftime("%Y-%m-%d %H:%M:%S")
        return dateTimeStr

    @Variables.undoLocal
    def setStartTime(self, dateTime):
        startYear = dateTime.year
        startDay = dateTime.timetuple().tm_yday
        startHour = dateTime.hour
        startMin = dateTime.minute
        startSec = dateTime.second
        self.__node_atmos.xmlSetData('start_year', startYear)
        self.__node_atmos.xmlSetData('start_day', startDay)
        self.__node_atmos.xmlSetData('start_hour', startHour)
        self.__node_atmos.xmlSetData('start_min', startMin)
        self.__node_atmos.xmlSetData('start_sec', startSec)

    #-------------------------------------------------------------------------
    #-----------------Activate chemistry block--------------------------------
    @Variables.noUndo
    def getChemistryStatus(self):
        """
        Return if Chemistry status is 'on' or 'off'.
        """
        node = self.__node_atmos.xmlInitChildNode(self.act_chemistry)
        if not node[self.status]:
            status = self.__default[self.act_chemistry]
            self.setChemistryStatus(status)
        return node[self.status]


    @Variables.undoLocal
    def setChemistryStatus(self, status):
        """
        Set Chemistry status to 'on' / 'off'.
        """
        self.isOnOff(status)
        self.__node_atmos.xmlInitChildNode(self.act_chemistry)[self.status] = status

        if status == 'off':
            for tag in ['activate_chemistry']:
                for node in self.case.xmlGetNodeList(tag):
                    node['status'] = "off"

    #-------------------------------------------------------------------------
    #-----------------Read meteo file------------------------------------------
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
        elif status == 'on':
            self.setLargeScaleMeteoStatus('off');


    @Variables.noUndo
    def getMeteoDataFileName(self):
        """
        Return the name of the meteo data file.
        """
        f = self.__node_atmos.xmlGetString('meteo_data')
        if f is None:
            f = self.__default['meteo_data']
            self.setMeteoDataFileName(f)
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
                self.__removeProperty(node, "real_temperature")
                FluidCharacteristicsModel(self.case).setPropertyMode('density', 'constant')

        else:
            self.__removeScalar(node, 'ym_water')
            self.__removeScalar(node, 'number_of_droplets')
            self.__removeProperty(node, 'liquid_water')
            self.__removeProperty(node, "real_temperature")


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

    def checkGetandSetLargeScaleMeteoStatus(self):
        """Check whether the AtmosphericFlowsModel class could be set and get the meteo data status"""
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.constant)
        mdl.setLargeScaleMeteoStatus('on')

        doc = """<atmospheric_flows model="constant">
                    <large_scale_meteo status="on"/>
                 </atmospheric_flows>"""

        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set large scale meteo status'
        assert mdl.getLargeScaleMeteoStatus() == 'on', \
            'Could not get large scale meteo status'

    def checkGetandSetChemistryStatus(self):
        """Check whether the AtmosphericFlowsModel class could be set and get the meteo data status"""
        mdl = AtmosphericFlowsModel(self.case)
        mdl.setAtmosphericFlowsModel(AtmosphericFlowsModel.constant)
        mdl.setChemistryStatus('on')

        doc = """<atmospheric_flows model="constant">
                    <activate_chemistry status="on"/>
                 </atmospheric_flows>"""

        assert mdl.atmosphericFlowsNode() == self.xmlNodeFromString(doc), \
            'Could not set activate chemistry status'
        assert mdl.getChemistryStatus() == 'on', \
            'Could not get activate chemistry status'


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
