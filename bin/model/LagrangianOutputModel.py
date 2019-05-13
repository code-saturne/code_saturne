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
This module defines the lagrangian two phase flow modelling management.

This module contains the following classes and function:
- LagrangianOutputModel
- LagrangianOutputTestCase
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import sys, unittest, logging


#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Model, Variables


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------


logging.basicConfig()
log = logging.getLogger("LagrangianOutputModel")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# lagrangian model class
#-------------------------------------------------------------------------------


class LagrangianOutputModel(Model):
    """
    Manage the input/output markups in the xml doc about Lagrangian module.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_lagr = self.case.root().xmlInitNode('lagrangian', 'model')
        self._setDefaultVariables()


    def _defaultLagrangianOutputValues(self):
        """
        Return a dictionnary which contains default values.
        """
        default = {}
        default['resident_time']              = "off"
        default['diameter']                   = "off"
        default['temperature']                = "off"
        default['velocity_particles']         = "off"
        default['velocity_fluid_seen']        = "off"
        default['mass']                       = "off"
        default['shrinking_core_diameter']    = "off"
        default['raw_coal_mass_fraction']     = "off"
        default['char_mass_fraction']         = "off"
        default['moisture_mass_fraction']     = "off"
        return default


    def _setDefaultVariables(self):
        """
        Set variables and properties if lagrangian model is on.
        """
        self.node_output = self.node_lagr.xmlInitChildNode('output')


    @Variables.undoLocal
    def setFluidVelocityStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'velocity of the locally undisturbed fluid' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_velocity = self.node_output.xmlInitChildNode('velocity_fluid_seen', 'status')
        node_velocity['status'] = status


    @Variables.noUndo
    def getFluidVelocityStatus(self):
        """
        Return status for association of the variable 'velocity of the locally
        undisturbed fluid' with the display.
        """
        node_velocity = self.node_output.xmlInitChildNode('velocity_fluid_seen', 'status')
        status = node_velocity['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['velocity_fluid_seen']
            self.setFluidVelocityStatus(status)
        return status


    @Variables.undoLocal
    def setParticlesVelocityStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'particle velocity' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_velocity = self.node_output.xmlInitChildNode('velocity_particles', 'status')
        node_velocity['status'] = status


    @Variables.noUndo
    def getParticlesVelocityStatus(self):
        """
        Return status for association of the variable 'particle velocity'
        with the display.
        """
        node_velocity = self.node_output.xmlInitChildNode('velocity_particles', 'status')
        status = node_velocity['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['velocity_particles']
            self.setParticlesVelocityStatus(status)
        return status


    @Variables.undoLocal
    def setResidentTimeStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'resident time' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_rtime = self.node_output.xmlInitChildNode('resident_time', 'status')
        node_rtime['status'] = status


    @Variables.noUndo
    def getResidentTimeStatus(self):
        """
        Return status for association of the variable 'resident time'
        with the display.
        """
        node_rtime = self.node_output.xmlInitChildNode('resident_time', 'status')
        status = node_rtime['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['resident_time']
            self.setResidentTimeStatus(status)
        return status


    @Variables.undoLocal
    def setParticleDiameterStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'particle diameter' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_diam = self.node_output.xmlInitChildNode('diameter', 'status')
        node_diam['status'] = status


    @Variables.noUndo
    def getParticleDiameterStatus(self):
        """
        Return status for association of the variable 'particle diameter'
        with the display.
        """
        node_diam = self.node_output.xmlInitChildNode('diameter', 'status')
        status = node_diam['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['diameter']
            self.setParticleDiameterStatus(status)
        return status


    @Variables.undoLocal
    def setParticleTemperatureStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'particle temperature' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_temp = self.node_output.xmlInitChildNode('temperature', 'status')
        node_temp['status'] = status


    @Variables.noUndo
    def getParticleTemperatureStatus(self):
        """
        Return status for association of the variable 'particle temperature'
        with the display.
        """
        node_temp = self.node_output.xmlInitChildNode('temperature', 'status')
        status = node_temp['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['temperature']
            self.setParticleTemperatureStatus(status)
        return status


    @Variables.undoLocal
    def setParticleMassStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'particle mass' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_mass = self.node_output.xmlInitChildNode('mass', 'status')
        node_mass['status'] = status


    @Variables.noUndo
    def getParticleMassStatus(self):
        """
        Return status for association of the variable 'particle mass'
        with the display.
        """
        node_mass = self.node_output.xmlInitChildNode('mass', 'status')
        status = node_mass['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['mass']
            self.setParticleMassStatus(status)
        return status


    @Variables.undoLocal
    def setCoalParticleDiameterStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'shrinking core diameter of the coal particles' with the display
        (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_diam = self.node_output.xmlInitChildNode('shrinking_core_diameter', 'status')
        node_diam['status'] = status


    @Variables.noUndo
    def getCoalParticleDiameterStatus(self):
        """
        Return status for association of the variable
        'shrinking core diameter of the coal particles' with the display.
        """
        node_diam = self.node_output.xmlInitChildNode('shrinking_core_diameter', 'status')
        status = node_diam['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['shrinking_core_diameter']
            self.setCoalParticleDiameterStatus(status)
        return status


    @Variables.undoLocal
    def setCoalParticleMassStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'mass of reactive coal of the coal particles' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_mass = self.node_output.xmlInitChildNode('raw_coal_mass_fraction', 'status')
        node_mass['status'] = status


    @Variables.noUndo
    def getCoalParticleMassStatus(self):
        """
        Return status for association of the variable 'mass of reactive coal of the coal particles'
        with the display.
        """
        node_mass = self.node_output.xmlInitChildNode('raw_coal_mass_fraction', 'status')
        status = node_mass['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['raw_coal_mass_fraction']
            self.setCoalParticleMassStatus(status)
        return status


    @Variables.undoLocal
    def setCokeParticleMassStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'mass of char of the coal particles' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_mass = self.node_output.xmlInitChildNode('char_mass_fraction', 'status')
        node_mass['status'] = status


    @Variables.noUndo
    def getCokeParticleMassStatus(self):
        """
        Return status for association of the variable 'mass of char of the coal particles'
        with the display.
        """
        node_mass = self.node_output.xmlInitChildNode('char_mass_fraction', 'status')
        status = node_mass['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['char_mass_fraction']
            self.setCokeParticleMassStatus(status)
        return status


    @Variables.undoLocal
    def setMoistureMassStatus(self, status):
        """
        Update the status markup from the XML document to associate the variable
        'mass of moisture' with the display (trajectory or particles) mode.
        """
        self.isOnOff(status)
        node_mass = self.node_output.xmlInitChildNode('moisture_mass_fraction', 'status')
        node_mass['status'] = status


    @Variables.noUndo
    def getMoistureMassStatus(self):
        """
        Return status for association of the variable 'mass of moisture'
        with the display.
        """
        node_mass = self.node_output.xmlInitChildNode('moisture_mass_fraction', 'status')
        status = node_mass['status']
        if not status:
            status = self._defaultLagrangianOutputValues()['moisture_mass_fraction']
            self.setMoistureMassStatus(status)
        return status


#-------------------------------------------------------------------------------
# LagrangianOutput test case
#-------------------------------------------------------------------------------


class LagrangianOutputTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        self.case = Case()
        XMLinit(self.case).initialize()


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case


    def checkLagrangianOutputInstantiation(self):
        """
        Check whether the LagrangianOutputModel class could be instantiated
        """
        model = None
        model = LagrangianOutputModel(self.case)

        assert model != None, 'Could not instantiate LagrangianOutputModel'


    def checkLagrangianOutputDefaultValues(self):
        """
        Check the default values
        """
        model = LagrangianOutputModel(self.case)
        doc = """
        <output>
        <resident_time status="off"/>
        <diameter status="off"/>
        <velocity_particles status="off"/>
        <velocity_fluid_seen status="off"/>
        <mass status="off"/>
        <shrinking_core_diameter status="off"/>
        <raw_coal_mass_fraction status="off"/>
        <char_mass_fraction status="off"/>
        </output>"""

        assert model.node_output == self.xmlNodeFromString(doc),\
               'Could not get default values for LagrangianOutputModel model '


    def checkSetandGetListingFrequency(self):
        """
        Check whether the listing frequency method could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        value = mdl.getListingFrequency()
        assert value == 1 ,\
        'Could not get default value for '
        mdl.setListingFrequency(1234)
        doc = """
        <listing_printing_frequency>
        1234
        </listing_printing_frequency>
        """

        assert mdl.node_output.xmlInitChildNode('listing_printing_frequency') == self.xmlNodeFromString(doc) ,\
            'Could not set values for'


    def checkSetandGetParticlesStatus(self):
        """
        Check whether the particles mode method could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getParticlesStatus()
        assert status == 'off' ,\
        'Could not get default values for particles status'
        mdl.setParticlesStatus('on')
        doc = """
        <particles status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('particles') == self.xmlNodeFromString(doc) ,\
            'Could not set values for particles mode status'


    def checkSetandGetFluidVelocityStatus(self):
        """
        Check whether the method for 'velocity of the locally undisturbed fluid' association
        with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getFluidVelocityStatus()
        assert status == 'off' ,\
        'Could not get default values for velocity_fluid_seen status'
        mdl.setFluidVelocityStatus('on')
        doc = """
        <velocity_fluid_seen status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('velocity_fluid_seen') == self.xmlNodeFromString(doc) ,\
            'Could not set values for velocity_fluid_seen status '


    def checkSetandGetParticlesVelocityStatus(self):
        """
        Check whether the method for 'particle velocity' association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getParticlesVelocityStatus()
        assert status == 'off' ,\
        'Could not get default values for velocity_particles status'
        mdl.setParticlesVelocityStatus('on')
        doc = """
        <velocity_particles status="on"/>
       """

        assert mdl.node_output.xmlInitChildNode('velocity_particles') == self.xmlNodeFromString(doc) ,\
            'Could not set values for velocity_particles status'


    def checkSetandGetResidentTimeStatus(self):
        """
        Check whether the method for 'resident time' association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getResidentTimeStatus()
        assert status == 'off' ,\
        'Could not get default values for resident time status'
        mdl.setResidentTimeStatus('on')
        doc = """
        <resident_time status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('resident_time') == self.xmlNodeFromString(doc) ,\
            'Could not set values for resident time status'


    def checkSetandGetParticleDiameterStatus(self):
        """
        Check whether the method for 'particle diameter' association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getParticleDiameterStatus()
        assert status == 'off' ,\
        'Could not get default values for diameter status'
        mdl.setParticleDiameterStatus('on')
        doc = """
        <diameter status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('diameter') == self.xmlNodeFromString(doc) ,\
            'Could not set values for diameter status'


    def checkSetandGetParticleTemperatureStatus(self):
        """
        Check whether the method for 'particle temperature' association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getParticleTemperatureStatus()
        assert status == 'off' ,\
        'Could not get default values for temperature status'
        mdl.setParticleTemperatureStatus('on')
        doc = """
        <temperature status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('temperature') == self.xmlNodeFromString(doc) ,\
            'Could not set values for temperature status'


    def checkSetandGetParticleMassStatus(self):
        """
        Check whether the method for 'particle mass' association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getParticleMassStatus()
        assert status == 'off' ,\
        'Could not get default values for mass status'
        mdl.setParticleMassStatus('on')
        doc = """
        <mass status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('mass') == self.xmlNodeFromString(doc) ,\
            'Could not set values for mass status'


    def checkSetandGetCoalParticleDiameterStatus(self):
        """
        Check whether the method for 'shrinking core diameter of the coal particles'
        association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getCoalParticleDiameterStatus()
        assert status == 'off' ,\
        'Could not get default values for shrinking_core_diameter status'
        mdl.setCoalParticleDiameterStatus('on')
        doc = """
        <shrinking_core_diameter status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('shrinking_core_diameter') == self.xmlNodeFromString(doc) ,\
            'Could not set values for shrinking_core_diameter status'


    def checkSetandGetCoalParticleMassStatus(self):
        """
        Check whether the method for 'mass of reactive coal of the coal particles'
        association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getCoalParticleMassStatus()
        assert status == 'off' ,\
        'Could not get default values for raw_coal_mass_fraction status'
        mdl.setCoalParticleMassStatus('on')
        doc = """
        <raw_coal_mass_fraction status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('raw_coal_mass_fraction') == self.xmlNodeFromString(doc) ,\
            'Could not set values for raw_coal_mass_fraction status'


    def checkSetandGetCokeParticleMassStatus(self):
        """
        Check whether the method for 'mass of char of the coal particles'
        association with display could be set and get
        """
        mdl = LagrangianOutputModel(self.case)
        status = mdl.getCokeParticleMassStatus()
        assert status == 'off' ,\
        'Could not get default values for status'
        mdl.setCokeParticleMassStatus('on')
        doc = """
        <char_mass_fraction status="on"/>
        """

        assert mdl.node_output.xmlInitChildNode('char_mass_fraction') == self.xmlNodeFromString(doc) ,\
            'Could not set values for status'



def suite():
    testSuite = unittest.makeSuite(LagrangianOutputTestCase, "check")
    return testSuite


def runTest():
    print("LagrangianOutputTestCase TODO*********.")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
