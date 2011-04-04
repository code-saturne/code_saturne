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
This module initialize model dynamics variables and model scalars

This module contents the following classes:
- InitializationModel
- InitializationTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, string, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import XMLmodel, ModelTest
from Base.XMLvariables import Model
from Pages.TurbulenceModel import TurbulenceModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.CoalCombustionModel import CoalCombustionModel
from Pages.ElectricalModelsModel import ElectricalModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.LocalizationModel import LocalizationModel

#-------------------------------------------------------------------------------
# Variables and Scalar model initialization modelling class
#-------------------------------------------------------------------------------

class InitializationModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        models = self.case.xmlGetNode('thermophysical_models')
        self.node_veloce   = models.xmlGetNode('velocity_pressure')
        self.node_turb     = models.xmlGetNode('turbulence', 'model')
        self.node_therm    = models.xmlGetNode('thermal_scalar', 'model')

        self.VelocityList = ('velocity_U', 'velocity_V', 'velocity_W')
        self.Turb_var_List = ('turb_k', 'turb_eps',
                              'component_R11', 'component_R22', 'component_R33',
                              'component_R12', 'component_R13', 'component_R23',
                              'turb_phi', 'turb_fb', 'turb_omega')

        self.turb = TurbulenceModel(self.case)
        self.turbulenceModes = ('values',
                                'reference_velocity',
                                'reference_velocity_length')


    def __defaultValues(self):
        """
        Private method.
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['init_mode']     = "reference_velocity"
        default['velocity']      = 0.0
        default['turb_velocity'] = 1.0
        default['length']        = 1

        # Initialization of k, eps, phi, fb, omega, R11, R22, R33,
        # with standard formulae

        Uref  = float(1.)
        Almax = float(1.)
        Cmu   = 0.09

        default['turb_k']        = 1.5 * pow((0.02*Uref), 2)
        default['turb_eps']      = pow(default['turb_k'], 1.5) / Cmu / Almax
        default['turb_phi']      = float(2./3.)
        default['turb_fb']       = 0
        default['turb_omega']    = default['turb_eps'] / default['turb_k'] / Cmu
        default['component_R11'] = pow((0.02*Uref), 2)
        default['component_R22'] = pow((0.02*Uref), 2)
        default['component_R33'] = pow((0.02*Uref), 2)
        default['component_R12'] = 0.0
        default['component_R23'] = 0.0
        default['component_R13'] = 0.0
        default['temperature']   = 0.0
        default['pressure']      = 0.0

        return default


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())


    def __initTurbulenceInitialValues(self, zone):
        """
        Private method.
        Initialize values for turbulence variables.
        """
        self.__verifyZone(zone)
        turb_model = self.turb.getTurbulenceModel()

        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            for txt in ('turb_k', 'turb_eps'):
                self.getTurbulenceInitialValue(zone, txt)

        elif turb_model in ('Rij-epsilon', 'Rij-epsilon-SSG'):
            for txt in ('component_R11',
                        'component_R22',
                        'component_R33',
                        'component_R12',
                        'component_R13',
                        'component_R23',
                        'turb_eps'):
                self.getTurbulenceInitialValue(zone, txt)

        elif turb_model == 'v2f-phi':
            for txt in ('turb_k',
                        'turb_eps',
                        'turb_phi',
                        'turb_fb'):
                self.getTurbulenceInitialValue(zone, txt)

        elif turb_model == 'k-omega-SST':
            for txt in ('turb_k', 'turb_omega'):
                self.getTurbulenceInitialValue(zone, txt)


    def __setDefaultTurbulenceInitialValues(self, zone, choice):
        """
        Private method.
        Set the values of initialization by default according
        to the choice of the initialization of turbulence
        """
        self.__verifyZone(zone)
        self.isInList(choice, self.turbulenceModes)

        if choice == 'reference_velocity':
            self.getReferenceVelocity()
        elif choice == 'reference_velocity_length':
            self.getReferenceVelocityAndLength()
        elif choice == 'values':
            self.__initTurbulenceInitialValues(zone)


    def setInitialVelocity(self, zone, var_name, var_init):
        """
        Public method.
        Set values for initialization of velocity.
        """
        self.__verifyZone(zone)
        self.isFloat(var_init)
        self.isInList(var_name, self.VelocityList)

        node = self.node_veloce.xmlGetNode('variable', name=var_name)
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        if var_init == 0:
            n.xmlRemoveNode()
        else:
            n.xmlSetTextNode(var_init)


    def getInitialVelocity(self, zone):
        """
        Public method.
        Return values for initialization of velocity in a list.
        """
        self.__verifyZone(zone)
        velocity_list = []

        for var_name in self.VelocityList:
            node = self.node_veloce.xmlGetNode('variable', name=var_name)
            if not node:
                msg = "There is an error: this node " + str(node) + "should be existed"
                raise ValueError(msg)
            v = node.xmlGetDouble('initial_value', zone_id=zone)
            if v == None:
                v = self.__defaultValues()['velocity']
            velocity_list.append(v)

        return velocity_list


    def setInitialTurbulenceChoice(self, zone, init_mode):
        """
        Public method.
        Set the initialization mode in the attribute choice.
        """
        self.__verifyZone(zone)
        self.isInList(init_mode, self.turbulenceModes)

        node_init = self.node_turb.xmlInitNode('initialization')
        node_init['choice'] = init_mode
        self.__setDefaultTurbulenceInitialValues(zone, init_mode)


    def getInitialTurbulenceChoice(self, zone):
        """
        Public method.
        Return the turbulence initialization choice.
        """
        self.__verifyZone(zone)
        node_init = self.node_turb.xmlInitNode('initialization')
        choice = node_init['choice']

        if choice not in self.turbulenceModes:
            choice = self.__defaultValues()['init_mode']
            node_init['choice'] = choice
        self.__setDefaultTurbulenceInitialValues(zone, choice)

        return choice


    def setReferenceVelocity(self, velocity):
        """
        Public method.
        Set the reference velocity for initialization of turbulence (UREF).
        """
        self.isFloat(velocity)

        node_init = self.node_turb.xmlGetNode('initialization')
        if not node_init:
            msg = "There is an error: this node " + str(node_init) + "should be existed"
            raise ValueError(msg)

        self.isInList(node_init['choice'], ('reference_velocity','reference_velocity_length'))

        node_init.xmlSetData('reference_velocity', velocity)


    def getReferenceVelocity(self):
        """
        Public method.
        Return the reference velocity for initialization of turbulence(UREF).
        """
        node_init = self.node_turb.xmlGetNode('initialization',
                                              choice='reference_velocity')
        if not node_init:
            msg = "There is an error: this node " + str(node_init) + "should be existed"
            raise ValueError(msg)

        v = node_init.xmlGetDouble('reference_velocity')
        if v == None:
            v = self.__defaultValues()['turb_velocity']
            self.setReferenceVelocity(v)

        return v


    def setReferenceLength(self, length):
        """
        Public method.
        Set the reference length (ALMAX).
        """
        self.isStrictPositiveFloat(length)

        node_init = self.node_turb.xmlGetNode('initialization',
                                               choice='reference_velocity_length')
        if not node_init:
            msg = "There is an error: this node " + str(node_init) + "should be existed"
            raise ValueError(msg)

        node_init.xmlSetData('reference_length', length)


    def getReferenceVelocityAndLength(self):
        """
        Public method.
        Return the reference velocity and length.
        """
        node_init = self.node_turb.xmlGetNode('initialization',
                                               choice='reference_velocity_length')
        if not node_init:
            msg = "There is an error: this node " + str(node_init) + "should be existed"
            raise ValueError(msg)

        v = node_init.xmlGetDouble('reference_velocity')
        if v == None:
            v = self.__defaultValues()['turb_velocity']
            self.setReferenceVelocity(v)

        l = node_init.xmlGetDouble('reference_length')
        if l == None:
            l = self.__defaultValues()['length']
            self.setReferenceLength(l)

        return v, l


    def setTurbulenceInitialValue(self, zone, var_name, var_init):
        """
        Public method.
        Set values for initialization of turbulence variables.
        """
        self.__verifyZone(zone)
        self.isInList(var_name, self.Turb_var_List)
        self.isFloat(var_init)

        node = self.node_turb.xmlGetNode('variable', name=var_name)
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        node.xmlInitNode('initial_value', zone_id=zone)
        node.xmlSetData('initial_value', var_init, zone_id=zone)


    def getTurbulenceInitialValue(self, zone, var_name):
        """
        Public method.
        Return values for initialization of turbulence variable.
        """
        self.__verifyZone(zone)
        self.isInList(var_name, self.Turb_var_List)

        node = self.node_turb.xmlGetNode('variable', name=var_name)
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        init = node.xmlGetDouble('initial_value', zone_id=zone)
        if init == None:
            init = self.__defaultValues()[var_name]
            self.setTurbulenceInitialValue(zone, var_name, init)

        return init


    def getTurbulenceVariableLabel(self, var_name):
        """
        Public method.
        Return the label of the turbulence variables. Only for the view.
        """
        self.isInList(var_name, self.Turb_var_List)

        label = Tool.dicoLabel(var_name)
        node = self.node_turb.xmlGetNode('variable', name=var_name)
        if node:
            label = node['label']

        return label


    def getVelocityLabel(self):
        """
        Public method.
        Return the list of users labels for velocity. Only for the view
        """
        label_list = []

        for var_name in self.VelocityList:
            node = self.node_veloce.xmlGetNode('variable', name=var_name)
            if not node:
                msg = "There is an error: this node " + str(node) + "should be existed"
                raise ValueError(msg)
            label_list.append(node['label'])

        return label_list


#-------------------------------------------------------------------------------
# InitializationModel test case
#-------------------------------------------------------------------------------


class InitializationTestCase(ModelTest):
    """
    Unittest.
    """
    def checkInitializationInstantiation(self):
        """Check whether the InitializationModel class could be instantiated."""
        model = None
        model = InitializationModel(self.case)
        assert model != None, 'Could not instantiate InitializationModel'


    def checkSetAndGetInitialVelocity(self):
        """Check whether the velocity initial value could be set and get"""
        mdl = InitializationModel(self.case)
        zone = '1'
        mdl.setInitialVelocity(zone, 'velocity_U', 123)
        mdl.setInitialVelocity(zone, 'velocity_V', 1.e+2)
        mdl.setInitialVelocity(zone, 'velocity_W', 0)

        doc = '''<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U">
                            <initial_value zone_id="1">123</initial_value>
                    </variable>
                    <variable label="VelocitV" name="velocity_V">
                            <initial_value zone_id="1">100.0</initial_value>
                    </variable>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>'''

        assert mdl.node_veloce == self.xmlNodeFromString(doc), \
          'Could not set the initial values of velocity'

        assert mdl.getInitialVelocity(zone) == [123.0, 100.0, 0.0], \
          'Could not get the initial velocity values'


    def checkSetAndGetInitialTurbulenceChoice(self):
        """Check whether the velocity initial value choice could be set."""
        mdl = InitializationModel(self.case)
        zone = '1'
        mdl.turb.setTurbulenceModel('k-epsilon')
        mdl.setInitialTurbulenceChoice('1', 'reference_velocity')

        doc = '''<turbulence model="k-epsilon">
                    <variable label="TurbEner" name="turb_k"/>
                    <variable label="Dissip" name="turb_eps"/>
                    <property label="TurbVisc" name="turb_viscosity"/>
                    <initialization choice="reference_velocity">
                            <reference_velocity>1</reference_velocity>
                    </initialization>
                 </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set choice for initialization of the turbulence'

        assert mdl.getInitialTurbulenceChoice('1') == 'reference_velocity',\
           'Could not get choice of initialization of the turbulence'


    def checkSetAndGetReferenceVelocity(self):
        """Check whether the reference velocity value could be set and get"""
        mdl = InitializationModel(self.case)
        zone = '1'
        mdl.turb.setTurbulenceModel('k-epsilon')
        mdl.setReferenceVelocity(5)
        doc = '''<turbulence model="k-epsilon">
                    <variable label="TurbEner" name="turb_k"/>
                    <variable label="Dissip" name="turb_eps"/>
                    <property label="TurbVisc" name="turb_viscosity"/>
                    <initialization choice="reference_velocity">
                            <reference_velocity>5</reference_velocity>
                    </initialization>
                 </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set values of reference velocity for initialization of the turbulence'

        assert mdl.getReferenceVelocity() == 5,\
           'Could not get values of reference velocity for initialization of the turbulence'


    def checkSetAndGetReferenceLength(self):
        """Check whether the reference lenght value could be set and get"""
        mdl = InitializationModel(self.case)
        zone = '1'
        mdl.turb.setTurbulenceModel('k-epsilon')
        mdl.setInitialTurbulenceChoice('1', 'reference_velocity_length')
        mdl.setReferenceVelocity(5)
        mdl.setReferenceLength(155.8)

        doc = '''<turbulence model="k-epsilon">
                    <variable label="TurbEner" name="turb_k"/>
                    <variable label="Dissip" name="turb_eps"/>
                    <property label="TurbVisc" name="turb_viscosity"/>
                    <initialization choice="reference_velocity_length">
                            <reference_velocity>5</reference_velocity>
                            <reference_length>155.8</reference_length>
                    </initialization>
                 </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set values of reference length for initialization of the turbulence'

        assert mdl.getReferenceVelocityAndLength() == (5, 155.8),\
           'Could not get values of reference length for initialization of the turbulence'


    def checkSetAndGetTurbulenceInitialValue(self):
        """Check whether the initial values of turbulence could be set.and get"""
        mdl = InitializationModel(self.case)
        zone = '1'
        mdl.turb.setTurbulenceModel('Rij-epsilon')
        mdl.setTurbulenceInitialValue(zone, 'component_R11', 0.0011)
        mdl.setTurbulenceInitialValue(zone, 'component_R22', 0.0022)
        mdl.setTurbulenceInitialValue(zone, 'component_R33', 0.0033)
        doc = '''<turbulence model="Rij-epsilon">
                    <variable label="TurbEner" name="turb_k"/>
                    <variable label="Dissip" name="turb_eps"/>
                    <property label="TurbVisc" name="turb_viscosity"/>
                    <initialization choice="reference_velocity">
                            <reference_velocity>1</reference_velocity>
                    </initialization>
                    <variable label="R11" name="component_R11">
                            <initial_value zone_id="1">0.0011</initial_value>
                    </variable>
                    <variable label="R22" name="component_R22">
                            <initial_value zone_id="1">0.0022</initial_value>
                    </variable>
                    <variable label="R33" name="component_R33">
                            <initial_value zone_id="1">0.0033</initial_value>
                    </variable>
                    <variable label="R12" name="component_R12"/>
                    <variable label="R13" name="component_R13"/>
                    <variable label="R23" name="component_R23"/>
                 </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set initial values of the turbulence'

        assert mdl.getTurbulenceInitialValue(zone, 'component_R22') == 0.0022,\
           'Could not get initial values of the turbulence'


def suite():
    testSuite = unittest.makeSuite(InitializationTestCase, "check")
    return testSuite


def runTest():
    print("InitializationTestCase - OK !!!!")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
