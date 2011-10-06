# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
This module defines the values of reference.

This module contains the following classes and function:
- MobileMeshModel
- MobileMeshTestCase
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
from Base.XMLmodel import ModelTest
from Pages.OutputControlModel import OutputControlModel

#-------------------------------------------------------------------------------
# Mobil Mesh model class
#-------------------------------------------------------------------------------

class MobileMeshModel(Model):
    """
    Manage the input/output markups in the xml doc about mobil mesh
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_ale    = self.node_models.xmlInitChildNode('ale_method', 'status')

        self.out = OutputControlModel(self.case)


    def __defaultInitialValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['nalinf' ]  = 0
        default['iortvm' ]  = 'isotrop'
        default['mei'    ]  = 'user_subroutine'
        default['formula_isotrop']  = 'mesh_vi1 = 1;'
        default['formula_orthotrop'] = 'mesh_vi1 = 1;\nmesh_vi2 = 1;\nmesh_vi3 = 1;'
        default['ale_method']  = 'off'

        return default


    def __setVariablesandProperties(self):
        """
        Set variables and properties if ALE method is activated.
        """
        self.node_ale.xmlInitChildNode('variable', name='mesh_velocity_U', label='mesh_u')
        self.node_ale.xmlInitChildNode('variable', name='mesh_velocity_V', label='mesh_v')
        self.node_ale.xmlInitChildNode('variable', name='mesh_velocity_W', label='mesh_w')
        self.node_ale.xmlInitChildNode('property', name='mesh_viscosity_1', label='mesh_vi1')

        # find node for property / choice and set to default value if require
        node = self.node_ale.xmlGetChildNode('property', name='mesh_viscosity_1')
        if not node['choice']:
            node['choice'] = self.__defaultInitialValues()['mei']

        self.__updateNodeViscosity()


    def __updateNodeViscosity(self):
        """
        Update properties beyond mesh visosity is isotrope or not.
        """
        if self.getViscosity() == 'orthotrop':
            self.node_ale.xmlInitChildNode('property', name='mesh_viscosity_2', label='mesh_vi2')
            self.node_ale.xmlInitChildNode('property', name='mesh_viscosity_3', label='mesh_vi3')

            # find choice for mesh_viscosity_1, create default value if require
            node = self.node_ale.xmlGetChildNode('property', name='mesh_viscosity_1')
            if not node['choice']:
                node['choice'] = self.__defaultInitialValues()['mei']
            mei = node['choice']

            # Syncrhonize other properties
            for n in ('mesh_viscosity_2', 'mesh_viscosity_3'):
                node = self.node_ale.xmlGetChildNode('property', name=n)
                node['choice'] = mei

        else:
            node1 = self.node_ale.xmlGetChildNode('property', name='mesh_viscosity_2')
            node2 = self.node_ale.xmlGetChildNode('property', name='mesh_viscosity_3')

            if node1:
                node1.xmlRemoveNode()
            if node2:
                node2.xmlRemoveNode()


    def __removeVariablesandProperties(self):
        """
        Remove variables and properties if ALE method is disabled.
        """
        self.node_ale.xmlRemoveChild('variable')
        self.node_ale.xmlRemoveChild('property')

##    def setMethod(self, answer):
##        """
##        Set method of activation of ALE into xml file.
##        """
##        self.isOnOff(answer)
##        typ = ''
##        typ = self.out.getTypePostMeshes()
##        self.node_ale['status'] = answer
##       if answer == 'on':
##            if typ not in ('10', '11', '12'):
##                typ = '10'
##            self.__setVariablesandProperties()
##        else:
##            if typ not in ('0', '1', '2'):
##                typ = '0'
##        self.out.setTypePostMeshes(typ)
##
##
##    def getMethod(self):
##        """
##        Get method of activation of ALE from xml file.
##        """
##        if self.node_ale['status'] == '':
##            status = 'off'
##            self.setMethod(status)
##        else:
##            status = self.node_ale['status']
##        return status

    def getMethod(self):
        """
        Get status on balise "ALE" from xml file
        """
        status = self.node_ale['status']
        if not status:
            status = self.__defaultInitialValues()['ale_method']
            self.setMethod(status)
        return status


    def setMethod(self, status):
        """
        Put status on balise "ALE" in xml file
        """
        self.isOnOff(status)
        typ = ''
        typ = self.out.getWriterTimeDependency("-1")
        self.node_ale['status'] = status
        if status == 'on':
            if typ not in ('10', '11', '12'):
                typ = '10'
            self.__setVariablesandProperties()
        else:
            if typ not in ('0', '1', '2'):
                typ = '0'
            self.__removeVariablesandProperties()
        self.out.setWriterTimeDependency("-1", 'fixed_mesh')


    def setSubIterations(self, value):
        """
        Set value of fluid initialization sub iterations into xml file.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_ale.xmlSetData('fluid_initialization_sub_iterations', value)


    def getSubIterations(self):
        """
        Get value of fluid initialization sub iterations from xml file.
        """
        nalinf = self.node_ale.xmlGetInt('fluid_initialization_sub_iterations')
        if not nalinf:
            nalinf = self.__defaultInitialValues()['nalinf']
            self.setSubIterations(nalinf)
        return nalinf


    def setViscosity(self, value):
        """
        Set value of mesh viscosity into xml file.
        """
        self.isInList(value, ['isotrop', 'orthotrop'])
        node = self.node_ale.xmlInitChildNode('mesh_viscosity')
        node['type'] = value
        self.__updateNodeViscosity()

        # update formula
        if self.getMEI() == 'user_function':
            self.setFormula(self.getDefaultFormula())


    def getViscosity(self):
        """
        Get value of mesh viscosity from xml file.
        """
        iortvm = self.__defaultInitialValues()['iortvm']
        node = self.node_ale.xmlGetChildNode('mesh_viscosity')
        if node :
            iortvm = node['type']
        else:
            self.setViscosity(iortvm)

        return iortvm


    def setMEI(self, value):
        """
        Set value of spatial distribution of the viscosity of the mesh.
        """
        self.isInList(value, ['user_subroutine', 'user_function'] )

        # do something only if mei has changed
        if self.getMEIWithoutDefaultValue() != value:
            for node in self.node_ale.xmlGetNodeList('property'):
                node['choice'] = value
            if value == 'user_subroutine':
                self.node_ale.xmlRemoveChild('formula')
            else:
                self.setFormula(self.getDefaultFormula())


    def getMEI(self):
        """
        Get value of spatial distribution of the viscosity of the mesh.
        """
        # Get the first node
        mei  = self.getMEIWithoutDefaultValue()

        if not mei:
            mei = self.__defaultInitialValues()['mei']
            self.setMEI(mei)
        return mei


    def getMEIWithoutDefaultValue(self):
        """
        Get value of spatial distribution of the viscosity of the mesh.
        Return null if no value is set
        """
        # Get the first node
        node = self.node_ale.xmlGetNode('property', label='mesh_vi1')
        mei  = None

        if node:
            mei = node['choice']
        return mei


    def setFormula(self, value):
        """
        Set the formula for the viscosity of mesh
        """
        self.node_ale.xmlSetData('formula', value)


    def getFormula(self):
        """
        Get the formula for the viscosity of mesh
        """
        formula = self.node_ale.xmlGetString('formula')
        if not formula:
            formula = self.getDefaultFormula()
            self.setFormula(formula)
        return formula


    def getDefaultFormula(self):
        """
        Get the default formula base on viscosity type
        """
        viscosity = self.getViscosity()
        return self.__defaultInitialValues()['formula_'+ viscosity ]


#-------------------------------------------------------------------------------
# MobileMesh Model test case
#-------------------------------------------------------------------------------


class MobileMeshTestCase(ModelTest):
    """
    """
    def checkMobileMeshInstantiation(self):
        """Check whether the TurbulenceModel class could be instantiated"""
        model = None
        model = MobileMeshModel(self.case)
        assert model != None, 'Could not instantiate MobileMeshModel'

    def checkGetandSetMethod(self):
        """Check whether the MobileMeshModel class could be set and get method"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        doc = """<ale_method status="on">
                    <variable label="mesh_u" name="mesh_velocity_U"/>
                    <variable label="mesh_v" name="mesh_velocity_V"/>
                    <variable label="mesh_w" name="mesh_velocity_W"/>
                    <property choice="user_subroutine" label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="isotrop"/>
                 </ale_method>"""
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set mobil mesh model method'
        assert mdl.getMethod() == 'on',\
            'Could not get mobil mesh model method'

    def checkGetandSetSubIterations(self):
        """Check whether the MobileMeshModel class could be set and get sub iterations"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        mdl.setSubIterations(12)
##
        doc = """<ale_method status="on">
                    <variable label="mesh_u" name="mesh_velocity_U"/>
                    <variable label="mesh_v" name="mesh_velocity_V"/>
                    <variable label="mesh_w" name="mesh_velocity_W"/>
                    <property choice="user_subroutine" label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="isotrop"/>
                    <fluid_initialization_sub_iterations>12</fluid_initialization_sub_iterations>
                </ale_method>"""

        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set mobil mesh model sub iterations'
        assert mdl.getSubIterations() == 12,\
            'Could not get mobil mesh model sub iteration'


    def checkGetandSetViscosity(self):
        """Check whether the MobileMeshModel class could be set and get viscosity"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        mdl.setViscosity('orthotrop')

        doc = """<ale_method status="on">
                    <variable label="mesh_u" name="mesh_velocity_U"/>
                    <variable label="mesh_v" name="mesh_velocity_V"/>
                    <variable label="mesh_w" name="mesh_velocity_W"/>
                    <property choice="user_subroutine" label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="orthotrop"/>
                    <property choice="user_subroutine" label="mesh_vi2" name="mesh_viscosity_2"/>
                    <property choice="user_subroutine" label="mesh_vi3" name="mesh_viscosity_3"/>
                </ale_method>"""
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set mobil mesh model visocity type'
        assert mdl.getViscosity() == 'orthotrop',\
            'Could not get mobil mesh model viscosity type'

    def checkGetAndSetMEI(self):
        """Check whether the MobileMeshModel class could be set and get mei"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        mdl.setMEI('user_subroutine')

        doc = """<ale_method status="on">
                    <variable label="mesh_u" name="mesh_velocity_U"/>
                    <variable label="mesh_v" name="mesh_velocity_V"/>
                    <variable label="mesh_w" name="mesh_velocity_W"/>
                    <property choice="user_subroutine" label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="isotrop"/>
                    </ale_method> """
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set mei'
        assert mdl.getMEI() == 'user_subroutine',\
            'Could not get mei'

    def checkGetAndSetFormula(self):
        """Check whether the MobileMeshModel class could be set and get formula"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        mdl.setFormula('mesh_vi1 = 1000;')

        doc = """<ale_method status="on">
                    <variable label="mesh_u" name="mesh_velocity_U"/>
                    <variable label="mesh_v" name="mesh_velocity_V"/>
                    <variable label="mesh_w" name="mesh_velocity_W"/>
                    <property choice="user_subroutine" label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="isotrop"/>
                    <formula>mesh_vi1 = 1000;</formula>
                    </ale_method> """
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set formula'
        assert mdl.getFormula() == 'mesh_vi1 = 1000;',\
            'Could not get formula'


def suite():
    testSuite = unittest.makeSuite(MobileMeshTestCase, "check")
    return testSuite

def runTest():
    print("MobileMeshTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
