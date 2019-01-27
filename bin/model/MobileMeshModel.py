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

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.OutputControlModel import OutputControlModel
from code_saturne.model.Boundary import Boundary
from code_saturne.model.LocalizationModel import LocalizationModel

#-------------------------------------------------------------------------------
# Mobil Mesh model class
#-------------------------------------------------------------------------------

class MobileMeshModel(Model):
    """
    Manage the input/output markups in the xml doc about mobile mesh
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case


    def __defaultInitialValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['nalinf' ]  = 0
        default['iortvm' ]  = 'isotrop'
        default['formula_isotrop']  = 'mesh_viscosity_1 = 1;'
        default['formula_orthotrop'] = 'mesh_viscosity_1 = 1;\nmesh_viscosity_2 = 1;\nmesh_viscosity_3 = 1;'
        default['ale_method']  = 'off'

        return default


    def __getMobileMeshNode(self):
        """
        Return ALE node, creating it if needed.
        """
        node_models = self.case.xmlInitNode('thermophysical_models')
        node_ale = node_models.xmlInitChildNode('ale_method')

        return node_ale


    def __getMobileMeshNodeTry(self):
        """
        Return ALE node if present
        """
        node_ale = None
        node_models = self.case.xmlGetNode('thermophysical_models')
        if node_models:
            node_ale = node_models.xmlGetChildNode('ale_method')

        return node_ale


    def isMobileMeshCompatible(self):
        """
        Indicate if the ALE method may be used
        """
        compat = True
        node_pm = self.case.xmlGetNode('thermophysical_models')
        if node_pm:
            node = node_pm.xmlGetNode('groundwater_model',  'model')
            if node and node['model'] != 'off':
                compat = False
        if self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            compat = False

        return compat


    def __setVariablesandProperties(self):
        """
        Set variables and properties if ALE method is activated.
        """
        node_ale = self.__getMobileMeshNode()
        node_ale.xmlInitChildNode('variable', name='mesh_velocity', label='Mesh Velocity', dimension=3)
        node_ale.xmlInitChildNode('property', name='mesh_viscosity_1', label='mesh_vi1')

        # find node for property / choice and set to default value if require
        node = node_ale.xmlGetChildNode('property', name='mesh_viscosity_1')

        self.__updateNodeViscosity()


    def __updateNodeViscosity(self):
        """
        Update properties beyond mesh visosity is isotrope or not.
        """
        node_ale = self.__getMobileMeshNode()
        if self.getViscosity() == 'orthotrop':
            node_ale.xmlInitChildNode('property', name='mesh_viscosity_2', label='mesh_vi2')
            node_ale.xmlInitChildNode('property', name='mesh_viscosity_3', label='mesh_vi3')

            # find choice for mesh_viscosity_1, create default value if require
            node = node_ale.xmlGetChildNode('property', name='mesh_viscosity_1')

            # Syncrhonize other properties
            for n in ('mesh_viscosity_2', 'mesh_viscosity_3'):
                node = node_ale.xmlGetChildNode('property', name=n)

        else:
            node1 = node_ale.xmlGetChildNode('property', name='mesh_viscosity_2')
            node2 = node_ale.xmlGetChildNode('property', name='mesh_viscosity_3')

            if node1:
                node1.xmlRemoveNode()
            if node2:
                node2.xmlRemoveNode()


    def __removeVariablesandProperties(self):
        """
        Remove variables and properties if ALE method is disabled.
        """
        node_ale = self.__getMobileMeshNodeTry()
        if node_ale:
            node_ale.xmlRemoveChild('variable')
            node_ale.xmlRemoveChild('property')


    @Variables.noUndo
    def getMethod(self):
        """
        Get status on tag "ALE" from xml file
        """
        status = None
        node_ale = self.__getMobileMeshNodeTry()
        if node_ale:
            status = node_ale['status']
        if not status:
            status = self.__defaultInitialValues()['ale_method']
        return status


    @Variables.undoGlobal
    def setMethod(self, status):
        """
        Put status on balise "ALE" in xml file
        """
        self.isOnOff(status)

        node_out = OutputControlModel(self.case)
        typ = node_out.getWriterTimeDependency("-1")
        node_ale = self.__getMobileMeshNode()
        old_status = node_ale['status']
        node_ale['status'] = status
        if status == 'on':
            if typ == 'fixed_mesh':
                typ = 'transient_coordinates'
            self.__setVariablesandProperties()
            if status and old_status != status:
                node_out.setWriterTimeDependency("-1", typ)
        else:
            if typ == 'transient_coordinates':
                typ = 'fixed_mesh'
            self.__removeVariablesandProperties()

            for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                if zone.getNature() == "free_surface":
                    Boundary("free_surface", zone.getLabel(), self.case).deleteFreeSurface()
                    LocalizationModel('BoundaryZone', self.case).deleteZone(zone.getLabel())
            if old_status and old_status != status:
                node_out.setWriterTimeDependency("-1", typ)


    @Variables.undoLocal
    def setSubIterations(self, value):
        """
        Set value of fluid initialization sub iterations into xml file.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)

        node_ale = self.__getMobileMeshNode()
        node_ale.xmlSetData('fluid_initialization_sub_iterations', value)


    @Variables.noUndo
    def getSubIterations(self):
        """
        Get value of fluid initialization sub iterations from xml file.
        """
        nalinf = None
        node_ale = self.__getMobileMeshNodeTry()
        if node_ale:
            nalinf = node_ale.xmlGetInt('fluid_initialization_sub_iterations')
        if not nalinf:
            nalinf = self.__defaultInitialValues()['nalinf']
        return nalinf


    @Variables.undoGlobal
    def setViscosity(self, value):
        """
        Set value of mesh viscosity into xml file.
        """
        self.isInList(value, ['isotrop', 'orthotrop'])
        node_ale = self.__getMobileMeshNode()
        node = node_ale.xmlInitChildNode('mesh_viscosity')
        node['type'] = value
        self.__updateNodeViscosity()


    @Variables.noUndo
    def getViscosity(self):
        """
        Get value of mesh viscosity from xml file.
        """
        iortvm = self.__defaultInitialValues()['iortvm']
        node_ale = self.__getMobileMeshNodeTry()
        if node_ale:
            node = node_ale.xmlGetChildNode('mesh_viscosity')
        if node :
            iortvm = node['type']

        return iortvm


    @Variables.undoLocal
    def setFormula(self, value):
        """
        Set the formula for the viscosity of mesh
        """
        node_ale = self.__getMobileMeshNode()
        node_ale.xmlSetData('formula', value)


    @Variables.noUndo
    def getFormula(self):
        """
        Get the formula for the viscosity of mesh
        """
        formula = None
        node_ale = self.__getMobileMeshNodeTry()
        if node_ale:
            formula = node_ale.xmlGetString('formula')
        if not formula:
            formula = self.getDefaultFormula()
        return formula


    @Variables.noUndo
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
                    <variable label="Mesh Velocity" name="mesh_velocity"/>
                    <property label="mesh_vi1" name="mesh_viscosity_1"/>
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
                    <variable label="Mesh Velocity" name="mesh_velocity"/>
                    <property label="mesh_vi1" name="mesh_viscosity_1"/>
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
                    <variable label="Mesh Velocity" name="mesh_velocity"/>
                    <property label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="orthotrop"/>
                    <property label="mesh_vi2" name="mesh_viscosity_2"/>
                    <property label="mesh_vi3" name="mesh_viscosity_3"/>
                </ale_method>"""
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set mobil mesh model visocity type'
        assert mdl.getViscosity() == 'orthotrop',\
            'Could not get mobil mesh model viscosity type'


    def checkGetAndSetFormula(self):
        """Check whether the MobileMeshModel class could be set and get formula"""
        mdl = MobileMeshModel(self.case)
        mdl.setMethod('on')
        mdl.setFormula('mesh_viscosity_1 = 1000;')

        doc = """<ale_method status="on">
                    <variable label="Mesh Velocity" name="mesh_velocity"/>
                    <property label="mesh_vi1" name="mesh_viscosity_1"/>
                    <mesh_viscosity type="isotrop"/>
                    <formula>mesh_viscosity_1 = 1000;</formula>
                    </ale_method> """
        assert mdl.node_ale == self.xmlNodeFromString(doc),\
            'Could not set formula'
        assert mdl.getFormula() == 'mesh_viscosity_1 = 1000;',\
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
