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
This module defines the XML calls for preprocessor execution
This module contains the following classes and function:
- FansModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import XMLmodel, ModelTest

#-------------------------------------------------------------------------------
# Class Fans Model
#-------------------------------------------------------------------------------

class FansModel(Variables, Model):

    """
    This class manages the fan objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case = case
        node_model = self.case.xmlInitNode('thermophysical_models')
        self.node_model_fans = node_model.xmlInitNode('fans')


    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['inlet_axis_x']   = 0.
        defvalue['inlet_axis_y']   = 0.
        defvalue['inlet_axis_z']   = 0.
        defvalue['outlet_axis_x']  = 0.1
        defvalue['outlet_axis_y']  = 0.
        defvalue['outlet_axis_z']  = 0.
        defvalue['fan_radius']     = 0.7
        defvalue['blades_radius']  = 0.5
        defvalue['hub_radius']     = 0.1
        defvalue['axial_torque']   = 0.01
        defvalue['curve_coeffs_x'] = 0.6
        defvalue['curve_coeffs_y'] = -0.1
        defvalue['curve_coeffs_z'] = -0.05
        defvalue['mesh_dimension'] = 3

        return defvalue


    @Variables.undoGlobal
    def addFan(self):
        """
        """
        nb = len(self.getFanList())
        node = self.node_model_fans.xmlAddChild('fan', fan_id = nb)


    @Variables.undoGlobal
    def delFan(self, idx):
        """
        """
        self.isInt(idx)
        lst = self.getFanList()
        n = lst[idx]
        for node in lst:
            nodeid = int(node['fan_id'])
            if nodeid > int(idx):
                node['fan_id'] = str(nodeid-1)
        n.xmlRemoveNode()


    @Variables.noUndo
    def getFanList(self):
        """
        """
        listNode = self.node_model_fans.xmlGetNodeList('fan')

        return listNode


    @Variables.noUndo
    def getFanProperty(self, idx, prop):
        """
        """
        self.isInt(idx)
        node = self.node_model_fans.xmlGetNode('fan', fan_id = idx)
        n = node.xmlGetNode(prop)
        if not n or n == "None":
            value = self.defaultValues()[prop]
            self.setFanProperty(idx, prop, value)
        value = node.xmlGetDouble(prop)

        return value


    @Variables.undoLocal
    def setFanProperty(self, idx, prop, val):
        """
        """
        self.isInt(idx)
        node = self.node_model_fans.xmlGetNode('fan', fan_id = idx)
        node.xmlSetData(prop, val)


    @Variables.noUndo
    def getFanMeshDimension(self, idx):
        """
        """
        self.isInt(idx)
        node = self.node_model_fans.xmlGetNode('fan', fan_id = idx)
        n = node.xmlGetNode('mesh_dimension')
        if not n or n == "None":
            value = self.defaultValues()['mesh_dimension']
            self.setFanMeshDimension(idx, value)
        value = node.xmlGetInt('mesh_dimension')

        return value


    @Variables.undoLocal
    def setFanMeshDimension(self, idx, val):
        """
        """
        self.isInt(idx)
        node = self.node_model_fans.xmlGetNode('fan', fan_id = idx)
        node.xmlSetData('mesh_dimension', val)


class FansStatus(Variables, Model):

    """
    This class checks for fan objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case


    def getFanCount(self):
        """
        Returns info on fans: -1 if node is not present, count if present
        """
        count = -1

        node_model = self.case.xmlGetNode('thermophysical_models')
        if node_model:
            node_fans = node_model.xmlGetNode('fans')
            if node_fans:
                count =len(node_fans.xmlGetNodeList('fan'))

        return count


    def initFans(self):
        """
        Ensures a fans node is present
        """
        node_model = self.case.xmlInitNode('thermophysical_models')
        node_model.xmlInitNode('fans')


    def cleanFans(self):
        """
        Removes fan node if empty
        """
        count = -1
        node_model = self.case.xmlGetNode('thermophysical_models')
        if node_model:
            node_fans = node_model.xmlGetNode('fans')
            if node_fans:
                count =len(node_fans.xmlGetNodeList('fan'))

        if count == 0:
            node_fans.xmlRemoveNode()


#-------------------------------------------------------------------------------
# Fans test case
#-------------------------------------------------------------------------------

class FansTestCase(ModelTest):
    """
    """

def suite():
    testSuite = unittest.makeSuite(FansTestCase, "check")
    return testSuite


def runTest():
    print("FansTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
