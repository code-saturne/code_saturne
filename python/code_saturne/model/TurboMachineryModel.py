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
- TurboMachineryModel
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
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Class TurboMachinery Model
#-------------------------------------------------------------------------------

class TurboMachineryModel(Variables, Model):

    """
    This class manages the turbomachinery objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case = case
        self.node_therm   = self.case.xmlGetNode('thermophysical_models')
        self.node_turbo   = self.node_therm.xmlInitNode('turbomachinery')
        self.node_join    = self.node_turbo.xmlInitNode('joining')

    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['selector']       = "all[]"
        defvalue['fraction']       = 0.1
        defvalue['plane']          = 25.0
        defvalue['verbosity']      = 0
        defvalue['visualization']  = 0
        defvalue['angle']          = 0.01
        defvalue['model']          = 'off'
        defvalue['velocity']       = 1.
        defvalue['invariant']      = 0.
        defvalue['axis']           = 0.

        return defvalue


    def _getJoinNode(self, join_id):
        """
        Get node for a given joining
        """
        node = None
        listNode = self.node_join.xmlGetNodeList('face_joining')
        if join_id < len(listNode):
            node = listNode[join_id]

        return node

    def _updateJoinSelectionNumbers(self):
        """
        Update names of join selection
        """
        listNode = self.node_join.xmlGetNodeList('face_joining')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    def _addJoinSelect(self, node, select):
        """
        Private method: Add faces to node join with dictionary select.
        """
        for sel, txt in [ (select['selector'],  'selector'),
                          (select['fraction'],  'fraction'),
                          (select['plane'],     'plane'),
                          (select['verbosity'], 'verbosity'),
                          (select['visualization'], 'visualization')]:
            if sel:
                node.xmlSetData(txt, sel)
            else:
                node.xmlRemoveChild(txt)


    def _getFaces(self, node):
        """
        Private method: Return values found for joining for a given node
        """
        default = {}
        default['selector'] =""
        default['fraction'] = ""
        default['plane'] = ""
        default['verbosity'] = ""
        default['visualization'] = ""

        if node:
            default['selector']  = node.xmlGetString('selector')
            default['fraction']  = node.xmlGetString('fraction')
            default['plane']     = node.xmlGetString('plane')
            default['verbosity'] = node.xmlGetString('verbosity')
            default['visualization'] = node.xmlGetString('visualization')
            if not default['selector']:
                default['selector'] = self.defaultValues()['selector']
            if not default['fraction']:
                default['fraction'] = self.defaultValues()['fraction']
            if not default['plane']:
                default['plane'] = self.defaultValues()['plane']
            if not default['verbosity']:
                default['verbosity'] = self.defaultValues()['verbosity']
            if not default['visualization']:
                default['visualization'] = self.defaultValues()['visualization']

        else:
            default = {}

        return default


    def _removeJoinChildren(self, node):
        """
        Private method: Remove all child nodes of node for one selection
        """
        for tag in ('selector',
                    'fraction',
                    'plane',
                    'verbosity',
                    'visualization',):
            node.xmlRemoveChild(tag)


    @Variables.noUndo
    def getJoinSelectionsCount(self):
        """
        Public method.

        @return: number of join faces selections
        @rtype: C{int}
        """
        return len(self.node_join.xmlGetNodeList('face_joining'))


    @Variables.undoGlobal
    def addJoinFaces(self, select):
        """
        Add faces selection for face joining.
        Select is a dictionary with 'selector', 'fraction', 'plane', 'verbosity', 'visualization'
        """
        nb = self.getJoinSelectionsCount()
        name = str(nb +1)
        node = self.node_join.xmlAddChild('face_joining', name=name)
        self._addJoinSelect(node, select)


    @Variables.noUndo
    def getJoinFaces(self, join_id):
        """
        Return faces selection named 'number' for face joining .
        """
        node = self._getJoinNode(join_id)
        return self._getFaces(node)


    @Variables.undoGlobal
    def replaceJoinFaces(self, join_id, select):
        """
        Replace values of faces selection named 'number' for face joining, by select
        """
        node = self._getJoinNode(join_id)
        self._removeJoinChildren(node)
        self._addJoinSelect(node, select)


    @Variables.undoGlobal
    def deleteJoinFaces(self, join_id):
        """
        Delete faces selection named 'number' for face joining
        """
        node = self._getJoinNode(join_id)
        node.xmlRemoveNode()
        if join_id < self.getJoinSelectionsCount():
            self._updateJoinSelectionNumbers()


    @Variables.undoGlobal
    def setTurboMachineryModel(self, model) :
        """
        Put model for turbo machinery
        """
        self.isInList(model, ('off', 'frozen', 'transient',
                              'frozen_coupled', 'transient_coupled'))
        self.node_turbo.xmlSetAttribute(model = model)

        if model == 'off':
            for node in self.getRotorList():
                node.xmlRemoveNode()
        else :
            if len(self.getRotorList()) == 0:
                self.addRotor()


    @Variables.noUndo
    def getTurboMachineryModel(self) :
        """
        Get model for turbo machinery
        """
        model = self.node_turbo['model']
        if model is None :
            model = self.defaultValues()['model']
            self.setTurboMachineryModel(model)

        return model


    @Variables.undoGlobal
    def addRotor(self):
        """
        Add a rotor
        """
        rotor_id = 0
        if self.getRotorList() != None:
            rotor_id = len(self.getRotorList())

        self.node_turbo.xmlInitNode('rotor', rotor_id = rotor_id)
        self.getRotorVelocity(rotor_id)
        self.getRotorCriteria(rotor_id)


    @Variables.undoGlobal
    def delRotor(self, rotor_id):
        """
        Suppress a rotor
        """
        lst = self.getRotorList()
        n = lst[rotor_id]
        for node in lst:
            nodeid = int(node['rotor_id'])
            if nodeid > int(rotor_id):
                node['rotor_id'] = str(nodeid-1)
        n.xmlRemoveNode()


    def _getRotorNode(self, rotor_id):
        """
        Get node for a given rotor
        """
        if rotor_id < len(self.getRotorList()):
            node = self.getRotorList()[rotor_id]

        return node


    def getRotorList(self):
        """
        Get rotor list
        """
        listNode = self.node_turbo.xmlGetNodeList('rotor')

        return listNode


    @Variables.noUndo
    def getRotorVelocity(self, rotor_id):
        """
        return values for rotor velocity
        """
        node = self.node_turbo.xmlGetNode('rotor', rotor_id = rotor_id)
        n = node.xmlGetNode('velocity')
        if not n or n == "None":
            vel = self.defaultValues()['velocity']
            self.setRotorVelocity(rotor_id, vel)
            n = node.xmlGetNode('velocity')
        vel = n.xmlGetDouble('value')

        return vel


    @Variables.undoLocal
    def setRotorVelocity(self, rotor_id, velocity):
        """
        Put values for rotor velocity
        """
        self.isFloat(velocity)
        node = self.node_turbo.xmlGetNode('rotor', rotor_id = rotor_id)
        n = node.xmlInitChildNode('velocity', status="constant")
        n.xmlSetData('value', velocity)


    @Variables.noUndo
    def getRotorCriteria(self, rotor_id):
        """
        return values for rotor criteria
        """
        node = self.node_turbo.xmlGetNode('rotor', rotor_id = rotor_id)
        n = node.xmlGetNode('criteria')
        if not n or n == "None":
            crit = self.defaultValues()['selector']
            self.setRotorCriteria(rotor_id, crit)
            n = node.xmlGetNode('criteria')
        criteria = str(n.xmlGetTextNode())

        return criteria


    @Variables.undoLocal
    def setRotorCriteria(self, rotor_id, criteria):
        """
        Put values for rotor criteria
        """
        node = self.node_turbo.xmlGetNode('rotor', rotor_id = rotor_id)
        node.xmlSetData('criteria', criteria)


    @Variables.noUndo
    def getRotationDirection(self, rotor_id):
        """
        Get values for director vector rotation for periodic translation
        """
        node = self._getRotorNode(rotor_id)

        n = node.xmlInitChildNode('rotation')
        rx = n.xmlGetString('axis_x')
        if not rx:
            rx = self.defaultValues()['axis']
        ry = n.xmlGetString('axis_y')
        if not ry:
            ry = self.defaultValues()['axis']
        rz = n.xmlGetString('axis_z')
        if not rz:
            rz = self.defaultValues()['axis']

        return rx, ry, rz


    @Variables.undoLocal
    def setRotationVector(self, rotor_id, dir, valcoor):
        """
        Put values for director vector rotation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(dir, ("axis_x", "axis_y", "axis_z"))

        node = self._getRotorNode(rotor_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData(dir,valcoor)


    @Variables.noUndo
    def getRotationCenter(self, rotor_id):
        """
        Get coordinates of center of rotation for periodic transformation
        """
        node = self._getRotorNode(rotor_id)

        n = node.xmlInitChildNode('rotation')
        px = n.xmlGetString('invariant_x')
        if not px:
            px = self.defaultValues()['invariant']
        py = n.xmlGetString('invariant_y')
        if not py:
            py = self.defaultValues()['invariant']
        pz = n.xmlGetString('invariant_z')
        if not pz:
            pz = self.defaultValues()['invariant']

        return px, py, pz


    @Variables.undoGlobal
    def setRotationCenter(self, rotor_id, pos, val):
        """
        Put coordinates of center of rotation for periodic transformation
        """
        self.isFloat(val)
        self.isInList(pos, ('invariant_x', 'invariant_y', 'invariant_z'))

        node = self._getRotorNode(rotor_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData(pos, val)


#-------------------------------------------------------------------------------
# TurboMachinery test case
#-------------------------------------------------------------------------------
class TurboMachineryTestCase(ModelTest):
    """
    """


def suite():
    testSuite = unittest.makeSuite(TurboMachineryTestCase, "check")
    return testSuite


def runTest():
    print("TurboMachineryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
