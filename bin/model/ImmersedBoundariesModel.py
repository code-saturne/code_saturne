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
This module defines the Immersed boundaries model data management.

This module contains the following classes and function:
- ImmersedBoundariesModel
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Cathare coupling model class
#-------------------------------------------------------------------------------

class ImmersedBoundariesModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc
    """

    # ----------------------------------
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.__node_models = self.case.xmlGetNode('thermophysical_models')
        self.__node_ibm    = self.__node_models.xmlInitNode('immersed_boundaries')
    # ----------------------------------

    # ----------------------------------
    def defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['OnOff']            = 'off'
        default['fsi_object_name']  = "fsi_object"
        default['fsi_moving']       = "non_moving"
        default['fsi_interaction']  = "off"
        default['method']           = "explicit"
        default['object_density']   = 0.0
        default['object_stiffness'] = 0.0
        default['object_damping']   = 0.0
        default['object_init']      = "0.0"
        return default
    # ----------------------------------

    # ----------------------------------
    def getNumberOfFSIObjects(self):

        return len(self.__node_ibm.xmlGetNodeList('ibm_object'))
    # ----------------------------------

    # ----------------------------------
    def setOnOff(self, state):

        self.__node_ibm.xmlSetData('ibm_state', state)


    def getOnOff(self):

        state = self.__node_ibm.xmlGetString('ibm_state')
        if state == None or state == '':
            state = self.defaultValues()['OnOff']

        return state
    # ----------------------------------

    # ----------------------------------
    def setMethod(self, method):

        self.__node_ibm.xmlSetData('ibm_method', method)


    def getMethod(self):

        method = self.__node_ibm.xmlGetString('ibm_method')
        if method == None or method == '':
            method = self.defaultValues()['method']

        return method
    # ----------------------------------

    #------------------------------------------------------------------
    # Helper function
    #------------------------------------------------------------------
    def __getStringData(self, index, name, setFunction):
        """
        Get string value from xml file.
        """
        self.isLowerOrEqual(index+1, self.getNumberOfFSIObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[index]
        value = node.xmlGetString(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getIntData(self, index, name, setFunction):
        """
        Get int value from xml file.
        """
        self.isLowerOrEqual(index+1, self.getNumberOfFSIObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[index]
        value = node.xmlGetInt(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getDefaultDataIfNone(self, index, value, name, setFunction):
        """
        Get default value if value is none.
        """
        if value == None or value == "":
            value = self.defaultValues()[name]
            setFunction(index+1, value)
        return value

    # ----------------------------------
    @Variables.undoGlobal
    def addFSIObject(self, name, motion, interaction):

        num = self.getNumberOfFSIObjects()

        node_new = self.__node_ibm.xmlAddChild('ibm_object')

        num += 1

        self.setObjectName(num, name)
        self.setObjectMotion(num, motion)
        self.setObjectInteraction(num, interaction)

        return num


    @Variables.undoLocal
    def deleteFSIObject(self, num):

        self.isLowerOrEqual(num, self.getNumberOfFSIObjects())
        node_list = self.__node_ibm.xmlGetNodeList('ibm_object')
        node = node_list[num-1]
        node.xmlRemoveNode()
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectName(self, num, name):

        self.isLowerOrEqual(num, self.getNumberOfFSIObjects())
        self.isStr(name)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_name', name)


    @Variables.noUndo
    def getObjectName(self, num):

        return self.__getStringData(num-1, 'object_name',
                                    self.setObjectName)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectMotion(self, num, motion):

        self.isLowerOrEqual(num, self.getNumberOfFSIObjects())
        self.isStr(motion)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_motion', motion)

        if motion == 'non_moving':
            self.setObjectInteraction(num, "Off")


    @Variables.noUndo
    def getObjectMotion(self, num):

        return self.__getStringData(num-1, 'object_motion',
                                    self.setObjectMotion)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInteraction(self, num, interaction):

        self.isLowerOrEqual(num, self.getNumberOfFSIObjects())
        self.isStr(interaction)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_interaction', interaction)


    @Variables.noUndo
    def getObjectInteraction(self, num):

        return self.__getStringData(num-1, 'object_interaction',
                                    self.setObjectInteraction)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectDensity(self, num, density):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_density', density)


    @Variables.noUndo
    def getObjectDensity(self, num):

        return self.__getStringData(num-1, 'object_density',
                                    self.setObjectDensity)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectStiffness(self, num, stiffness):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_stiffness', stiffness)


    @Variables.noUndo
    def getObjectStiffness(self, num):

        return self.__getStringData(num-1, 'object_density',
                                    self.setObjectStiffness)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectDamping(self, num, damping):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_damping', damping)


    @Variables.noUndo
    def getObjectDamping(self, num):

        return self.__getStringData(num-1, 'object_damping',
                                    self.setObjectDamping)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectEqPosition(self, num, xeq=None, yeq=None, zeq=None):

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_eq = node.xmlInitNode('equilibrium_position')
        if xeq:
            node_eq.xmlSetData('xeq',xeq)
        if yeq:
            node_eq.xmlSetData('yeq',yeq)
        if zeq:
            node_eq.xmlSetData('zeq',zeq)


    @Variables.noUndo
    def getObjectEqPosition(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_eq = node.xmlGetChildNode('equilibrium_position')

        if node_eq == None:
            self.setObjectEqPosition(num,
                                     self.defaultValues()['object_init'],
                                     self.defaultValues()['object_init'],
                                     self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']

        else:
            return node_eq.xmlGetString('xeq'), \
                   node_eq.xmlGetString('yeq'), \
                   node_eq.xmlGetString('zeq')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitPosition(self, num, xini=None, yini=None, zini=None):

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ini = node.xmlInitNode('initial_position')
        if xini:
            node_ini.xmlSetData('xini', xini)
        if yini:
            node_ini.xmlSetData('yini', yini)
        if zini:
            node_ini.xmlSetData('zini', zini)


    @Variables.noUndo
    def getObjectInitPosition(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ini = node.xmlGetChildNode('initial_position')

        if node_ini == None:
            self.setObjectInitPosition(num,
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_ini.xmlGetString('xini'), \
                   node_ini.xmlGetString('yini'), \
                   node_ini.xmlGetString('zini')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitVel(self, num, vx=None, vy=None, vz=None):

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlInitNode('initial_velocity')
        if vx:
            node_vel.xmlSetData('vx', vx)
        if vy:
            node_vel.xmlSetData('vy', vy)
        if vz:
            node_vel.xmlSetData('vz', vz)


    @Variables.noUndo
    def getObjectInitVel(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlGetChildNode('initial_velocity')

        if node_vel == None:
            self.setObjectInitVel(num,
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_vel.xmlGetString('vx'), \
                   node_vel.xmlGetString('vy'), \
                   node_vel.xmlGetString('vz')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitAcc(self, num, ax=None, ay=None, az=None):

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_acc = node.xmlInitNode('initial_acceleration')
        if ax:
            node_acc.xmlSetData('ax', ax)
        if ay:
            node_acc.xmlSetData('ay', ay)
        if az:
            node_acc.xmlSetData('az', az)


    @Variables.noUndo
    def getObjectInitAcc(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlGetChildNode('initial_acceleration')

        if node_vel == None:
            self.setObjectInitAcc(num,
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']

        else:
            return node_vel.xmlGetString('ax'), \
                   node_vel.xmlGetString('ay'), \
                   node_vel.xmlGetString('az')
    # ----------------------------------

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
