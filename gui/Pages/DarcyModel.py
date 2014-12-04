# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
- DarcyModel
- DarcyTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common                     import *
from code_saturne.Base.XMLvariables               import Variables, Model
from code_saturne.Base.XMLmodel                   import ModelTest
from code_saturne.Pages.TurbulenceModel           import TurbulenceModel
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel

#-------------------------------------------------------------------------------
# Mobil Mesh model class
#-------------------------------------------------------------------------------

class DarcyModel(Model):
    """
    Manage the input/output markups in the xml doc about mobil mesh
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_darcy = self.node_models.xmlInitChildNode('darcy_model', 'model')


    def __defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['permeability' ]    = 'isotropic'
        default['diffusion' ]       = 'isotropic'
        default['flow' ]            = 'steady'
        default['criterion' ]       = 'pressure'
        default['darcy_model']      = 'off'
        default['gravity']          = 'off'
        default['axis']             = '0.'

        return default


    @Variables.noUndo
    def getDarcyModel(self):
        """
        Get the darcy model
        """
        mdl = self.node_darcy['model']
        if mdl == "":
            mdl = self.__defaultValues()['darcy_model']
            self.setDarcyModel(mdl)
        return mdl


    @Variables.undoLocal
    def setDarcyModel(self, choice):
        """
        Put the darcy model
        """
        self.isInList(choice, ['off', 'darcy'])
        self.node_darcy['model'] = choice

        ### TODO remove node when set off
        if choice != "off":
            TurbulenceModel(self.case).setTurbulenceModel('off')
            FluidCharacteristicsModel(self.case).setPropertyMode('density', 'constant')
            FluidCharacteristicsModel(self.case).setInitialValue('density', 1.)


    @Variables.noUndo
    def getPermeabilityType(self):
        """
        Get the permeability model
        """
        node = self.node_darcy.xmlInitChildNode('permeability')
        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['permeability']
            self.setPermeabilityType(mdl)
        return mdl


    @Variables.undoLocal
    def setPermeabilityType(self, choice):
        """
        Put the permeability model
        """
        self.isInList(choice, ['isotropic', 'anisotropic'])
        node = self.node_darcy.xmlInitChildNode('permeability')
        oldchoice = node['model']

        node['model'] = choice

        if oldchoice != None and oldchoice != choice:
            node.xmlRemoveChild('formula')


    @Variables.noUndo
    def getDiffusionType(self):
        """
        Get the diffusion model
        """
        node = self.node_darcy.xmlInitChildNode('diffusion')
        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['diffusion']
            self.setDiffusionType(mdl)
        return mdl


    @Variables.undoLocal
    def setDiffusionType(self, choice):
        """
        Put the diffusion model
        """
        self.isInList(choice, ['isotropic', 'anisotropic'])
        node = self.node_darcy.xmlInitChildNode('diffusion')
        node['model'] = choice


    @Variables.noUndo
    def getFlowType(self):
        """
        Get flow type : steady or unsteady
        """
        node = self.node_darcy.xmlInitChildNode('flowType')
        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['flow']
            self.setFlowType(mdl)
        return mdl


    @Variables.undoLocal
    def setFlowType(self, choice):
        """
        Put flow type : steady or unsteady
        """
        self.isInList(choice, ['steady', 'unsteady'])
        node = self.node_darcy.xmlInitChildNode('flowType')
        node['model'] = choice


    @Variables.noUndo
    def getCriterion(self):
        """
        Get convergence criterion of Newton scheme
        """
        node = self.node_darcy.xmlInitChildNode('criterion')
        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['criterion']
            self.setCriterion(mdl)
        return mdl


    @Variables.undoLocal
    def setCriterion(self, choice):
        """
        Put convergence criterion of Newton scheme
        """
        self.isInList(choice, ['pressure', 'velocity'])
        node = self.node_darcy.xmlInitChildNode('criterion')
        node['model'] = choice


    @Variables.noUndo
    def getGravity(self):
        """
        Return if gravity is taken into account
        """
        node = self.node_darcy.xmlInitNode('gravity', 'status')
        status = node['status']
        if not status:
            v = self.__defaultValues()['gravity']
            self.setGravity(v)
        return status


    @Variables.undoLocal
    def setGravity(self, v):
        """
        Put if gravity is taken into account
        """
        self.isOnOff(v)
        node = self.node_darcy.xmlInitNode('gravity', 'status')
        node['status'] = v


    @Variables.noUndo
    def getGravityVector(self):
        """
        Return gravity vector
        """
        node = self.node_darcy.xmlInitNode('gravity')
        n = node.xmlInitChildNode('vector')
        rx = n.xmlGetString('axis_x')
        if not rx:
            rx = self.__defaultValues()['axis']
        ry = n.xmlGetString('axis_y')
        if not ry:
            ry = self.__defaultValues()['axis']
        rz = n.xmlGetString('axis_z')
        if not rz:
            rz = self.__defaultValues()['axis']

        return rx, ry, rz


    @Variables.noUndo
    def setGravityVector(self, direction, valcoor):
        """
        Put values for director gravity vector
        """
        self.isFloat(valcoor)
        self.isInList(direction, ("axis_x", "axis_y", "axis_z"))


        node = self.node_darcy.xmlInitNode('gravity')
        n = node.xmlGetChildNode('vector')
        n.xmlSetData(direction,valcoor)

