# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2015 EDF S.A.
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
- GroundwaterModel
- GroundwaterTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common                     import *
from code_saturne.model.XMLvariables               import Variables, Model
from code_saturne.model.XMLmodel                   import ModelTest
from code_saturne.model.TurbulenceModel           import TurbulenceModel
from code_saturne.model.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# Mobil Mesh model class
#-------------------------------------------------------------------------------

class GroundwaterModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about mobil mesh
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_darcy = self.node_models.xmlInitChildNode('groundwater_model', 'model')

        self.sca_mo       = DefineUserScalarsModel(self.case)


    def __defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['permeability' ]     = 'isotropic'
        default['flow' ]             = 'steady'
        default['groundwater_model'] = 'off'
        default['gravity']           = 'off'
        default['unsaturated']       = 'true'
        default['chemistry_model']   = 'Kd'
        default['fo_decay_rate']     = 0.0
        return default


    @Variables.noUndo
    def getGroundwaterModel(self):
        """
        Get the Groundwater model
        """
        mdl = self.node_darcy['model']
        if mdl == "":
            mdl = self.__defaultValues()['groundwater_model']
            self.setGroundwaterModel(mdl)
        return mdl


    @Variables.undoLocal
    def setGroundwaterModel(self, choice):
        """
        Put the Groundwater model
        """
        self.isInList(choice, ['off', 'groundwater'])
        old_choice = self.node_darcy['model']
        self.node_darcy['model'] = choice

        self.node_models = self.case.xmlInitNode('thermophysical_models')

        if choice != "off":
            TurbulenceModel(self.case).setTurbulenceModel('off')
            FluidCharacteristicsModel(self.case).setPropertyMode('density', 'constant')
            FluidCharacteristicsModel(self.case).setInitialValue('density', 1.)

            node = self.node_models.xmlInitNode('velocity_pressure')
            nn =  node.xmlGetNode('variable', name='pressure')
            if nn:
                nn.xmlRemoveNode()
            self.setNewVariable(node, 'hydraulic_head')
            nn =  node.xmlGetNode('property', name='total_pressure')
            if nn:
                nn.xmlRemoveNode()
            nn = node.xmlGetNode('property', name='yplus')
            if nn:
                nn.xmlRemoveNode()
            nn =  node.xmlGetNode('property', name='stress')
            if nn:
                nn.xmlRemoveNode()
            nn =  node.xmlGetNode('property', name='stress_normal')
            if nn:
                nn.xmlRemoveNode()
            nn =  node.xmlGetNode('property', name='stress_tangential')
            if nn:
                nn.xmlRemoveNode()
            node_control = self.case.xmlGetNode('analysis_control')
            node_time    = node_control.xmlInitNode('time_parameters')
            nn = node_time.xmlGetNode('property', name='courant_number')
            if nn:
                nn.xmlRemoveNode()
            nn = node_time.xmlGetNode('property', name='fourier_number')
            if nn:
                nn.xmlRemoveNode()

        elif old_choice and old_choice != "off":
            TurbulenceModel(self.case).setTurbulenceModel("k-epsilon-PL")
            node = self.node_models.xmlInitNode('velocity_pressure')
            nn =  node.xmlGetNode('variable', name='hydraulic_head')
            if nn:
                nn.xmlRemoveNode()
            self.setNewVariable(node, 'pressure')
            self.setNewProperty(node, 'total_pressure')
            self.setNewProperty(node, 'yplus')
            n = self.setNewProperty(node, 'stress')
            n['support'] = 'boundary'
            n['label'] = 'Stress'
            if not node.xmlGetChildNode('property', name='stress_tangential'):
                n = self.setNewProperty(node, 'stress_tangential')
                n['label'] = 'Stress, tangential'
                n['support'] = 'boundary'
                n.xmlInitNode('postprocessing_recording')['status']= "off"
            if not node.xmlGetChildNode('property', name='stress_normal'):
                n = self.setNewProperty(node, 'stress_normal')
                n['label'] = 'Stress, normal'
                n['support'] = 'boundary'
                n.xmlInitNode('postprocessing_recording')['status']= "off"


    @Variables.noUndo
    def getPermeabilityType(self):
        """
        Get the permeability model
        """
        node = self.node_darcy.xmlInitChildNode('permeability')
        mdl = node['model']
        if mdl is None:
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
    def getFlowType(self):
        """
        Get flow type : steady or unsteady
        """
        node = self.node_darcy.xmlInitChildNode('flowType')
        mdl = node['model']
        if mdl is None:
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
    def getUnsaturatedZone(self):
        """
        Get unsaturated zone status : True or False
        """
        node = self.node_darcy.xmlInitChildNode('unsaturatedZone')
        mdl = node['model']
        if mdl is None:
            mdl = self.__defaultValues()['unsaturated']
            self.setUnsaturatedZone(mdl)
        return mdl


    @Variables.undoLocal
    def setUnsaturatedZone(self, choice):
        """
        Get unsaturated zone status : True or False
        """
        self.isInList(choice, ['true', 'false'])
        node = self.node_darcy.xmlInitChildNode('unsaturatedZone')
        node['model'] = choice


    @Variables.noUndo
    def getGravity(self):
        """
        Return if gravity is taken into account
        """
        node = self.node_darcy.xmlInitChildNode('gravity', 'status')
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
        node = self.node_darcy.xmlInitChildNode('gravity', 'status')
        node['status'] = v


    @Variables.noUndo
    def getDecayRate(self, scalar_name):
        """
        Get value for the first-order decay rate of one scalar
        """
        prop = 'fo_decay_rate'

        self.isNotInList(scalar_name, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_name, self.sca_mo.getUserScalarNameList())

        nodeScalar = self.node_darcy.xmlInitChildNode('scalar', name=scalar_name)
        value = nodeScalar.xmlGetDouble(prop)

        if value is None:
            value = self.__defaultValues()[prop]
            self.setDecayRate(scalar_name, value)

        return value


    @Variables.undoLocal
    def setDecayRate(self, scalar_name, value):
        """
        Set value for the first-order decay rate of one scalar
        """
        self.isNotInList(scalar_name, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_name, self.sca_mo.getUserScalarNameList())

        nodeScalar = self.node_darcy.xmlInitChildNode('scalar', name=scalar_name)
        nodeScalar.xmlSetData('fo_decay_rate', value)


    @Variables.noUndo
    def getChemistryModel(self, scalar_name):
        """
        Get choice of the chemistry model (representing soil-water partition) of one scalar
        """
        prop = 'chemistry_model'

        self.isNotInList(scalar_name, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_name, self.sca_mo.getUserScalarNameList())

        nodeScalar = self.node_darcy.xmlInitChildNode('scalar', name=scalar_name)
        choice = nodeScalar['chemistry_model']

        if choice is None:
            choice = self.__defaultValues()[prop]
            self.setChemistryModel(scalar_name, choice)

        return choice


    @Variables.undoLocal
    def setChemistryModel(self, scalar_name, choice):
        """
        Set choice of the chemistry model of one scalar: Kd or EK
        """
        self.isNotInList(scalar_name, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_name, self.sca_mo.getUserScalarNameList())
        self.isInList(choice, ['Kd', 'EK'])

        nodeScalar = self.node_darcy.xmlInitChildNode('scalar', name=scalar_name)
        nodeScalar['chemistry_model'] = choice
