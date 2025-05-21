# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the solid heat transfer model management.

This module contains the following classes and function:
- HTSModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest
import sys

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
#from code_saturne.model.ThermalScalarModel import ThermalScalarModel

#-------------------------------------------------------------------------------
# Homogeneous model class
#-------------------------------------------------------------------------------

class HTSModel(Variables, Model):
    """
    Enable solid heat transfer model
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_thermo = self.case.xmlGetNode('thermophysical_models')
        self.node_vp     = self.node_thermo.xmlGetNode('velocity_pressure')
        self.node_hts    = self.node_thermo.xmlGetNode('hts_model')
        self.node_np     = self.case.xmlInitNode('numerical_parameters')
        self.node_prop   = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = self.node_prop.xmlInitNode('fluid_properties')

        self.hts_choice = ['off', 'hts_conduction']
        self.var_list   = [] #CK : ['temperature']


    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['activation'] = "off"

        return default


    @Variables.undoGlobal
    def setHTSModel(self, model):
        """
        Enable or disable homogeneous model
        Add and remove the variables and properties associated
        """
        self.isInList(model, self.hts_choice)

        if self.node_hts:
            oldModel = self.node_hts['model']
        else:
            oldModel = 'off'

        if model != 'off':
            self.node_hts  = self.node_thermo.xmlInitNode('hts_model')
            self.node_hts['model'] = model
            if self.node_vp:
                self.node_vp.xmlRemoveChild('variable')

                self.node_vp.xmlRemoveChild('property')
            #ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
            for v in self.var_list:
                self.setNewVariable(self.node_hts, v, tpe="model", label=v)
        elif oldModel and oldModel != "off":
            self.__removeVariablesAndProperties()
            #ThermalScalarModel(self.case).setThermalModel('off')
            self.node_np.xmlRemoveChild('hydrostatic_pressure')
            self.node_hts.xmlRemoveNode()
            self.node_hts = None

    @Variables.noUndo
    def getHTSModel(self):
        """
        Return model
        """
        node = self.node_thermo.xmlGetNode('hts_model')
        if node:
            status = node['model']
            if status not in self.hts_choice:
                status = self._defaultValues()['activation']
                self.setHTSModel(status)
        else:
            status = 'off'
        return status


    @Variables.noUndo
    def getHTSName(self):
        """
        Get name for thermal scalar
        """
        name = 'off'
        if self.node_hts:
            node = self.node_hts.xmlGetNode('variable', type='model')
            if node:
                name = node['name']

        return name


    def __removeVariablesAndProperties(self):
        """
        Remove variables and properties associated to current model.
        """
        for v in self.var_list:
            self.node_hts.xmlRemoveChild('variable', name=v)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
