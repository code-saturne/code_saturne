# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This module defines the cooling towers model management.

This module contains the following classes and function:
- CoolingTowersModel
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
from code_saturne.model.ThermalScalarModel import ThermalScalarModel

#-------------------------------------------------------------------------------
# Cooling Towers model class
#-------------------------------------------------------------------------------

class CoolingTowersModel(Variables, Model):
    """
    Enable homogeneous mutliphase model
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_pmodel = self.case.xmlGetNode('thermophysical_models')
        self.node_ctwr   = self.node_pmodel.xmlGetNode('cooling_towers')

        self.node_np     = self.case.xmlInitNode('numerical_parameters')
        node_prop        = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = node_prop.xmlInitNode('fluid_properties')

        self.cooling_towers_choice = ['off', 'poppe', 'merkel']
        self.var_list   = [('y_p', 'Yp rain'),
                           ('y_p_t_l', 'Yp.Tp rain'),
                           ('y_l_packing', 'Yl packing'),
                           ('enthalpy_liquid', 'Enthalpy liq packing'),
                           ('ym_water', 'Ym water bulk')]
        self.prop_list   = [('humidity', 'Humidity'),
                            ('x_s', 'Humidity sat'),
                            ('enthalpy', 'Enthalpy humid air'),
                            ('temperature_liquid', 'Temperature liq packing'),
                            ('vertvel_l', 'Velocity liq packing'),
                            ('t_rain', 'Temperature rain'),
                            ('x_c', 'Gas mass fraction')]

    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['activation'] = "off"

        return default

    @Variables.undoGlobal
    def setCoolingTowersModel(self, model):
        """
        Enable or disable homogeneous model
        Add and remove the variables and properties associated
        """
        self.isInList(model, self.cooling_towers_choice)
        self.node_ctwr = self.node_pmodel.xmlInitNode('cooling_towers')
        oldModel = self.node_ctwr['model']
        self.node_ctwr['model'] = model

        if model != 'off':
            ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
            for v in self.var_list:
                self.setNewVariable(self.node_ctwr, v[0], tpe="model", label=v[1])
            for v in self.prop_list:
                self.setNewProperty(self.node_ctwr, v[0], label=v[1])
        elif oldModel and oldModel != "off":
            self.__removeVariablesAndProperties()
            self.node_np.xmlRemoveChild('ym_water')
            ThermalScalarModel(self.case).setThermalModel('off')

    @Variables.noUndo
    def getCoolingTowersModel(self):
        """
        Return model
        """
        node = self.node_pmodel.xmlGetNode('cooling_towers')
        if not node:
            return self._defaultValues()['activation']
        status = node['model']
        if status not in self.cooling_towers_choice:
            status = self._defaultValues()['activation']
            self.setCoolingTowersModel(status)
        return status

    def __removeVariablesAndProperties(self):
        """
        Remove variables and properties associated to current model.
        """
        for v in self.var_list:
            self.node_ctwr.xmlRemoveChild('variable', name=v[0])
        for v in self.prop_list:
            self.node_ctwr.xmlRemoveChild('property', name=v[0])

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
