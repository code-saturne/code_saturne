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
This module defines the homogeneous model management.

This module contains the following classes and function:
- HgnModel
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

#-------------------------------------------------------------------------------
# Homogeneous model class
#-------------------------------------------------------------------------------

class HgnModel(Variables, Model):
    """
    Enable homogeneous mutliphase model
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case

        self.node_thermo = self.case.xmlGetNode('thermophysical_models')
        self.node_hgn   = self.node_thermo.xmlInitNode('hgn_model')
        self.node_np     = self.case.xmlInitNode('numerical_parameters')
        self.node_prop   = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = self.node_prop.xmlInitNode('fluid_properties')
        self.node_ref    = self.node_thermo.xmlInitNode('reference_values')

        self.hgn_choice = ['off', 'no_mass_transfer', 'merkle_model']
        self.var_list   = ['void_fraction']


    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['activation'] = "off"

        return default


    @Variables.undoGlobal
    def setHgnModel(self, model):
        """
        Enable or disable homogeneous model
        Add and remove the variables and properties associated
        """
        self.isInList(model, self.hgn_choice)
        oldModel = self.node_hgn['model']
        self.node_hgn['model'] = model

        if model == 'off':
            self.__removeVariablesAndProperties()
            self.node_np.xmlRemoveChild('hydrostatic_pressure')
        elif oldModel and oldModel != "off":
            for v in self.var_list:
                self.setNewVariable(self.node_hgn, v, tpe="model", label=v)

    @Variables.noUndo
    def getHgnModel(self):
        """
        Return model
        """
        node = self.node_thermo.xmlInitNode('hgn_model')
        status = node['model']
        if status not in self.hgn_choice:
            status = self._defaultValues()['activation']
            self.setHgnModel(status)
        return status


    def __removeVariablesAndProperties(self):
        """
        Remove variables and properties associated to current model.
        """
        for v in self.var_list:
            self.node_hgn.xmlRemoveChild('variable', name=v)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
