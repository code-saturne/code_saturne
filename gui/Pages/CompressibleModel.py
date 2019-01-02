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
This module defines the compressible model management.

This module contains the following classes and function:
- CompressibleModel
- CompressibleTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest
import sys

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Pages.LocalizationModel import LocalizationModel
from code_saturne.Pages.Boundary import Boundary
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel

#-------------------------------------------------------------------------------
# Compressible model class
#-------------------------------------------------------------------------------

class CompressibleModel(Variables, Model):
    """
    Active de compressible model
    There is three different thermodynamic law
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case

        self.node_thermo = self.case.xmlGetNode('thermophysical_models')
        self.node_comp   = self.node_thermo.xmlInitNode('compressible_model')
        self.node_np     = self.case.xmlInitNode('numerical_parameters')
        self.node_prop   = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = self.node_prop.xmlInitNode('fluid_properties')
        self.node_ref    = self.node_thermo.xmlInitNode('reference_values')

        self.comp_choice = ['off', 'constant_gamma', 'variable_gamma', 'van_der_waals']
        self.var_list   = ['temperature']


    def _defaultCompressibleValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['activation'] = "off"

        return default


    @Variables.undoGlobal
    def setCompressibleModel(self, model):
        """
        Active or desactive the compressible model
        Add and remove the variables and properties associated
        """
        self.isInList(model, self.comp_choice)
        oldModel = self.node_comp['model']
        if oldModel != model:
            self.node_comp['model'] = model
            if model == 'off':
                for zone in LocalizationModel('BoundaryZone', self.case).getZones():
                    if zone.getNature() == "outlet":
                        Boundary("compressible_outlet", zone.getLabel(), self.case).deleteCompressibleOutlet()
                    if zone.getNature() == "inlet":
                        Boundary("inlet", zone.getLabel(), self.case).deleteCompressibleInlet()
                self.__removeVariablesAndProperties()
                self.node_np.xmlRemoveChild('hydrostatic_equilibrium')
                self.node_fluid.xmlRemoveChild('property', name = 'volume_viscosity')
                self.node_ref.xmlRemoveChild('mass_molar')
                self.node_ref.xmlRemoveChild('temperature')
                ThermalScalarModel(self.case).setThermalModel('off')
            else :
                ThermalScalarModel(self.case).setThermalModel('total_energy')
                for v in self.var_list:
                    self.setNewVariable(self.node_comp, v, tpe="model", label=v)
                from code_saturne.Pages.TurbulenceModel import TurbulenceModel
                TurbulenceModel(self.case).setTurbulenceModel('off')
                del TurbulenceModel


    @Variables.noUndo
    def getCompressibleModel(self):
        """
        Return the model of the compressible
        """
        node = self.node_thermo.xmlInitNode('compressible_model')
        status = node['model']
        if status not in self.comp_choice:
            status = self._defaultCompressibleValues()['activation']
            self.setCompressibleModel(status)
        return status


    def __removeVariablesAndProperties(self):
        """
        Delete variables and property that are useless accordingly to the model.
        """
        for v in self.var_list:
            self.node_comp.xmlRemoveChild('variable', name=v)


#-------------------------------------------------------------------------------
# TurbulenceModel test case
#-------------------------------------------------------------------------------

class CompressibleModelTestCase(ModelTest):
    """
    """
    def checkCompressibleInstantiation(self):
        """Check whether the Compressible Model class could be instantiated"""
        model = None
        model = CompressibleModel(self.case)
        assert model != None, 'Could not instantiate Compressible Model'


def suite():
    testSuite = unittest.makeSuite(CompressibleModelTestCase, "check")
    return testSuite


def runTest():
    print("CompressibleModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
