# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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
This module defines the thermal scalar management.

This module contains the following classes and function:
- ThermalScalarModel
- ThermalScalarTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest
import sys

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import ModelTest
from Base.XMLvariables import Variables, Model
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.ThermalRadiationModel import ThermalRadiationModel
from Pages.ConjugateHeatTransferModel import ConjugateHeatTransferModel

#-------------------------------------------------------------------------------
# Thermal scalar model class
#-------------------------------------------------------------------------------

class ThermalScalarModel(DefineUserScalarsModel, Variables, Model):
    """
    The thermal scalar could be temperature
    (Celsius degrees or Kelvin) either enthalpy.
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        DefineUserScalarsModel.__init__(self, case)

        self.node_models   = self.case.xmlGetNode('thermophysical_models')
        self.node_therm    = self.node_models.xmlInitChildNode('thermal_scalar', 'model')
        self.node_prop     = self.case.xmlGetNode('physical_properties').xmlGetNode('fluid_properties')

        self.node_coal  = self.node_models.xmlGetChildNode('solid_fuels',        'model')
        self.node_joule = self.node_models.xmlGetChildNode('joule_effect',       'model')
        self.node_gas   = self.node_models.xmlGetChildNode('gas_combustion',     'model')
        self.node_ray   = self.node_models.xmlGetChildNode('radiative_transfer', 'model')
        self.node_atmo  = self.node_models.xmlGetChildNode('atmospheric_flows',  'model')
        self.node_comp  = self.node_models.xmlGetChildNode('compressible_model',  'model')

        self.old_scaTh = "off"

        self.thermalModel = ('off',
                             'temperature_celsius',
                             'temperature_kelvin',
                             'enthalpy',
                             'potential_temperature',
                             'liquid_potential_temperature',
                             'total_energy')


    def _defaultThermalScalarValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['thermal_scalar'] = "off"

        return default


    def _setNewThermalScalar(self, thermal_scalar):
        """
        Put a new thermal scalar name on to the XML document and its dependances
        """
        self.isInList(thermal_scalar, self.thermalModel)
        node = self.node_therm.xmlInitNode('scalar', name=thermal_scalar, type='thermal')

        if not node['label']:
            lab_def = Tool.dicoLabel(node['name'])
            node['label'] = lab_def
        else:
            lab_def = node['label']

        self.setNewFluidProperty(self.node_prop, 'thermal_conductivity')

        self.setScalarBoundaries()


    def _removeThermalTimeStep(self):
        """
        Private method : remove node 'thermal_time_step' in time_parameters
        we call function from TimeStepModel
        """
        from Pages.TimeStepModel import TimeStepModel
        TimeStepModel(self.case).RemoveThermalTimeStepNode()
        del TimeStepModel


    def thermalScalarModelsList(self):
        """
        Create a tuple with the thermal scalar allowed by the calculation
        features (multi-phases model, and reactive flow models).
        """
        thermalScalarList = self.thermalModel

        for node in (self.node_gas, self.node_coal, self.node_joule, self.node_atmo, self.node_comp):
            if node['model'] != 'off':
                thermalScalarList = ('off',)

        return thermalScalarList


    def isSpecificPhysicActiv(self):
        """
        """
        spec = False
        for node in (self.node_gas, self.node_coal, self.node_joule, self.node_atmo, self.node_comp):
            if node != None:
                if node['model'] != 'off':
                    spec = True

        return spec


    @Variables.undoGlobal
    def setThermalModel(self, thermal_scalar):
        """
        Update the thermal model and create the thermal scalar markup from the XML document.
        """
        self.isInList(thermal_scalar, self.thermalModel)

        self.node_therm['model'] = thermal_scalar

        if thermal_scalar != 'off':
            node = self.scalar_node.xmlGetNode('scalar', type='thermal')
            if node:
                if node['name'] != thermal_scalar:
                    self.deleteScalar(node['label'])
            self._setNewThermalScalar(thermal_scalar)

        else:
            node = self.scalar_node.xmlGetNode('scalar', type='thermal')
            if node:
                self.deleteScalar(node['label'])
            self._removeThermalTimeStep()
            ThermalRadiationModel(self.case).setRadiativeModel('off')
            ConjugateHeatTransferModel(self.case).deleteConjugateHeatTransfer()

        t_outputs = (("tplus", "Tplus", False),
                     ("input_thermal_flux", "Thermal flux", True),
                     ("boundary_temperature", "Boundary temperature", True),
                     ("boundary_layer_nusselt", "Boundary layer Nusselt", False))

        if thermal_scalar != 'off' or self.isSpecificPhysicActiv():
            for v in t_outputs:
                n = self.node_therm.xmlInitChildNode('property',
                                                     name=v[0],
                                                     support="boundary")
                if not n['label']:
                    n['label'] = v[1]
                if not v[2]:
                    n.xmlInitNode('postprocessing_recording')['status']= "off"

        else:
            for v in t_outputs:
                self.node_therm.xmlRemoveChild('property',
                                               name=v[0],
                                               support="boundary")

    @Variables.noUndo
    def getThermalScalarModel(self):
        """
        Get name of thermal scalar (not label)
        """
        model = self.node_therm['model']
        if not model:
            model = self._defaultThermalScalarValues()['thermal_scalar']
            self.setThermalModel(model)
        return model


    @Variables.noUndo
    def getThermalScalarLabel(self):
        """
        Get label for thermal scalar
        """
        label = ""
        node = self.node_therm.xmlGetNode('scalar', type='thermal')
        if node:
            label = node['label']

        return label


#-------------------------------------------------------------------------------
# ThermalScalar Model test case
#-------------------------------------------------------------------------------


class ThermalScalarTestCase(ModelTest):
    """
    """
    def checkThermalInstantiation(self):
        """
        Check whether the ThermalScalarModel class could be instantiated
        """
        model = None
        model = ThermalScalarModel(self.case)
        assert model != None, 'Could not instantiate ThermalScalarModel'

    def checkThermalScalarModelsList(self):
        """
        Check whether the ThermalScalarModelsList could be get
        """
        mdl = ThermalScalarModel(self.case)
        mdl.node_gas['model'] = 'on'
        assert mdl.getThermalScalarModel() == ('off'),\
            'Could not use the thermalScalarModelsList method'

    def checkUpdateandGetThermalScalarModel(self):
        """
        Check whether a new thermal scalar could be set and get
        """
        mdl = ThermalScalarModel(self.case)
        mdl.setThermalModel('temperature_kelvin')
        doc = '''<thermal_scalar>
                    <scalar label="TempK" name="temperature_kelvin" type="thermal">
                        <initial_value zone_id="1">293.15</initial_value>
                        <min_value>-1e+12</min_value>
                        <max_value>1e+12</max_value>
                    </scalar>
                 </thermal_scalar>'''

        assert mdl.scalar_node == self.xmlNodeFromString(doc), \
            'Could not update thermal scalar in ThermalScalarModel'
        assert mdl.getThermalScalarModel() == 'temperature_kelvin', \
            'Could not get the thermal scalar'


def suite():
    testSuite = unittest.makeSuite(ThermalScalarTestCase, "check")
    return testSuite


def runTest():
    print("ThermalScalarTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
