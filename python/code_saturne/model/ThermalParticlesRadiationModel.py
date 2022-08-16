
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
This module defines the Thermal Radiation model

This module contains the following classes and function:
- ThermalParticlesRadiationModel
- ThermalParticlesRadiationTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.XMLvariables import Model, Variables

#-------------------------------------------------------------------------------
# ThermalRadiation model class
#-------------------------------------------------------------------------------

class ThermalParticlesRadiationModel(Variables, Model):

    def __init__(self, case):
        """
        Constructor
        """
        self.case = case

        # Default values
        self._isActivated = "off"
        self._elasticity = "0.9"
        self._emissivity = "1.0"

        tm_node = self.case.xmlGetNode('thermophysical_models')
        self.xml_node = tm_node.xmlInitNode("interparticles_radiative_transfer")

        self.initialize_xml()
        self.read_xml()

    def initialize_xml(self):
        for childNode, value in zip(["status", "elasticity", "emissivity"],
                [self._isActivated, self._elasticity, self._emissivity]):
            if self.xml_node.xmlGetChildNode(childNode) == None:
                self.xml_node.xmlInitChildNode(childNode)
                self.xml_node.xmlSetData(childNode, value)
        return

    def read_xml(self):
        values = []
        self._isActivated = self.xml_node.xmlGetChildString("status")
        self._elasticity= self.xml_node.xmlGetChildString("elasticity")
        self._emissivity= self.xml_node.xmlGetChildString("emissivity")

    @property
    def isActivated(self):
        return self._isActivated

    @property
    def elasticity(self):
        return self._elasticity

    @property
    def emissivity(self):
        return self._emissivity

    @isActivated.setter
    def isActivated(self, status):
        self.xml_node.xmlSetData("status", status)
        self._isActivated = status

    @elasticity.setter
    def elasticity(self, value):
        self.xml_node.xmlSetData("elasticity", value)
        self._elasticity = value

    @emissivity.setter
    def emissivity(self, value):
        self.xml_node.xmlSetData("emissivity", value)
        self._emissivity = value
