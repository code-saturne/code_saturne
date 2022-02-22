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
This modules defines properties for internal coupling.
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import XMLmodel, ModelTest

from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.NotebookModel import NotebookModel

class InternalCouplingModel(Variables, Model):
    """
    Class to manipulate internal coupling parameters in xml file.
    """

    def __init__(self, case):
        """Constructor"""
        self.case = case

        n = self.case.xmlGetNode("thermophysical_models")
        self.node = n.xmlGetChildNode("internal_coupling")
        if not self.node:
            self.node = n.xmlInitChildNode("internal_coupling")

        self.node_zones = self.node.xmlInitChildNode("solid_zones")
        self.node_scals = self.node.xmlInitChildNode("coupled_scalars")


    def addZone(self, zone_name):
        """Add a coupled zone"""

        self.node_zones.xmlInitChildNode("zone", name=zone_name)


    def removeZone(self, zone_name):
        """Remove a zone from the coupled list"""

        n = self.node_zones.xmlGetChildNode("zone", name=zone_name)
        if n:
            n.xmlRemoveNode()


    def getZonesList(self):
        """Get list of coupled zones."""

        lst = []
        for n in self.node_zones.xmlGetChildNodeList("zone"):
            lst.append(n['name'])

        return lst


    def addScalar(self, scalar_name):
        """Add a coupled scalar"""

        self.node_scals.xmlInitChildNode("scalar", name=scalar_name)


    def removeScalar(self, scalar_name):
        """Remove a scalar from the coupled list"""
        n = self.node_scals.xmlGetChildNode("scalar", name=scalar_name)
        if n:
            n.xmlRemoveNode()


    def getScalarsList(self):
        """Get list of coupled scalars."""

        lst = []
        for n in self.node_scals.xmlGetChildNodeList("scalar"):
            lst.append(n['name'])

        return lst


    def isScalarCoupled(self, scalar_name):
        """Check if a scalar is coupled."""

        retval = False
        if scalar_name in self.getScalarsList():
            retval = True

        return retval

    def getListOfAvailableScalars(self):
        """Get list of scalars which may be coupled"""

        lst = []

        thm = ThermalScalarModel(self.case)
        thermal = thm.getThermalScalarName()
        if thermal and thermal != "":
            lst.append(thermal)

        scm = DefineUserScalarsModel(self.case)
        for scalar in scm.getUserScalarNameList():
            lst.append(scalar)

        return lst
