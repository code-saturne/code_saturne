# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
This module defines the XML calls for preprocessor execution
This module contains the following classes and function:
- BalanceModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables
from code_saturne.model.XMLvariables import Model
from code_saturne.model.BalanceModel import *

#-------------------------------------------------------------------------------
# Class Balance Model
#-------------------------------------------------------------------------------

class BalanceModelNeptune(BalanceModel, Variables, Model):

    """
    This class manages the scalar balance objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        BalanceModel.__init__(self, case)


    @Variables.noUndo
    def __getVariables__(self):

        # Get fields:
        fd = []
        fd.append('none')
        thermo = self.case.xmlGetNode('thermophysical_models')
        fields = thermo.xmlGetNode('fields')
        for node in fields.xmlInitChildNodeList('field'):
            field_id = node.xmlGetAttribute('field_id')
            fd.append(field_id)

        l = []
        for variableType in ('variable', 'property', 'scalar'):
            for field in fd:
                for node in self.case.xmlGetNodeList(variableType, field_id = field):
                    l.append(node)

        return l


    @Variables.noUndo
    def getScalarVariables(self):
        """
        Creates a dictionnary to connect name and label of
        scalar variables (different pressure).
        """
        self.dicoLabel2Name = {}

        for nodeList in [self.getVariables()]:

            for node in nodeList:

                name = node['name']
                label = node['label']
                if not label:
                    raise ValueError("Node has no label")

                dim = node['dimension']
                if not dim or int(dim) == 1:
                    if name != "pressure":
                        self.dicoLabel2Name[label] = (name, str(0))

        return list(self.dicoLabel2Name.keys())
