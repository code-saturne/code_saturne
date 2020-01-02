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
This module initialize model dynamics variables and model scalars

This module contents the following classes:
- MainFieldsSourceTermsModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Variables and Scalar model initialization modelling class
#-------------------------------------------------------------------------------

class MainFieldsSourceTermsModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.models     = self.case.xmlGetNode('thermophysical_models')
        self.node_sterm = self.models.xmlInitNode('source_terms')
        self.mfm        = MainFieldsModel(self.case)

        self.notebook = NotebookModel(self.case)


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())

    def getKnownFields(self, fieldId):
        field_name = self.mfm.getFieldLabelsList()[int(fieldId)-1]

        known_fields = [('enthalpy_'+field_name, 'enthalpy_'+str(fieldId))]

        return known_fields

    def getThermalFormulaComponents(self, zone, fieldId, scalar):

        exp = self.getThermalFormula(zone, fieldId, scalar)
        if not exp:
            exp = self.getDefaultThermalFormula(scalar)
        req = [('S', 'thermal source term'),
               ('dS', 'thermal source term derivative')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('t', 'current time'),
               ('volume', 'Source terms zone volume')]

        for knf in self.getKnownFields(fieldId):
            sym.append(knf)

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoGlobal
    def setThermalFormula(self, zone, fieldId, scalar, formula):
        """
        Public method.
        Set the formula for tharmal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, ['enthalpy'])
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)
        n = node.xmlInitChildNode('thermal_formula',
                                  name = scalar,
                                  zone_id=zone,
                                  field_id=fieldId)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getThermalFormula(self, zone, fieldId, scalar):
        """
        Public method.
        Return the formula for thermal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, ['enthalpy'])
        node = self.node_sterm
        if not node:
            msg = "There is an error: this node " + str(node) + "should exist"
            raise ValueError(msg)

        formula = node.xmlGetString('thermal_formula',
                                    name = scalar,
                                    zone_id=zone,
                                    field_id=fieldId)

        return formula


    def getDefaultThermalFormula(self, scalar):
        """
        Public method.
        Return the default formula for thermal scalars.
        """
        formula = """S = 0;\ndS = 0;"""

        return formula


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
