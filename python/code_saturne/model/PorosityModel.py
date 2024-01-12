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
This module defines the Porosity model data management.

This module contains the following classes and function:
- PorosityModel
- PorosityTestCase
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.LocalizationModel import LocalizationModel, VolumicLocalizationModel, Zone
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Porosity model class
#-------------------------------------------------------------------------------

class PorosityModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about Porosity
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_domain  = self.case.xmlGetNode('solution_domain')
        self.node_volzone = self.node_domain.xmlGetNode('volumic_conditions')
        self.node_porosit = self.node_models.xmlInitNode('porosities')


        self.notebook = NotebookModel(self.case)
        self.choicevalue = ('choice')

        self.getNameAndLocalizationZone()


    def __defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['kxx']     = 0.0
        default['choice'] = 'isotropic'
        return default


    @Variables.noUndo
    def getNameAndLocalizationZone(self):
        """
        Return name and localization zone from volume regions definitions.
        """
        zoneDico = {}
        zonesList = LocalizationModel('VolumicZone', self.case).getZones()
        for zone in zonesList:
            if zone.isNatureActivated("porosity"):
                label = zone.getLabel()
                zoneid = zone.getCodeNumber()
                localization = zone.getLocalization()
                zoneDico[label] = (zoneid, localization)
                self.setNameAndLabelZone(zoneid)

        return zoneDico


    @Variables.undoGlobal
    def setNameAndLabelZone(self, zoneid):
        """
        Set name and label zone for porosity markups.
        """
        self.node_porosit.xmlInitChildNode('porosity', zone_id=zoneid)

        # Initialize mdl at None for sanity sake
        mdl = None

        # Use existing model or default if none is defined
        try:
            n = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)
            mdl = n['model']
        except:
            pass

        # If model is None set to default.
        # Could be the result of 'try' failure or un-initialized zone
        if not mdl:
            mdl = self.__defaultValues()['choice']

        self.setPorosityModel(zoneid, mdl)


    @Variables.noUndo
    def getPorosityModel(self, zoneid):
        """
        Get the Transfo Matrix choice
        """
        self.isInt(int(zoneid))
        node = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)

        try:
            return node["model"]
        except TypeError as e:
            message = "Porosity not activated for zone : {0}".format(zoneid)
            raise Exception(message).with_traceback(e.__traceback__)

    @Variables.undoLocal
    def setPorosityModel(self, zoneid, choice):
        """
        Get the Transfo Matrix choice
        """
        self.isInt(int(zoneid))
        self.isInList(choice, ['isotropic', 'anisotropic', 'integral'])
        node = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)

        oldchoice = node['model']

        node['model'] = choice

        if oldchoice != None and oldchoice != choice:
            node.xmlRemoveChild('formula')


    @Variables.undoLocal
    def setPorosityFormula(self, zoneid, formula):
        """
        Public method.
        Set the formula for the porosity.
        """
        self.isInt(int(zoneid))
        node = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)
        n = node.xmlInitChildNode('formula')
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getPorosityFormula(self, zoneid):
        """
        Public method.
        Return the formula for the porosity.
        """
        self.isInt(int(zoneid))
        node = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)

        formula = "porosity = 1.0;"
        if node:
            formula = node.xmlGetString('formula')

        return formula


    @Variables.noUndo
    def getDefaultPorosityFormula(self, choice):
        """
        Public method.
        Return the default formula for the porosity.
        """
        self.isInList(choice, ['isotropic', 'anisotropic', 'integral'])
        if choice == 'anisotropic':
            formula = """porosity = 1.;
tensorial_porosity[XX] = 1.;
tensorial_porosity[YY] = 1.;
tensorial_porosity[ZZ] = 1.;
tensorial_porosity[XY] = 0.;
tensorial_porosity[XZ] = 0.;
tensorial_porosity[YZ] = 0.;"""
        else:
            formula = """porosity = 1.;"""

        return formula


    @Variables.noUndo
    def getPorosityFormulaComponents(self, zoneid):

        exp = self.getPorosityFormula(zoneid)

        if self.getPorosityModel(zoneid) in ('isotropic', 'integral'):
            req = [('porosity', 'Porosity')]
        else:
            req = [('porosity', 'Porosity'),
                   ('tensorial_porosity[XX]', 'Tensor Porosity'),
                   ('tensorial_porosity[YY]', 'Tensor Porosity'),
                   ('tensorial_porosity[ZZ]', 'Tensor Porosity'),
                   ('tensorial_porosity[XY]', 'Tensor Porosity'),
                   ('tensorial_porosity[XZ]', 'Tensor Porosity'),
                   ('tensorial_porosity[YZ]', 'Tensor Porosity')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, [], sym

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
