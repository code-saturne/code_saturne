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
            if zone.getNature()['porosity'] == 'on':
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
        self.getPorosityModel(zoneid)


    @Variables.noUndo
    def getPorosityModel(self, zoneid):
        """
        Get the Transfo Matrix choice
        """
        self.isInt(int(zoneid))
        node = self.node_porosit.xmlGetNode('porosity', zone_id=zoneid)

        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['choice']
            self.setPorosityModel(zoneid, mdl)
        return mdl


    @Variables.undoLocal
    def setPorosityModel(self, zoneid, choice):
        """
        Get the Transfo Matrix choice
        """
        self.isInt(int(zoneid))
        self.isInList(choice, ['isotropic', 'anisotropic'])
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

        formula = node.xmlGetString('formula')
        return formula


    @Variables.noUndo
    def getDefaultPorosityFormula(self, choice):
        """
        Public method.
        Return the default formula for the porosity.
        """
        self.isInList(choice, ['isotropic', 'anisotropic'])
        if choice == 'isotropic':
            formula = """porosity=1.;"""
        else:
            formula = """porosity=1.;
porosity[XX]=1.;
porosity[YY]=1.;
porosity[ZZ]=1.;
porosity[XY]=0.;
porosity[XZ]=0.;
porosity[YZ]=0.;"""

        return formula


    @Variables.noUndo
    def getPorosityFormulaComponents(self, zoneid):

        exp = self.getPorosityFormula(zoneid)

        if self.getPorosityModel(zoneid) == 'isotropic':
            req = [('porosity', 'Porosity')]
        else:
            req = [('porosity', 'Porosity'),
                   ('porosity[XX]', 'Porosity'),
                   ('porosity[YY]', 'Porosity'),
                   ('porosity[ZZ]', 'Porosity'),
                   ('porosity[XY]', 'Porosity'),
                   ('porosity[XZ]', 'Porosity'),
                   ('porosity[YZ]', 'Porosity')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, [], sym

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
