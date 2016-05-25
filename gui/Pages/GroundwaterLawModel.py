# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2015 EDF S.A.
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
This module defines the Groundwater Law model data management.

This module contains the following classes and function:
- GroundwaterLawModel
- GroundwaterLawTestCase
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Pages.LocalizationModel import LocalizationModel, VolumicLocalizationModel, Zone
from code_saturne.Pages.GroundwaterModel import GroundwaterModel
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# GroundwaterLaw model class
#-------------------------------------------------------------------------------

class GroundwaterLawModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about Groundwater Law
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_domain  = self.case.xmlGetNode('solution_domain')
        self.node_volzone = self.node_domain.xmlGetNode('volumic_conditions')
        self.node_darcy   = self.node_models.xmlInitNode('groundwater')

        self.sca_mo       = DefineUserScalarsModel(self.case)

        self.choicevalue = ('choice')

        self.getNameAndLocalizationZone()


    def __defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['choice']                = 'VanGenuchten'
        default['ks']                    = 0.3
        default['ks_xx']                 = 0.3
        default['ks_yy']                 = 0.3
        default['ks_zz']                 = 0.3
        default['ks_xy']                 = 0.
        default['ks_xz']                 = 0.
        default['ks_yz']                 = 0.
        default['thetas']                = 0.3
        default['thetar']                = 0.078
        default['n']                     = 1.56
        default['l']                     = 0.5
        default['alpha']                 = 0.036
        default['longitudinal']          = 1.0
        default['transverse']            = 0.0
        default['isotropic']             = 0.0
        default['diffusion_choice']      = 'variable'
        return default


    @Variables.noUndo
    def getNameAndLocalizationZone(self):
        """
        Return name and localization zone from volume regions definitions.
        """
        zoneDico = {}
        zonesList = LocalizationModel('VolumicZone', self.case).getZones()
        for zone in zonesList:
            if zone.getNature()['groundwater_law'] == 'on':
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
        self.node_darcy.xmlInitChildNode('groundwater_law', zone_id=zoneid)
        self.getGroundwaterLawModel(zoneid)


    @Variables.noUndo
    def getGroundwaterLawModel(self, zoneid):
        """
        Get the darcy law choice
        """
        self.isInt(int(zoneid))
        node = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)

        mdl = node['model']
        if mdl == None:
            mdl = self.__defaultValues()['choice']
            self.setGroundwaterLawModel(zoneid, mdl)
        return mdl


    @Variables.undoLocal
    def setGroundwaterLawModel(self, zoneid, choice):
        """
        Get the darcy law choice
        """
        self.isInt(int(zoneid))
        self.isInList(choice, ['user', 'VanGenuchten'])
        node = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)

        oldchoice = node['model']

        node['model'] = choice

        # TODO
        if oldchoice != None and oldchoice != choice:
            if choice == "user":
                node.xmlRemoveChild('VanGenuchten_parameters')
            else:
                node.xmlRemoveChild('formula')


    @Variables.undoLocal
    def setValue(self, zoneid, variable, value):
        """
        Input value for variable
        """
        self.isInt(int(zoneid))
        self.isFloat(value)
        self.isInList(variable, ['ks', 'ks_xx',  'ks_yy',  'ks_zz',
                                 'ks_xy',  'ks_xz',  'ks_yz',
                                 'thetas', 'thetar','n','l','alpha'])

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        node = nodeZone.xmlInitChildNode('VanGenuchten_parameters')

        node.xmlSetData(variable, value)


    @Variables.noUndo
    def getValue(self, zoneid, variable):
        """
        Return value for variable
        """
        self.isInt(int(zoneid))
        self.isInList(variable, ['ks', 'ks_xx',  'ks_yy',  'ks_zz',
                                 'ks_xy',  'ks_xz',  'ks_yz',
                                 'thetas', 'thetar','n','l','alpha'])

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        node = nodeZone.xmlInitChildNode('VanGenuchten_parameters')

        value = node.xmlGetDouble(variable)

        if value == None:
            value = self.__defaultValues()[variable]
            self.setValue(zoneid, variable, value)
        return value


    @Variables.undoLocal
    def setDispersionCoefficient(self, zoneid, variable, value):
        """
        Input value for variable
        """
        self.isInt(int(zoneid))
        self.isFloat(value)
        self.isInList(variable, ['longitudinal', 'transverse', 'isotropic'])

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        node = nodeZone.xmlInitChildNode('diffusion_coefficient')

        node.xmlSetData(variable, value)


    @Variables.noUndo
    def getDispersionCoefficient(self, zoneid, variable):
        """
        Return value for variable
        """
        self.isInt(int(zoneid))
        self.isInList(variable, ['longitudinal', 'transverse', 'isotropic'])

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        node = nodeZone.xmlInitChildNode('diffusion_coefficient')

        value = node.xmlGetDouble(variable)

        if value == None:
            value = self.__defaultValues()[variable]
            self.setDispersionCoefficient(zoneid, variable, value)
        return value


    @Variables.undoLocal
    def setGroundwaterLawFormula(self, zoneid, formula):
        """
        Public method.
        Set the formula for Groundwater
        """
        self.isInt(int(zoneid))
        node = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        n = node.xmlInitChildNode('formula')
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getGroundwaterLawFormula(self, zoneid):
        """
        Public method.
        Return the formula for Groundwater
        """
        self.isInt(int(zoneid))
        node = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)

        formula = node.xmlGetString('formula')
        return formula


    @Variables.noUndo
    def getDefaultGroundwaterLawFormula(self):
        """
        Public method.
        Return the default formula for Groundwater.
        """

        if GroundwaterModel(self.case).getPermeabilityType() == 'anisotropic':
            formula = """capacity = 0.;
saturation = 1.;
permeability[XX]=1.;
permeability[YY]=1.;
permeability[ZZ]=1.;
permeability[XY]=0.;
permeability[XZ]=0.;
permeability[YZ]=0.;"""
        else:
            formula = """capacity = 0.;
saturation = 1.;
permeability=1.;"""

        return formula


    @Variables.undoLocal
    def setScalarDiffusivityChoice(self, scalar_label, zoneid, choice):
        """
        Set choice of diffusivity's property for an additional_scalar
        with label scalar_label
        """
        self.isInt(int(zoneid))
        self.isNotInList(scalar_label, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_label, self.sca_mo.getUserScalarNameList())
        self.isInList(choice, ('constant', 'variable'))

        name = self.sca_mo.getScalarDiffusivityName(scalar_label)

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        n = nodeZone.xmlInitChildNode('variable', name=name)
        n_diff = n.xmlInitChildNode('property')
        n_diff['choice'] = choice


    @Variables.noUndo
    def getScalarDiffusivityChoice(self, scalar_label, zoneid):
        """
        Get choice of diffusivity's property for an additional_scalar
        with label scalar_label
        """
        self.isInt(int(zoneid))
        self.isNotInList(scalar_label, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_label, self.sca_mo.getUserScalarNameList())

        name = self.sca_mo.getScalarDiffusivityName(scalar_label)

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        n = nodeZone.xmlInitChildNode('variable', name=name)
        choice = n.xmlInitChildNode('property')['choice']
        if not choice:
            choice = self.__defaultValues()['diffusion_choice']
            self.setScalarDiffusivityChoice(scalar_label, zoneid, choice)

        return choice


    @Variables.noUndo
    def getDiffFormula(self, scalar_label, zoneid):
        """
        Return a formula for I{tag} 'density', 'molecular_viscosity',
        'specific_heat' or 'thermal_conductivity'
        """
        self.isInt(int(zoneid))
        self.isNotInList(scalar_label, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_label, self.sca_mo.getUserScalarNameList())

        name = self.sca_mo.getScalarDiffusivityName(scalar_label)

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        n = nodeZone.xmlInitChildNode('variable', name=name)
        node = n.xmlGetNode('property')
        formula = node.xmlGetString('formula')
        if not formula:
            formula = self.getDefaultFormula(scalar_label)
            self.setDiffFormula(scalar_label, zoneid, formula)
        return formula


    @Variables.noUndo
    def getDefaultFormula(self, scalar_label):
        """
        Return default formula
        """
        self.isNotInList(scalar_label, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_label, self.sca_mo.getUserScalarNameList())

        name = self.sca_mo.getScalarDiffusivityName(scalar_label)

        formula = str(name) + " = 0.;\n" + str(scalar_label) + "_delay = 1.;"

        return formula


    @Variables.undoLocal
    def setDiffFormula(self, scalar_label, zoneid, str):
        """
        Gives a formula for 'density', 'molecular_viscosity',
        'specific_heat'or 'thermal_conductivity'
        """
        self.isInt(int(zoneid))
        self.isNotInList(scalar_label, self.sca_mo.getScalarsVarianceList())
        self.isInList(scalar_label, self.sca_mo.getUserScalarNameList())

        name = self.sca_mo.getScalarDiffusivityName(scalar_label)

        nodeZone = self.node_darcy.xmlGetNode('groundwater_law', zone_id=zoneid)
        n = nodeZone.xmlInitChildNode('variable', name=name)
        node = n.xmlGetNode('property')
        node.xmlSetData('formula', str)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
