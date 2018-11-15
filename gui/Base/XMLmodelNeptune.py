# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
This module defines the XML data model in which the user defines the physical
options of the treated case. This module defines also a very usefull
function for the NavigatorTree display updating, taking into account
the current options selected by the user.

This module contains the following classes and function:
- XMLmodel
- XMLmodelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.XMLvariablesNeptune import Variables
from code_saturne.Base import Toolbox
from code_saturne.Base.XMLmodel import XMLmodel

#-------------------------------------------------------------------------------
# class XMLmodel
#-------------------------------------------------------------------------------

class TurbulenceModelsDescribing:
    """
    """
    continuousTurbulenceModels = ['none', 'mixing_length',
                                  'k-epsilon', 'k-epsilon_linear_production',
                                  'rij-epsilon_ssg', 'rij-epsilon_ebrsm',
                                  'les_smagorinsky', 'les_wale']

    dispersedTurbulenceModels  = ['none','tchen','q2-q12', 'r2-q12', 'r2-r12-tchen']

    continuousCouplingModels = ['none','separate_phase','separate_phase_cond']
    dispersedCouplingModels  = ['none','small_inclusions','large_inclusions']

    ThermalTurbFluxModels = ['sgdh', 'ggdh']

    turbulenceVariables = {}
    turbulenceVariables['none'] = []
    turbulenceVariables['mixing_length'] = []
    turbulenceVariables['k-epsilon'] = ['TurbKineEner_k', 'TurbDissip']
    turbulenceVariables['k-epsilon_linear_production'] = ['TurbKineEner_k', 'TurbDissip']
    turbulenceVariables['rij-epsilon_ssg'] = ['ReynoldsStressXX', 'ReynoldsStressXY', 'ReynoldsStressXZ',
                        'ReynoldsStressYY', 'ReynoldsStressYZ', 'ReynoldsStressZZ', 'TurbDissip']
    turbulenceVariables['rij-epsilon_ebrsm'] = ['ReynoldsStressXX', 'ReynoldsStressXY', 'ReynoldsStressXZ',
                        'ReynoldsStressYY', 'ReynoldsStressYZ', 'ReynoldsStressZZ', 'TurbDissip']
    turbulenceVariables['les_smagorinsky'] = []
    turbulenceVariables['les_wale'] = []
    turbulenceVariables['tchen'] = []
    turbulenceVariables['q2-q12'] = ['TurbKineEner_q2', 'Covariance_q12']
    turbulenceVariables['r2-q12'] = ['ReynoldsStressXX', 'ReynoldsStressXY', 'ReynoldsStressXZ',
                                     'ReynoldsStressYY', 'ReynoldsStressYZ', 'ReynoldsStressZZ','Covariance_q12']
    turbulenceVariables['r2-r12-tchen'] = ['ReynoldsStressXX', 'ReynoldsStressXY', 'ReynoldsStressXZ',
                                           'ReynoldsStressYY', 'ReynoldsStressYZ', 'ReynoldsStressZZ',
                                           'R12XX','R12XY','R12XZ','R12YY','R12YZ','R12ZZ']

    turbulenceVariables['all'] = turbulenceVariables['k-epsilon'] \
                               + turbulenceVariables['k-epsilon_linear_production'] \
                               + turbulenceVariables['rij-epsilon_ssg'] \
                               + turbulenceVariables['rij-epsilon_ebrsm'] \
                               + turbulenceVariables['les_smagorinsky'] \
                               + turbulenceVariables['les_wale'] \
                               + turbulenceVariables['q2-q12'] \
                               + turbulenceVariables['r2-q12'] \
                               + turbulenceVariables['r2-r12-tchen']

    turbulenceProperties = {}
    turbulenceProperties['none'] = []
    turbulenceProperties['mixing_length'] = ["turb_viscosity"]
    turbulenceProperties['k-epsilon'] = ["turb_viscosity"]
    turbulenceProperties['k-epsilon_linear_production'] = ["turb_viscosity"]
    turbulenceProperties['rij-epsilon_ssg'] = ["turb_viscosity"]
    turbulenceProperties['rij-epsilon_ebrsm'] = ["turb_viscosity"]
    turbulenceProperties['les_smagorinsky'] = ["turb_viscosity"]
    turbulenceProperties['les_wale'] = ["turb_viscosity"]
    turbulenceProperties['tchen'] = ["TurbKineEner_q2", "Covariance_q12", "turb_viscosity"]
    turbulenceProperties['q2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-r12-tchen'] = ["turb_viscosity"]



class FieldAttributesDescribing:
    """
    """
    typeChoiceValues = ['continuous', 'dispersed', 'auto']
    phaseValues = ['liquid', 'gas', 'particle']


class XMLmodel(XMLmodel, Variables):
    #TODO revoir toute la classe  : a supprimer ????? (cf output)
    """
    This class initialize the XML contents of the case.
    """
    def __init__(self, case):
        """
        """
        self.case = case
        self.root = self.case.root()
        self.node_models = self.case.xmlGetNode('closure_modeling')


#-------------------------------------------------------------------------------
# XMLmodel test case
#-------------------------------------------------------------------------------

class ModelTest(unittest.TestCase):
    """
    Class beginning class test case of all pages
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from code_saturne.Base.XMLengine import Case, XMLDocument
        from code_saturne.Base.XMLinitializeNeptune import XMLinit
        Toolbox.GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case)
        self.doc = XMLDocument()

    def tearDown(self):
        """This method is executed after all "check" methods."""
        del self.case
        del self.doc

    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        n = self.doc.parseString(string).root()
        self.doc.xmlCleanAllBlank(n)
        return self.doc.root()
