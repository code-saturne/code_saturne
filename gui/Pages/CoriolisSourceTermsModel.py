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
This module contains the following classes and function:
- CoriolisSourceTermsModel
- CoriolisSourceTermsView
- PageText
- CoriolisSourceTermsTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from math import sqrt
import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Coriolis source terms model class
#-------------------------------------------------------------------------------

class CoriolisSourceTermsModel(Model):

    def __init__(self, case):

        """
        Constuctor.
        """
        self.case = case

        self.node_pp      = self.case.xmlGetNode('physical_properties')
        self.node_omega = self.node_pp.xmlInitNode('omega')

        self.nodes = ['omega_x', 'omega_y', 'omega_z']


    def defaultCoriolisSourceTermsValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['omega'] = 0.0

        return default


    @Variables.noUndo
    def getOmega(self, var):
        """
        Return value of omega for var
        """
        self.isInList(var, self.nodes)
        omega = self.node_omega.xmlGetDouble(var)
        if omega == None:
            omega = self.defaultCoriolisSourceTermsValues()['omega']
            self.setOmega(var, omega)

        return omega


    @Variables.undoLocal
    def setOmega(self, txml, value):
        """
        Put value of omega for txml balise
        """
        self.isInList(txml, self.nodes)
        self.isFloat(value)
        self.node_omega.xmlSetData(txml, value)


#-------------------------------------------------------------------------------
# CoriolisSourceTerms test case
#-------------------------------------------------------------------------------


class CoriolisSourceTermsTestCase(ModelTest):
    """
    """
    def checkCoriolisSourceTermsinstantiation(self):
        """
        Check whether the CoriolisSourceTermsModel class could be instantiated
        """
        model = None
        model = CoriolisSourceTermsModel(self.case)
        assert model != None, 'Could not instantiate CoriolisSourceTerms'

    def checkGetandSetOmega(self):
        """
        Check whether the omega terms could be set and get
        """
        mdl = CoriolisSourceTermsModel(self.case)
        mdl.setOmega('omega_x', 1.1)
        mdl.setOmega('omega_y', 2.2)
        mdl.setOmega('omega_z', 3.3)
        doc = '''<omega>
                         <omega_x>
                                 1.1
                         </omega_x>
                         <omega_y>
                                 2.2
                         </omega_y>
                         <omega_z>
                                 3.3
                         </omega_z>
                 </omega>'''
        assert mdl.node_omega == self.xmlNodeFromString(doc), \
              'Could not set omega'
        assert mdl.getOmega('omega_y') == 2.2, \
              'Could not get omega'


def suite():
    testSuite = unittest.makeSuite(CoriolisSourceTermsTestCase, "check")
    return testSuite


def runTest():
    print("CoriolisSourceTermsTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
