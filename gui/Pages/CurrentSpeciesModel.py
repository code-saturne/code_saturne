# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, re, sys
import string
import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
import Pages.CoalCombustionModel as CoalCombustion
import Pages.CommonCombustion as CommonCombustion

#-------------------------------------------------------------------------------
# Class Model Main
#-------------------------------------------------------------------------------

class CurrentSpeciesModel:

    def __init__(self, case):
        """ constructor """
        self.case = case
        model = CoalCombustion.CoalCombustionModel(self.case)
        if model.getCoalCombustionModel() != "off":
            import Pages.CoalThermoChemistry as CoalThermoChemistry
            model = CoalThermoChemistry.CoalThermoChemistryModel('dp_FCP', self.case)
            self.species = model.getSpecies()

    def getSpecies(self):
        return self.species

#-------------------------------------------------------------------------------
# Class CurrentSpeciesModelTestCase
#-------------------------------------------------------------------------------

class CurrentSpeciesModelTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from Base.XMLengine import Case, XMLDocument
        from Base.XMLinitialize import XMLinit
        Tool.GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case)
        self.doc = XMLDocument()


    def tearDown(self):
        """This method is executed after all "check" methods."""
        del self.case
        del self.doc


    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()


    def checkCurrentSpeciesModelInstantiation(self):
        """Check whether the NOMModel class could be instantiated"""
        model = None
        model = CurrentSpeciesModel(self.case)
        assert model != None, 'Could not instantiate CurrentSpeciesModel'


def suite():
    testSuite = unittest.makeSuite(CurrentSpeciesModelTestCase, "check")
    return testSuite

def runTest():
    print("CurrentSpeciesModelTestCase - A FAIRE**************")
    runner = unittest.TextTestRunner()
    runner.run(suite())
