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

"""
This module defines the storage system type.

This module contains the following classes and function:
- MatisseTypeModel
- MatisseTypeModelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool

#-------------------------------------------------------------------------------
# Turbulence model class
#-------------------------------------------------------------------------------

class MatisseTypeModel:
    def __init__(self, case):
        self.case        = case
        self.node_matisse= self.case.root().xmlInitChildNode('matisse')
        self.node_compute= self.node_matisse.xmlInitChildNode('compute')
        self.node_geom   = self.node_compute.xmlInitChildNode('geometry')
        self.node_phys   = self.node_compute.xmlInitChildNode('physical_model')

        self.node_type   = self.node_geom.xmlInitChildNode('typent', 'label')
        self.node_alveo  = self.node_phys.xmlInitChildNode('ialveo', 'status')

        self.matisseTypes = ('vault', 'emm','djw')


    def defaultMatisseTypeValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['typent'] = 'vault'
        return default


    def setMatisseType(self, matisse_type):
        """
        Input storage matisse type
        """
        if matisse_type not in self.matisseTypes :
            print("MatisseTypeModel: error with the storage concept list")
            sys.exit(1)

        if matisse_type == 'djw' :
            self.node_type['label'] = 'vault'
            self.node_alveo['status'] = 'on'
        else :
            self.node_alveo['status'] = 'off'
            self.node_type['label'] = matisse_type


    def getMatisseType(self):
        """
        Return the storage matisse type
        """
        if self.node_alveo['status'] == 'on' :
            if self.node_type['label'] != 'vault':
                print("MatisseTypeModel: error with the storage concept input")
                sys.exit(1)
            matisse_type = 'djw'
        else :
            matisse_type = self.node_type['label']

        if matisse_type not in self.matisseTypes :
            matisse_type = self.defaultMatisseTypeValues()['typent']
            self.setMatisseType(matisse_type)
        return matisse_type



#-------------------------------------------------------------------------------
# MatisseTypeTestCase test case
#-------------------------------------------------------------------------------


class MatisseTypeTestCase(unittest.TestCase):
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



    def checkMatisseTypeInstantiation(self):
        """
        Check whether the TurbulenceModel class could be instantiated
        """
        model = None
        model = MatisseTypeModel(self.case)
        assert model != None, 'Could not instantiate MatisseTypeModel'


def suite():
    testSuite = unittest.makeSuite(MatisseTypeTestCase, "check")
    return testSuite


def runTest():
    print("MatisseTypeTestCase - TODO*****")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
