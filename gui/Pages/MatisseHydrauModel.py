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
- MatisseHydrauModel
- MatisseHydrauModelTestCase
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
import Pages.MatisseTypeModel as MatisseType

#-------------------------------------------------------------------------------
# Matisse Hydraulique model class
#-------------------------------------------------------------------------------

class MatisseHydrauModel:
    """
    Manage the input/output markups in the xml doc about matisse hydraulic load
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case               = case
        self.node_matisse       = self.case.root().xmlInitChildNode('matisse')
        self.node_compute       = self.node_matisse.xmlInitChildNode('compute')
        self.node_phymodel      = self.node_compute.xmlInitChildNode('physical_model')

        self.node_cofor         = self.node_phymodel.xmlInitChildNode('icofor','status')
        self.node_conlg         = self.node_phymodel.xmlInitChildNode('iconlg','status')
        
        self.status = ('on', 'off')


    def defaultMatisseHydrauValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        #
        # bool
        default['icofor'] = 'off'
        default['iconlg'] = 'off'

        #
        # double
        default['debmas'] = 0.
        default['pdccha'] = 1.7
        default['pdcfch'] = 0.
        default['dhchea'] = 1.7
        default['sdchea'] = 3.2
        default['pdcche'] = 2.
        default['pdccch'] = 0.
        default['dhches'] = 1.2
        default['sdches'] = 3.82
        default['pdcalg'] = 0.25
        default['pdcatv'] = 0.25
        default['argamt'] = 0.
        default['pdcslg'] = 0.25
        default['pdcstv'] = 0.25
        default['argavl'] = 0.
        default['amppdc'] = 1.
        default['dhalve'] = 0.13
        default['dpvent'] = 0.
        
        return default


    def setMatisseHydrauVar(self, tag, val):
        """
        """
        self.node_phymodel.xmlSetData(tag, val)


    def getMatisseHydrauDoubleVar(self, tag):
        """
        """
        val = self.node_phymodel.xmlGetDouble(tag)

        if val == "" or val == None:
            self.node_phymodel.xmlInitChildNode(tag)
            val = self.defaultMatisseHydrauValues()[tag]
            self.setMatisseHydrauVar(tag, val)

        return val


    def setConstrainedConvStatus(self, stat):
        """
        Input constrained convection status
        """
        if stat not in self.status :
            sys.exit(1)

        self.node_cofor['status'] = stat


    def getConstrainedConvStatus(self):
        """
        Return constrained convection status
        """
        stat = self.node_cofor['status']
            
        if stat not in self.status :
            stat = self.defaultMatisseHydrauValues()['icofor']
            self.setConstrainedConvStatus(stat)
        return stat


    def setInlineContainerNetworkStatus(self, stat):
        """
        Input containers network in line status
        """
        if stat not in self.status :
            sys.exit(1)

        self.node_conlg['status'] = stat


    def getInlineContainerNetworkStatus(self):
        """
        Return containers network in line status
        """
        stat = self.node_conlg['status']
            
        if stat not in self.status :
            stat = self.defaultMatisseHydrauValues()['iconlg']
            self.setInlineContainerNetworkStatus(stat)
        return stat

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
class MatisseHydrauModelTestCase(unittest.TestCase):
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

    def checkMatisseHydrauModelInstantiation(self):
        """Check whether the NOMModel class could be instantiated"""
        model = None
        model = MatisseHydrauModel(self.case)
        assert model != None, 'Could not instantiate MatisseHydrauModel'
        
        
        
def suite():
    testSuite = unittest.makeSuite(MatisseHydrauModelTestCase, "check")
    return testSuite

def runTest():
    print "MatisseHydrauModelTestCase - A FAIRE**************"
    runner = unittest.TextTestRunner()
    runner.run(suite())



#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------