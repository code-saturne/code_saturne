# -*- coding: utf-8 -*-
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
This module defines the values of reference.

This module contains the following classes and function:
- MemoryAllocationModel
- MemoryAllocationTestCase
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
from Base.XMLvariables import Variables, Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Memory allocation model class
#-------------------------------------------------------------------------------

class MemoryAllocationModel(Model):
    """
    Manage the input/output markups in the xml doc about Integer and real memory allocation
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_calcul = self.case.xmlGetNode('calculation_management')
        self.node_int    = self.node_calcul.xmlInitNode('integer_work_array')
        self.node_real   = self.node_calcul.xmlInitNode('real_work_array')


    def defaultValues(self):
        """
        Return reference values by default
        """
        default = {}
        default['icel'] = 0
        default['ifac'] = 0
        default['ifab'] = 0
        default['idls'] = 0
        default['rcel'] = 0
        default['rfac'] = 0
        default['rfab'] = 0
        default['rdls'] = 0
        return default


    def setIntegerNcelet(self, value):
        """
        Set number of cells with halo into integer work's array.
        """
        self.isPositiveInt(value)
        self.node_int.xmlSetData('ncelet', value)


    def setIntegerNfac(self, value):
        """
        Set number of internal faces into integer work's array.
        """
        self.isPositiveInt(value)
        self.node_int.xmlSetData('nfac', value)


    def setIntegerNfabor(self, value):
        """
        Set number of boundary faces into integer work's array.
        """
        self.isPositiveInt(value)
        self.node_int.xmlSetData('nfabor', value)


    def setIntegerDimless(self, value):
        """
        Set integer value for integer work's array.
        """
        self.isPositiveInt(value)
        self.node_int.xmlSetData('dimless', value)


    def setRealNcelet(self, value):
        """
        Set number of cells with halo into real work's array.
        """
        self.isPositiveInt(value)
        self.node_real.xmlSetData('ncelet', value)


    def setRealNfac(self, value):
        """
        Set number of internal faces into real work's array.
        """
        self.isPositiveInt(value)
        self.node_real.xmlSetData('nfac', value)


    def setRealNfabor(self, value):
        """
        Set number of boundary faces into real work's array.
        """
        self.isPositiveInt(value)
        self.node_real.xmlSetData('nfabor', value)


    def setRealDimless(self, value):
        """
        Set integer value into real work's array.
        """
        self.isPositiveInt(value)
        self.node_real.xmlSetData('dimless', value)


    def getIntegerNcelet(self):
        """
        Return number of cells with halo from integer work's array.
        """
        value = self.node_int.xmlGetInt('ncelet')
        if value == None:
            value = self.defaultValues()['icel']
            self.setIntegerNcelet(value)
        return value


    def getIntegerNfac(self):
        """
        Return number of cells with halo from integer work's array.
        """
        value = self.node_int.xmlGetInt('nfac')
        if value == None:
            value = self.defaultValues()['ifac']
        return value


    def getIntegerNfabor(self):
        """
        Return number of cells with halo from integer work's array.
        """
        value = self.node_int.xmlGetInt('nfabor')
        if value == None:
            value = self.defaultValues()['ifab']
        return value


    def getIntegerDimless(self):
        """
        Return number of cells with halo from integer work's array.
        """
        value = self.node_int.xmlGetInt('dimless')
        if value == None:
            value = self.defaultValues()['idls']
        return value


    def getRealNcelet(self):
        """
        Return number of cells with halo from real work's array.
        """
        value = self.node_real.xmlGetInt('ncelet')
        if value == None:
            value = self.defaultValues()['rcel']
        return value


    def getRealNfac(self):
        """
        Return number of cells with halo from real work's array.
        """
        value = self.node_real.xmlGetInt('nfac')
        if value == None:
            value = self.defaultValues()['rfac']
        return value


    def getRealNfabor(self):
        """
        Return number of cells with halo from real work's array.
        """
        value = self.node_real.xmlGetInt('nfabor')
        if value == None:
            value = self.defaultValues()['rfab']
        return value


    def getRealDimless(self):
        """
        Return number of cells with halo from real work's array.
        """
        value = self.node_real.xmlGetInt('dimless')
        if value == None:
            value = self.defaultValues()['rdls']
        return value


#-------------------------------------------------------------------------------
# MemoryAllocationModel test case
#-------------------------------------------------------------------------------

class MemoryAllocationTestCase(ModelTest):
    """
    """
    def checkMemoryAllocationInstantiation(self):
        """Check whether the MemoryAllocationModel class could be instantiated"""
        model = None
        model = MemoryAllocationModel(self.case)
        assert model != None, 'Could not instantiate MemoryAllocationModel'

    def checkGetandSetIntegerNcelet(self):
        """Check whether the integer ncelet memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setIntegerNcelet(111)
        doc = '''<integer_work_array>
                    <ncelet>111</ncelet>
                 </integer_work_array>'''
        assert mdl.node_int == self.xmlNodeFromString(doc),\
                'Could not set integer ncelet memory allocation model'
        assert mdl.getIntegerNcelet() == 111,\
                'Could not get integer ncelet memory allocation model'

    def checkGetandSetIntegerNfac(self):
        """Check whether the integer nfac memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setIntegerNfac(211)
        doc = '''<integer_work_array>
                    <nfac>211</nfac>
                 </integer_work_array>'''
        assert mdl.node_int == self.xmlNodeFromString(doc),\
                'Could not set integer nfac memory allocation model'
        assert mdl.getIntegerNfac() == 211,\
                'Could not get integer nfac memory allocation model'

    def checkGetandSetIntegerNfabor(self):
        """Check whether the integer nfabor memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setIntegerNfabor(33)
        doc = '''<integer_work_array>
                    <nfabor>33</nfabor>
                 </integer_work_array>'''
        assert mdl.node_int == self.xmlNodeFromString(doc),\
                'Could not set integer nfabor memory allocation model'
        assert mdl.getIntegerNfabor() == 33,\
                'Could not get integer nfabor memory allocation model'

    def checkGetandSetIntegerDimless(self):
        """Check whether the integer dimless memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setIntegerDimless(10)
        doc = '''<integer_work_array>
                    <dimless>10</dimless>
                 </integer_work_array>'''
        assert mdl.node_int == self.xmlNodeFromString(doc),\
                'Could not set integer dimless memory allocation model'
        assert mdl.getIntegerDimless() == 10,\
                'Could not get integer dimless memory allocation model'

    def checkGetandSetRealNcelet(self):
        """Check whether the real ncelet memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setRealNcelet(999)
        doc = '''<real_work_array>
                    <ncelet>999</ncelet>
                 </real_work_array>'''
        assert mdl.node_real == self.xmlNodeFromString(doc),\
                'Could not set real ncelet memory allocation model'
        assert mdl.getRealNcelet() == 999,\
                'Could not get real ncelet memory allocation model'

    def checkGetandSetRealNfac(self):
        """Check whether the real nfac memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setRealNfac(988)
        doc = '''<real_work_array>
                    <nfac>988</nfac>
                 </real_work_array>'''
        assert mdl.node_real == self.xmlNodeFromString(doc),\
                'Could not set real nfac memory allocation model'
        assert mdl.getRealNfac() == 988,\
                'Could not get real nfac memory allocation model'

    def checkGetandSetRealNfabor(self):
        """Check whether the real nfabor memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setRealNfabor(977)
        doc = '''<real_work_array>
                    <nfabor>977</nfabor>
                 </real_work_array>'''
        assert mdl.node_real == self.xmlNodeFromString(doc),\
                'Could not set real nfabor memory allocation model'
        assert mdl.getRealNfabor() == 977,\
                'Could not get real nfabor memory allocation model'

    def checkGetandSetRealDimless(self):
        """Check whether the real dimless memory allocation model could be set"""
        mdl = MemoryAllocationModel(self.case)
        mdl.setRealDimless(966)
        doc = '''<real_work_array>
                    <dimless>966</dimless>
                 </real_work_array>'''
        assert mdl.node_real == self.xmlNodeFromString(doc),\
                'Could not set real dimless memory allocation model'
        assert mdl.getRealDimless() == 966,\
                'Could not get real dimless memory allocation model'





def suite():
    testSuite = unittest.makeSuite(MemoryAllocationTestCase, "check")
    return testSuite

def runTest():
    print("MemoryAllocationTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
