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
This module defines the conjugate heat transfer model data management.

This module contains the following classes and function:
- ConjugateHeatTransferModel
- ConjugateHeatTransferModelTestCase
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

#-------------------------------------------------------------------------------
# Conjugate Heat Transfer model class
#-------------------------------------------------------------------------------

class ConjugateHeatTransferModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about HeadLosses
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.__node_models = self.case.xmlGetNode('thermophysical_models')
        self.__node_cht    = self.__node_models.xmlInitNode('conjugate_heat_transfer')
        self.__node_syr    = self.__node_cht.xmlInitNode('external_coupling')

    def defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['syrthes_name']       = "SYRTHES"
        default['verbosity']  = 0
        default['visualization']  = 1
        default['projection_axis']    = "off"
        default['selection_criteria'] = "all[]"
        return default


    def __getNumberOfSyrthesCoupling(self):
        return len(self.__node_syr.xmlGetNodeList('syrthes'))


    @Variables.undoLocal
    def deleteConjugateHeatTransfer(self):
        """
        Update the 'Conjugate heat transfer' status.
        @type status: C{String}
        @param status: set to "on" or "off" the conjugate heat transfer
        """
        self.__node_syr.xmlRemoveChild('syrthes')


    @Variables.noUndo
    def getSyrthesCouplingList(self):
        """
        @return: list of Syrthes coupling description.
        @rtype: C{List}
        """
        node_list = self.__node_syr.xmlGetNodeList('syrthes')
        lst = []
        for index in range(len(node_list)):
            num = index + 1
            syrthes_name  = self.getSyrthesInstanceName(num)
            verbosity     = self.getSyrthesVerbosity(num)
            visualization = self.getSyrthesVisualization(num)
            proj_axis     = self.getSyrthesProjectionAxis(num)
            location      = self.getSelectionCriteria(num)
            lst.append([syrthes_name, verbosity, visualization,
                         proj_axis, location])

        return lst


    @Variables.undoGlobal
    def addSyrthesCoupling(self, syrthes_name,
                           verbosity, visualization, proj_axis, location):
        """
        Add a new definition of a Syrthes coupling.

        @type syrthes_name: C{String}
        @param syrthes_name: Syrthes instance name
        @type verbosity: C{Int}
        @param verbosity: Syrthes verbosity
        @param visualization: Syrthes visualization output
        @type proj_axis: C{String}
        @param proj_axis: Syrthes projection axis
        @type location: C{String}
        @param location: selection criteria
        @return: new number of Syrthes coupling
        @rtype: C{Int}
        """
        num = len(self.__node_syr.xmlGetNodeList('syrthes'))
        node_new = self.__node_syr.xmlAddChild('syrthes')

        num = num + 1
        self.setSyrthesInstanceName(num, syrthes_name)
        self.setSyrthesVerbosity(num, verbosity)
        self.setSyrthesVisualization(num, visualization)
        self.setSyrthesProjectionAxis(num, proj_axis)
        self.setSelectionCriteria(num, location)

        return num


    @Variables.undoLocal
    def deleteSyrthesCoupling(self, num):
        """
        Delete a definition of a Syrthes coupling.

        @type num: C{Int}
        @param num: Syrthes coupling number
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        node_list = self.__node_syr.xmlGetNodeList('syrthes')
        node = node_list[num-1]
        node.xmlRemoveNode()

    #------------------------------------------------------------------
    # Helper function
    #------------------------------------------------------------------
    def __getStringData(self, index, name, setFunction):
        """
        Get string value from xml file.
        """
        self.isLowerOrEqual(index+1, self.__getNumberOfSyrthesCoupling())
        node = self.__node_syr.xmlGetNodeList('syrthes')[index]
        value = node.xmlGetString(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getIntData(self, index, name, setFunction):
        """
        Get int value from xml file.
        """
        self.isLowerOrEqual(index+1, self.__getNumberOfSyrthesCoupling())
        node = self.__node_syr.xmlGetNodeList('syrthes')[index]
        value = node.xmlGetInt(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getDefaultDataIfNone(self, index, value, name, setFunction):
        """
        Get default value if value is none.
        """
        if value == None or value == "":
            value = self.defaultValues()[name]
            setFunction(index+1, value)
        return value

    #------------------------------------------------------------------
    # Syrthes instance name
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setSyrthesInstanceName(self, num, value):
        """
        Set value of Syrthes instance name.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{String}
        @param value: Syrthes instance name
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isStr(value)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('syrthes_name', value)


    @Variables.noUndo
    def getSyrthesInstanceName(self, num):
        """
        Get value of Syrthes instance name.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes verbosity
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'syrthes_name',
                                    self.setSyrthesInstanceName)

    #------------------------------------------------------------------
    # Syrthes verbosity
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setSyrthesVerbosity(self, num, value):
        """
        Set value of Syrthes verbosity.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{Int}
        @param value: Syrthes verbosity
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('verbosity', value)


    @Variables.noUndo
    def getSyrthesVerbosity(self, num):
        """
        Get value of Syrthes verbosity.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes verbosity
        @rtype: C{Int}
        """
        return self.__getIntData(num-1,
                                 'verbosity',
                                 self.setSyrthesVerbosity)

    #------------------------------------------------------------------
    # Syrthes visualization output
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setSyrthesVisualization(self, num, value):
        """
        Set value of Syrthes visualization.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{Int}
        @param value: Syrthes visualization
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('visualization', value)


    @Variables.noUndo
    def getSyrthesVisualization(self, num):
        """
        Get value of Syrthes visualization.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes visualization
        @rtype: C{Int}
        """
        return self.__getIntData(num-1,
                                 'visualization',
                                 self.setSyrthesVisualization)

    #------------------------------------------------------------------
    # Projection axis
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setSyrthesProjectionAxis(self, num, value):
        """
        Set value of Syrthes projection axis.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{String}
        @param value: Syrthes projection axis
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isStr(value)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('projection_axis', value)


    @Variables.noUndo
    def getSyrthesProjectionAxis(self, num):
        """
        Get value of Syrthes projection axis.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes projection axis
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'projection_axis',
                                    self.setSyrthesVerbosity)

    #------------------------------------------------------------------
    # Selection criteria
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setSelectionCriteria(self, num, value):
        """
        Set value of selection criteria.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{String}
        @param value: selection criteria
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isStr(value)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('selection_criteria', value)


    @Variables.noUndo
    def getSelectionCriteria(self, num):
        """
        Get value of selection criteria.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: selection criteria
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'selection_criteria',
                                    self.setSelectionCriteria)

#-------------------------------------------------------------------------------
# ConjugateHeatTransferModel test case
#-------------------------------------------------------------------------------


class ConjugateHeatTransferModelTestCase(ModelTest):
    """
    """
    def checkConjugateHeatTransferInstantiation(self):
        """Check whether the ConjugateHeatTransfer class could be instantiated"""
        model = None
        model = ConjugateHeatTransferModel(self.case)
        assert model != None, 'Could not instantiate HeadLossesModel'


def suite():
    testSuite = unittest.makeSuite(ConjugateHeatTransferModelTestCase, "check")
    return testSuite

def runTest():
    print("ConjugateHeatTransferModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
