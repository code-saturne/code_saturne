# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2008 EDF S.A., France
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

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Variables, Model
from Base.XMLmodel import ModelTest
from Pages.OutputControlModel import OutputControlModel

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
        self.__case = case

        self.__node_models = self.__case.xmlGetNode('thermophysical_models')
        self.__node_cht    = self.__node_models.xmlInitNode('conjugate_heat_transfer')
        self.__node_syr    = self.__node_cht.xmlInitNode('external_coupling', status="on")


    def defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['syrthes_name']       = "SYRTHES"
        default['syrthes_app_num']    = 0
        default['projection_axis']    = "off"
        default['selection_criteria'] = "all[]"
        return default


    def __getNumberOfSyrthesCoupling(self):
        return len(self.__node_syr.xmlGetNodeList('syrthes'))


    def setConjugateHeatTransferStatus(self, status):
        """
        Update the 'Conjugate heat transfer' status.
        @type status: C{String}
        @param status: set to "on" or "off" the conjugate heat transfer
        """
        self.isOnOff(status)
        self.__node_syr['status'] = status

        if status == "off":
            self.__node_syr.xmlRemoveChild('syrthes')
        OutputControlModel(self.__case).setSyrthesBoundaryPostProStatus(status)


    def getSyrthesCouplingList(self):
        """
        @return: list of Syrthes coupling description.
        @rtype: C{List}
        """
        node_list = self.__node_syr.xmlGetNodeList('syrthes')
        list = []
        for index in range(len(node_list)):
            num = index + 1
            syrthes_name = self.getSyrthesInstanceName(num)
            app_num      = self.getSyrthesAppNumber(num)
            proj_axis    = self.getSyrthesProjectionAxis(num)
            location     = self.getSelectionCriteria(num)
            list.append([syrthes_name, app_num, proj_axis, location])

        return list


    def addSyrthesCoupling(self, syrthes_name, app_num, proj_axis, location):
        """
        Add a new definition of a Syrthes coupling.

        @type syrthes_name: C{String}
        @param syrthes_name: Syrthes instance name
        @type app_num: C{Int}
        @param app_num: Syrthes Application number
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
        self.setSyrthesAppNumber(num, app_num)
        self.setSyrthesProjectionAxis(num, proj_axis)
        self.setSelectionCriteria(num, location)
        self.setConjugateHeatTransferStatus("on")

        return num


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


    def getSyrthesInstanceName(self, num):
        """
        Get value of Syrthes instance name.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes Application number
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'syrthes_name',
                                    self.setSyrthesInstanceName)

    #------------------------------------------------------------------
    # Syrthes application number
    #------------------------------------------------------------------
    def setSyrthesAppNumber(self, num, value):
        """
        Set value of Syrthes Application number.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @type value: C{Int}
        @param value: Syrthes Application number
        """
        self.isLowerOrEqual(num, self.__getNumberOfSyrthesCoupling())
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node = self.__node_syr.xmlGetNodeList('syrthes')[num-1]
        node.xmlSetData('syrthes_app_num', value)


    def getSyrthesAppNumber(self, num):
        """
        Get value of Syrthes Application number.

        @type num: C{Int}
        @param num: Syrthes coupling number
        @return: Syrthes Application number
        @rtype: C{Int}
        """
        return self.__getIntData(num-1,
                                 'syrthes_app_num',
                                 self.setSyrthesAppNumber)

    #------------------------------------------------------------------
    # Projection axis
    #------------------------------------------------------------------
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
                                    self.setSyrthesAppNumber)

    #------------------------------------------------------------------
    # Selection criteria
    #------------------------------------------------------------------
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
