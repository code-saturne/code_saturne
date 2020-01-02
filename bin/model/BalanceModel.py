# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
This module defines the XML calls for preprocessor execution
This module contains the following classes and function:
- BalanceModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.OutputVolumicVariablesModel import OutputVolumicVariablesModel

#-------------------------------------------------------------------------------
# Class Balance Model
#-------------------------------------------------------------------------------

class BalanceModel(Variables, Model):

    """
    This class manages the scalar balance objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case = case
        self.node_ana        = self.case.xmlGetNode('analysis_control')
        self.node_scal_bal   = self.node_ana.xmlInitNode('scalar_balances')
        self.node_model      = self.case.xmlInitNode('thermophysical_models')
        self.node_model_vp   = self.node_model.xmlInitNode('velocity_pressure')
        self.node_var_vp     = self.node_model_vp.xmlGetNodeList('variable')
        self.__var_prop_list = self.getScalarVariables()


    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['selector']       = "all[]"

        return defvalue


    @Variables.noUndo
    def getScalarVariables(self):
        """
        Creates a dictionnary to connect name and label of
        scalar variables (different pressure).
        """
        self.dicoLabel2Name = {}
        model = XMLmodel(self.case)
        output = OutputVolumicVariablesModel(self.case)
        for nodeList in [self.node_var_vp,
                         model.getTurbVariable(),
                         output.getModelVariables('solid_fuels'),
                         output.getModelVariables('gas_combustion'),
                         output.getModelVariables('atmospheric_flows'),
                         output.getModelVariables('joule_effect'),
                         output.getThermalScalar(),
                         output.getAdditionalScalar()]:

            for node in nodeList:
                name = node['name']
                label = node['label']
                if not label:
                    raise ValueError("Node has no label")

                dim = node['dimension']
                if not dim or int(dim) == 1:
                    if name != "pressure":
                        self.dicoLabel2Name[label] = (name, str(0))

        return list(self.dicoLabel2Name.keys())


    @Variables.undoGlobal
    def addPressureDrop(self):
        """
        """
        nb = len(self.getPressureDropList())
        node = self.node_scal_bal.xmlAddChild('pressure_drop', pres_id = nb)


    @Variables.undoGlobal
    def delPressureDrop(self, idx):
        """
        """
        lst = self.getPressureDropList()
        n = lst[idx]
        for node in lst:
            nodeid = int(node['pres_id'])
            if nodeid > int(idx):
                node['pres_id'] = str(nodeid-1)
        n.xmlRemoveNode()


    @Variables.noUndo
    def getPressureDropList(self):
        """
        """
        listNode = self.node_scal_bal.xmlGetNodeList('pressure_drop')

        return listNode


    @Variables.noUndo
    def getPressureDropCriteria(self, idx):
        """
        """
        node = self.node_scal_bal.xmlGetNode('pressure_drop', pres_id = idx)
        n = node.xmlGetNode('criteria')
        if not n or n == "None":
            crit = self.defaultValues()['selector']
            self.setPressureDropCriteria(idx, crit)
            n = node.xmlGetNode('criteria')
        criteria = str(n.xmlGetTextNode())

        return criteria


    @Variables.undoLocal
    def setPressureDropCriteria(self, idx, criteria):
        """
        """
        node = self.node_scal_bal.xmlGetNode('pressure_drop', pres_id = idx)
        node.xmlSetData('criteria', criteria)


    @Variables.undoGlobal
    def addScalarBalance(self):
        """
        """
        nb = len(self.getScalarBalanceList())
        node = self.node_scal_bal.xmlAddChild('scalar_balance', balance_id = nb)


    @Variables.undoGlobal
    def deleteScalarBalance(self, idx):
        """
        """
        lst = self.getScalarBalanceList()
        n = lst[idx]
        for node in lst:
            nodeid = int(node['balance_id'])
            if nodeid > int(idx):
                node['balance_id'] = str(nodeid-1)
        n.xmlRemoveNode()


    @Variables.noUndo
    def getScalarBalanceList(self):
        """
        """
        listNode = self.node_scal_bal.xmlGetNodeList('scalar_balance')

        return listNode


    @Variables.noUndo
    def getScalarBalanceCriteria(self, idx):
        """
        """
        node = self.node_scal_bal.xmlGetNode('scalar_balance', balance_id = idx)
        n = node.xmlGetNode('criteria')
        if not n or n == "None":
            crit = self.defaultValues()['selector']
            self.setScalarBalanceCriteria(idx, crit)
            n = node.xmlGetNode('criteria')
        criteria = str(n.xmlGetTextNode())

        return criteria


    @Variables.undoLocal
    def setScalarBalanceCriteria(self, idx, criteria):
        """
        """
        node = self.node_scal_bal.xmlGetNode('scalar_balance', balance_id = idx)
        node.xmlSetData('criteria', criteria)


    @Variables.undoLocal
    def setVariable(self, idx, lst):
        """
        Public method.
        """
        node = self.node_scal_bal.xmlGetNode('scalar_balance', balance_id = idx)
        node.xmlRemoveChild('var_prop')
        for var in lst:
            self.isInList(var, self.__var_prop_list)
            (name, comp) = self.dicoLabel2Name[var]
            node.xmlAddChild('var_prop', name=name, component=comp)


    @Variables.noUndo
    def getVariable(self, idx):
        """
        Public method.
        """
        node = self.node_scal_bal.xmlGetNode('scalar_balance', balance_id = idx)

        lst = []
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['component']) :
                    lst.append(name)
        return lst


#-------------------------------------------------------------------------------
# Balance test case
#-------------------------------------------------------------------------------
class BalanceTestCase(ModelTest):
    """
    """


def suite():
    testSuite = unittest.makeSuite(BalanceTestCase, "check")
    return testSuite


def runTest():
    print("BalanceTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
