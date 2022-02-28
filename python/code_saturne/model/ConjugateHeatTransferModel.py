# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest

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
        self.__node_syr = self.__node_cht.xmlInitNode('external_coupling')
        self.__node_inst = self.__node_syr.xmlInitNode('syrthes_instances')

    @Variables.undoLocal
    def deleteConjugateHeatTransfer(self):
        """
        Update the 'Conjugate heat transfer' status.
        @type status: C{String}
        @param status: set to "on" or "off" the conjugate heat transfer
        """
        self.__node_inst.xmlRemoveNode()

    def getNumberOfSyrthesCoupling(self):
        return len(self.__node_inst.xmlGetNodeList('instance'))

    @Variables.noUndo
    def getSyrthesInstancesList(self):
        """
        @return: list of syrthes instances
        """

        lst = []
        node_list = self.__node_inst.xmlGetNodeList('instance')

        for node in node_list:
            syrthes_name = node.xmlGetAttribute("name")
            if syrthes_name and syrthes_name not in lst:
                lst.append(syrthes_name)

        return lst

    @Variables.noUndo
    def getSyrthesCouplingList(self):
        """
        @return: list of Syrthes coupling description.
        @rtype: C{List}
        """
        node_list = self.__node_inst.xmlGetNodeList('instance')
        lst = []
        for node in node_list:
            syrthes_name = node.xmlGetAttribute("name")
            boundary_labels = self.getBoundaryLabelList(syrthes_name)
            lst.append([syrthes_name, boundary_labels])
        return lst

    @Variables.undoGlobal
    def addSyrthesCoupling(self, syrthes_name):
        """
        Add a new definition of a Syrthes coupling.

        @type syrthes_name: C{String}
        @param syrthes_name: Syrthes instance name
        """
        node_new = self.__node_inst.xmlInitChildNode('instance', name=syrthes_name)
        return


    @Variables.undoGlobal
    def deleteSyrthesCoupling(self, syrthes_name, boundary_label):
        """
        Delete a syrthes coupling. If the syrthes instance has no
        other coupled boundaries, remove it from the instances list.
        """

        node = self.__node_inst.xmlGetChildNode('instance', name=syrthes_name)
        if node:
            # First delete coupling
            for node_b in node.xmlGetChildNodeList("coupled_boundary"):
                if node_b.xmlGetAttribute("label") == boundary_label:
                    node_b.xmlRemoveNode()
                    break

            # If instance has no other coupling, delete it
            if not node.xmlGetChildNodeList("coupled_boundary"):
                node.xmlRemoveNode()



    def addBoundaryLabel(self, syrthes_name, boundary_label):

        # If boundary allready added, exit
        if boundary_label in self.getBoundaryLabelList(syrthes_name):
            return

        node_instance = self.__node_inst.xmlGetChildNode('instance', name=syrthes_name)
        if node_instance is None:
            raise ValueError("SYRTHES instance not found : ", syrthes_name)
        node_zone = node_instance.xmlInitChildNode('coupled_boundary', label=boundary_label)

    def getBoundaryLabelList(self, syrthes_name):
        node_instance = self.__node_inst.xmlGetChildNode('instance', name=syrthes_name)
        if node_instance is None:
            raise ValueError("SYRTHES instance not found : ", syrthes_name)
        return [n.xmlGetAttribute("label") for n in node_instance.xmlGetChildNodeList("coupled_boundary")]

    def flushBoundaryLabels(self):
        for node_instance in self.__node_inst.xmlGetChildNodeList('instance'):
            node_instance.xmlRemoveChildren()

    @Variables.undoLocal
    def setSyrthesVerbosity(self, value):
        self.__node_syr.xmlSetData('verbosity', value)

    @Variables.noUndo
    def getSyrthesVerbosity(self):
        return self.__node_syr.xmlGetChildString('verbosity')

    @Variables.undoLocal
    def setSyrthesVisualization(self, value):
        self.__node_syr.xmlSetData('visualization', value)

    @Variables.noUndo
    def getSyrthesVisualization(self):
        return self.__node_syr.xmlGetChildString('visualization')

    @Variables.undoLocal
    def setSyrthesProjectionAxis(self, value):
        self.__node_syr.xmlSetData('projection_axis', value)

    @Variables.noUndo
    def getSyrthesProjectionAxis(self):
        return self.__node_syr.xmlGetChildString('projection_axis')

    @Variables.undoLocal
    def setSyrthesTolerance(self, value):
        self.__node_syr.xmlSetData('tolerance', value)

    @Variables.noUndo
    def getSyrthesTolerance(self):

        val = self.__node_syr.xmlGetChildString('tolerance')
        if not val or val == "":
            val = "0.1"
            self.setSyrthesTolerance(val)

        return val

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
