# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
This module defines the 'PerformanceTuning' page.

This module defines the following classes:
- PerformanceTuningModel
- PerformanceTuningTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# PerformanceTuning model class
#-------------------------------------------------------------------------------

class PerformanceTuningModel(Model):
    """
    Manage the input/output markups in the xml doc about PerformanceTuninging
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case
        node_mgt = self.case.xmlInitNode('calculation_management')
        self.node_part = node_mgt.xmlInitNode('partitioning')


    def _defaultPartitionValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['partition_input'] = "off"
        default['ignore_periodicity'] = "off"
        return default


    def getPartitionInputPath(self):
        """
        Return restart path if applicable
        """
        node = self.node_part.xmlInitNode('partition_input', 'path')
        partition_input = node['path']
        if not partition_input:
            partition_input = None
            self.setPartitionInputPath(partition_input)
        return partition_input


    def setPartitionInputPath(self, v):
        """
        Set partition path if applicable
        """
        node = self.node_part.xmlInitNode('partition_input', 'path')
        if v:
            node['path'] = v
        else:
            node.xmlRemoveNode()


    def getPartitionInputPath(self):
        """
        Return restart path if applicable
        """
        node = self.node_part.xmlInitNode('partition_input', 'path')
        partition_input = node['path']
        if not partition_input:
            partition_input = None
            self.setPartitionInputPath(partition_input)
        return partition_input


    def getPartitionType(self):
        """
        Get partition type.
        """
        val = self.node_part.xmlGetString('type')
        if not val:
            val = 'default'

        return val


    def setPartitionType(self, p):
        """
        Set partition type.
        """
        self.isInList(p, ('default', 'scotch', 'metis',
                          'morton sfc', 'morton sfc cube',
                          'hilbert sfc', 'hilbert sfc cube', 'block'))
        if p == 'default':
            node = self.node_part.xmlGetNode('type')
            if node:
                node.xmlRemoveNode()
        else:
            self.node_part.xmlSetData('type', p)


    def getPartitionOut(self):
        """
        Get partition type.
        """
        val = self.node_part.xmlGetString('output')
        if not val:
            val = 'default'

        return val


    def setPartitionOut(self, p):
        """
        Set partition type.
        """
        self.isInList(p, ('no', 'default', 'yes'))
        if p == 'default':
            node = self.node_part.xmlGetNode('output')
            if node:
                node.xmlRemoveNode()
        else:
            self.node_part.xmlSetData('output', p)


    def getPartitionList(self):
        """
        Get partitions list.
        """
        val = self.node_part.xmlGetString('partition_list')
        if not val:
            val = ''

        return val


    def setPartitionList(self, parts):
        """
        Set partitions list.
        """
        if not parts:
            node = self.node_part.xmlGetNode('partition_list')
            if node:
                node.xmlRemoveNode()
        else:
            self.node_part.xmlSetData('partition_list', parts)


    def getPartitionRankStep(self):
        """
        Get partitions list.
        """
        val = self.node_part.xmlGetString('rank_step')
        if not val:
            val = 1

        return val


    def setPartitionRankStep(self, rank_step):
        """
        Set partitions list.
        """
        if rank_step < 2:
            node = self.node_part.xmlGetNode('rank_step')
            if node:
                node.xmlRemoveNode()
        else:
            self.node_part.xmlSetData('rank_step', rank_step)


    def getIgnorePerio(self):
        """
        """
        node = self.node_part.xmlInitNode('ignore_periodicity', 'status')
        status = node['status']
        if not status:
            v = self._defaultPartitionValues()['ignore_periodicity']
            self.setIgnorePerio(v)
        return status


    def setIgnorePerio(self, v):
        """
        """
        self.isOnOff(v)
        node = self.node_part.xmlInitNode('ignore_periodicity', 'status')
        if v == 'on':
            node['status'] = v
        else:
            node.xmlRemoveNode()

#-------------------------------------------------------------------------------
# PartitionModel test case
#-------------------------------------------------------------------------------


class PerformanceTuningTestCase(ModelTest):
    """
    """
    def checkPartitionInstantiation(self):
        """
        Check whether the PartitionModel class could be instantiated
        """
        model = None
        model = PartitionModel(self.case)
        assert model != None, 'Could not instantiate PartitionModel'

    def checkSetandGetPartInput(self):
        """
        Check whether the partition method could be set and get
        """
        model = PartitionModel(self.case)
        model.setPartitionInputPath("RESU/test/partition_output")
        doc= '''<partitioning>
                    <partition path="RESU/test/partition"/>
                </partitioning>'''

        assert model.node_part == self.xmlNodeFromString(doc),\
                    'Could not set partition in Partition model'
        assert model.getPartitionInputPath() == 'RESU/test/partition_output',\
                    'Could not get partition in Partition model'

def runTest():
    print("PerformanceTuningTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
