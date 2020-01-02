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
This module modifies the parameters used by the script file
- ScriptRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, os.path, shutil, sys, types, re

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from code_saturne.model.SolutionDomainModel import MeshModel, SolutionDomainModel
from code_saturne.model.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.model.XMLvariables import Variables, Model

#-------------------------------------------------------------------------------
# Class ScriptRunningModel
#-------------------------------------------------------------------------------

class ScriptRunningModel(Model):
    """
    This class modifies some running model options
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_mgt = self.case.xmlInitNode('calculation_management')

        self.parameters = None


    @Variables.noUndo
    def getPreproMode(self):
        """
        Get run type.
        """
        val = self.node_mgt.xmlGetString('run_type')

        if not val or val == 'standard':
            return False
        else:
            return True


    @Variables.noUndo
    def getTrace(self):
        """
        Get logging options.
        """
        trace = False
        trace_type = 'no'

        node = self.node_mgt.xmlGetNode('logging')
        if node:
            trace_type = node['main']

        if not trace_type:
            if self.case['salome']:
                trace_type = 'stdout'

        if trace_type == 'stdout':
            trace = True

        return trace


    @Variables.undoLocal
    def setTrace(self, trace):
        """
        Set tracing options.
        """

        node = self.node_mgt.xmlInitNode('logging')

        if not trace and not self.getLogParallel():
            if node:
                node.xmlRemoveNode()
        else:
            if trace:
                node['main'] = 'stdout'
            else:
                del node['main']


    @Variables.noUndo
    def getLogParallel(self):
        """
        Get logging options.
        """
        logp = False

        node = self.node_mgt.xmlGetNode('logging')
        if node:
            if node['parallel'] and node['parallel'] != 'null':
                logp = True

        return logp


    @Variables.undoLocal
    def setLogParallel(self, logp):
        """
        Set logging options.
        """
        node = self.node_mgt.xmlInitNode('logging')

        if not logp and not self.getTrace():
            if node:
                node.xmlRemoveNode()
        else:
            if logp:
                node['parallel'] = 'listing'
            else:
                del node['parallel']


    @Variables.noUndo
    def getString(self, key):
        """
        Get entry by named string.
        """
        val = self.node_mgt.xmlGetString(key)
        if not val:
            val = ''

        return val


    @Variables.undoLocal
    def setString(self, key, string):
        """
        Set entry by named string.
        """
        if not string:
            node = self.node_mgt.xmlGetNode(key)
            if node:
                node.xmlRemoveNode()
        else:
            self.node_mgt.xmlSetData(key, string)


#-------------------------------------------------------------------------------
# ScriptRunningModel test class
#-------------------------------------------------------------------------------

class ScriptRunningModelTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all 'check' methods.
        """
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        from code_saturne.model.Common import GuiParam
        GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case).initialize()

        domain = SolutionDomainModel(self.case)
        domain.addMesh('mail1.des', 'des')
        domain.addMesh('mail2.des', 'des')
        domain.addMesh('mail3.des', 'des')

        self.case['xmlfile'] = 'NEW.xml'

def suite():
    testSuite = unittest.makeSuite(ScriptRunningModelTestCase, "check")
    return testSuite


def runTest():
    print("ScriptRunningModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of ScriptRunningModel
#-------------------------------------------------------------------------------
