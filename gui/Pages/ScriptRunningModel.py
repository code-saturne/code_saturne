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
This module modifies the parameters used by the script file
- ScriptRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, os.path, shutil, sys, string, types, re

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import code_saturne.Base.Toolbox as Tool
from code_saturne.Pages.SolutionDomainModel import MeshModel, SolutionDomainModel
from code_saturne.Pages.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.Base.XMLvariables import Variables, Model

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
    def getRunType(self, preproview = None):
        """
        Get run type.
        """
        val = self.node_mgt.xmlGetString('run_type')

        if not preproview:
            val = 'standard'

        if not val:
            if not preproview:
                val = 'standard'
            else:
                val = 'mesh quality'
            self.setRunType(val)
        elif (preproview and val == 'standard'):
            val = 'mesh quality'
            self.setRunType(val)

        return val


    @Variables.undoLocal
    def setRunType(self, run):
        """
        Set run type.
        """
        self.isInList(run, ('none', 'mesh preprocess', 'mesh quality', 'standard'))
        self.node_mgt.xmlSetData('run_type', run)


    @Variables.noUndo
    def getLogType(self):
        """
        Get logging options.
        """
        log_type = ['listing', 'null']

        node = self.node_mgt.xmlGetNode('logging')
        if node:
            log_type[0] = node['main']
            log_type[1] = node['parallel']

        if not log_type[0]:
            if self.case['salome']:
                log_type[0] = 'stdout'
            else:
                log_type[0] = 'listing'

        if not log_type[1]:
            log_type[1] = 'null'

        return log_type


    @Variables.undoLocal
    def setLogType(self, log_type):
        """
        Set logging options.
        """
        if log_type[0] == 'listing' and not self.case['salome']:
            log_type[0] = None
        if log_type[1] == 'null':
            log_type[1] = None

        node = self.node_mgt.xmlInitNode('logging')

        if not log_type[0] and not log_type[1]:
            if node:
                node.xmlRemoveNode()
        else:
            if log_type[0]:
                node['main'] = log_type[0]
            else:
                del node['main']
            if log_type[1]:
                node['parallel'] = log_type[1]
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
        from code_saturne.Base.XMLengine import Case
        from code_saturne.Base.XMLinitialize import XMLinit
        from code_saturne.Base.Toolbox import GuiParam
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
