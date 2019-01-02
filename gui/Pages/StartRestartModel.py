# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the 'Start/Restart' page.

This module defines the following classes:
- StartRestartModel
- StartRestartTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types
import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Start-Restart model class
#-------------------------------------------------------------------------------

class StartRestartModel(Model):
    """
    Manage the input/output markups in the xml doc about Start and Restart
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case
        node_magt = self.case.xmlInitNode('calculation_management')
        self.node_start = node_magt.xmlInitNode('start_restart')


    def _defaultStartRestartValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['restart']                = "off"
        default['frozen_field']           = "off"
        default['restart_with_auxiliary'] = "on"
        default['restart_rescue']         = 0
        default['period_rescue']          = "4 output"
        return default


    @Variables.noUndo
    def getRestartPath(self):
        """
        Return restart path if applicable; use '*' for automatic mode.
        """
        node = self.node_start.xmlInitNode('restart', 'path')
        restart = node['path']
        if not restart:
            restart = None
            self.setRestartPath(restart)
        return restart


    @Variables.undoLocal
    def setRestartPath(self, v):
        """
        Set restart path if applicable; use '*' for automatic mode.
        """
        node = self.node_start.xmlInitNode('restart', 'path')
        if v:
            node['path'] = v
        else:
            node.xmlRemoveNode()
            for n in self.case.xmlGetNodeList('time_average'):
                n.xmlRemoveChild('restart_from_time_average')


    @Variables.noUndo
    def getRestartMeshPath(self):
        """
        Return restart mesh path if applicable
        """
        node = self.node_start.xmlInitNode('restart_mesh', 'path')
        restart_mesh = node['path']
        if not restart_mesh:
            restart_mesh = None
            self.setRestartMeshPath(restart_mesh)
        return restart_mesh


    @Variables.undoLocal
    def setRestartMeshPath(self, v):
        """
        Set restart mesh path if applicable
        """
        node = self.node_start.xmlInitNode('restart_mesh', 'path')
        if v:
            node['path'] = v
        else:
            node.xmlRemoveNode()


    @Variables.noUndo
    def getFrozenField(self):
        """
        Return if the velocity and the pressure are solved
        """
        node = self.node_start.xmlInitNode('frozen_field', 'status')
        status = node['status']
        if not status:
            v = self._defaultStartRestartValues()['frozen_field']
            self.setFrozenField(v)
        return status


    @Variables.undoLocal
    def setFrozenField(self, v):
        """
        """
        self.isOnOff(v)
        node = self.node_start.xmlInitNode('frozen_field', 'status')
        node['status'] = v


    @Variables.noUndo
    def getRestartWithAuxiliaryStatus(self):
        """
        Return status of reading auxiliary restart file for advanced options.
        """
        node = self.node_start.xmlInitNode('restart_with_auxiliary', 'status')
        status = node['status']
        if not status:
            status = self._defaultStartRestartValues()['restart_with_auxiliary']
            self.setRestartWithAuxiliaryStatus(status)
        return status


    @Variables.noUndo
    def getRestartRescue(self):
        """
        Return frequency for restart checkpoints from advanced options.
        """
        val = self.node_start.xmlGetInt('restart_rescue')
        if val == None or val == 0:
            period = self._defaultStartRestartValues()['period_rescue']
            val = self._defaultStartRestartValues()['restart_rescue']
            self.setRestartRescue(val)
        else:
            if val == -2:
                period = "Never"
            elif val == -1:
                period = "At the end"
            else:
                period = "Frequency"
        return val, period


    @Variables.undoLocal
    def setRestartWithAuxiliaryStatus(self, status):
        """
        Input status of reading auxiliary restart file for advanced options.
        """
        self.isOnOff(status)
        node = self.node_start.xmlInitNode('restart_with_auxiliary', 'status')
        node['status'] = status


    @Variables.undoLocal
    def setRestartRescue(self, freq):
        """
        Inputfrequency for restart checkpoints from advanced options.
        """
        self.isInt(freq)
        self.node_start.xmlSetData('restart_rescue', freq)


#-------------------------------------------------------------------------------
# StartRestartModel test case
#-------------------------------------------------------------------------------


class StartRestartTestCase(ModelTest):
    """
    """
    def checkStartRestartInstantiation(self):
        """
        Check whether the StartRestartModel class could be instantiated
        """
        model = None
        model = StartRestartModel(self.case)
        assert model != None, 'Could not instantiate StartRestartModel'

    def checkSetandGetRestart(self):
        """
        Check whether the restart method could be set and get
        """
        model = StartRestartModel(self.case)
        model.setRestartPath("RESU/test/restart")
        doc= '''<start_restart>
                    <restart path="RESU/test/restart"/>
                </start_restart>'''

        assert model.node_start == self.xmlNodeFromString(doc),\
                    'Could not set restart in StartRestart model'
        assert model.getRestartPath() == 'RESU/test/restart',\
                    'Could not get restart in StartRestart model'

    def checkSetandGetFrozenStatus(self):
        """
        Check whether the Frozen status method could be set and get
        """
        model = StartRestartModel(self.case)
        model.setRestart("on")
        model.setFrozenField('on')
        doc = '''<start_restart>
                    <restart status="on"/>
                    <frozen_field status="on"/>
                 </start_restart>'''
        assert model.node_start == self.xmlNodeFromString(doc),\
                'Could not set frozen status in StartRestart model'
        assert model.getFrozenField() == "on",\
                'Could not get frozen status in StartRestart model'

    def checkSetAuxiliaryRestartStatus(self):
        """
        Check whether the  Auxiliary Restart Status method
        could be set and get
        """
        model = StartRestartModel(self.case)
        model.setRestart("on")
        model.setRestartWithAuxiliaryStatus('on')
        doc= '''<start_restart>
                    <restart status="on"/>
                    <restart_with_auxiliary status="on"/>
                </start_restart>'''
        assert model.node_start == self.xmlNodeFromString(doc),\
                'Could not set auxiliary restart status in StartRestart model'
        assert model.getRestartWithAuxiliaryStatus() == "on",\
                'Could not get auxiliary restart status in StartRestart model'

    def checkSetandGetRestartRescue(self):
        """
        Check whether the  Restart rescue could be set and get
        """
        model = StartRestartModel(self.case)
        model.setRestart("on")
        model.setRestartRescue(15)
        doc = '''<start_restart>
                    <restart status="on"/>
                    <restart_rescue>15</restart_rescue>
                 </start_restart>'''
        assert model.node_start == self.xmlNodeFromString(doc),\
                'Could not set restart rescue value in StartRestart model'
        freq, period = model.getRestartRescue()
        assert freq == 15,\
                'Could not get restart rescue value in StartRestart model'
        assert period == "Frequency",\
                'Could not get restart rescue period in StartRestart model'

        model.setRestartRescue(-1)
        freq, period = model.getRestartRescue()
        assert freq == -1,\
                'Could not get restart rescue value in StartRestart model'
        assert period == "At the end",\
                'Could not get restart rescue period in StartRestart model'


def suite():
    testSuite = unittest.makeSuite(StartRestartTestCase, "check")
    return testSuite


def runTest():
    print("StartRestartTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
