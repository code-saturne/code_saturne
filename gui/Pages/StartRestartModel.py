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
This module defines the 'Start/Restart' page.

This module defines the following classes:
- StartRestartModel
- StartRestartTestCase
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
        node_magt = self.case.xmlInitNode('calcul_management')
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


    def getRestart(self):
        """
        Return if this calcul is restarted with a restart file
        """
        node = self.node_start.xmlInitNode('restart', 'status')
        restart = node['status']
        if not restart:
            restart = self._defaultStartRestartValues()['restart']
            self.setRestart(restart)
        return restart


    def setRestart(self, v):
        """
        Set status if we continue calculation or not
        """
        self.isOnOff(v)
        node = self.node_start.xmlInitNode('restart', 'status')
        node['status'] = v
        if v == 'off':
            self.node_start.xmlRemoveChild('current_restart')
            for n in self.case.xmlGetNodeList('time_average'):
                n.xmlRemoveChild('restart_from_time_average')



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


    def setFrozenField(self, v):
        """
        """
        self.isOnOff(v)
        node = self.node_start.xmlInitNode('frozen_field', 'status')
        node['status'] = v


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
            if val == -1:
                period = "At the end"
            else:
                period = "Frequency"
        return val, period


    def setRestartWithAuxiliaryStatus(self, status):
        """
        Input status of reading auxiliary restart file for advanced options.
        """
        self.isOnOff(status)
        node = self.node_start.xmlInitNode('restart_with_auxiliary', 'status')
        node['status'] = status


    def setRestartRescue(self, freq):
        """
        Inputfrequency for restart checkpoints from advanced options.
        """
        self.isInt(freq)
        self.node_start.xmlSetData('restart_rescue', freq)


    def getRestartDirectory(self):
        """ Convenient method only for the View """
        return self.node_start.xmlGetString('current_restart')


    def setRestartDirectory(self, dir):
        """ Convenient method only for the View """
#        if not os.path.isdir(self.case['resu_path'] + "/" + dir):
#            raise ValueError,  "Invalid restart directory %s" % dir
        self.node_start.xmlSetData('current_restart', dir)


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
        model.setRestart("on")
        doc= '''<start_restart>
                    <restart status="on"/>
                </start_restart>'''

        assert model.node_start == self.xmlNodeFromString(doc),\
                    'Could not set restart in StartRestart model'
        assert model.getRestart() == 'on',\
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

##    def checkSetandGetMainandAxiliaryRestartFormat(self):
##        """
##        Check whether the Main and Auxiliary Restart format method
##        could be set and get
##        """
##        model = StartRestartModel(self.case)
##        model.setRestart("on")
##        model.setRestartWithAuxiliaryStatus('on')
##        model.setMainRestartFormat('binary')
##        model.setAuxiliaryRestartFormat('ascii')
##        doc= '''<start_restart>
##                    <restart status="on"/>
##                    <restart_with_auxiliary status="on"/>
##                    <main_restart format="binary"/>
##                    <auxiliary_restart format="ascii"/>
##                </start_restart>'''
##        assert model.node_start == self.xmlNodeFromString(doc),\
##                'Could not set main and auxiliary restart formats in StartRestart model'
##        assert model.getMainRestartFormat() == "binary",\
##                'Could not get main restart format in StartRestart model'
##        assert model.getAuxiliaryRestartFormat() == "ascii",\
##                'Could not get auxiliary restart format in StartRestart model'

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

    def checksetRestartDirectory(self):
        """
        Check whether the setRestartDirectory method could be set
        """
        model = StartRestartModel(self.case)
        model.setRestart("on")
        model.setRestartDirectory('RESTART.11221504')
        doc = '''<start_restart>
                    <restart status="on"/>
                    <current_restart>RESTART.11221504</current_restart>
                 </start_restart>'''
        assert model.node_start == self.xmlNodeFromString(doc),\
                'Could not update restart directory in StartRestart model'

        model.setRestart("off")
        doc = '''<start_restart>
                    <restart status="off"/>
                 </start_restart>'''
        assert model.node_start == self.xmlNodeFromString(doc),\
                'Could not remove the restart directory in StartRestart model'


def suite():
    testSuite = unittest.makeSuite(StartRestartTestCase, "check")
    return testSuite


def runTest():
    print "StartRestartTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
