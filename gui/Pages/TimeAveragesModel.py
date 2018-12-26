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
This module defines the Time averages page.

This module defines the following classes:
- TimeAveragesModel
- TimeAveragesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
from code_saturne.Pages.StartRestartModel import StartRestartModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class TimeAveragesModel(Model):
    """
    Class to open TimeAverages Page.
    """
    def __init__(self, case):
        """
        Simple constructor.
        """
        self.case = case
        self.node_anal     = self.case.xmlInitNode('analysis_control')
        self.node_mean     = self.node_anal.xmlInitNode('time_averages')

        self.__var_prop_list = self.__updateDicoLabel2Name()


    def defaultValues(self):
        """
        Public method.
        @return: dictionary with default values.
        @rtype: C{Dict}
        """
        value = {}
        value['start']     = 1
        value['timestart'] = 0.
        value['restart']   = -2
        value['name']     = "TimeAverage"

        return value


    def __updateDicoLabel2Name(self):
        """
        Private method.
        Gets a dictionnary to connect name and label from
        variables and properties.
        """
        mdl = OutputVolumicVariablesModel(self.case)

        self.dicoLabel2Name = mdl.getVolumeFieldsLabel2Name(time_averages=False)

        return self.dicoLabel2Name


    @Variables.undoGlobal
    def addTimeAverage(self):
        """
        Public method.
        Add a new time average and return a default name
        """
        name = self.defaultValues()['name']
        def_name = name + str(len(self.getTimeAverageNames()) + 1)

        # define default name
        if def_name in self.getTimeAverageNames():
            i = 2
            while def_name in self.getTimeAverageNames():
                def_name = name + str(len(self.getTimeAverageNames()) + i)
                i = i + 1

        node = self.node_mean.xmlInitNode('time_average', name = def_name, label = def_name)
        node['id'] = len(self.getTimeAverageNames())
        ntdmom = self.defaultValues()['start']
        ttdmom = -1.0
        imoold = self.defaultValues()['restart']
        self.setTimeStepStart(def_name, ntdmom)
        self.setTimeStart(def_name, ttdmom)
        self.setRestart(def_name, imoold)

        return def_name, ntdmom, ttdmom, imoold


    @Variables.undoGlobal
    def deleteTimeAverage(self, name):
        """
        Public method.
        @type name: C{String}
        @param name: name of the time average to delete.
        """
        node = self.node_mean.xmlGetNode('time_average', name=name)
        if node:
            nb = node['id']
            node.xmlRemoveNode()

            # renumbering of all time averages
            for p in range(int(nb)+1, self.getNumberOfTimeAverage()+2):
                t = self.node_mean.xmlGetNode('time_average', id=p)
                t['id'] = p - 1


    @Variables.noUndo
    def getTimeAverageData(self, imom):
        """
        Public method.
        @return: data for time average number I{imom}.
        @rtype: C{Tuple}
        """
        self.isInt(imom)
        lst = []
        restart = self.defaultValues()['restart']
        node = self.node_mean.xmlGetNode('time_average', id=imom)
        start = node.xmlGetInt('time_step_start')
        timestart = node.xmlGetDouble('time_start')
        restart = node.xmlGetInt('restart_from_time_average')

        for var in node.xmlGetChildNodeList('var_prop'):
            for name in list(self.dicoLabel2Name.keys()):
                if self.dicoLabel2Name[name] == (var['name'], var['component']):
                    lst.append(name)
        return node['name'], start, timestart, restart, lst


    @Variables.undoLocal
    def setName(self, old_name, name):
        """
        Public method.
        """
        self.isInList(old_name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=old_name)
        node['name'] = name
        node['label'] = name


    @Variables.undoLocal
    def setTimeStart(self, name, start):
        """
        Public method.
        """
        self.isFloat(start)
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        node.xmlSetData('time_start', start)


    @Variables.noUndo
    def getTimeStart(self, name):
        """
        Public method.
        """
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        return node.xmlGetDouble('time_start')


    @Variables.undoLocal
    def setTimeStepStart(self, name, start):
        """
        Public method.
        """
        self.isInt(start)
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        node.xmlSetData('time_step_start', start)


    @Variables.noUndo
    def getTimeStepStart(self, name):
        """
        Public method.
        """
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        return node.xmlGetInt('time_step_start')


    @Variables.undoLocal
    def setRestart(self, name, restart):
        """
        Public method.
        """
        self.isInt(restart)
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        if restart != -2:
            node.xmlSetData('restart_from_time_average', restart)
        else:
            node.xmlRemoveChild('restart_from_time_average')


    @Variables.noUndo
    def getRestart(self, name):
        """
        Public method.
        """
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        return node.xmlGetInt('restart_from_time_average')


    @Variables.undoLocal
    def setVariable(self, name, lst):
        """
        Public method.
        """
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)
        node.xmlRemoveChild('var_prop')
        for var in lst:
            self.isInList(var, self.__var_prop_list)
            (name, comp) = self.dicoLabel2Name[var]
            node.xmlAddChild('var_prop', name=name, component=comp)


    @Variables.noUndo
    def getVariable(self, name):
        """
        Public method.
        """
        self.isInList(name, self.getTimeAverageNames())
        node = self.node_mean.xmlInitNode('time_average', name=name)

        lst = []
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['component']) :
                    lst.append(name)
        return lst


    @Variables.noUndo
    def getTimeAverageNames(self):
        """
        Public method.
        @return: list of time averages names.
        @rtype: C{List} of C{String}
        """
        lst = []
        for node in self.node_mean.xmlGetNodeList('time_average'):
            name = node['name']
            lst.append(name)
        return lst


    @Variables.noUndo
    def getNumberOfTimeAverage(self):
        """
        Public method.
        @return: number of time averages already defined.
        @rtype: C{Int}
        """
        return len(self.node_mean.xmlGetNodeList('time_average'))

#-------------------------------------------------------------------------------
# TimeAveragesModel test case
#-------------------------------------------------------------------------------

class TimeAveragesTestCase(ModelTest):
    """
    unittest
    """
    def checkTimeAveragesInstantiation(self):
        """Check whether the TimeAveragesModel class could be instantiated"""
        model = None
        model = TimeAveragesModel(self.case)
        assert model != None, 'Could not instantiate TimeAveragesModel'


    def checksetTimeAverage(self):
        """Check whether the TimeAveragesModel class could be set a average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setTimeAverage(2, 'deux', 20, 1, ['Pressure'])

        doc = '''<time_averages>
                    <time_average id="1" name="moyenne">
                            <var_prop name="velocity" component="0"/>
                            <var_prop name="velocity" component="1"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" name="deux">
                            <var_prop name="pressure"/>
                            <time_step_start>20</time_step_start>
                    </time_average>
                 </time_averages>'''

        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not set some averages in TimeAveragesModel'


    def checkreplaceTimeAverage(self):
        """Check whether the TimeAveragesModel class could be replaced one average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setTimeAverage(2, 'deux', 20, 1, ['Pressure'])
        mdl.setTimeAverage(3, 'trois', 33, 1, ['Pressure', 'VelocitU'])
        mdl.replaceTimeAverage(2, 'SECOND', 12, 1, ['VelocitW', 'VelocitV'])
        mdl.replaceTimeAverage(3, 'trois', 33, 1, ['Pressure', 'VelocitW'])
        doc = '''<time_averages>
                    <time_average id="1" name="moyenne">
                            <var_prop name="velocity" component="0"/>
                            <var_prop name="velocity" component="1"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" name="SECOND">
                            <var_prop name="velocity" component="2"/>
                            <var_prop name="velocity" component="1"/>
                            <time_step_start>12</time_step_start>
                    </time_average>
                    <time_average id="3" name="trois">
                            <var_prop name="pressure"/>
                            <var_prop name="velocity" component="2"/>
                            <time_step_start>33</time_step_start>
                    </time_average>
                 </time_averages>'''
        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not replace one average in TimeAveragesModel'


    def checkdeleteTimeAverage(self):
        """Check whether the TimeAveragesModel class could be deleted one average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setTimeAverage(2, 'deux', 20, 1, ['VelocitU'])
        mdl.setTimeAverage(3, 'trois', 33, 1, ['VelocitW', 'VelocitU'])
        mdl.deleteTimeAverage(2)
        doc = '''<time_averages>
                    <time_average id="1" name="moyenne">
                            <var_prop name="velocity" component="0"/>
                            <var_prop name="velocity" component="1"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" name="trois">
                            <var_prop name="velocity" component="2"/>
                            <var_prop name="velocity" component="0"/>
                            <time_step_start>33</time_step_start>
                    </time_average>
                 </time_averages>'''

        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not delete one average in TimeAveragesModel'


def suite():
    testSuite = unittest.makeSuite(TimeAveragesTestCase, "check")
    return testSuite


def runTest():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
