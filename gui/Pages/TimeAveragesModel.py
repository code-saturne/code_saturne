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
This module defines the Time averages page.

This module defines the following classes:
- TimeAveragesModel
- TimeAveragesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import XMLmodel, ModelTest
from Base.XMLvariables import Model
from Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
from Pages.StartRestartModel import StartRestartModel

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
        self.node_anal      = self.case.xmlInitNode('analysis_control')
        self.node_mean      = self.node_anal.xmlInitNode('time_averages')
        self.node_model     = self.case.xmlInitNode('thermophysical_models')
        self.node_model_vp  = self.node_model.xmlInitNode('velocity_pressure')
        self.node_var_vp    = self.node_model_vp.xmlGetNodeList('variable')
        self.node_pro_vp    = self.node_model_vp.xmlGetNodeList('property')

        self.__updateDicoLabel2Name()


    def __defaultValues(self):
        """
        Private method.
        Return a dictionnary with default values
        """
        value = {}
        if StartRestartModel(self.case).getRestart() == 'off':
            value['restart']    = 0
        else:
            value['restart']    = -1

        return value


    def __updateDicoLabel2Name(self):
        """
        Private method.
        Gets a dictionnary to connect name and label from
        variables and properties.
        """
        # FIXME: merge this method with the same one in ProfilesView
        self.dicoLabel2Name = {}
        model = XMLmodel(self.case)
        output = OutputVolumicVariablesModel(self.case)
        for nodeList in [self.node_var_vp,
                         self.node_pro_vp,
                         model.getTurbNodeList(),
                         output.getFluidProperty(),
                         output.getAdditionalScalarProperty(),
                         output.getTimeProperty(),
                         output.getListOfTimeMeans(),
                         output.getPuCoalScalProper(),
                         output.getThermalScalar(),
                         output.getAdditionalScalar()]:

            for node in nodeList:
                name = node['name']
                label = node['label']
                if not label:
                    raise ValueError, "Node has no label"

                if not (node['support'] and node['support'] == "boundary"):
                    if name != 'local_time_step':
                        self.dicoLabel2Name[label] = name

        return self.dicoLabel2Name.keys()


    def ___updateAverageNumbers(self, imom):
        """
        Private method.
        Update order of average.
        """
        for index in self.getAverageList():
            lab, start, restart, list = self.getAverageInformations(index)
            if index >= imom:
                nb = str(index - 1)
                if lab == 'Moy' + str(index):
                    lab = 'Moy' + nb
                node = self.node_mean.xmlGetNode('time_average', id=index)
                node['id'] = nb
                node['label'] = lab


    def __updateAverage(self, nb, label, start, restart, list):
        """
        Private method.
        Update data for average I{llabel}.
        """
        node = self.node_mean.xmlInitNode('time_average', id=nb)
        node['label'] = label
        for var in list:
            self.isInList(var, self.dicoLabel2Name.keys())
            node.xmlAddChild('var_prop', name=self.dicoLabel2Name[var])
        node.xmlSetData('time_step_start', start)
        if StartRestartModel(self.case).getRestart() != 'off' and restart != -1:
            node.xmlSetData('restart_from_time_average', restart)
        else:
            if node.xmlGetInt('restart_from_time_average'):
                node.xmlRemoveChild('restart_from_time_average')


    def setAverage(self, nb, label, start, restart, list):
        """
        Public method.
        Sets list of variables or properties used for calculation of average
        for average number nb, and time's start and restart values.
        """
        self.isNotInList(nb, self.getAverageList())
        for i in (nb, start):
            self.isInt(i)
        if StartRestartModel(self.case).getRestart() != 'off':
            self.isInt(restart)
        self.isNotInList(label, self.getAverageLabelsList())

        self.__updateAverage(nb, label, start, restart, list)


    def replaceAverage(self, nb, label, start, restart, list):
        """
        Public method.
        Replaces list of variables or properties used for calculation of mean
        for average number nb, or label, or time's start or restart values.
        """
        self.isInList(nb, self.getAverageList())
        for i in (nb, start):
            self.isInt(i)
        if StartRestartModel(self.case).getRestart() != 'off':
            self.isInt(restart)

        node = self.node_mean.xmlGetNode('time_average', id=nb)
        if node:
            node.xmlRemoveChild('var_prop')
            node.xmlRemoveChild('time_step_start')
            node.xmlRemoveChild('restart_from_time_average')
        self.__updateAverage(nb, label, start, restart, list)


    def deleteAverage(self, nb):
        """
        Public method.
        Sets list of variables or properties used for calculation of mean
        for average number nb, and time's start and restart values.
        """
        self.isInt(nb)
        node = self.node_mean.xmlGetNode('time_average', id=nb)
        if node:
            node.xmlRemoveNode()
##            node.xmlRemoveChild('var_prop')
##            node.xmlRemoveChild('time_step_start')
##            if node.xmlGetInt('restart_from_time_average'):
##                node.xmlRemoveChild('restart_from_time_average')
        self.___updateAverageNumbers(nb)


    def getAverageList(self):
        """
        Public method.
        Gets list of averages and return list of imom number.
        """
        list = []
        for node in self.node_mean.xmlGetNodeList('time_average'):
            list.append(int(node['id']))

        return list


    def getAverageInformations(self, imom):
        """
        Public method.
        Gets average with imom and return time_step_start, restart_from_time_average,
        and list of variables or properties ....
        """
        self.isInt(imom)
        list = []
        restart = self.__defaultValues()['restart']
        node = self.node_mean.xmlGetNode('time_average', id=imom)
        start = node.xmlGetInt('time_step_start')
        if StartRestartModel(self.case).getRestart() != 'off':
            restart = node.xmlGetInt('restart_from_time_average')

        for var in node.xmlGetChildNodeList('var_prop'):
            for label in self.dicoLabel2Name.keys():
                if self.dicoLabel2Name[label] == var['name']:
                    list.append(label)
        return node['label'], start, restart, list


    def getAverageLabelsList(self):
        """
        Public method.
        Gets list of averages's labels.
        """
        list = []
        for node in self.node_mean.xmlGetNodeList('time_average'):
            label = node['label']
            list.append(label)
        return list


    def getAverageRestart(self, nb):
        """
        Public method. Only for GUI (StartRestartModel class).
        Get restart to know if balise exists.
        """
        self.isInt(nb)
        restart = ''
        node = self.node_mean.xmlGetNode('time_average', id=nb)
        if node:
            restart = node.xmlGetInt('restart_from_time_average')
        return restart

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


    def checkSetAverage(self):
        """Check whether the TimeAveragesModel class could be set a average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setAverage(2, 'deux', 20, 1, ['Pressure'])

        doc = '''<time_averages>
                    <time_average id="1" label="moyenne">
                            <var_prop name="velocity_U"/>
                            <var_prop name="velocity_V"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" label="deux">
                            <var_prop name="pressure"/>
                            <time_step_start>20</time_step_start>
                    </time_average>
                 </time_averages>'''

        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not set some averages in TimeAveragesModel'


    def checkReplaceAverage(self):
        """Check whether the TimeAveragesModel class could be replaced one average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setAverage(2, 'deux', 20, 1, ['Pressure'])
        mdl.setAverage(3, 'trois', 33, 1, ['Pressure', 'VelocitU'])
        mdl.replaceAverage(2, 'SECOND', 12, 1, ['VelocitW', 'VelocitV'])
        mdl.replaceAverage(3, 'trois', 33, 1, ['Pressure', 'VelocitW'])
        doc = '''<time_averages>
                    <time_average id="1" label="moyenne">
                            <var_prop name="velocity_U"/>
                            <var_prop name="velocity_V"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" label="SECOND">
                            <var_prop name="velocity_W"/>
                            <var_prop name="velocity_V"/>
                            <time_step_start>12</time_step_start>
                    </time_average>
                    <time_average id="3" label="trois">
                            <var_prop name="pressure"/>
                            <var_prop name="velocity_W"/>
                            <time_step_start>33</time_step_start>
                    </time_average>
                 </time_averages>'''
        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not replace one average in TimeAveragesModel'


    def checkDeleteAverage(self):
        """Check whether the TimeAveragesModel class could be deleted one average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setAverage(2, 'deux', 20, 1, ['VelocitU'])
        mdl.setAverage(3, 'trois', 33, 1, ['VelocitW', 'VelocitU'])
        mdl.deleteAverage(2)
        doc = '''<time_averages>
                    <time_average id="1" label="moyenne">
                            <var_prop name="velocity_U"/>
                            <var_prop name="velocity_V"/>
                            <time_step_start>10</time_step_start>
                    </time_average>
                    <time_average id="2" label="trois">
                            <var_prop name="velocity_W"/>
                            <var_prop name="velocity_U"/>
                            <time_step_start>33</time_step_start>
                    </time_average>
                 </time_averages>'''

        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not delete one average in TimeAveragesModel'


    def checkDeleteAverage(self):
        """Check whether the TimeAveragesModel class could be get restart average"""
        mdl = TimeAveragesModel(self.case)
        mdl.setAverage(1, 'moyenne', 10, 1, ['VelocitU', 'VelocitV'])
        mdl.setAverage(2, 'deux', 20, 1, ['VelocitU'])
        mdl.setAverage(3, 'trois', 33, 1, ['VelocitW', 'VelocitU'])

        assert mdl.getAverageRestart(1) == None,\
            'Could not get restart average in TimeAveragesModel'

        from Pages.StartRestartModel import StartRestartModel
        StartRestartModel(self.case).setRestart('on')
        mdl.replaceAverage(2, 'deux', 20, 3, ['VelocitU'])

        assert mdl.getAverageRestart(2) == 3,\
            'Could not get restart average in TimeAveragesModel'


def suite():
    testSuite = unittest.makeSuite(TimeAveragesTestCase, "check")
    return testSuite


def runTest():
    print __file__
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------