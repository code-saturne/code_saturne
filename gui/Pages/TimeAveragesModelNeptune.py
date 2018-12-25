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

import sys, unittest
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Base.XMLengine import *
from code_saturne.Base.Common import LABEL_LENGTH_MAX


from code_saturne.Pages.StartRestartModel import StartRestartModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class TimeAveragesModel(Variables, Model):
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
        value['label']     = "TimeAverage"

        return value


    def __updateDicoLabel2Name(self):
        """
        Private method.
        Gets a dictionnary to connect name and label from
        variables and properties.
        """
        self.dicoLabel2Name = {}
        # get all variable nodes
        nodeList = self.getVariablesPropertiesList('no', 'no', 'no')

        for node in nodeList:
            name = node['name']
            label = node['label']
            fieldId = node['field_id']
            if not label:
                raise ValueError("Node has no label")

            dim = node['dimension']
            if int(dim) > 1:
                for ii in range(int(dim)):
                    label1 = label + "[" + str(ii) + "]"
                    if not (node['support'] and node['support'] == "boundary"):
                        self.dicoLabel2Name[label1] = (name, fieldId, str(ii))
            else:
                if not (node['support'] and node['support'] == "boundary"):
                    if name != 'local_time_step':
                        self.dicoLabel2Name[label] = (name, fieldId, str(0))

        return list(self.dicoLabel2Name.keys())


    @Variables.undoGlobal
    def addTimeAverage(self):
        """
        Public method.
        Add a new time average and return a default label
        """
        label = self.defaultValues()['label']
        def_label = label + str(len(self.getTimeAverageLabels()) + 1)

        # define default label
        if def_label in self.getTimeAverageLabels():
            i = 2
            while def_label in self.getTimeAverageLabels():
                def_label = label + str(len(self.getTimeAverageLabels()) + i)
                i = i + 1

        node = self.node_mean.xmlInitNode('time_average', label = def_label, name = def_label)
        node.xmlInitChildNode('listing_printing', status='on')
        node.xmlInitChildNode('postprocessing_recording', status='on')
        node['dimension'] = "1"
        node['field_id']  = "none"
        node['id'] = len(self.getTimeAverageLabels())
        ntdmom = self.defaultValues()['start']
        ttdmom = -1.0
        imoold = self.defaultValues()['restart']
        self.setTimeStepStart(def_label, ntdmom)
        self.setTimeStart(def_label, ttdmom)
        self.setRestart(def_label, imoold)

        return def_label, ntdmom, ttdmom, imoold


    @Variables.undoGlobal
    def deleteTimeAverage(self, label):
        """
        Public method.
        @type label: C{String}
        @param label: label of the time average to delete.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlGetNode('time_average', label=label)
        if node:
            nb = node['id']
            node.xmlRemoveNode()

            # renumerotation of all time average
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
        if StartRestartModel(self.case).getRestartPath():
            restart = node.xmlGetInt('restart_from_time_average')

        for var in node.xmlGetChildNodeList('var_prop'):
            for label in list(self.dicoLabel2Name.keys()):
                if self.dicoLabel2Name[label] == (var['name'], var['field_id'], var['component']):
                    lst.append(label)
        return node['label'], start, timestart, restart, lst


    @Variables.undoLocal
    def setLabel(self, old_label, label):
        """
        Public method.
        """
        self.isInList(old_label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=old_label)
        node['label'] = label
        node['name'] = label


    @Variables.undoLocal
    def setTimeStart(self, label, start):
        """
        Public method.
        """
        self.isFloat(start)
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        node.xmlSetData('time_start', start)


    @Variables.noUndo
    def getTimeStart(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        return node.xmlGetDouble('time_start')


    @Variables.undoLocal
    def setTimeStepStart(self, label, start):
        """
        Public method.
        """
        self.isInt(start)
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        node.xmlSetData('time_step_start', start)


    @Variables.noUndo
    def getTimeStepStart(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        return node.xmlGetInt('time_step_start')


    @Variables.undoLocal
    def setRestart(self, label, restart):
        """
        Public method.
        """
        self.isInt(restart)
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        if StartRestartModel(self.case).getRestartPath():
            if restart != -2:
                node.xmlSetData('restart_from_time_average', restart)
            else:
                node.xmlRemoveChild('restart_from_time_average')
        else:
            node.xmlRemoveChild('restart_from_time_average')


    @Variables.noUndo
    def getRestart(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        return node.xmlGetInt('restart_from_time_average')


    @Variables.undoLocal
    def setVariable(self, label, lst):
        """
        Public method.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)
        node.xmlRemoveChild('var_prop')
        for var in lst:
            self.isInList(var, self.__var_prop_list)
            (name, field, comp) = self.dicoLabel2Name[var]
            node.xmlAddChild('var_prop', name=name, field_id=field, component=comp)


    @Variables.noUndo
    def getVariable(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getTimeAverageLabels())
        node = self.node_mean.xmlInitNode('time_average', label=label)

        lst = []
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['field_id'], var['component']) :
                    lst.append(name)
        return lst


    @Variables.noUndo
    def getTimeAverageLabels(self):
        """
        Public method.
        @return: list of time averages labels.
        @rtype: C{List} of C{String}
        """
        lst = []
        for node in self.node_mean.xmlGetNodeList('time_average'):
            label = node['label']
            lst.append(label)
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
        """Check whether the TimeAveragesModel class could be set an average"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        mdl.setTimeAverage('deux', 20, 1, ['Pressure'])
        doc = '''<time_averages>
                         <time_average id="1" label="moyenne">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <time_step_start>
                                         10
                                 </time_step_start>
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </time_average>
                         <time_average id="2" label="deux">
                                 <var_prop field_id="none" name="Pressure"/>
                                 <time_step_start>
                                         20
                                 </time_step_start>
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </time_average>
                 </time_averages>'''
        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not set an average'


    def checkreplaceTimeAverage(self):
        """Check whether the TimeAveragesModel class could be replaced an average"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        mdl.setTimeAverage('deux', 20, 1, ['Pressure'])
        mdl.setTimeAverage('trois', 33, 1, ['Pressure', 'U1'])
        mdl.replaceTimeAverage('deux','SECOND', 12, 1, ['W1', 'V1'])
        mdl.replaceTimeAverage('trois','DREI', 33, 1, ['Pressure', 'W1'])
        doc = '''<time_averages>
                         <time_average id="1" label="moyenne">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <time_step_start>
                                         10
                                 </time_step_start>
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </time_average>
                         <time_average id="2" label="SECOND">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <var_prop field_id="1" name="VelocityZ"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <time_step_start>
                                         12
                                 </time_step_start>
                         </time_average>
                         <time_average id="3" label="DREI">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <var_prop field_id="none" name="Pressure"/>
                                 <var_prop field_id="1" name="VelocityZ"/>
                                 <time_step_start>
                                         33
                                 </time_step_start>
                         </time_average>
                 </time_averages>'''
        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not replace one average in TimeAveragesModel'


    def checkdeleteTimeAverage(self):
        """Check whether the TimeAveragesModel class could be deleted time average"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        mdl.setTimeAverage('deux', 20, 1, ['U1'])
        mdl.setTimeAverage('trois', 33, 1, ['W1', 'U1'])
        mdl.deleteTimeAverage('deux')
        doc = '''<time_averages>
                         <time_average id="1" label="moyenne">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <time_step_start>
                                         10
                                 </time_step_start>
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </time_average>
                         <time_average id="2" label="trois">
                                 <var_prop field_id="1" name="VelocityZ"/>
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <time_step_start>
                                         33
                                 </time_step_start>
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </time_average>
                 </time_averages>'''

        assert mdl.node_mean == self.xmlNodeFromString(doc),\
            'Could not delete one average in TimeAveragesModel'


    def checkGetTimeAverageData(self):
        """Check whether the TimeAveragesModel class could get time average data"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        assert mdl.getTimeAverageData(1) == ('moyenne', 10, -2, ['U1', 'V1']),\
            'Could not get TimeAverageData'


    def checkGetTimeAverageLabels(self):
        """Check whether the TimeAveragesModel class could get time average Labels"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        mdl.setTimeAverage('deux', 20, 1, ['U1'])
        assert mdl.getTimeAverageLabels() == ['moyenne', 'deux'],\
            'Could not get TimeAverageLabels'


    def checkGetNumberOfTimeAverage(self):
        """Check whether the TimeAveragesModel class could get NumberOfTimeAverage"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = TimeAveragesModel(self.case)
        mdl.setTimeAverage('moyenne', 10, 1, ['U1', 'V1'])
        mdl.setTimeAverage('deux', 20, 1, ['U1'])
        assert mdl.getNumberOfTimeAverage() == 2,\
            'Could not get NumberOfTimeAverage'


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
