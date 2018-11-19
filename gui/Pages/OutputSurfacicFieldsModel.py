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
from code_saturne.Base.XMLmodelNeptune import *


class OutputSurfacicFieldsModel(Model):
    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor
        """
        self.case = case
        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeproperty = self.node_models.xmlInitNode('properties')
        self.listNodeSurface = self._getListOfSurfacicProperties()

        self.dicoLabelName = {} # AZ test
        self.list_name = []     #AZ
        self.list_label = []
        self._updateDicoLabelName()


    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['status']    = "on"

        return default


    def _getListOfSurfacicProperties(self):
        """
        Private method: return node of yplus, effort
        """
        nodeList = []
        for node in self.XMLNodeproperty.xmlGetNodeList('property'):
            if node['support']  and node['support'] == "boundary":
                nodeList.append(node)
        return nodeList


    def _updateDicoLabelName(self):
        """
        Private method: update dictionaries of labels for all properties .....
        """
        for node in self.listNodeSurface:
            name = node['name']
            if not name:
                name = node['label']
            if not node['label']:
                msg = "xml node named "+ name +" has no label"
                raise ValueError(msg)
            self.list_label.append(node['label'])

            field_id = node['field_id']
            if field_id != "none":
                name += str(field_id)
            self.dicoLabelName[name] = node['label']
            self.list_name.append(name)


    def getLabelsList(self):
        """
        Return list of labels for all properties .....Only for the View
        """
        lst = []
        for node in self.listNodeSurface:
            lst.append(node['label'])
        return lst


    @Variables.undoLocal
    def setPropertyLabel(self, old_label, new_label):
        """
        Replace old_label by new_label for node with name and old_label. Only for the View
        """
        self.isInList(old_label, self.getLabelsList())
        if old_label != new_label:
            self.isNotInList(new_label, self.getLabelsList())
        for node in self.listNodeSurface:
            if node['label'] == old_label:
                node['label'] = new_label
        self._updateDicoLabelName()


    @Variables.noUndo
    def getPostProcessing(self, label):
        """ Return status of post processing for node withn label 'label'"""
        self.isInList(label, self.getLabelsList())
        status = self._defaultValues()['status']
        for node in self.listNodeSurface:
            if node['label'] == label:
                node_post = node.xmlGetChildNode('postprocessing_recording', 'status')
                if node_post:
                    status = node_post['status']
        return status


    @Variables.undoLocal
    def setPostProcessing(self, label, status):
        """ Put status of post processing for node with label 'label'"""
        self.isOnOff(status)
        self.isInList(label, self.getLabelsList())
        for node in self.listNodeSurface:
            if node['label'] == label:
                node.xmlInitChildNode('postprocessing_recording')['status'] = status

#-------------------------------------------------------------------------------
# OutputVariableModel test case
#-------------------------------------------------------------------------------

class OutputSurfacicVariablesTestCase(ModelTest):
    """
    """
    def checkOutputSurfacicVariablesInstantiation(self):
        """
        Check whether the OutputSurfacicVariables Model class could be instantiated
        """
        model = None
        model = OutputSurfacicVariablesModel(self.case)
        assert model != None, 'Could not instantiate OutputSurfacicVariablesModel'

    def checkSetPropertyLabel(self):
        """
        Check whether the OutputSurfacicVariablesModel class could be set a label
        of property
        """
        model = OutputSurfacicVariablesModel(self.case)
        model.setPropertyLabel('Yplus', 'toto')
        node = model.node_models.xmlInitNode('velocity_pressure')
        doc = '''<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="toto" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                 </velocity_pressure>'''
        assert node == self.xmlNodeFromString(doc),\
            'Could not set label of property in output surfacic variables model'

    def checkSetandGetPostProcessing(self):
        """
        Check whether the OutputSurfacicVariablesModel class could be set and
        get status for post processing of named property
        """
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')
        del ThermalRadiationModel
        model = OutputSurfacicVariablesModel(self.case)

        model.setPostProcessing('Flux_convectif','off')
        doc = '''<property label="Flux_convectif" name="flux_convectif" support="boundary">
                    <postprocessing_recording status="off"/>
                 </property>'''

        assert model.getPostProcessing('Flux_convectif') == 'off',\
                'Could not get status for post processing of named property'


def suite():
    testSuite = unittest.makeSuite(OutputSurfacicVariablesTestCase, "check")
    return testSuite


def runTest():
    print("OutputSurfacicVariablesTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


