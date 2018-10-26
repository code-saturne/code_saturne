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
Output of volume variables for NEPTUNE_CFD module
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool

from code_saturne.Base.XMLmodelNeptune import XMLmodel, ModelTest
from code_saturne.Base.XMLvariablesNeptune import Variables
from code_saturne.Base.XMLvariables import Model

from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel

from code_saturne.Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
from code_saturne.Pages.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesModelNeptune(OutputVolumicVariablesModel, MainFieldsModel, Variables, Model):

    def __init__(self, case):
        """
        Constuctor.
        """
        MainFieldsModel.__init__(self, case)
        self.case = case
        self.node_models    = self.case.root().xmlInitNode('closure_modeling')
        self.analysis_ctrl  = self.case.root().xmlInitNode('analysis_control')
        self.fluid_prop = self.case.root().xmlInitNode('physical_properties')
        self.node_model_vp  = self.node_models.xmlInitNode('velocity_pressure')
        self.node_ale       = self.node_models.xmlInitNode('ale_method')

        self.node_output    = self.analysis_ctrl.xmlInitNode('output')
        self.node_probe     = self.node_output.xmlGetNodeList('probe','name')
        self.node_means     = self.analysis_ctrl.xmlInitNode('time_averages')
        self.node_error     = self.analysis_ctrl.xmlInitNode('error_estimator')

        model = XMLmodel(self.case)

        self.updateList()

    def defaultValues(self):
        default = {}
        default['listing'] = 'on'
        default['writer']  = 'on'
        return default

    def updateList(self):
        """
        """

        model = XMLmodel(self.case)

        self.listNode = []

        fieldId = 'none'
        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            for node in self.case.xmlGetNodeList(variableType, field_id = "none"):
                if not node['name'].startswith("User_"):
                    self.listNode.append([node, fieldId])


        #On recupere les fields :
        fd = []
        fd.append('none')
        thermo = self.case.xmlGetNode('thermophysical_models')
        fields = thermo.xmlGetNode('fields')
        for node in fields.xmlInitChildNodeList('field'):
            field_id = node.xmlGetAttribute('field_id')
            fd.append(field_id)

        for variableType in ('variable', 'property', 'scalar'):
            for field in fd :
                for node in self.case.xmlGetNodeList(variableType, field_id = field):
                    self.listNode.append([node, field])

        self.dicoLabelName = {}
        self.list_name = []
        self._updateDictLabelName()


    #This method was written based on neptune_cfd.gui.OutputFieldsModel.OutputFieldsModel.getListingStatus method
    @Variables.noUndo
    def getPrintingStatus(self, name, fieldId):
        """
        return status for listing output for variable name on fieldId
        """
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            value = self.defaultValues()['listing']
            n = node.xmlGetNode('listing_printing')
            if n :
                value = n['status']
            return value
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    #This method was written based on neptune_cfd.gui.OutputFieldsModel.OutputFieldsModel.getPostProcessingStatus method
    @Variables.noUndo
    def getPostStatus(self, name, fieldId):
        """
        return status for post processing for variable name on fieldId
        """

        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            value = self.defaultValues()['writer']
            n = node.xmlGetNode('postprocessing_recording')
            if n :
                value = n['status']
            return value
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    @Variables.noUndo
    def getMonitorStatus(self, name, fieldId):
        """
        Return status of markup monitoring from node with name. Only for the View
        """

        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)
        status = self._defaultValues()['status']

        for variableType in ('variable', 'property', 'scalar', 'time_average'):
            for parent in self.case.xmlGetNodeList(variableType, field_id = str(fieldId)):
                if parent['name'] == name :
                    node_post = parent.xmlGetChildNode('probes_recording', 'status')
                    if node_post :
                        status = node_post['status']
        return status


    #This method was written based on neptune_cfd.gui.OutputFieldsModel.OutputFieldsModel.setListingStatus method
    @Variables.undoLocal
    def setPrintingStatus(self, name, status, fieldId):
        """
        return status for listing output for variable name on fieldId
        """

        self.isOnOff(status)
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            n = node.xmlInitNode('listing_printing')
            n['status'] = status
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    #This method was written based on neptune_cfd.gui.OutputFieldsModel.OutputFieldsModel.setPostProcessingStatus method
    @Variables.undoLocal
    def setPostStatus(self, name, status, fieldId):
        """
        return status for post processing for variable name on fieldId
        """

        self.isOnOff(status)
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            n = node.xmlInitNode('postprocessing_recording')
            n['status'] = status
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    #This method was written based on neptune_cfd.gui.OutputFieldsModel.OutputFieldsModel.setPostProcessingStatus method
    @Variables.undoLocal
    def setMonitorStatus(self, name, status, fieldId):
        """
        return status for post monitoring for variable name on fieldId
        """

        self.isOnOff(status)
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            if status == 'off':
                node.xmlInitChildNode('probes_recording')['status'] = status
            else:
                if node.xmlGetChildNode('probes_recording'):
                    node.xmlRemoveChild('probes_recording')
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    def getVariables(self):

        #On recupere les fields :
        fd = []
        fd.append('none')
        thermo = self.case.xmlGetNode('thermophysical_models')
        fields = thermo.xmlGetNode('fields')
        for node in fields.xmlInitChildNodeList('field'):
            field_id = node.xmlGetAttribute('field_id')
            fd.append(field_id)

        l = []
        for variableType in ('variable', 'property', 'scalar'):
            for field in fd :
                for node in self.case.xmlGetNodeList(variableType, field_id = field):
                    l.append(node)

        return l

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
