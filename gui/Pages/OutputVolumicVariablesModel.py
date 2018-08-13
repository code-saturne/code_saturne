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
"""

#-------------------------------------------------------------------------------
# Library odules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesModel(Variables, Model):

    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case
        self.node_models    = self.case.xmlInitNode('thermophysical_models')
        self.analysis_ctrl  = self.case.xmlInitNode('analysis_control')
        self.fluid_prop     = self.case.xmlInitNode('physical_properties')
        self.node_model_vp  = self.node_models.xmlInitNode('velocity_pressure')
        self.node_ale       = self.node_models.xmlGetChildNode('ale_method')
        self.node_output    = self.analysis_ctrl.xmlInitNode('output')
        self.node_probe     = self.node_output.xmlGetNodeList('probe','name')
        self.node_means     = self.analysis_ctrl.xmlInitNode('time_averages')
        self.node_error     = self.analysis_ctrl.xmlInitNode('error_estimator')

        model = XMLmodel(self.case)

        self.updateList()


# Following private methods: (to see for gathering eventually)

    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['status']    = "on"
        default['estimator'] = "0"

        return default


    def updateList(self):
        """
        """
        model = XMLmodel(self.case)

        self.listNodeVolum = (self._getListOfVelocityPressureVariables(),
                              model.getTurbNodeList(),
                              self.getThermalScalar(),
                              self.getAdditionalScalar(),
                              self.getAdditionalScalarProperty(),
                              self.getFluidProperty(),
                              self.getTimeProperty(),
                              self.getMeteoScalProper(),
                              self.getElecScalProper(),
                              self.getPuCoalScalProper(),
                              self.getGasCombScalProper(),
                              self._getWeightMatrixProperty(),
                              self.getListOfTimeAverage(),
                              self._getListOfAleMethod(),
                              self._getThermalRadiativeProperties(),
                              self._getListOfEstimator())

        self.listNode = []
        for part in self._getListOfVelocityPressureVariables():
            self.listNode.append([part, 'base'])
        for part in model.getTurbNodeList():
            self.listNode.append([part, 'turbulence'])
        for part in self.getThermalScalar():
            self.listNode.append([part, 'thermal'])
        for part in self._getThermalRadiativeProperties():
            self.listNode.append([part, 'thermal'])
        for part in self.getPuCoalScalProper():
            self.listNode.append([part, 'coal'])
        for part in self.getGasCombScalProper():
            self.listNode.append([part, 'gas'])
        for part in self.getMeteoScalProper():
            self.listNode.append([part, 'atmospheric'])
        for part in self.getElecScalProper():
            self.listNode.append([part, 'electric'])
        for part in self.getAdditionalScalar():
            self.listNode.append([part, 'other'])
        for part in self.getAdditionalScalarProperty():
            self.listNode.append([part, 'other'])
        for part in self.getFluidProperty():
            self.listNode.append([part, 'physical_properties'])
        for part in self.getTimeProperty():
            self.listNode.append([part, 'other'])
        for part in self.getListOfTimeAverage():
            self.listNode.append([part, 'other'])
        for part in self._getListOfAleMethod():
            self.listNode.append([part, 'other'])
        for part in self._getListOfEstimator():
            self.listNode.append([part, 'estimator'])

        self.dicoLabelName = {}
        self.list_name = []
        self._updateDictLabelName()


    def _updateDictLabelName(self):
        """
        Update dictionaries of labels for all variables, properties .....
        """
        for nodeList in self.listNode:
            node = nodeList[0]
            tpe  = nodeList[1]

            name = node['name']
            if not name:
                name = node['label']
            if not node['label']:
                msg = "xml node named "+ name +" has no label"
                raise ValueError(msg)
            self.dicoLabelName[name] = node['label']
            self.list_name.append([name, tpe])


    def _getListOfVelocityPressureVariables(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        for tag in ('variable', 'property'):
            for node in self.node_model_vp.xmlGetNodeList(tag):
                if not node['support']:
                    nodeList.append(node)
        return nodeList


    def _getListOfEstimator(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        for node in self.node_error.xmlGetNodeList('property'):
            if not node['support']:
                nodeList.append(node)
        return nodeList


    def _getWeightMatrixProperty(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        node0 = self.case.xmlGetNode('numerical_parameters')
        node1 = node0.xmlGetNode('velocity_pressure_coupling', 'status')
        if node1:
            if node1['status'] == 'on':
                nodeList = node0.xmlGetNodeList('property')
        return nodeList


    def _getListOfAleMethod(self):
        """
        Private method: return list of variables and properties for ale method if it's activated
        """
        nodeList = []
        if self.node_ale['status'] == 'on':
            for tag in ('variable', 'property'):
                for node in self.node_ale.xmlGetChildNodeList(tag):
                    nodeList.append(node)

        return nodeList


    def _getThermalRadiativeProperties(self):
        """
        Private method: return list of volumic properties for thermal radiation
        """
        nodeList = []
        if ThermalRadiationModel(self.case).getRadiativeModel() != "off":
            self.node_ray = self.node_models.xmlGetNode('radiative_transfer')
            for node in self.node_ray.xmlGetChildNodeList('property'):
                if not node['support']:
                    nodeList.append(node)
        return nodeList


# Following methods also called by ProfilesModel and TimeAveragesModel

    @Variables.noUndo
    def getThermalScalar(self):
        """
        Return node of thermal scalar (idem ds NumericalParamEquationModel)
        """
        node_models = self.case.xmlGetNode('thermophysical_models')
        node = node_models.xmlGetNode('thermal_scalar')
        return node.xmlGetNodeList('variable', type='thermal')


    @Variables.noUndo
    def getPuCoalScalProper(self):
        """
        Return list fo nodes of pulverized coal.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodList = []
        node = self.node_models.xmlGetNode('solid_fuels', 'model')
        model = node['model']
        varList = []
        if model != 'off':
            for var in ('variable', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getGasCombScalProper(self):
        """
        Return list of nodes of gas combustion.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodList = []
        node = self.node_models.xmlGetNode('gas_combustion', 'model')
        model = node['model']
        varList = []
        if model != 'off':
            for var in ('variable', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getMeteoScalProper(self):
        """
        Return list fo nodes of atmospheric flows.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodList = []
        node = self.node_models.xmlGetNode('atmospheric_flows', 'model')
        if not node: return []
        model = node['model']
        varList = []
        if model != 'off':
            for var in ('variable', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getElecScalProper(self):
        """
        Return list fo nodes of electric flows.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodList = []
        node = self.node_models.xmlGetNode('joule_effect', 'model')
        if not node: return []
        model = node['model']
        varList = []
        if model != 'off':
            for var in ('variable', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getModelVariables(self, model_name):
        """
        Return list of variable nodes for a given model.
        """
        nodList = []
        node = self.node_models.xmlGetNode(model_name, 'model')
        if not node: return []
        model = node['model']
        varList = []
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
            for nodvar in nodList:
                varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getAdditionalScalar(self):
        """
        Return list of nodes of user scalars
        Also called by ProfilesModel and TimeAveragesModel
        (idem ds NumericalParamEquationModel named getAdditionalScalarNodes)
        """
        node = self.case.xmlGetNode('additional_scalars')
        return node.xmlGetNodeList('variable', type='user')


    @Variables.noUndo
    def getAdditionalScalarProperty(self):
        """
        Return list of nodes of properties of user scalars
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodeList = []
        for node in self.getAdditionalScalar():
            L = node.xmlGetNode('property', choice='variable')
            if L:
                nodeList.append(L)
        return nodeList


    @Variables.noUndo
    def getFluidProperty(self):
        """
        Return list of nodes of fluid properties
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodeList = []
        model = self.getThermalScalar()

        node = self.fluid_prop.xmlGetNode('fluid_properties')
        if node:
            for prop in ('density',
                         'molecular_viscosity',
                         'specific_heat',
                         'thermal_conductivity'):
                L = node.xmlGetNode('property', name=prop, choice='variable')
                if L:
                    nodeList.append(L)

        return nodeList


    @Variables.noUndo
    def getTimeProperty(self):
        """
        Return list fo nodes of properties of time_parameters.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodeList = []

        node1 = self.analysis_ctrl.xmlGetNode('time_parameters')

        if node1:
            if node1.xmlGetInt('time_passing'):
                node2 = node1.xmlGetNode('property', name='local_time_step')
                if node2:
                    nodeList.append(node2)

            for prop in ('courant_number', 'fourier_number'):
                L = node1.xmlGetNode('property', name=prop)
                if L: nodeList.append(L)

        return nodeList


    @Variables.noUndo
    def getListOfTimeAverage(self):
        """
        Return list of time averages variables
        Also called by ProfilesModel
        """
        nodeList = []
        for node in self.node_means.xmlGetNodeList('time_average'):
            nodeList.append(node)

        return nodeList


#Following methods only called by the View
    @Variables.noUndo
    def getLabelsList(self):
        """
        Return list of labels for all variables, properties .....Only for the View
        """
        lst = []
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                lst.append(node['label'])
        return lst


    @Variables.noUndo
    def getNamesList(self):
        """
        Return list of names for all variables, properties .....Only for the View
        """
        lst = []
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                lst.append(node['name'])
        return lst


    @Variables.noUndo
    def getProbeList(self):
        """ Return list of node for probes """
        probeList = []
        for node in self.node_probe:
            probeList.append(node['name'])
        return probeList


    @Variables.noUndo
    def getPrintingStatus(self, name):
        """
        Return status of markup printing from node with name. Only for the View
        """
        self.isInList(name, self.getNamesList())
        status = self._defaultValues()['status']
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    node_printing = node.xmlGetChildNode('listing_printing', 'status')
                    if node_printing:
                        status = node_printing['status']
        return status


    @Variables.noUndo
    def getPostStatus(self, name):
        """
        Return status of markup  post processing from node with name. Only for the View
        """
        self.isInList(name, self.getNamesList())
        status = self._defaultValues()['status']
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    node_post = node.xmlGetChildNode('postprocessing_recording', 'status')
                    if node_post:
                        status = node_post['status']
        return status


    @Variables.noUndo
    def getMonitorStatus(self, name):
        """
        Return status of markup monitoring from node with name. Only for the View
        """
        self.isInList(name, self.getNamesList())
        status = self._defaultValues()['status']
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    node_post = node.xmlGetChildNode('probes_recording', 'status')
                    if node_post:
                        status = node_post['status']
        return status


    @Variables.undoLocal
    def setVariableLabel(self, old_label, new_label):
        """
        Replace old_label by new_label for node with name and old_label. Only for the View
        """
        # fusion de cette methode avec DefineUserScalarsModel.renameScalarLabel
        self.isInList(old_label, self.getLabelsList())
        self.isNotInList(new_label, [""])

        if old_label != new_label:
            self.isNotInList(new_label, self.getLabelsList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == old_label:
                    node['label'] = new_label

        self._updateDictLabelName()
        self._updateBoundariesNodes(old_label, new_label)

        for node in self.case.xmlGetNodeList('formula'):
            f = node.xmlGetTextNode()
            if f:
                f.replace(old_label, new_label)
                node.xmlSetTextNode(f)


    def _updateBoundariesNodes(self, old_label, new_label):
        """
        Update good label for boundaries nodes with name and label. Only for the View
        """
        self.node_bc  = self.case.xmlInitNode('boundary_conditions')
        self.node_var = self.node_bc.xmlInitNodeList('variable')

        for node in self.node_var:
            if node['label'] == old_label:
                node['label'] = new_label


    @Variables.undoLocal
    def setPrintingStatus(self, name, status):
        """
        Put status for printing from node with name and label
        """
        self.isOnOff(status)
        self.isInList(name, self.getNamesList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    if status == 'off':
                        node.xmlInitChildNode('listing_printing')['status'] = status
                    else:
                        if node.xmlGetChildNode('listing_printing'):
                            node.xmlRemoveChild('listing_printing')


    @Variables.noUndo
    def getVariableLabel(self, name) :
        """
        return label of name variable
        """
        for variableType in ('variable', 'property') :
            node = self.case.xmlGetNode(variableType, name = name)
            if node != None:
                break

        if node != None:
            label = node['label']
            return label
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    @Variables.undoLocal
    def setPostStatus(self, name, status):
        """
        Put status for postprocessing from node with name and label
        """
        self.isOnOff(status)
        self.isInList(name, self.getNamesList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    if status == 'off':
                        node.xmlInitChildNode('postprocessing_recording')['status'] = status
                    else:
                        if node.xmlGetChildNode('postprocessing_recording'):
                            node.xmlRemoveChild('postprocessing_recording')


    @Variables.undoLocal
    def setMonitorStatus(self, name, status):
        """
        Put status for monitoring from node with name and label
        """
        self.isOnOff(status)
        self.isInList(name, self.getNamesList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['name'] == name:
                    if status == 'off':
                        node.xmlInitChildNode('probes_recording')['status'] = status
                    else:
                        if node.xmlGetChildNode('probes_recording'):
                            node.xmlRemoveChild('probes_recording')


    @Variables.noUndo
    def getEstimatorModel(self, name):
        """
        Return model for an error estimator
        """
        self.isInList(name, ["Correction", "Drift", "Prediction", "Total"])
        status = self._defaultValues()['estimator']

        nn = self.node_error.xmlGetChildNode(name, 'model')
        if nn:
            status = nn['model']

        return status


    @Variables.undoLocal
    def setEstimatorModel(self, name, model):
        """
        Put model for an error estimator
        """
        self.isInList(model, ['0', '1', '2'])
        self.isInList(name, ["Correction", "Drift", "Prediction", "Total"])
        status = self._defaultValues()['estimator']

        if model != status:
            if self.node_error.xmlGetChildNode(name):
                self.node_error.xmlRemoveChild(name)
            self.node_error.xmlInitChildNode(name)['model'] = model
        else:
            if self.node_error.xmlGetChildNode(name):
                self.node_error.xmlRemoveChild(name)

        # add and remove field associated
        if model != status:
            nn = self.node_error.xmlGetNode(name)
            if name == "Correction":
                self.setNewProperty(nn, "est_error_cor_" + model)
            elif name == "Drift":
                self.setNewProperty(nn, "est_error_der_" + model)
            elif name == "Prediction":
                self.setNewProperty(nn, "est_error_pre_" + model)
            elif name == "Total":
                self.setNewProperty(nn, "est_error_tot_" + model)

        self.updateList()


#-------------------------------------------------------------------------------
# OutputVolumicVariablesModel Test Class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesModelTestCase(ModelTest):
    """
    Unittest
    """
    def checkOutputVolumicVariablesModelInstantiation(self):
        """Check whether the OutputVolumicVariablesModel class could be instantiated"""
        mdl = None
        mdl = OutputVolumicVariablesModel(self.case)
        assert mdl != None, 'Could not instantiate OutputVolumicVariablesModel'


    def checkSetVariableLabel(self):
        """
        Check whether the OutputVolumicVariablesModel class could be set a label
        of property
        """
        model = OutputVolumicVariablesModel(self.case)
        model.setVariableLabel('VelocitV', 'vitV')
        node = model.node_models.xmlInitNode('velocity_pressure')
        doc = '''<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="vitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>'''
        assert node == self.xmlNodeFromString(doc),\
            'Could not set label of property in output volumic variables model'

    def checkSetAndGetPrintingStatus(self):
        """
        Check whether the OutputVolumicVariablesModel class could be
        set and get status for printing listing
        """
        from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

        mdl = OutputVolumicVariablesModel(self.case)
        mdl.setPrintingStatus('TempC', 'off')
        node_out = mdl.case.xmlGetNode('additional_scalars')
        doc = '''<additional_scalars>
                    <variable label="TempC" name="temperature_celsius" type="thermal">
                        <initial_value zone_id="1">20.0</initial_value>
                        <min_value>-1e+12</min_value>
                        <max_value>1e+12</max_value>
                        <listing_printing status="off"/>
                    </variable>
                 </additional_scalars>'''

        assert node_out == self.xmlNodeFromString(doc),\
            'Could not set status of listing printing in output volumic variables model'
        assert mdl.getPrintingStatus('TempC') == 'off',\
            'Could not get status of listing printing in output volumic variables model'

    def checkSetAndGetPostStatus(self):
        """
        Check whether the OutputVolumicVariablesModel class could be
        set and get status for printing
        """
        from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

        mdl = OutputVolumicVariablesModel(self.case)
        mdl.setPostStatus('TempC', 'off')
        node_out = mdl.case.xmlGetNode('additional_scalars')
        doc = '''<additional_scalars>
                    <variable label="TempC" name="temperature_celsius" type="thermal">
                        <initial_value zone_id="1">20.0</initial_value>
                        <min_value>-1e+12</min_value>
                        <max_value>1e+12</max_value>
                        <postprocessing_recording status="off"/>
                    </variable>
                 </additional_scalars>'''

        assert node_out == self.xmlNodeFromString(doc),\
            'Could not set status of post processing in output volumic variables model'
        assert mdl.getPostStatus('TempC') == 'off',\
            'Could not get status of post processing in output volumic variables model'

def suite():
    testSuite = unittest.makeSuite(OutputVolumicVariablesModelTestCase, "check")
    return testSuite

def runTest():
    print("OutputVolumicVariablesModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
