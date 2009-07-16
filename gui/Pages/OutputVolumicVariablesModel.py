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
"""

#-------------------------------------------------------------------------------
# Library odules import
#-------------------------------------------------------------------------------

import string, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import XMLmodel, ModelTest
from Base.XMLvariables import Model
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.ThermalRadiationModel import ThermalRadiationModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesModel(Model):

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

        model = XMLmodel(self.case)

        self.listNodeVolum = (self._getListOfVelocityPressureVariables(),
                              model.getTurbNodeList(),
                              self.getThermalScalar(),
                              self.getAdditionalScalar(),
                              self.getAdditionalScalarProperty(),
                              self.getFluidProperty(),
                              self.getTimeProperty(),
                              self.getMeteoScalProper(),
                              self.getPuCoalScalProper(),
                              self._getWeightMatrixProperty(),
                              self.getListOfTimeMeans(),
                              self._getListOfAleMethod(),
                              self._getThermalRadiativeProperties())

        self.dicoLabelName = {}
        self.list_name = []
        self._updateDicoLabelName()


# Following private methods: (to see for gathering eventually)

    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['status']    = "on"

        return default


    def _updateDicoLabelName(self):
        """
        Update dictionaries of labels for all variables, properties .....
        """
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                name = node['name']
                if not name: name = node['label']
                if not node['label']:
                    msg = "xml node named "+ name +" has no label"
                    raise ValueError, msg
                self.dicoLabelName[name] = node['label']
                self.list_name.append(name)


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


    def getThermalScalar(self):
        """
        Return node of thermal scalar (idem ds NumericalParamEquationModel)
        """
        node = self.case.xmlGetNode('additional_scalars')
        return node.xmlGetNodeList('scalar', type='thermal')



    def getPuCoalScalProper(self):
        """
        Return list fo nodes of pulverized coal.
        Also called by ProfilesModel and TimeAveragesModel
        """
        nodList = []
        node = self.node_models.xmlGetNode('pulverized_coal', 'model')
        model = node['model']
        varList = []
        if model != 'off':
            for var in ('scalar', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList


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
            for var in ('scalar', 'property'):
                nodList = node.xmlGetNodeList(var)
                for nodvar in nodList:
                    varList.append(nodvar)
        return varList



    def getAdditionalScalar(self):
        """
        Return list of nodes of user scalars
        Also called by ProfilesModel and TimeAveragesModel
        (idem ds NumericalParamEquationModel named getAdditionalScalarNodes)
        """
        node = self.case.xmlGetNode('additional_scalars')
        return node.xmlGetNodeList('scalar', type='user')


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


    def getListOfTimeMeans(self):
        """
        Return list of time averages variables
        Also called by ProfilesModel
        """
        nodeList = []
        for node in self.node_means.xmlGetNodeList('time_average'):
            nodeList.append(node)

        return nodeList


#Following methods only called by the View
    def getLabelsList(self):
        """
        Return list of labels for all variables, properties .....Only for the View
        """
        list = []
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                list.append(node['label'])
        return list


    def getVariableProbeList(self):
        """ Return list of node for probes """
        probeList = []
        for node in self.node_probe:
            probeList.append(node['name'])
        return probeList


    def getProbesList(self, label):
        """
        Return list of probes if it exists for node['name'] = name. Only for the View
        """
        self.isInList(label, self.getLabelsList())
        list = self.getVariableProbeList()
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == label:
                    node_probes = node.xmlGetChildNode('probes')
                    if node_probes:
                        nb_probes = node_probes['choice']
                        if nb_probes == '0':
                            list = []
                        elif nb_probes > '0':
                            list = []
                            for n in node_probes.xmlGetChildNodeList('probe_recording'):
                                list.append(n['name'])
        return list


    def getPrintingStatus(self, label):
        """
        Return status of markup printing from node with label. Only for the View
        """
        self.isInList(label, self.getLabelsList())
        status = self._defaultValues()['status']
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == label:
                    node_printing = node.xmlGetChildNode('listing_printing', 'status')
                    if node_printing:
                        status = node_printing['status']
        return status


    def getPostStatus(self, label):
        """
        Return status of markup  post processing from node with label. Only for the View
        """
        self.isInList(label, self.getLabelsList())
        status = self._defaultValues()['status']
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == label:
                    node_post = node.xmlGetChildNode('postprocessing_recording', 'status')
                    if node_post:
                        status = node_post['status']
        return status


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

        self._updateDicoLabelName()
        self._updateBoundariesNodes(old_label, new_label)

        for node in self.case.xmlGetNodeList('formula'):
            f = node.xmlGetTextNode().replace(old_label, new_label)
            node.xmlSetTextNode(f)


    def _updateBoundariesNodes(self, old_label, new_label):
        """
        Update good label for boundaries nodes with name and label. Only for the View
        """
        self.node_bc  = self.case.xmlInitNode('boundary_conditions')
        self.node_var = self.node_bc.xmlInitNodeList('variable')
        self.node_sca = self.node_bc.xmlInitNodeList('scalar')

        for node in [self.node_var, self.node_sca]:
            for nodebc in node:
                if nodebc['label'] == old_label:
                    nodebc['label'] = new_label


    def setPrintingStatus(self, label, status):
        """
        Put status for balise printing from node with name and label
        """
        self.isOnOff(status)
        self.isInList(label, self.getLabelsList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == label:
                    if status == 'off':
                        node.xmlInitChildNode('listing_printing')['status'] = status
                    else:
                        if node.xmlGetChildNode('listing_printing'):
                            node.xmlRemoveChild('listing_printing')


    def setPostStatus(self, label, status):
        """
        Put status for balise postprocessing from node with name and label
        """
        self.isOnOff(status)
        self.isInList(label, self.getLabelsList())
        for nodeList in self.listNodeVolum:
            for node in nodeList:
                if node['label'] == label:
                    if status == 'off':
                        node.xmlInitChildNode('postprocessing_recording')['status'] = status
                    else:
                        if node.xmlGetChildNode('postprocessing_recording'):
                            node.xmlRemoveChild('postprocessing_recording')


    def updateProbes(self, label, list):
        """
        Update probe_recording markups if it exists
        """
        self.isInList(label, self.getLabelsList())
        nb = len(string.split(list))
        if nb == len(self.getVariableProbeList()):
            return
        else:
            for nodeList in self.listNodeVolum:
                for node in nodeList:
                    if node['label'] == label:
                        try:
                            node.xmlRemoveChild('probes')
                        except:
                            pass
                        n = node.xmlInitNode('probes', choice=str(nb))
                        if nb > 0:
                            for i in string.split(list):
                                n.xmlInitChildNodeList('probe_recording',name=i)

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

##    def checkGetPuCoalScalProper(self):
##        """
##        Check whether the OutputVolumicVariablesModel class could be get
##        properties of pulverized coal
##        """
##        mdl = OutputVolumicVariablesModel(self.case)
##        from CoalCombustionModel import CoalCombustionModel
##        CoalCombustionModel(self.case).setCoalCombustionModel('coal_homo')
##        del CoalCombustionModel
##        mdl.getPuCoalScalProper()
##        node = mdl.node_models.xmlGetNode('pulverized_coal', 'model')
##        print node
##        print mdl.getPuCoalScalProper()

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
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>'''
        assert node == self.xmlNodeFromString(doc),\
            'Could not set label of property in output volumic variables model'

    def checkSetAndGetPrintingStatus(self):
        """
        Check whether the OutputVolumicVariablesModel class could be
        set and get status for printing listing
        """
        from ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

        mdl = OutputVolumicVariablesModel(self.case)
        mdl.setPrintingStatus('Temp.C', 'off')
        node_out = mdl.case.xmlGetNode('additional_scalars')
        doc = '''<additional_scalars>
                    <scalar label="Temp.C" name="temperature_celsius" type="thermal">
                        <initial_value zone="1">20.0</initial_value>
                        <min_value>-1e+12</min_value>
                        <max_value>1e+12</max_value>
                        <listing_printing status="off"/>
                    </scalar>
                 </additional_scalars>'''

        assert node_out == self.xmlNodeFromString(doc),\
            'Could not set status of listing printing in output volumic variables model'
        assert mdl.getPrintingStatus('Temp.C') == 'off',\
            'Could not get status of listing printing in output volumic variables model'

    def checkSetAndGetPostStatus(self):
        """
        Check whether the OutputVolumicVariablesModel class could be
        set and get status for printing
        """
        from ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

        mdl = OutputVolumicVariablesModel(self.case)
        mdl.setPostStatus('Temp.C', 'off')
        node_out = mdl.case.xmlGetNode('additional_scalars')
        doc = '''<additional_scalars>
                    <scalar label="Temp.C" name="temperature_celsius" type="thermal">
                        <initial_value zone="1">20.0</initial_value>
                        <min_value>-1e+12</min_value>
                        <max_value>1e+12</max_value>
                        <postprocessing_recording status="off"/>
                    </scalar>
                 </additional_scalars>'''

        assert node_out == self.xmlNodeFromString(doc),\
            'Could not set status of post processing in output volumic variables model'
        assert mdl.getPostStatus('Temp.C') == 'off',\
            'Could not get status of post processing in output volumic variables model'

    def checkSetAndGetPostStatusForRadiativeProperties(self):
        """
        Check whether the OutputVolumicVariablesModel class could be
        set and get status for post processing of radaitive property
        """
        from Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')
        del ThermalRadiationModel

        mdl = OutputVolumicVariablesModel(self.case)
        mdl.setPostStatus('Srad', 'off')
        node_out = mdl.case.xmlGetNode('radiative_transfer')

        doc = '''<radiative_transfer model="dom">
                    <property label="Srad" name="srad">
                        <postprocessing_recording status="off"/>
                    </property><property label="Qrad" name="qrad"/>
                    <property label="Absorp" name="absorp"/>
                    <property label="Emiss" name="emiss"/>
                    <property label="CoefAb" name="coefAb"/>
                    <property label="Wall_temp" name="wall_temp" support="boundary"/>
                    <property label="Flux_incident" name="flux_incident" support="boundary"/>
                    <property label="Th_conductivity" name="thermal_conductivity" support="boundary"/>
                    <property label="Thickness" name="thickness" support="boundary"/>
                    <property label="Emissivity" name="emissivity" support="boundary"/>
                    <property label="Flux_net" name="flux_net" support="boundary"/>
                    <property label="Flux_convectif" name="flux_convectif" support="boundary"/>
                    <property label="Coeff_ech_conv" name="coeff_ech_conv" support="boundary"/>
                    <restart status="off"/>
                    <directions_number>32</directions_number>
                    <absorption_coefficient type="constant">0</absorption_coefficient>
                 </radiative_transfer>'''

        assert node_out == self.xmlNodeFromString(doc),\
        'Could not set status of post processing for radiative property \
                   in output volumic variables model'
        assert mdl.getPostStatus('Srad') == 'off',\
        'Could not get status of post processing for radiative property \
                   in output volumic variables model'


def suite():
    testSuite = unittest.makeSuite(OutputVolumicVariablesModelTestCase, "check")
    return testSuite

def runTest():
    print "OutputVolumicVariablesModelTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
