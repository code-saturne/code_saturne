# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
Output of volume variables
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.ThermalRadiationModel import ThermalRadiationModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesModel(Variables, Model):

    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_models    = self.case.xmlGetNode('thermophysical_models')
        self.analysis_ctrl  = self.case.xmlGetNode('analysis_control')
        self.fluid_prop     = self.case.xmlGetNode('physical_properties')
        self.have_multizone_phys_prop = self.__has_multiple_phys_prop_nodes__()

        if self.node_models:
            self.node_model_vp  = self.node_models.xmlGetNode('velocity_pressure')
            self.node_ale       = self.node_models.xmlGetChildNode('ale_method')
        else:
            self.node_model_vp  = None
            self.node_ale       = None
        if self.analysis_ctrl:
            self.node_output    = self.analysis_ctrl.xmlInitNode('output')
            self.node_means     = self.analysis_ctrl.xmlGetNode('time_averages')
            self.node_error     = self.analysis_ctrl.xmlGetNode('error_estimator')
        else:
            self.node_output    = None
            self.node_means     = None
            self.node_error     = None
        if self.node_output:
            self.node_probe     = self.node_output.xmlGetNodeList('probe','name')
        else:
            elf.node_probe     = None

        self.node_lagr_vstats = None
        node_lagr  = self.case.xmlGetNode('lagrangian', 'model')
        if node_lagr:
            node_stat = node_lagr.xmlGetChildNode('statistics')
            if node_stat:
                self.node_lagr_vstats = node_stat.xmlGetChildNode('volume')

        self.updateList()


    # Following private methods: (to see for gathering eventually)


    def __has_multiple_phys_prop_nodes__(self):
        """
        Determine if we have mutiple physical properties nodes
        """
        # We do not use the LocalizationModel here, as it would lead to
        # circular dependencies, so directly use XML definitions.
        m_pp_zones = 0
        n = self.case.xmlGetNode('solution_domain')
        if n:
            n = n.xmlGetNode('volumic_conditions')
            if n:
                nl = n.xmlGetChildNodeList('zone', 'label', 'id')
                for n in nl:
                    if n['physical_properties'] == 'on':
                        m_pp_zones += 1
                        if m_pp_zones > 1:
                            return True
        return False


    def _defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['status']    = "on"
        default['estimator'] = "0"

        return default


    def __getVolFieldsNodeList(self, root, constant, time_averages):
        """
        Returns the list of active volume fields (variables, propeties,
        and possibly postprocessing fields).
        """

        tags = ['variable', 'property']
        if time_averages:
            tags.append('time_average')

        l = []
        for tag in tags:
            for node in root.xmlGetNodeList(tag):
                if node['support']:
                    if node['support'] == 'boundary':
                        continue
                if tag == 'property':
                    if not (self.have_multizone_phys_prop or constant):
                        choice = node['choice']
                        if choice and choice == 'constant':
                            continue
                if node.xmlGetAttribute('status', default='on') == 'on':
                    l.append(node)

        return l


    def getVolumeFieldsNodeList(self, constant=False, time_averages=True):
        """
        Returns the list of active volume fields (variables, properties,
        and possibly postprocessing fields).
        """

        l = self.__getVolFieldsNodeList(self.case, constant, time_averages)

        n_bc = self.case.xmlGetNode('boundary_conditions')
        if n_bc:
            l_c = self.__getVolFieldsNodeList(n_bc, constant, time_averages)
            for n in l_c:
                l.remove(n)

        return l


    def getVolumeFieldsLabel2Name(self,
                                  constant=False,
                                  time_averages=True,
                                  get_components=False):
        """
        Returns a dictionnary to relate name and label to volume fields
        (variables, properties) for time averages and profiles.
        """

        dicoLabel2Name = {}

        for node in self.getVolumeFieldsNodeList(constant, time_averages):

            name = self.__nodeName__(node)
            label = node['label']
            if not label:
                label = name

            # Check that we wish to show the field
            if not (node['support'] and node['support'] == "boundary"):

                # Check if we want to also get components
                dim = node['dimension']

                # Add main field. "-1" is used as default since for
                # vectors/tensors it yields all components
                dicoLabel2Name[label] = (name, str(-1), str(dim))

                if dim and int(dim) > 1 and get_components:
                    # If we consider the Rij tensor, the user will see
                    # R11, R22, ... in the GUI instead of Rij[0], Rij[1], ...
                    # This choice was considered as the clearest.
                    if name == 'rij':
                        rij_lbls = ['R11', 'R22', 'R33', 'R12', 'R23', 'R13']
                        for ii in range(int(dim)):
                            label1 = rij_lbls[ii]
                            dicoLabel2Name[label1] = (name, str(ii))

                    elif 'reynolds_stress' in name:
                        comp = label.split('reynolds_stress')[-1]
                        rij_lbls = ['reynolds_stress11',
                                    'reynolds_stress22',
                                    'reynolds_stress33',
                                    'reynolds_stress12',
                                    'reynolds_stress23',
                                    'reynolds_stress13']
                        for ii in range(int(dim)):
                            label1 = rij_lbls[ii] + comp
                            dicoLabel2Name[label1] = (name, str(ii))

                    else:
                        for ii in range(int(dim)):
                            label1 = label + "[" + str(ii) + "]"
                            dicoLabel2Name[label1] = (name, str(ii))

                is_scalar = not( (dim != None) and (int(dim)>1) )
                if is_scalar:
                    dicoLabel2Name[label] = (name, None)

        return dicoLabel2Name

    def getVariablesAtNode(self, node, time_averages=False, get_components=False):
        recognized_variables = []
        dicoLabel2Name =self.getVolumeFieldsLabel2Name(time_averages=time_averages, get_components=get_components)
        for xml_variable in node.xmlGetChildNodeList('var_prop'):
            for name, label in dicoLabel2Name.items():
                component = label[1]
                is_recognized = (label[0] == xml_variable['name'])
                if component:
                    is_recognized = (label == (xml_variable['name'], xml_variable['component']))
                if is_recognized:
                    recognized_variables.append(name)
        return recognized_variables

    def setVariablesAtNode(self, node, variables, time_averages=False, get_components=False):
        node.xmlRemoveChild('var_prop')
        dicoLabel2Name = self.getVolumeFieldsLabel2Name(time_averages=time_averages,
                                      get_components=get_components)
        authorized_variables = dicoLabel2Name.keys()
        for var in variables:
            self.isInList(var, authorized_variables)
            (name, comp) = dicoLabel2Name[var]
            if comp is None:
                comp = -1
            node.xmlAddChild('var_prop', name=name, component=comp)

    def __updateListFilter(self, node_list, category, remain_list=None):
        """
        Update name list using a location filter
        """
        self.listNodeVolum = self.getVolumeFieldsNodeList()

        list_nodes = self.getVolumeFieldsNodeList()

        # Main/known categories

        for node in node_list:
            if node['support']:
                if node['support'] == 'boundary':
                    continue
            elif node.xmlGetAttribute('status', default='on') == 'off':
                continue
            name = self.__nodeName__(node)
            label = node['label']
            if not label:
                label = name
            self.dicoLabelName[name] = label
            self.list_name.append([self.__nodeName__(node), category])

            if remain_list:
                for node in node_list:
                    try:
                        remain_list.remove(node)
                    except Exception:
                        pass


    def updateList(self):
        """
        Update list of volume nodes and matching names and category
        """
        self.listNodeVolum = self.getVolumeFieldsNodeList()

        list_remain = self.getVolumeFieldsNodeList()

        self.dicoLabelName = {}
        self.list_name = []

        fields_node = None
        fields_list = {}

        # Main/known categories

        self.__updateListFilter(self.__getListOfVelocityPressureVariables__(),
                                'Base', list_remain)
        self.__updateListFilter(self.__getListOfTurbulenceVariables__(),
                                'Turbulence', list_remain)
        self.__updateListFilter(self.getThermalScalar(),
                                'Thermal', list_remain)
        self.__updateListFilter(self.__getThermalRadiativeProperties__(),
                                'Thermal', list_remain)
        self.__updateListFilter(self.__getPuCoalScalProper__(),
                                'Coal', list_remain)
        self.__updateListFilter(self.__getGasCombScalProper__(),
                                'Gas', list_remain)
        self.__updateListFilter(self.__getMeteoScalProper__(),
                                'Atmospheric', list_remain)
        self.__updateListFilter(self.__getElecScalProper__(),
                                'Electric', list_remain)
        self.__updateListFilter(self.getAdditionalScalar(),
                                'Other', list_remain)
        self.__updateListFilter(self.__getAdditionalScalarProperty__(),
                                'Other', list_remain)
        self.__updateListFilter(self.__getFluidProperty__(),
                                'Physical properties', list_remain)
        self.__updateListFilter(self.__getTimeProperty__(),
                                'Other', list_remain)
        self.__updateListFilter(self.__getListOfTimeAverage__(),
                                'Time moments', list_remain)
        self.__updateListFilter(self.__getListOfLagrVolumeStats__(),
                                'Lagrangian statistics', list_remain)
        self.__updateListFilter(self.__getListOfAleMethod__(),
                                'Other', list_remain)
        self.__updateListFilter(self.__getListOfEstimator__(),
                                'Estimator', list_remain)

        # Remaining nodes

        for node in list_remain:

            # location filter

            if node['support']:
                if node['support'] == 'boundary':
                    continue

            # Name and label

            if node['name']:
                name = self.__nodeName__(node)
            else:
                name = node['label']
            if not name:
                continue
            if not node['label']:
                node['label'] = name

            self.dicoLabelName[name] = node['label']

            # Category based on parent node
            # (could also be based on dictionary)

            category = node.xmlGetParentName()
            if category == 'fluid_properties':
                category = 'Physical properties'
            else:
                category = category.replace('_', ' ').capitalize()

            # For NCFD multiphase, use the field_id as a parent category
            if node['field_id']:
                if node['field_id'] != 'none':
                    if fields_node is None:
                        fields_node = self.node_models.xmlGetNode('fields')
                        for nf in fields_node.xmlGetNodeList('field'):
                            fields_list[nf['field_id']] = nf['label']

                    category = fields_list[node['field_id']]



            self.list_name.append([self.__nodeName__(node), category])


    def __getNode__(self, name):
        """
        Return a node matching a given name
        """

        for node in self.listNodeVolum:
            n_name = node['name']
            if node['field_id']:
                if node['field_id'] != 'none':
                    n_name += '_' + str(node['field_id'])
            if n_name == name:
                return node

        return None


    def __nodeName__(self, node):
        """
        Return a node matching a given name
        """

        n_name = node['name']
        if not n_name:
            n_name = node['label']
        else:
            if node['field_id']:
                if node['field_id'] != 'none':
                    n_name += '_' + str(node['field_id'])

        return n_name


    def __getListOfVelocityPressureVariables__(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        if self.node_model_vp:
            for tag in ('variable', 'property'):
                for node in self.node_model_vp.xmlGetNodeList(tag):
                    nodeList.append(node)
        return nodeList


    def __getListOfTurbulenceVariables__(self):
        """
        Private method: return node of properties of weight matrix
        """
        node_models = self.case.xmlGetNode('thermophysical_models')
        if node_models:
            node_turb = node_models.xmlGetNode('turbulence', 'model')

        nodeList = []
        if node_turb:
            for tag in ('variable', 'property'):
                for node in node_turb.xmlGetNodeList(tag):
                    if node['name'] == 'rij':
                        node['dimension'] = 6
                    nodeList.append(node)
        return nodeList


    def __getListOfEstimator__(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        if self.node_error:
            for node in self.node_error.xmlGetNodeList('property'):
                nodeList.append(node)
        return nodeList


    def __getWeightMatrixProperty__(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []

        node0 = self.case.xmlGetNode('numerical_parameters')
        if node0:
            node1 = node0.xmlGetNode('velocity_pressure_coupling', 'status')
            if node1:
                if node1['status'] == 'on':
                    nodeList = node0.xmlGetNodeList('property')
        return nodeList


    def __getListOfAleMethod__(self):
        """
        Private method: return list of variables and properties for ALE method
        """
        nodeList = []
        if self.node_ale:
            if self.node_ale['status'] == 'on':
                for tag in ('variable', 'property'):
                    for node in self.node_ale.xmlGetChildNodeList(tag):
                        nodeList.append(node)

        return nodeList


    def __getThermalRadiativeProperties__(self):
        """
        Private method: return list of volumic properties for thermal radiation
        """
        nodeList = []
        if ThermalRadiationModel(self.case).getRadiativeModel() != "off":
            self.node_ray = self.node_models.xmlGetNode('radiative_transfer')
            for node in self.node_ray.xmlGetChildNodeList('property'):
                nodeList.append(node)
        return nodeList


    @Variables.noUndo
    def getThermalScalar(self):
        """
        Return node of thermal scalar (idem ds NumericalParamEquationModel)
        """
        node_models = self.case.xmlGetNode('thermophysical_models')
        node = node_models.xmlGetNode('thermal_scalar')
        if node:
            return node.xmlGetNodeList('variable', type='thermal')
        else:
            return []


    @Variables.noUndo
    def __getPuCoalScalProper__(self):
        """
        Return list fo nodes of pulverized coal.
        """
        varList = []
        if self.node_models:
            node = self.node_models.xmlGetNode('solid_fuels', 'model')
            if node:
                model = node['model']
                if model != 'off':
                    for var in ('variable', 'property'):
                        nodList = node.xmlGetNodeList(var)
                        for nodvar in nodList:
                            varList.append(nodvar)
        return varList


    @Variables.noUndo
    def __getGasCombScalProper__(self):
        """
        Return list of nodes of gas combustion.
        """
        varList = []
        if self.node_models:
            node = self.node_models.xmlGetNode('gas_combustion', 'model')
            if node:
                model = node['model']
                if model != 'off':
                    for var in ('variable', 'property'):
                        nodList = node.xmlGetNodeList(var)
                        for nodvar in nodList:
                            varList.append(nodvar)
        return varList


    @Variables.noUndo
    def __getMeteoScalProper__(self):
        """
        Return list fo nodes of atmospheric flows.
        """
        varList = []
        if self.node_models:
            node = self.node_models.xmlGetNode('atmospheric_flows', 'model')
            if node:
                model = node['model']
                if model != 'off':
                    for var in ('variable', 'property'):
                        nodList = node.xmlGetNodeList(var)
                        for nodvar in nodList:
                            varList.append(nodvar)
        return varList


    @Variables.noUndo
    def __getElecScalProper__(self):
        """
        Return list fo nodes of electric flows.
        """
        varList = []
        if self.node_models:
            node = self.node_models.xmlGetNode('joule_effect', 'model')
            if not node: return []
            model = node['model']
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
        varList = []
        if self.node_models:
            node = self.node_models.xmlGetNode(model_name, 'model')
            if node:
                model = node['model']
                if model != 'off':
                    nodList = node.xmlGetNodeList('variable')
                    for nodvar in nodList:
                        varList.append(nodvar)
        return varList


    @Variables.noUndo
    def getAdditionalScalar(self):
        """
        Return list of nodes of user scalars
        """
        node = self.case.xmlGetNode('additional_scalars')
        if node:
            return node.xmlGetNodeList('variable', type='user')
        else:
            return []


    @Variables.noUndo
    def __getAdditionalScalarProperty__(self):
        """
        Return list of nodes of properties of user scalars
        """
        nodeList = []
        for node in self.getAdditionalScalar():
            L = node.xmlGetNode('property', choice='variable')
            if L:
                nodeList.append(L)
        return nodeList


    @Variables.noUndo
    def __getFluidProperty__(self):
        """
        Return list of nodes of fluid properties
        """
        nodeList = []

        if self.fluid_prop:
            node = self.fluid_prop.xmlGetNode('fluid_properties')

            if node:
                for prop in ('density',
                             'specific_heat',
                             'thermal_conductivity'):
                    L = node.xmlGetNode('property', name=prop, choice='variable')
                    if not L:
                        L = node.xmlGetNode('property', name=prop, choice='user_law')
                    if self.have_multizone_phys_prop and not L:
                        L = node.xmlGetNode('property', name=prop, choice='constant')
                    if L:
                        nodeList.append(L)
                for prop in ('molecular_viscosity',
                             'dynamic_diffusion'):
                    L = node.xmlGetNode('property', name=prop, choice='variable')
                    if not L:
                        L = node.xmlGetNode('property', name=prop, choice='user_law')
                    if L:
                        nodeList.append(L)

        return nodeList


    @Variables.noUndo
    def __getTimeProperty__(self):
        """
        Return list fo nodes of properties of time_parameters.
        """
        nodeList = []

        if self.analysis_ctrl:
            node1 = self.analysis_ctrl.xmlGetNode('time_parameters')
        else:
            node1 = None

        if node1:
            if node1.xmlGetInt('time_passing'):
                node2 = node1.xmlGetNode('property', name='local_time_step')
                if node2:
                    nodeList.append(node2)

            # We are looping on the list since for NCFD we have one of each
            # for each phase. If not courant/fourier is available, the list
            # is empty, hence the second for loop does nothing.
            for prop in ('courant_number', 'fourier_number'):
                for _n in node1.xmlGetChildNodeList('property', name=prop):
                    nodeList.append(_n)

        return nodeList


    @Variables.noUndo
    def __getListOfTimeAverage__(self):
        """
        Return list of time averages variables
        """
        nodeList = []
        if self.node_means:
            for node in self.node_means.xmlGetNodeList('time_average'):
                nodeList.append(node)

        return nodeList


    def __getListOfLagrVolumeStats__(self):
        """
        Private method: return node of properties of weight matrix
        """
        nodeList = []
        if self.node_lagr_vstats:
            for node in self.node_lagr_vstats.xmlGetNodeList('property'):
                nodeList.append(node)

        return nodeList


    #Following methods only called by the View
    @Variables.noUndo
    def getLabelsList(self):
        """
        Return list of labels for all variables, properties .....Only for the View
        """
        lst = []
        for node in self.listNodeVolum:
            lst.append(node['label'])
        return lst


    @Variables.noUndo
    def getNamesList(self):
        """
        Return list of names for all variables, properties .....Only for the View
        """
        lst = []
        for node in self.listNodeVolum:
            lst.append(self.__nodeName__(node))
        return lst


    @Variables.noUndo
    def getProbeList(self):
        """ Return list of node for probes """
        probeList = []
        for node in self.node_probe:
            probeList.append(self.__nodeName__(node))
        return probeList


    @Variables.noUndo
    def getPrintingStatus(self, name):
        """
        Return status of markup printing from node with name. Only for the View
        """
        self.isInList(name, self.getNamesList())
        status = self._defaultValues()['status']
        node = self.__getNode__(name)
        if node:
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
        node = self.__getNode__(name)
        if node:
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
        node = self.__getNode__(name)
        if node:
            node_post = node.xmlGetChildNode('probes_recording', 'status')
            if node_post:
                status = node_post['status']
        return status


    @Variables.undoLocal
    def setVariableLabel(self, old_label, new_label):
        """
        Replace old_label by new_label for node with name and old_label.
        Only for the View
        """
        # fusion de cette methode avec DefineUserScalarsModel.renameScalarLabel
        self.isInList(old_label, self.getLabelsList())
        self.isNotInList(new_label, [""])

        if old_label != new_label:
            self.isNotInList(new_label, self.getLabelsList())
        for node in self.listNodeVolum:
            if node['label'] == old_label:
                node['label'] = new_label

        self.updateList()
        self._updateBoundariesNodes(old_label, new_label)

        for node in self.case.xmlGetNodeList('formula'):
            f = node.xmlGetTextNode()
            if f:
                f.replace(old_label, new_label)
                node.xmlSetTextNode(f)


    def _updateBoundariesNodes(self, old_label, new_label):
        """
        Update good label for boundaries nodes with name and label.
        Only for the View
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
        node = self.__getNode__(name)
        if node:
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
        self.isInList(name, self.getNamesList())
        node = self.__getNode__(name)
        if node:
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
        node = self.__getNode__(name)
        if node:
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

        if self.node_error:
            nn = self.node_error.xmlGetChildNode(name, 'model')
            if nn:
                status = nn['model']
                if status == "1":
                    status = "2"

        return status


    @Variables.undoLocal
    def setEstimatorModel(self, name, model):
        """
        Put model for an error estimator
        """
        self.isInList(model, ['0', '1', '2'])
        self.isInList(name, ["Correction", "Drift", "Prediction", "Total"])
        status = self._defaultValues()['estimator']

        if not self.analysis_ctrl:
            self.analysis_ctrl = self.case.xmlInitNode('analysis_control')
        if not self.node_error:
            self.node_error = self.analysis_ctrl.xmlInitNode('error_estimator')

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
                self.setNewProperty(nn, "est_error_cor_2")
            elif name == "Drift":
                self.setNewProperty(nn, "est_error_der_2")
            elif name == "Prediction":
                self.setNewProperty(nn, "est_error_pre_2")
            elif name == "Total":
                self.setNewProperty(nn, "est_error_tot_2")

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
        from code_saturne.model.ThermalScalarModel import ThermalScalarModel
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
        from code_saturne.model.ThermalScalarModel import ThermalScalarModel
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
