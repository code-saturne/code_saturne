# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module defines the XML data model in which the user defines the physical
options of the treated case.

This module contains the following classe:
- XMLinit
- XMLinitTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest, re

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLinitialize import BaseXmlInit
from code_saturne.model.LocalizationModel import Zone, LocalizationModel
from code_saturne.model.OutputControlModel import OutputControlModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel
from code_saturne.model.NeptuneWallTransferModel import NeptuneWallTransferModel
from code_saturne.model.NucleateBoilingModel import NucleateBoilingModel

from code_saturne.model.EosWrapper import eosWrapper

#-------------------------------------------------------------------------------
# class XMLinitNeptune for neptune_cfd solver
#-------------------------------------------------------------------------------

class XMLinitNeptune(BaseXmlInit):
    """
    This class initializes the XML parameter file for the neptune_cfd solver.
    """
    def initialize(self, prepro = False):
        """
        Verify that all Headings exist only once in the XMLDocument and
        create the missing heading.
        """
        msg = self.__initHeading(prepro)
        if msg:
            return msg

        OutputControlModel(self.case).addDefaultWriter()
        OutputControlModel(self.case).addDefaultMesh()

        if not prepro:
            self._backwardCompatibility()

            # Initialization (order is important)

            self.XMLNodeAna      = self.case.xmlGetNode('analysis_control')
            self.XMLNodeAverage  = self.XMLNodeAna.xmlInitNode('time_averages')
            self.XMLThermo       = self.case.xmlGetNode('thermophysical_models')
            self.__XMLThermo     = self.XMLThermo.xmlInitNode('thermodynamics')
            self.XMLNodeVariable = self.XMLThermo.xmlInitNode('variables')
            OutputControlModel(self.case).addDefaultWriter()
            OutputControlModel(self.case).addDefaultMesh()
            self.XMLUserScalar   = self.case.xmlGetNode('additional_scalars')
            self.XMLUser         = self.XMLUserScalar.xmlInitNode('users')

            # Add surfacic properties to post-process
            # TODO : check if NCFD version allows it !!!
            self.XMLNodeProperty = self.XMLThermo.xmlInitNode('properties')
            node = self.XMLNodeProperty
            n = self.setNewProperty(node, 'stress')
            n['support'] = 'boundary'
            n['label'] = 'Stress'
            if not node.xmlGetChildNode('property', name='stress_tangential'):
                n = self.setNewProperty(node, 'stress_tangential')
                n['label'] = 'Stress, tangential'
                n['support'] = 'boundary'
                n.xmlInitNode('postprocessing_recording')['status']= "off"
            if not node.xmlGetChildNode('property', name='stress_normal'):
                n = self.setNewProperty(node, 'stress_normal')
                n['label'] = 'Stress, normal'
                n['support'] = 'boundary'
                n.xmlInitNode('postprocessing_recording')['status']= "off"


            # First Volume Zone definition for all cells -> initialization

            zones = LocalizationModel("VolumicZone", self.case).getZones()
            iok = 0
            for zone in zones:
                if zone.getLabel() == 'all_cells':
                    iok = 1
            if iok == 0:
                zone = Zone("VolumicZone",
                            case=self.case,
                            label='all_cells',
                            localization='all[]',
                            nature="initialization:physical_properties")
                LocalizationModel("VolumicZone", self.case).addZone(zone)
                zone = LocalizationModel("VolumicZone", self.case).getCodeNumberOfZoneLabel('all_cells')

            # If EOS is not avalaible, check if EOS is needed by the current set up.

            if eosWrapper().isActive() == 0:
                from code_saturne.model.MainFieldsModel import MainFieldsModel
                from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel

                for fieldId in MainFieldsModel(self.case).getFieldIdList():
                    if ThermodynamicsModel(self.case).getMaterials(fieldId) != "user_material":
                        msg = "The current GUI does not found EOS, but this file of parameters has" \
                              " been generated with EOS. \n\n Please check the availability of "  \
                              "the prerequisite EOS."

            # Initialize fields
            if not self.XMLThermo.xmlGetNode('fields'):
                from code_saturne.model.MainFieldsModel import MainFieldsModel
                MainFieldsModel(self.case).setPredefinedFlow("None")
                del MainFieldsModel

        return msg


    def __initHeading(self, prepro):
        """
        Create if necessary headings from the root element of the case.
        """
        msg = ""
        tagList = ('solution_domain',
                   'analysis_control',
                   'calculation_management')
        if not prepro :
            tagList += ('thermophysical_models',
                        'additional_scalars',
                        'closure_modeling',
                        'boundary_conditions',
                        'numerical_parameters')

        for tag in tagList:
            nodeList = self.case.root().xmlInitChildNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. " \
                      "There is more than one occurence of the tag: \n\n" + tag +  \
                      "\n\nThe application will finish. Sorry."

        for tag in tagList:
            nodeList = self.case.xmlInitNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. " \
                      "There is more than one occurence of the tag: \n\n" + tag +  \
                      "\n\nThe application will finish. Sorry."

        return msg


    def _backwardCompatibilityOldVersion(self, from_vers):
        """
        Change XML in order to ensure backward compatibility for old version
        there is nothing to do for 2.1 to 2.2
        """
        if (from_vers <= "-1.0") or (from_vers[:3] < "3.0"):
            self.__backwardCompatibilityFrom_2_0()
            self.__backwardCompatibilityFrom_2_2()

        if from_vers[:3] < "5.0":
            if from_vers[:3] < "4.3":
                self.__backwardCompatibilityFrom_4_2()

            self.__backwardCompatibilityFrom_4_4()

        if from_vers[:3] < "7.0":
            if from_vers[:3] < "6.1":
                self.__backwardCompatibilityFrom_6_0()
            if from_vers[:3] < "6.2":
                self.__backwardCompatibilityFrom_6_1()
            if from_vers[:3] < "6.4":
                self.__backwardCompatibilityFrom_6_3()

        if from_vers[:3] < "8.0":
            if from_vers[:3] < "7.1":
                self.__backwardCompatibilityFrom_7_0()
            if from_vers[:3] < "7.2":
                self.__backwardCompatibilityFrom_7_1()


    def __backwardCompatibilityFrom_2_0(self):
        """
        Change XML in order to ensure backward compatibility from 2.0 to 2.1
        """
        # time averages
        idTA = 0
        for Node in self.case.xmlGetNodeList('time_average'):
            idTA = idTA + 1
            if Node['name'] is None:
                Node['name'] = "TimeAverage_" + str(idTA)

        # Profiles
        for node in self.case.xmlGetNodeList('profile'):
            if node:
                n = node.xmlGetNode("output_type")
                if n is None:
                    freq = node.xmlGetInt("output_frequency")
                    if freq == -1:
                        node.xmlSetData('output_type', "end")
                    else:
                        node.xmlSetData('output_type', "frequency")

        # Users
        nodeInit = self.case.xmlGetNode('users')
        if nodeInit:
            nodeList = nodeInit.xmlGetNodeList('variable')
            for i in range(len(nodeList)):
                varNode = nodeList[i]
                if varNode['support'] is None:
                    varNode['support'] = "cells"

        # <enthalpy_model field_id="2" model="user_function"/>
        for node in self.case.xmlGetNodeList('enthalpy_model'):
            if node['model'] == "user_function":
                node['model'] = "no_source_term"

        nfield = self.case.xmlGetNode('fields')
        if nfield != None:
            nfield1 = nfield.xmlGetNode('field', field_id="1")
            nh = nfield1.xmlGetNode('hresolution')
            if nh['status'] == "off":
                nodeSurf = self.case.xmlGetNode('property', name='surface_tension')
                if nodeSurf is None:
                    XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
                    XMLNodeproperty = XMLNodethermo.xmlInitNode('properties')
                    Variables(self.case).setNewVariableProperty("property", "constant", XMLNodeproperty, "none", "surface_tension", "Surf_tens")

        XMLNodeNonCondens = self.case.xmlGetNode('non_condensable_list')
        if XMLNodeNonCondens != None:
            nodeList = XMLNodeNonCondens.xmlGetNodeList('variable')
            for i in range(len(nodeList)):
                oldName = "non_condensable" + str(i+1)
                newName = "MassFractionNonCondensableGas_" + str(i+1)
                for node in self.case.xmlGetNodeList('variable', name = oldName):
                    if node != None:
                        node['name'] = newName

        # update name for node property
        dicoO2N = {"drho_dh"                             : "d_rho_d_h",
                   "drho_dP"                             : "d_rho_d_P",
                   "Hsat1"                               : "SaturationEnthalpyLiquid",
                   "Hsat2"                               : "SaturationEnthalpyGas",
                   "TsatK"                               : "SaturationTemperature",
                   "DTSDPDerivative"                     : "d_Tsat_d_P",
                   "DHSDPDerivativeLiquid"               : "d_Hsat_d_P_Liquid",
                   "DHSDPDerivativeGas"                  : "d_Hsat_d_P_Gas"}

        for node in self.case.xmlGetNodeList('property'):
            if node['name'] in dicoO2N.keys():
                old_name = node['name']
                node['name'] = dicoO2N[old_name]


    def __backwardCompatibilityFrom_2_2(self):
        """
        Change XML in order to ensure backward compatibility from 2.2 to 2.4-alpha.
        """
        for node in self.case.xmlGetNodeList('hresolution'):
            if node['model'] is None:
                if node['status'] == 'off':
                    node['model'] = 'off'
                else:
                    node['model'] = 'total_enthalpy'


    def __backwardCompatibilityFrom_4_2(self):
        """
        Change XML in order to ensure backward compatibility from versions prior to 4.3
        Reason: Renaming of wall_temperature as boundary_temperature
        """
        # For versions prior to 5.0,renaming of wall_temperature as boundary_temperature
        for node in self.case.xmlGetNodeList('property'):
            if node['name'] == 'wall_temperature':
                node['name']  = 'boundary_temperature'

    def __backwardCompatibilityFrom_4_4(self):

        # For versions prior to 5.0,renaming of wall_temperature as boundary_temperature
        for node in self.case.xmlGetNodeList('property'):
            if node['name'] == 'wall_temperature':
                node['name']  = 'boundary_temperature'

            if node['name'] == 'wall_friction_velocity':
                self.case.xmlRemoveChild('property',
                                         name='wall_friction_velocity',
                                         field_id='none')


        # Add the choice between SGDH and GGDH turbulent thermal flux models
        cnode = self.case.xmlGetNode('closure_modeling')
        tnode = cnode.xmlGetNode('turbulence')

        thermo_node = self.case.xmlGetNode('thermophysical_models')
        fnode = thermo_node.xmlGetNode('fields')

        from code_saturne.model.MainFieldsModel import MainFieldsModel
        field_names = MainFieldsModel(self.case).getFieldLabelsList()

        turb_dico = {'TurbDissip':'epsilon',
                     'epsilon':'epsilon',
                     'TurbKineEner_k':'k',
                     'k':'k',
                     'turb_viscosity':'turb_viscosity',
                     'ReynoldsStress':'reynolds_stress'}

        rij_comp = {'ReynoldsStressXX':'0',
                    'ReynoldsStressYY':'1',
                    'ReynoldsStressZZ':'2',
                    'ReynoldsStressXY':'3',
                    'ReynoldsStressYZ':'4',
                    'ReynoldsStressXZ':'5'}

        if tnode != None:
            tvn   = tnode.xmlGetNode('variables')
            for node in tnode.xmlGetNodeList('field'):
                if node['turb_flux'] is None:
                    node['turb_flux'] = 'sgdh'

            # Check for missing alpha in the EBRSM model
            for node in tnode.xmlGetNodeList('field'):
                if node['model'] == 'rij-epsilon_ebrsm':
                    fid = node['field_id']
                    na = tvn.xmlGetNode('variable',name='alpha', field_id=fid)
                    if na is None:
                        Variables(self.case).setNewVariableProperty("variable", "",
                                                                    tvn,
                                                                    fid,
                                                                    'alpha',
                                                                    'alpha_'+field_names[int(fid)-1])


            # Renaming of Rij tensor
            for node in fnode.xmlGetNodeList('field'):
                fieldId = node['field_id']
                field_name = field_names[int(fieldId)-1]

                rn = tvn.xmlGetNode('variable',
                                    name="ReynoldsStressXX",
                                    field_id=fieldId)
                if rn != None:
                    rn['name']  = "reynolds_stress"
                    rn['label'] = "reynolds_stress_"+field_name
                    rn['dimension']   = 6

                    for comp in ["XY", "XZ", "YY", "YZ", "ZZ"]:
                        tvn.xmlRemoveChild('variable',
                                           name="ReynoldsStress"+comp,
                                           field_id=fieldId)

            # Renaming k and espilon
            for node in tnode.xmlGetNodeList("variable")+tnode.xmlGetNodeList("property"):
                fieldId = node['field_id']
                field_name = field_names[int(fieldId)-1]
                for tv in turb_dico.keys():
                    if tv in node['name']:
                        node['name']  = turb_dico[tv]
                        node['label'] = turb_dico[tv]+"_"+field_name



        # Modify the rad transfer xml node name for particles to allow a correct
        # workflow with the RTE SOLVER
        tpnode = self.case.xmlGetNode('thermophysical_models')
        fnode  = tpnode.xmlGetNode('fields')
        if fnode != None:
            for node in fnode.xmlGetNodeList('field'):
                rn = node.xmlGetNode('radiative_transfer')
                if rn != None:
                    st = rn['status']
                    node.xmlRemoveChild('radiative_transfer')
                    node.xmlInitChildNode('particles_radiative_transfer', status=st)

        # Renaming of Pressure
        vnode = thermo_node.xmlGetNode('variables')
        pnode = thermo_node.xmlGetNode('properties')
        ncnode = thermo_node.xmlGetNode('non_condensable_list')

        rdico = {'Enthalpy':'enthalpy',
                 'enthalpy':'enthalpy',
                 'Pressure':'pressure',
                 'pressure':'pressure',
                 'Velocity':'velocity',
                 'velocity':'velocity',
                 'VolumeFraction':'volume_fraction',
                 'volume_fraction':'volume_fraction',
                 'Temperature':'temperature',
                 'temperature':'temperature',
                 'SaturationEnthalpyLiquid':'SaturationEnthalpyLiquid',
                 'SaturationEnthalpyGas':'SaturationEnthalpyGas',
                 'mass_trans':'mass_trans',
                 'molecular_viscosity':'molecular_viscosity',
                 'specific_heat':'specific_heat',
                 'thermal_conductivity':'thermal_conductivity',
                 'drag_coefficient':'drag_coefficient',
                 'density':'density',
                 'Diameter':'diameter',
                 'diameter':'diameter',
                 'DriftComponent':'drift_component',
                 'drift_component':'drift_component',
                 'emissivity':'emissivity',
                 'elasticity':'elasticity',
                 'Xd':'Xd'}

        ldico = {'Enthalpy':'Enthalpy',
                 'enthalpy':'Enthalpy',
                 'Pressure':'Pressure',
                 'pressure':'Pressure',
                 'Velocity':'Velocity',
                 'velocity':'Velocity',
                 'VolumeFraction':'vol_f',
                 'volume_fraction':'vol_f',
                 'Temperature':'temp',
                 'temperature':'temp',
                 'SaturationEnthalpyLiquid':'HsatLiquid',
                 'SaturationEnthalpyGas':'HsatGas',
                 'mass_trans':'mass_trans',
                 'molecular_viscosity':'molecular_viscosity',
                 'specific_heat':'specific_heat',
                 'thermal_conductivity':'thermal_conductivity',
                 'drag_coefficient':'drag_coef',
                 'density':'density',
                 'Diameter':'diameter',
                 'diameter':'diameter',
                 'DriftComponent':'drift_component',
                 'drift_component':'drift_component',
                 'emissivity':'emissivity',
                 'elasticity':'elasticity',
                 'Xd':'Xd'}

        for ii in range(20):
            rdico['MassFractionNonCondensableGas_'+str(ii)]='mass_fraction_non_condensable_gas_'+str(ii)
            rdico['mass_fraction_non_condensable_gas_'+str(ii)]='mass_fraction_non_condensable_gas_'+str(ii)
            ldico['MassFractionNonCondensableGas_'+str(ii)]='mass_fraction_non_condensable_gas_'+str(ii)
            ldico['mass_fraction_non_condensable_gas_'+str(ii)]='mass_fraction_non_condensable_gas_'+str(ii)


        old_mei_names = {'VolumeFraction':['alpha','vol_f'],
                         'volume_fraction':['alpha','vol_f']}

        if ncnode != None:
            for node in ncnode.xmlGetNodeList('variable'):
                ncname = node['name']
                if 'MassFractionNonCondensableGas' in ncname:
                    node['name'] = ncname.replace('MassFractionNonCondensableGas',
                                                  'mass_fraction_non_condensable_gas')

        if vnode != None:
            for node in vnode.xmlGetNodeList('variable'):
                vname = node['name']
                if vname in rdico.keys():
                    node['name'] = rdico[vname]

                    field_id = node['field_id']
                    if field_id == 'none':
                        label = ldico[vname]
                    else:
                        label = ldico[vname] + '_' + field_names[int(field_id)-1]
                    node['label'] = label

                    for nzi in node.xmlGetNodeList('initial_value'):
                        nf = nzi.xmlGetNode('formula')
                        if nf != None:
                            f  = nzi.xmlGetString('formula')
                            if vname in old_mei_names.keys():
                                n2r = old_mei_names[vname][0]
                                n2a = old_mei_names[vname][1]
                            else:
                                n2r = vname
                                n2a = rdico[vname]
                            nf.xmlSetTextNode(f.replace(n2r, n2a))

            bcnode = self.case.xmlGetNode('boundary_conditions')
            bc_list = ['inlet', 'wall', 'outlet']
            for bc_type in bc_list:
                for nb in bcnode.xmlGetNodeList(bc_type):
                    for nv in nb.xmlGetNodeList('variable'):
                        if nv['name'] in rdico.keys():
                            nv['name'] = rdico[nv['name']]

        if pnode != None:
            for node in pnode.xmlGetNodeList('property'):
                pname = node['name']
                if pname in rdico.keys():
                    node['name'] = rdico[pname]

                    field_id = node['field_id']
                    if field_id == 'none':
                        label = ldico[pname]
                    else:
                        label = ldico[pname] + '_' + field_names[int(field_id)-1]
                    node['label'] = label

        # User arrays should be properties, not scalars
        XMLUserScalar = self.case.xmlGetNode('additional_scalars')
        if XMLUserScalar:
           XMLUser = XMLUserScalar.xmlInitNode('users')
           if XMLUser:
              for node in XMLUser.xmlGetNodeList('variable'):
                 newnode = XMLUser.xmlInitNode('property', name=node['name'])
                 newnode.xmlChildsCopy(node)
                 for tag in ('dimension', 'field_id', 'label', 'support'):
                    if node[tag]:
                       newnode[tag] = node[tag]
                 node.xmlRemoveNode()


        # Rename variable names in time averages and profile
        for pp_type in ('time_average', 'profile'):
            for node in self.case.xmlGetNodeList(pp_type):
                for vn in node.xmlGetNodeList('var_prop'):
                    vn_split = vn['name'].split('_')
                    old_name = vn['name'].split('_')[0]
                    if len(vn_split) > 1:
                        field_id = vn['name'].split('_')[1]
                    else:
                        field_id = None

                    if old_name in rdico.keys():
                        vn['name'] = rdico[old_name]
                        if field_id:
                            vn['name'] += '_' + field_id

                    elif old_name in rij_comp.keys():
                        vn['component'] = rij_comp[old_name]
                        vn['name'] = 'reynolds_stress_'
                        if field_id:
                            vn['name'] += field_id

                    elif old_name in turb_dico.keys():
                        vn['name'] = turb_dico[old_name]
                        if field_id:
                            vn['name'] += '_' + field_id

                    elif bool(re.search('MassFractionNonCondensableGas', old_name)):
                        vn['name'] = vn['name'].replace('MassFractionNonCondensableGas',
                                                        'mass_fraction_non_condensable_gas')



    def __backwardCompatibilityFrom_6_0(self):
        """
        Change XML in order to ensure backward compatibility from v6.0 to v6.1.
        """

        # Retrocompatibility: the use of probes in neptune is now the same as for code_saturne
        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            for parent in self.case.xmlGetNodeList(variableType):
                parent.xmlRemoveChild('probe_recording')
                if parent.xmlGetNode('no_probe') != None :
                    parent.xmlRemoveChild('no_probe')
                    probes_recording = parent.xmlInitNode('probes_recording', 'status')
                    probes_recording['status'] = 'off'

        timelst = []
        timetup = {}
        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                time_node = node.xmlGetNode("time_start")
                if not time_node:
                    node.xmlSetData('time_start', -1.)
                # construc map for profiles
                if node['name']:
                    timetup[node['name']] = node['label']
                    timelst.append(node['name'])
                    # now name = label
                    node['name'] = node['label']
                if node['field_id'] == "None":
                    node['field_id'] = "none"
                if not node['name']:
                    node['name'] = node['label']

                for nodevar in node.xmlGetNodeList('var_prop'):
                    name = nodevar['name']
                    if nodevar['field_id'] and nodevar['field_id'] != "none":
                        name += '_' + nodevar['field_id']
                    component = nodevar['component']
                    nodevar.xmlRemoveNode()
                    if component:
                        newnode = node.xmlInitNode('var_prop', name=name, component=component)
                    else:
                        newnode = node.xmlInitNode('var_prop', name=name)

        for node in self.case.xmlGetNodeList('profile'):
            if node:
                for nodevar in node.xmlGetNodeList('var_prop'):
                    if nodevar['name'] in timelst:
                        nodevar['name'] = timetup[nodevar['name']]
                    name = nodevar['name']
                    if nodevar['field_id'] and nodevar['field_id'] != "none":
                        name += '_' + nodevar['field_id']
                    component = nodevar['component']
                    nodevar.xmlRemoveNode()
                    if component:
                        newnode = node.xmlInitNode('var_prop',
                                                   name=name,
                                                   component=component)
                    else:
                        newnode = node.xmlInitNode('var_prop', name=name)

        # suppress gradient and flux reconstruction if needed
        for node in self.case.xmlGetNodeList('variable'):
            n = node.xmlGetNode('flux_reconstruction')
            if n:
                node.xmlRemoveChild('flux_reconstruction')
            n = node.xmlGetNode('gradient_reconstruction')
            if n:
                node.xmlRemoveChild('gradient_reconstruction')

        # update for cdudn
        self.XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
        self.__XMLNodefields = self.XMLNodethermo.xmlInitNode('fields')
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            fieldId = node['field_id']
            XMLWallNode = node.xmlGetNode('wall_model')

            if XMLWallNode:
                mdl = XMLWallNode['model']
                node.xmlRemoveChild('wall_model')

                for n in self.case.xmlGetNodeList('wall'):
                    n.xmlInitChildNode('wall_model', field_id = fieldId, model = mdl)


    def __backwardCompatibilityFrom_6_1(self):
        """
        Change XML in order to ensure backward compatibility from v6.1 to v6.2.
        """

        # Probes update to allow renaming
        XMLAnaControl = self.case.xmlGetNode('analysis_control')
        XMLOutput = XMLAnaControl.xmlGetNode('output')

        if XMLOutput.xmlGetNodeList('probe') != None:
            for node in XMLOutput.xmlGetNodeList('probe'):
                if node['id'] is None:
                    node['id'] = node['name']

    def __backwardCompatibilityFrom_6_3(self):
        """
        Change XML to ensure backward compatibility from v6.2/v6.3 to v7.0
        """

        main_xml_model = MainFieldsModel(self.case)
        forces_xml_model = InterfacialForcesModel(self.case)

        # Obtain heat and mass transfer status by checking existence of interfacial enthalpy transfer couple
        interfacial_enthalpy_node = self.case.xmlGetNode("interfacial_enthalpy")
        heat_mass_transfer_status = "off"
        if (interfacial_enthalpy_node != None):
            if interfacial_enthalpy_node.xmlGetChildNodeList("enthalpy") != []:
                heat_mass_transfer_status = "on"
        main_xml_model.setHeatMassTransferStatus(heat_mass_transfer_status)

        flow_node = self.case.xmlGetNode("predefined_flow")
        flow_choice = "None"

        if flow_node != None:
            flow_choice = flow_node.xmlGetAttribute("choice")
            # Replace (Free surface flow + LIM + bubbles for LIM) by (Multiregime)
            if flow_choice == "free_surface":
                force_node = self.case.xmlGetNode('closure_modeling').xmlGetNode('interfacial_forces').xmlGetNode(
                    'continuous_field_momentum_transfer')
                if force_node != None:
                    if (force_node.xmlGetAttribute("model") == "Large_Interface_Model") and (
                            force_node.xmlGetChildString("BubblesForLIM") == "on"):
                        flow_choice = "multiregime"
            main_xml_model.setPredefinedFlow(flow_choice,
                                             overwriteEnergy=False)

        # Continuous interfacial forces

        force_node = self.case.xmlGetNode('closure_modeling').xmlGetNode('interfacial_forces').xmlGetNode(
            'continuous_field_momentum_transfer')
        if force_node != None:
            try:
                force_model = force_node.xmlGetAttribute("model")
                surface_tension = force_node.xmlGetChildString("ITMSurfaceTension")
                interface_sharpening = force_node.xmlGetChildString("InterfaceSharpening")
                unsharpened_cells = force_node.xmlGetChildString("UnsharpenedCells")

                # For previous cases, force a drag model
                if force_model == "none":
                    force_model = forces_xml_model.defaultValuesContinuous()['continuousdragmodel']
                if flow_choice == "multiregime":
                    force_model = "G_Large_Interface_Model"

                forces_xml_model.setContinuousCouplingModel("1", "2", force_model)
                model = "none"
                if surface_tension == "on":
                    model = "Brackbill"
                forces_xml_model.setSurfaceTensionModel("1", "2", model)

                model = "none"
                if interface_sharpening == "on":
                    if unsharpened_cells == "on":
                        model = "Olsson_Partial_Interface_Sharpening"
                    else:
                        model = "Olsson_Interface_Sharpening"
                forces_xml_model.setInterfaceSharpeningModel("1", "2", model)
                force_node.xmlRemoveNode()
            except ValueError:  # Remnants of continuous interfacial forces but no complete information, so do nothing
                pass

        # Activate heat and mass transfer if separate_phase_cond model is detected
        cfc_node = self.case.xmlGetNode('closure_modeling').xmlGetNode('turbulence').xmlGetNode(
            'continuous_field_coupling')
        if cfc_node != None:
            cfc_model = cfc_node.xmlGetAttribute("model")
            if cfc_model == "separate_phase_cond":
                main_xml_model.setHeatMassTransferStatus("on")

        # Set wall transfer model type
        mass_transfer_node = self.case.xmlGetNode("mass_transfer_model")
        if mass_transfer_node != None:
            walltm_node = mass_transfer_node.xmlGetNode("wall_transfer_type")
            if (walltm_node is None) and (heat_mass_transfer_status == "on"):
                if flow_choice in ["free_surface", "boiling_flow", "multiregime"]:
                    walltm = "nucleate_boiling"
                    NucleateBoilingModel(self.case).resetToDefaultValues()
                elif flow_choice == "droplet_flow":
                    walltm = "droplet_evaporation_condensation"
                else:
                    walltm = "none"
                NeptuneWallTransferModel(self.case).wall_transfer_type = walltm

    def __backwardCompatibilityFrom_7_0(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        # Check that all zones have the 'physical_properties' flag
        XMLVolumicNode = self.case.xmlGetNode('volumic_conditions')
        for node in XMLVolumicNode.xmlGetChildNodeList('zone', 'label', 'id'):
            if not node['physical_properties']:
                # The 'all_cells' zone is always used for physical properties
                # for the v7.0
                if node['label'] == "all_cells":
                    node['physical_properties'] = 'on'
                else:
                    node['physical_properties'] = 'off'


    def __backwardCompatibilityFrom_7_1(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        # Update pressure boundary condition for outlets
        XMLBCNode = self.case.xmlGetNode('boundary_conditions')
        for node in XMLBCNode.xmlGetChildNodeList('outlet', field_id='none'):
            p_node = node.xmlGetChildNode('dirichlet', name='pressure')
            if p_node:
                # Get previous value of pressure
                p_val = node.xmlGetDouble('dirichlet', name='pressure')

                # Remove old node
                p_node.xmlRemoveNode()

                # Create new node witch correct choice value. Note that previous
                # versions indicated 'dirichlet' while the bnd condition used was
                # the dpdndtau (mean pressure value)
                p_node = node.xmlInitChildNode('variable', choice='dpdndtau', name='pressure')
                p_node.xmlSetData('value', str(p_val))


    def _backwardCompatibilityCurrentVersion(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        self.__backwardCompatibilityFrom_7_1()
#-------------------------------------------------------------------------------
# XMLinit test case
#-------------------------------------------------------------------------------

class XMLinitTestCaseNeptune(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from code_saturne.gui import XMLengine
        GuiParam.lang = 'en'
        self.doc = XMLengine.XMLDocument("")
        self.case = XMLengine.Case(None)


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case
        del self.doc


    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()


    def checkXMLinitInstantiation(self):
        """
        Check whether the Case class could be instantiated
        """
        xmldoc = None
        xmldoc = XMLinit(self.case)
        assert xmldoc != None, 'Could not instantiate XMLinit'


    def checkInitHeading(self):
        """
        Check whether the headings markups could be initialized
        """
        XMLinit(self.case)

        assert self.case.root() == self.xmlNodeFromString(doc), \
               'Could not use the constructor of the XMLinit class'


def suite():
    testSuite = unittest.makeSuite(XMLinitTestCase, "check")
    return testSuite


def runTest():
    print("XMLinitTestCase to be completed...")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of XMLinit
#-------------------------------------------------------------------------------
