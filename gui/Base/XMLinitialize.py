# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.XMLvariables import Variables
from code_saturne.Base import Toolbox

from code_saturne.Pages.LocalizationModel import Zone, LocalizationModel
from code_saturne.Pages.OutputControlModel import OutputControlModel
from code_saturne.Pages.MobileMeshModel import MobileMeshModel
from code_saturne.Pages.TurbulenceModel import TurbulenceModel
from code_saturne.Pages.InitializationModel import InitializationModel
from code_saturne.Pages.TimeStepModel import TimeStepModel
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.ElectricalModel import ElectricalModel
from code_saturne.Pages.GasCombustionModel import GasCombustionModel
from code_saturne.Pages.GroundwaterModel import GroundwaterModel
from code_saturne.Pages.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.Pages.LagrangianModel import LagrangianModel
from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel

#-------------------------------------------------------------------------------
# class XMLinit
#-------------------------------------------------------------------------------

class XMLinit(Variables):
    """
    This class initializes the XML contents of the case.
    """
    def __init__(self, case):
        """
        """
        self.case = case


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
            self.__backwardCompatibility()

            # Initialization (order is important)

            grdflow = GroundwaterModel(self.case).getGroundwaterModel()

            self.node_models = self.case.xmlInitNode('thermophysical_models')
            node = self.node_models.xmlInitNode('velocity_pressure')
            if grdflow == 'groundwater':
                self.setNewVariable(node, 'hydraulic_head')
            else:
                self.setNewVariable(node, 'pressure')
            self.setNewVariable(node, 'velocity', dim = '3')
            self.setNewProperty(node, 'total_pressure')

            if grdflow != 'groundwater':
                n = self.setNewProperty(node, 'yplus')
                n['support'] = 'boundary'
                n['label'] = 'Yplus'
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

            MobileMeshModel(self.case).getMethod()
            TurbulenceModel(self.case).getTurbulenceModel()

            # First Volume Zone definition for all cells -> initialization

            zones = LocalizationModel("VolumicZone", self.case).getZones()
            iok = 0
            for zone in zones:
                if zone.getLabel() == 'all_cells':
                    iok = 1
            if iok == 0:
                zone = Zone("VolumicZone", case=self.case, label='all_cells', localization='all[]')
                LocalizationModel("VolumicZone", self.case).addZone(zone)
                zone = LocalizationModel("VolumicZone", self.case).getCodeNumberOfZoneLabel('all_cells')
                InitializationModel(self.case).getInitialTurbulenceChoice(zone)

            # Time settings

            TimeStepModel(self.case).getTimeStep()
            TimeStepModel(self.case).getTimePassing()

            # Thermodynamics definitinon

            m = FluidCharacteristicsModel(self.case)
            for tag in ('density',
                        'molecular_viscosity',
                        'specific_heat',
                        'thermal_conductivity'):
                m.getInitialValue(tag)

            # Calculation features

            ThermalScalarModel(self.case).getThermalScalarModel()
            CoalCombustionModel(self.case).getCoalCombustionModel()
            GasCombustionModel(self.case).getGasCombustionModel()
            ElectricalModel(self.case).getElectricalModel()
            ThermalRadiationModel(self.case).getRadiativeModel()
            AtmosphericFlowsModel(self.case).getAtmosphericFlowsModel()
            LagrangianModel(self.case).getLagrangianModel()

        return msg


    def __initHeading(self, prepro):
        """
        Create if necessary headings from the root element of the case.
        """
        msg = ""
        tagList = ('solution_domain',
                   'analysis_control',
                   'calculation_management')
        if not prepro:
            tagList += ('thermophysical_models',
                       'additional_scalars',
                       'physical_properties',
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


    def __renameSingle(self, parent_tag, old_tag, new_tag):
        """
        Rename some nodes in order to ensure backward compatibility.
        """

        # renames

        node = self.case.xmlGetNode(parent_tag)
        if node:
            oldnode = node.xmlGetNode(old_tag)
            if oldnode:
                newnode = node.xmlInitNode(new_tag)
                newnode.xmlChildsCopy(oldnode)
                oldnode.xmlRemoveNode()


    def __clean_version(self, vers):
        """
        Simplify version history number, replacing "alpha" or "beta"
        version with previous version  to force ensuring of backward
        compatibility.
        """
        known_versions = ["3.0", "3.1", "3.2", "3.3",
                          "4.0", "4.1", "4.2", "4.3",
                          "5.0", "5.1", "5.2", "5.3",
                          "6.0"]
        j = -2
        for i in range(0, len(known_versions)):
            if vers.find(known_versions[i]) == 0:
                j = i
                for e in ("-alpha", "-beta"):
                    if vers.find(e) > -1:
                        j = i-1
                break
        if j == -1:
            if not vers:
                vers = "-1.0"
        elif j > 0:
            vers = known_versions[j]
        return vers


    def __backwardCompatibility(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        cur_vers = self.case['package'].version

        if self.case.root()["solver_version"]:
            his_r = self.case.root()["solver_version"]
            history = his_r.split(";")
            last_vers = self.__clean_version(history[len(history) - 1])
            if last_vers == cur_vers:
                self.__backwardCompatibilityCurrentVersion()
            else:
                self.__backwardCompatibilityOldVersion(last_vers)
                self.__backwardCompatibilityCurrentVersion()
            his = ""
            vp = ""
            for v in history:
                vc = self.__clean_version(v)
                if vc != vp:
                    his += vc + ";"
                    vp = vc
            if cur_vers != vp:
                his += cur_vers + ";"
            his = his[:-1]
            if his != his_r:
                self.case.root().xmlSetAttribute(solver_version = his)

        else:
            vers = cur_vers
            self.case.root().xmlSetAttribute(solver_version = vers)

            # apply all backwardCompatibilities as we don't know
            # when it was created
            self.__backwardCompatibilityOldVersion("-1")
            self.__backwardCompatibilityCurrentVersion()


    def __backwardCompatibilityOldVersion(self, from_vers):
        """
        Change XML in order to ensure backward compatibility for old version
        """
        if from_vers <= "-1.0":
            self.__backwardCompatibilityBefore_3_0()

        if from_vers[:3] < "4.0.0":
            if from_vers[:3] < "3.1.0":
                self.__backwardCompatibilityFrom_3_0()
            if from_vers[:3] < "3.2.0":
                self.__backwardCompatibilityFrom_3_1()
            if from_vers[:3] < "3.3.0":
                self.__backwardCompatibilityFrom_3_2()
            if from_vers[:3] < "4.0.0":
                self.__backwardCompatibilityFrom_3_3()

        if from_vers[:3] < "5.0.0":
            if from_vers[:3] < "4.1.0":
                self.__backwardCompatibilityFrom_4_0()
            if from_vers[:3] < "4.2.0":
                self.__backwardCompatibilityFrom_4_1()
            if from_vers[:3] < "4.3.0":
                self.__backwardCompatibilityFrom_4_2()
            if from_vers[:3] < "5.0.0":
                self.__backwardCompatibilityFrom_4_3()

        if from_vers[:3] < "6.0.0":
            if from_vers[:3] < "5.1.0":
                self.__backwardCompatibilityFrom_5_0()
            if from_vers[:3] < "5.2.0":
                self.__backwardCompatibilityFrom_5_1()
            if from_vers[:3] < "5.3.0":
                self.__backwardCompatibilityFrom_5_2()
            if from_vers[:3] < "6.0.0":
                self.__backwardCompatibilityFrom_5_3()


    def __backwardCompatibilityBefore_3_0(self):
        """
        Change XML in order to ensure backward compatibility from 2.x to 3.0
        """
        for node in self.case.xmlGetNodeList('initial_value', 'zone'):
            node['zone_id'] = node['zone']

        for varNode in self.case.xmlGetNodeList('variable'):
            value = varNode.xmlGetDouble('solveur_precision')
            if value:
                varNode.xmlSetData('solver_precision', value)
                varNode.xmlRemoveChild('solveur_precision')

        XMLSolutionDomainNode = self.case.xmlInitNode('solution_domain')
        self.__XMLVolumicConditionsNode = XMLSolutionDomainNode.xmlInitNode('volumic_conditions')
        for node in self.__XMLVolumicConditionsNode.xmlGetNodeList('zone'):
            if node['id'] == None:
                node['id'] = node['name']

        oldnode = self.case.xmlGetNode('calcul_management')
        if oldnode:
            newnode = self.case.xmlInitNode('calculation_management')
            newnode.xmlChildsCopy(oldnode)
            oldnode.xmlRemoveNode()

        # Reference values
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        self.__XMLVelocityPressureNode = XMLThermoPhysicalNode.xmlInitNode('velocity_pressure')
        self.__RefValuesNode = XMLThermoPhysicalNode.xmlInitNode('reference_values')

        nodeP = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="pressure")
        if nodeP:
            value = nodeP.xmlGetDouble('reference_pressure')
            if value:
                self.__RefValuesNode.xmlSetData('pressure', value)
                nodeP.xmlRemoveChild('reference_pressure')

        nodeTurb = XMLThermoPhysicalNode.xmlInitNode('turbulence', 'model')

        for nodeInit in nodeTurb.xmlGetNodeList('initialization'):
            if nodeInit:
                value = nodeInit.xmlGetDouble('reference_velocity')
                if value:
                    self.__RefValuesNode.xmlSetData('velocity', value)
                    nodeInit.xmlRemoveChild('reference_velocity')

                value = nodeInit.xmlGetDouble('reference_length')
                if value:
                    self.__RefValuesNode.xmlSetData('length', value)
                    nodeInit.xmlRemoveChild('reference_length')

        for node in self.case.xmlGetNodeList('scalar'):
            value = node.xmlGetDouble('initial_value', zone_id="1")
            if value != None:
                formula = node['label'] + " = " + str(value) + ";"
                n = node.xmlInitChildNode('formula', zone_id="1")
                n.xmlSetTextNode(formula)
                node.xmlRemoveChild('initial_value', zone_id="1")

        # solver
        XMLNumParameterNode = self.case.xmlInitNode('numerical_parameters')
        node = XMLNumParameterNode.xmlGetNode('multigrid')
        if node:
            if node['status'] == "off":
                if nodeP:
                    nodeP.xmlInitNode('solver_choice', choice='conjugate_gradient')
            node.xmlRemoveNode()

        # hydrostatic pressure
        XMLPhysicalPropNode = self.case.xmlInitNode('physical_properties')
        node = XMLPhysicalPropNode.xmlGetNode('hydrostatic_pressure')
        if node:
            stat = node['status']
            XMLNumParameterNode.xmlInitNode('hydrostatic_pressure', status=stat)
            node.xmlRemoveNode()

        # Profiles
        compt = 0
        for node in self.case.xmlGetNodeList('profile'):
            nodeInit = node.xmlGetNode('x1')
            if nodeInit:
                node.xmlRemoveNode()
                compt = compt + 1
        if compt != 0:
            print("Profiles have been removed from your files due to  incompatibility")
            print("You must re-create them")

        # restart
        nodeR = self.case.xmlGetNode("start_restart")
        if nodeR:
            n = nodeR.xmlGetNode("restart", "status")
            if n:
                n.xmlRemoveNode()


    def __backwardCompatibilityFrom_3_0(self):
        """
        Change XML in order to ensure backward compatibility from 3.0 to 3.1
        """
        # Profiles
        for node in self.case.xmlGetNodeList('profile'):
            if node:
                n = node.xmlGetNode("output_type")
                if n == None:
                    freq = node.xmlGetInt("output_frequency")
                    if freq == -1:
                        node.xmlSetData('output_type', "end")
                    else:
                        node.xmlSetData('output_type', "frequency")


    def __backwardCompatibilityFrom_3_1(self):
        """
        Change XML in order to ensure backward compatibility from 3.1 to 3.2
        """
        # thermal scalar
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        for phys in ['solid_fuels', 'gas_combustion', 'joule_effect', 'atmospheric_flows']:
            node = XMLThermoPhysicalNode.xmlInitNode(phys, 'model')
            mdl = node['model']
            if mdl and mdl != 'off':
                if phys != 'atmospheric_flows':
                    n = node.xmlGetNode('scalar', name="Enthalpy")
                    if n:
                        n.xmlRemoveNode()
                    ThermalScalarModel(self.case).setThermalModel('enthalpy')
                else:
                    if mdl == "dry":
                        n = node.xmlGetNode('scalar', name="potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('potential_temperature')
                    elif mdl == "constant":
                        n = node.xmlGetNode('scalar', name="potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        n = node.xmlGetNode('scalar', name="liquid_potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('off')
                    else:
                        n = node.xmlGetNode('scalar', name="liquid_potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('liquid_potential_temperature')
        node = self.case.xmlGetNode('additional_scalars')
        n = node.xmlGetNode('scalar', type='thermal')
        if n:
            nth = XMLThermoPhysicalNode.xmlGetNode('thermal_scalar')
            nthvar = nth.xmlInitNode('scalar', 'type')
            nthvar['type']  = "thermal"
            nthvar['name']  = n['name']
            nthvar['label'] = n['label']
            nthvar.xmlChildsCopy(n)
            n.xmlRemoveNode()


    def __backwardCompatibilityFrom_3_2(self):
        """
        Change XML in order to ensure backward compatibility from 3.2 to 3.3
        """
        # thermal scalar
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        for phys in ['solid_fuels', 'gas_combustion', 'joule_effect', 'atmospheric_flows', 'compressible_model']:
            node = XMLThermoPhysicalNode.xmlInitNode(phys, 'model')
            mdl = node['model']
            if mdl and mdl != 'off':
                if phys != 'atmospheric_flows' and phys != 'compressible_model':
                    n = node.xmlGetNode('scalar', name="Enthalpy")
                    if n:
                        n.xmlRemoveNode()
                    ThermalScalarModel(self.case).setThermalModel('enthalpy')
                elif phys == 'atmospheric_flows':
                    if (mdl == "dry"):
                        n = node.xmlGetNode('scalar', name="potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('potential_temperature')
                    if (mdl == "humid"):
                        n = node.xmlGetNode('scalar', name="liquid_potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('liquid_potential_temperature')
                else:
                    n = node.xmlGetNode('scalar', name="EnergieT")
                    if n:
                        n.xmlRemoveNode()
                    ThermalScalarModel(self.case).setThermalModel('total_energy')

        # properties
        XMLPhysicalPropNode = self.case.xmlInitNode('physical_properties')
        nodeF = XMLPhysicalPropNode.xmlInitNode('fluid_properties')
        for prop in ['density', 'molecular_viscosity', 'specific_heat',
                     'thermal_conductivity', 'volume_viscosity']:
            node = nodeF.xmlGetNode('property', name=prop)
            if node:
                if node['choice'] == 'user_law':
                    node['choice'] = 'variable'

        node = self.case.xmlGetNode('additional_scalars')
        n = node.xmlGetNode('scalar', type='thermal')
        if n:
            nth = XMLThermoPhysicalNode.xmlGetNode('thermal_scalar')
            nthvar = nth.xmlInitNode('variable', 'type')
            nthvar['type']  = "thermal"
            nthvar['name']  = n['name']
            nthvar['label'] = n['label']
            nthvar.xmlChildsCopy(n)
            n.xmlRemoveNode()

        # Replace scalar by variable in xml
        for phys in ['solid_fuels', 'gas_combustion', 'joule_effect', 'atmospheric_flows', 'compressible_model', 'thermal_scalar']:
            nodeP = XMLThermoPhysicalNode.xmlInitNode(phys, 'model')
            for node in nodeP.xmlGetNodeList('scalar'):
                name = node['name']
                label = node['label']
                dim = node['dimension']
                tpe = node['type']
                newnode = nodeP.xmlInitNode('variable', name=name)
                if label != None:
                    newnode['label'] = label
                if dim != None:
                    newnode['dimension'] = dim
                if tpe != None:
                    newnode['type'] = tpe
                newnode.xmlChildsCopy(node)
                node.xmlRemoveNode()

        self.scalar_node = self.case.xmlGetNode('additional_scalars')
        for node in self.scalar_node.xmlGetNodeList('scalar'):
            name = node['name']
            label = node['label']
            dim = node['dimension']
            tpe = node['type']
            newnode = self.scalar_node.xmlInitNode('variable', name=name)
            if label != None:
                newnode['label'] = label
            if dim != None:
                newnode['dimension'] = dim
            if tpe != None:
                newnode['type'] = tpe
            newnode.xmlChildsCopy(node)
            node.xmlRemoveNode()

        n = XMLThermoPhysicalNode.xmlGetNode('variable', type='thermal')
        if n:
            for nf in n.xmlGetNodeList('formula'):
                if nf:
                    status = nf["status"]
                    if not(status) or status == "on":
                        content = nf.xmlGetTextNode()
                        # Substitute only perfectly matching labels
                        pattern = '\\b' + n['label'] + '\\b'
                        content = re.sub(pattern, n['name'], content)
                        nf.xmlSetTextNode(content)

        # update velocity node
        nodeV = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_U")
        if nodeV:
            nodeV['name'] = 'velocity'
            nodeV['dimension'] = '3'

            nodeTmp = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_V")
            if nodeTmp:
                nodeTmp.xmlRemoveNode()
            nodeTmp = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_W")
            if nodeTmp:
                nodeTmp.xmlRemoveNode()
        for node in self.case.xmlGetNodeList('profile'):
            if node:
                for n in node.xmlGetNodeList('var_prop'):
                    name = n['name']
                    if name in ["velocity_U", "velocity_V", "velocity_W"]:
                        if name == 'velocity_U':
                            component = '0'
                        if name == 'velocity_V':
                            component = '1'
                        if name == 'velocity_W':
                            component = '2'
                        name = 'velocity'
                        n['name'] = name
                        n['component'] = component
                    elif name != "velocity":
                        n['component'] = "0"

        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                for n in node.xmlGetNodeList('var_prop'):
                    name = n['name']
                    if name in ["velocity_U", "velocity_V", "velocity_W"]:
                        if name == 'velocity_U':
                            component = '0'
                        if name == 'velocity_V':
                            component = '1'
                        if name == 'velocity_W':
                            component = '2'
                        name = 'velocity'
                        n['name'] = name
                        n['component'] = component
                    elif name != "velocity":
                        n['component'] = "0"

        for node in self.case.xmlGetNodeList('dirichlet'):
            if node:
                name = node['name']
                if name in ["velocity_U", "velocity_V", "velocity_W"]:
                    if name == 'velocity_U':
                        component = '0'
                    if name == 'velocity_V':
                        component = '1'
                    if name == 'velocity_W':
                        component = '2'
                    name = 'velocity'
                    node['name'] = name
                    node['component'] = component

        dicoName = [("NP_CP",                        "n_p_"),
                    ("XCH_CP",                       "x_p_coal_"),
                    ("XCK_CP",                       "x_p_char_"),
                    ("ENT_CP",                       "x_p_h_"),
                    ("XWT_CP",                       "x_p_wt_"),
                    ("Fr_MV1",                       "fr_mv1_"),
                    ("Fr_MV2",                       "fr_mv2_"),
                    ("Fr_HET_O2",                    "fr_het_o2"),
                    ("Fr_HET_CO2",                   "fr_het_co2"),
                    ("Fr_HET_H2O",                   "fr_het_h2o"),
                    ("FR_HCN",                       "x_c_hcn"),
                    ("FR_NO",                        "x_c_no"),
                    ("FR_NH3",                       "x_c_nh3"),
                    ("FR_CO2",                       "x_c_co2"),
                    ("Enth_Ox",                      "x_c_h_ox"),
                    ("FR_H20",                       "fr_h2o"),
                    ("FR_OXYD2",                     "fr_oxyd2"),
                    ("FR_OXYD3",                     "fr_oxyd3"),
                    ("Var_F1F2",                     "f1f2_variance"),
                    ("scalar",                       "user_"),
                    ("PotElecReal",                  "elec_pot_r"),
                    ("POT_EL_I",                     "elec_pot_i"),
                    ("YM_ESL",                       "esl_fraction_"),
                    ("POT_VEC",                      "vec_potential_"),
                    ("Fra_MEL",                      "mixture_fraction"),
                    ("Var_FMe",                      "mixture_fraction_variance"),
                    ("Fra_GF",                       "fresh_gas_fraction"),
                    ("Fra_Mas",                      "mass_fraction"),
                    ("COYF_PP4",                     "mass_fraction_covariance"),
                    ("Var_FMa",                      "mass_fraction_variance"),
                    ("temperature_celsius",          "temperature"),
                    ("temperature_kelvin",           "temperature"),
                    ("TempK",                        "temperature"),
                    ("potential_temperature",        "temperature"),
                    ("liquid_potential_temperature", "temperature"),
                    ("component_R11",                "r11"),
                    ("component_R22",                "r22"),
                    ("component_R33",                "r33"),
                    ("component_R12",                "r12"),
                    ("component_R13",                "r13"),
                    ("component_R23",                "r23"),
                    ("StressTensor",                 "rij"),
                    ("turb_k",                       "k"),
                    ("turb_eps",                     "epsilon"),
                    ("turb_phi",                     "phi"),
                    ("turb_alpha",                   "alpha"),
                    ("turb_omega",                   "omega"),
                    ("nusa",                         "nu_tilda"),
                    ("volumic_viscosity",            "volume_viscosity")]
        dico = {}
        for (u,v) in dicoName:
            dico[u] = v
        for node in self.case.xmlGetNodeList('variable'):
            name = node["name"]
            if name:
                for key in dico.keys():
                    if name.startswith(key):
                        idx = name.find(key) + len(key)
                        node["name"] = dico[key] + name[idx:]
                        break

        XMLBoundaryNode = self.case.xmlInitNode('boundary_conditions')
        for node in XMLBoundaryNode.xmlGetNodeList('scalar'):
            name = node["name"]
            if name:
                for key in dico.keys():
                    if name.startswith(key):
                        idx = name.find(key) + len(key)
                        node["name"] = dico[key] + name[idx:]
                        break

        for node in self.case.xmlGetNodeList('var_prop'):
            name = node["name"]
            if name:
                for key in dico.keys():
                    if name.startswith(key):
                        idx = name.find(key) + len(key)
                        node["name"] = dico[key] + name[idx:]
                        break


        # update formula
        nth = XMLThermoPhysicalNode.xmlGetNode('thermal_scalar')
        nvel = XMLThermoPhysicalNode.xmlGetNode('velocity_pressure')

        for node in nvel.xmlGetNodeList('formula'):
            status = node["status"]
            if not(status) or status == "on":
                content = node.xmlGetTextNode()
                # Substitute only perfectly matching labels
                pattern = '\\bu\\b'
                content = re.sub(pattern, 'velocity[0]', content)
                pattern = '\\bv\\b'
                content = re.sub(pattern, 'velocity[1]', content)
                pattern = '\\bw\\b'
                content = re.sub(pattern, 'velocity[2]', content)
                pattern = '\\bP\\b'
                content = re.sub(pattern, 'pressure', content)
                node.xmlSetTextNode(content)

        for node in nth.xmlGetNodeList('formula'):
            status = node["status"]
            if not(status) or status == "on":
                content = node.xmlGetTextNode()
                # Substitute only perfectly matching labels
                pattern = '\\bT\\b'
                content = re.sub(pattern, 'temperature', content)
                pattern = '\\btemperature_celsius\\b'
                content = re.sub(pattern, 'temperature', content)
                pattern = '\\btemperature_kelvin\\b'
                content = re.sub(pattern, 'temperature', content)
                node.xmlSetTextNode(content)

        for node in XMLThermoPhysicalNode.xmlGetNodeList('formula'):
            status = node["status"]
            if not(status) or status == "on":
                content = node.xmlGetTextNode()
                nodeas = self.case.xmlGetNode('additional_scalars')
                nth = nodeas.xmlGetNode('scalar', type='thermal')
                if nth:
                    # Substitute only perfectly matching labels
                    pattern = '\\b' + nth['label'] + '\\b'
                    content = re.sub(pattern, nth['name'], content)
                node.xmlSetTextNode(content)

        for node in XMLPhysicalPropNode.xmlGetNodeList('formula'):
            nodeas = self.case.xmlGetNode('additional_scalars')
            nth = nodeas.xmlGetNode('scalar', type='thermal')
            if nth:
                content = node.xmlGetTextNode()
                # Substitute only perfectly matching labels
                pattern = '\\b' + nth['label'] + '\\b'
                content = re.sub(pattern, nth['name'], content)
                node.xmlSetTextNode(content)

        XMLAddScalar = self.case.xmlGetNode('additional_scalars')
        for node in XMLAddScalar.xmlGetNodeList('variable'):
            nfor = node.xmlGetNode('formula')
            if nfor:
                content = nfor.xmlGetTextNode()
                # Substitute only perfectly matching labels
                pattern = '\\b' + node['label'] + '\\b'
                content = re.sub(pattern, node['name'], content)
                nfor.xmlSetTextNode(content)

        for node in XMLBoundaryNode.xmlGetNodeList('turbulence'):
            if node["choice"] == "formula":
                nf = node.xmlGetNode('formula')
                if nf:
                    content = nf.xmlGetTextNode()
                    # Substitute only perfectly matching labels
                    pattern = '\\beps\\b'
                    content = re.sub(pattern, 'epsilon', content)
                    nf.xmlSetTextNode(content)

        # TODO update formula BC for turbulence
        #for node in XMLBoundaryNode.xmlGetNodeList('turbulence'):
        #    if node["choice"] = "formula":
        #        nf = node.xmlGetNode('formula')
        #        for key in dico.keys():
        #            if name.startswith(key):
        #                idx = name.find(key) + len(key)
        #                node["name"] = dico[key] + name[idx:]
        #                break

        dicoProp = [("Rho",                          "density"),
                    ("turb_viscosity",               "turbulent_viscosity"),
                    ("smagorinsky_constant",         "smagorinsky_constant^2"),
                    ("Temperature",                  "temperature"),
                    ("YM_Fuel",                      "ym_fuel"),
                    ("YM_Oxyd",                      "ym_oxyd"),
                    ("YM_Prod",                      "ym_prod"),
                    ("Mas_Mol",                      "molar_mass"),
                    ("T.SOURCE",                     "source_term"),
                    ("RHOL0",                        "rho_local_"),
                    ("TEML0",                        "temperature_local_"),
                    ("FMEL0",                        "ym_local_"),
                    ("FMAL0",                        "w_local_"),
                    ("AMPL0",                        "amplitude_local_"),
                    ("TSCL0",                        "chemical_st_local_"),
                    ("MAML0",                        "molar_mass_local_"),
                    ("Temp_GAZ",                     "t_gas"),
                    ("ROM_GAZ",                      "rho_gas"),
                    ("YM_CHx1m",                     "ym_chx1m"),
                    ("YM_CHx2m",                     "ym_chx2m"),
                    ("YM_CO",                        "ym_co"),
                    ("YM_H2S",                       "ym_h2s"),
                    ("YM_H2",                        "ym_h2"),
                    ("YM_HCN",                       "ym_hcn"),
                    ("YM_NH3",                       "ym_nh3"),
                    ("YM_O2",                        "ym_o2"),
                    ("YM_CO2",                       "ym_co2"),
                    ("YM_H2O",                       "ym_h2o"),
                    ("YM_SO2",                       "ym_so2"),
                    ("YM_N2",                        "ym_n2"),
                    ("XM",                           "xm"),
                    ("EXP1",                         "exp1"),
                    ("EXP2",                         "exp2"),
                    ("EXP3",                         "exp3"),
                    ("EXP4",                         "exp4"),
                    ("EXP5",                         "exp5"),
                    ("F_HCN_DEV",                    "f_hcn_dev"),
                    ("F_HCN_HET",                    "f_hcn_het"),
                    ("F_NH3_DEV",                    "f_nh3_dev"),
                    ("F_NH3_HET",                    "f_nh3_het"),
                    ("F_NO_HCN",                     "f_no_hcn"),
                    ("F_NO_NH3",                     "f_no_nh3"),
                    ("F_NO_HET",                     "f_no_het"),
                    ("F_NO_THE",                     "f_no_the"),
                    ("C_NO_HCN",                     "c_no_hcn"),
                    ("C_NO_NH3",                     "c_no_nh3"),
                    ("F_HCN_RB",                     "f_hcn_rb"),
                    ("C_NO_RB",                      "c_no_rb"),
                    ("EXP_RB",                       "exp_rb"),
                    ("Temp_CP",                      "t_p_"),
                    ("Frm_CP",                       "x_p_"),
                    ("Rho_CP",                       "rho_p_"),
                    ("Dia_CK",                       "diam_p_"),
                    ("Ga_DCH",                       "dissapear_rate_p_"),
                    ("Ga_DV1",                       "m_transfer_v1_p_"),
                    ("Ga_DV2",                       "m_transfer_v2_p_"),
                    ("Ga_HET_O2",                    "het_ts_o2_p_"),
                    ("Ga_HET_CO2",                   "het_ts_co2_p"),
                    ("Ga_HET_H2O",                   "het_ts_h2o_p"),
                    ("Ga_HET",                       "het_ts_coal"),
                    ("Ga_SEC",                       "dry_ts_p"),
                    ("Bilan_C",                      "x_carbone"),
                    ("Bilan_O",                      "x_oxygen"),
                    ("Bilan_H",                      "x_hydrogen"),
                    ("PuisJoul",                     "joule_power"),
                    ("Cour_re",                      "current_re"),
                    ("CouImag",                      "current_im"),
                    ("For_Lap",                      "laplace_force"),
                    ("Coef_Abso",                    "absorption_coeff"),
                    ("c_NO_HCN",                     "radiation_source"),
                    ("Sigma",                        "elec_sigma"),
                    ("IntLuminance_4PI",             "intensity"),
                    ("volumic_viscosity",            "volume_viscosity")]

        dicoP = {}
        for (u,v) in dicoProp:
            dicoP[u] = v
        for node in self.case.xmlGetNodeList('property'):
            name = node["name"]
            if name:
                for key in dicoP.keys():
                    if name.startswith(key) and name != "smagorinsky_constant^2":
                        idx = name.find(key) + len(key)
                        node["name"] = dicoP[key] + name[idx:]
                        break

        for node in self.case.xmlGetNodeList('var_prop'):
            name = node["name"]
            if name:
                for key in dicoP.keys():
                    if name.startswith(key) and name != "smagorinsky_constant^2":
                        idx = name.find(key) + len(key)
                        node["name"] = dicoP[key] + name[idx:]
                        break

        nodeCompress = XMLThermoPhysicalNode.xmlGetNode('compressible_model')
        if nodeCompress['model'] and nodeCompress['model'] != "off":
            n = nodeCompress.xmlGetNode("property", name = "density")
            if n:
                ndens = XMLPhysicalPropNode.xmlGetNode('property', name='density')
                ndens.xmlChildsCopy(n)
                n.xmlRemoveNode()
        for node in XMLPhysicalPropNode.xmlGetNodeList('property'):
            n = node.xmlGetNode('formula')
            if n:
                f = n.xmlGetTextNode()
                if f != None:
                    # Substitute only perfectly matching labels
                    pattern = '\\brho\\b'
                    content = re.sub(pattern, 'density', content)
                    pattern = '\\bmu\\b'
                    content = re.sub(pattern, 'molecular_viscosity', content)
                    pattern = '\\bcp\\b'
                    content = re.sub(pattern, 'specific_heat', content)
                    pattern = '\\blambda\\b'
                    content = re.sub(pattern, 'thermal_conductivity', content)
                    pattern = '\\bviscv\\b'
                    content = re.sub(pattern, 'volume_viscosity', content)
                    n.xmlSetTextNode(f)


    def __backwardCompatibilityFrom_3_3(self):
        """
        Change XML in order to ensure backward compatibility from 3.3 to 4.0
        """
        XMLAnaControl = self.case.xmlGetNode('analysis_control')
        self.scalar_node = self.case.xmlGetNode('additional_scalars')
        for node in self.scalar_node.xmlGetNodeList('variable'):
            name = node['name']
            label = node['label']
            if name == None:
                node['name'] = label
            if label == None:
                node['label'] = name
            for n in XMLAnaControl.xmlGetNodeList('var_prop'):
                if n['name'] == name:
                    n['name'] = node['name']
            for n in node.xmlGetNodeList('formula'):
                if n:
                    content = n.xmlGetTextNode()
                    # Substitute only perfectly matching labels
                    pattern = '\\b' + name + '\\b'
                    content = re.sub(pattern, name, content)
                    n.xmlSetTextNode(content)

        XMLBoundaryNode = self.case.xmlInitNode('boundary_conditions')
        for node in XMLBoundaryNode.xmlGetNodeList('scalar'):
            name = node['name']
            label = node['label']
            if name == None:
                node['name'] = label
            if label == None:
                node['label'] = name

        XMLThermoPhysicalModel = self.case.xmlGetNode('thermophysical_models')
        XMLAleMethod = XMLThermoPhysicalModel.xmlInitChildNode('ale_method')
        if XMLAleMethod:
            for node in XMLAleMethod.xmlGetNodeList('formula'):
                if node:
                    content = node.xmlGetTextNode()
                    # Substitute only perfectly matching labels
                    pattern = '\\bmesh_vi1\\b'
                    content = re.sub(pattern, "mesh_viscosity_1", content)
                    pattern = '\\bmesh_vi2\\b'
                    content = re.sub(pattern, "mesh_viscosity_2", content)
                    pattern = '\\bmesh_vi3\\b'
                    content = re.sub(pattern, "mesh_viscosity_3", content)
                    node.xmlSetTextNode(content)

        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                time_node = node.xmlGetNode("time_start")
                if not time_node:
                    node.xmlSetData('time_start', -1.)

        # update velocity node
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        self.__XMLVelocityPressureNode = XMLThermoPhysicalNode.xmlInitNode('velocity_pressure')
        nodeV = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity")
        if nodeV:
            if nodeV['label'] == 'VelocityX':
                nodeV['label'] = 'Velocity'
        #update input_thermal_flux
        nth = XMLThermoPhysicalNode.xmlGetNode('thermal_scalar')
        if nth:
            node = nth.xmlGetNode('property', name="input_thermal_flux")
            node2 = nth.xmlGetNode('property', name="thermal_flux")
            if node2 and node:
                node.xmlRemoveNode()
            elif node:
                node['name'] = "thermal_flux"


    def __backwardCompatibilityFrom_4_0(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        XMLThermoPhysicalModelNode = self.case.xmlGetNode('thermophysical_models')
        n = XMLThermoPhysicalModelNode.xmlGetNode('variable', type='thermal')
        if n:
            # try to get turbulent_flux_model
            ntfm = n.xmlGetString("turbulent_flux_model")
            if not ntfm:
                n.xmlSetData('turbulent_flux_model', "SGDH")

        # replace label by name in each formula
        for node in self.case.xmlGetNodeList('formula'):
            if node:
                status = node["status"]
                if not(status) or status == "on":
                    content = node.xmlGetTextNode()
                    for n in self.case.xmlGetNodeList('variable', 'name', 'label'):
                        # Substitute only perfectly matching labels
                        pattern = '\\b' + n['label'] + '\\b'
                        content = re.sub(pattern, n['name'], content)
                    node.xmlSetTextNode(content)

        for node in self.case.xmlGetNodeList('variable'):
            variance = node.xmlGetString('variance')
            if variance:
                for n in self.case.xmlGetNodeList('variable', 'name', 'label'):
                    if variance == n['label']:
                        node.xmlSetData('variance', n['name'])
                        break

        # update mesh velocity node
        # Note: it is important to do this only after updating formulas and
        #       not before, to apply updates in the order if code changes.
        XMLThermoPhysicalModel = self.case.xmlGetNode('thermophysical_models')
        XMLAleMethod = XMLThermoPhysicalModel.xmlGetChildNode('ale_method')
        if XMLAleMethod:
            nodeV = XMLAleMethod.xmlGetNode('variable', name="mesh_velocity_U")
            if nodeV:
                nodeV['name'] = 'mesh_velocity'
                nodeV['label'] = 'Mesh Velocity'
                nodeV['dimension'] = '3'

                nodeTmp = XMLAleMethod.xmlGetNode('variable', name="mesh_velocity_V")
                if nodeTmp:
                    nodeTmp.xmlRemoveNode()
                nodeTmp = XMLAleMethod.xmlGetNode('variable', name="mesh_velocity_W")
                if nodeTmp:
                    nodeTmp.xmlRemoveNode()

        # replace bounce by part_symmetry for lagrangian model on
        # symmetry
        XMLBoundaryNode = self.case.xmlInitNode('boundary_conditions')
        for node in XMLBoundaryNode.xmlGetNodeList('symmetry'):
            nn = node.xmlGetNode("particles")
            if nn:
                if nn["choice"] == "bounce":
                    nn["choice"] = "part_symmetry"

        # add lagrangian writer if needed
        XMLLagrangianModel = self.case.xmlGetNode('lagrangian')
        if XMLLagrangianModel:
            mdl = XMLLagrangianModel["model"]
            if mdl != "off":
                XMLAnaControl = self.case.xmlGetNode('analysis_control')
                node_out = XMLAnaControl.xmlGetNode('output')
                nn = node_out.xmlGetNode('writer', 'label', id = "-3")
                if nn == None:
                    nodeL = node_out.xmlInitNode('writer', id = "-3", label = 'particles')
                    nodeL.xmlInitNode('frequency', period = 'none')
                    nodeL.xmlInitNode('output_at_end', status = 'on')
                    nodeL.xmlInitNode('format', name = 'ensight', options = 'binary')
                    nodeL.xmlInitNode('directory', name = 'postprocessing')
                    nodeL.xmlInitNode('time_dependency', choice = 'transient_connectivity')

                nn = node_out.xmlGetNode('writer', 'label', id = "-4")
                if nn == None:
                    nodeT = node_out.xmlInitNode('writer', id = "-4", label = 'trajectories')
                    nodeT.xmlInitNode('frequency', period = 'none')
                    nodeT.xmlInitNode('output_at_end', status = 'on')
                    nodeT.xmlInitNode('format', name = 'ensight', options = 'binary')
                    nodeT.xmlInitNode('directory', name = 'postprocessing')
                    nodeT.xmlInitNode('time_dependency', choice = 'fixed_mesh')

                nn = node_out.xmlGetNode('mesh', id = "-3")
                if nn == None:
                    node1 = node_out.xmlInitNode('mesh', id = "-3",
                                                 label = 'particles',
                                                 type = 'particles')
                    node1.xmlInitNode('all_variables', status = 'on')
                    node1.xmlInitNode('location')
                    node1.xmlSetData('location','all[]')
                    node1.xmlInitNode('density')
                    node1.xmlSetData('density', 1)
                    node1.xmlInitNode('writer', id = '-3')

        lst = self.case.xmlGetNodeList('external_coupling')
        if len(lst) > 1:
            for i in range(len(lst)):
                lst[i].xmlRemoveNode()


    def __backwardCompatibilityFrom_4_1(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        # update wall functions settings
        XMLThermoPhysicalModelNode = self.case.xmlGetNode('thermophysical_models')
        XMLTurbModelNode = XMLThermoPhysicalModelNode.xmlGetNode('turbulence')
        if XMLTurbModelNode:
            scaleModelNode = XMLTurbModelNode.xmlGetNode('scale_model')
            wallFunctionNode = XMLTurbModelNode.xmlGetNode('wall_function')

            if scaleModelNode and not wallFunctionNode:
                scale = XMLTurbModelNode.xmlGetInt('scale_model')

                if scale == 0:
                    wallFunction = 2
                elif scale == 1:
                    wallFunction = 3
                elif scale == 2:
                    wallFunction = 4
                else:
                    wallFunction = 0

                model = XMLTurbModelNode['model']
                if model == 'v2f-BL-v2/k' or \
                   model == 'Rij-EBRSM':
                    wallFunction = 0

                XMLTurbModelNode.xmlSetData('wall_function', wallFunction)
                scaleModelNode.xmlRemoveNode()

        node = XMLThermoPhysicalModelNode.xmlGetNode('velocity_pressure')
        if node:
            for f_name in ['effort', 'effort_tangential', 'effort_normal']:
                nn = node.xmlGetChildNode('property', name=f_name)
                if nn:
                    s = None
                    l = None
                    ns = nn.xmlGetChildNode('postprocessing_recording')
                    if ns:
                        s = ns['status']
                    l = nn['label']
                    nn = node.xmlRemoveChild('property', name=f_name)
                    l = 'Stress' + l[7:]
                    nn = self.setNewProperty(node, 'stress' + f_name[6:])
                    nn['label'] = l
                    nn['support'] = 'boundary'
                    if s:
                        nn.xmlInitNode('postprocessing_recording')['status']= s


    def __backwardCompatibilityFrom_4_2(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        # renames
        node = self.case.xmlGetNode('calculation_management')
        if node:
            cmd = node.xmlGetString('valgrind')
            if cmd:
                node.xmlSetData('debug', cmd)
                nn = node.xmlGetChildNode('valgrind')
                if nn:
                    nn.xmlRemoveNode()

        # darcy_model -> groundwater
        XMLThermoPhysicalModelNode = self.case.xmlInitNode('thermophysical_models')
        oldnode = XMLThermoPhysicalModelNode.xmlGetNode('darcy_model')
        if oldnode:
            mdl = oldnode['model']
            newnode = XMLThermoPhysicalModelNode.xmlInitNode('groundwater_model', 'model')
            newnode['model'] = mdl
            newnode.xmlChildsCopy(oldnode)
            oldnode.xmlRemoveNode()

        node = XMLThermoPhysicalModelNode.xmlGetNode('groundwater_model')
        if node:
            n = node.xmlGetNode('diffusion')
            if n:
                mdl = n['model']
                node.xmlInitNode('dispersion')['model']= mdl
                n.xmlRemoveNode()
            n = node.xmlGetNode('criterion')
            if n:
                n.xmlRemoveNode()

        XMLDarcy = XMLThermoPhysicalModelNode.xmlGetNode('darcy')
        if XMLDarcy:
            XMLGround = XMLThermoPhysicalModelNode.xmlInitNode('groundwater')
            XMLGround.xmlChildsCopy(XMLDarcy)
            XMLDarcy.xmlRemoveNode()

            nodelist = XMLGround.xmlGetNodeList('darcy_law')
            for oldnode in nodelist:
                mdl = oldnode['model']
                zid = oldnode['zone_id']
                newnode = XMLGround.xmlInitNode('groundwater_law', 'model')
                newnode['model'] = mdl
                newnode['zone_id'] = zid
                newnode.xmlChildsCopy(oldnode)
                oldnode.xmlRemoveNode()

        XMLSolutionDomainNode = self.case.xmlInitNode('solution_domain')
        self.__XMLVolumicConditionsNode = XMLSolutionDomainNode.xmlInitNode('volumic_conditions')
        for node in self.__XMLVolumicConditionsNode.xmlGetNodeList('zone'):
            law = node['darcy_law']
            if law:
                node['groundwater_law'] = law

        # Profile titles
        for node in self.case.xmlGetNodeList('profile'):
            node.xmlDelAttribute('title')

        # Lagrangian Model
        stats_node = None
        XMLLagrangianNode = self.case.xmlGetNode('lagrangian')
        if XMLLagrangianNode:
            stats_node = XMLLagrangianNode.xmlGetNode('statistics')
        if stats_node:
            v_node = stats_node.xmlGetNode('volume')
            b_node = stats_node.xmlGetNode('boundary')
            if v_node:
                it_start = v_node.xmlGetInt('iteration_start_volume')
                if it_start != None:
                    stats_node.xmlSetData('iteration_start', it_start)
                    v_node.xmlRemoveChild('iteration_start_volume')
                it_start_st = v_node.xmlGetInt('iteration_steady_start_volume')
                if it_start_st != None:
                    stats_node.xmlSetData('iteration_steady_start', it_start_st)
                    v_node.xmlRemoveChild('iteration_steady_start_volume')
                threshold = v_node.xmlGetDouble('threshold_volume')
                if threshold != None:
                    stats_node.xmlSetData('threshold', threshold)
                    v_node.xmlRemoveChild('threshold_volume')
            if b_node:
                for t in ['iteration_start_boundary', 'threshold_boundary']:
                    n = b_node.xmlGetNode(t)
                    if n:
                        n.xmlRemoveNode()
        if XMLLagrangianNode:
            if XMLLagrangianNode['model'] == "on":
                node = XMLLagrangianNode.xmlGetNode('coupling_mode')
                if node:
                    XMLLagrangianNode['model'] = node['model']
                    node.xmlRemoveNode()
        XMLBoundaryNode = self.case.xmlInitNode('boundary_conditions')
        for nature in ('inlet', 'outlet', 'wall', 'symmetry',
                       'free_inlet_outlet', 'groundwater'):
            for node in XMLBoundaryNode.xmlGetNodeList(nature):
                node['field_id'] = 'none'


    def __backwardCompatibilityFrom_4_3(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        for f_type in ['variable', 'property']:
            for node in self.case.xmlGetNodeList(f_type):
                n = node.xmlGetChildNode('probes')
                if n:
                    if node.xmlGetChildNode('probes_recording'):
                        node.xmlRemoveChild('probes_recording')
                    if n['choice'] == '0':
                        node.xmlInitChildNode('probes_recording')['status'] = "off"
                    n.xmlRemoveNode()

        # rename pressure into hydraulic_head for ground water module
        XMLThermoPhysicalModelNode = self.case.xmlGetNode('thermophysical_models')
        XMLPhysicalModelGWNode = XMLThermoPhysicalModelNode.xmlGetNode('groundwater')
        if XMLPhysicalModelGWNode:

            for XMLGWLaw in XMLPhysicalModelGWNode.xmlGetNodeList('groundwater_law'):
                if XMLGWLaw:
                    for XMLGWScalar in XMLGWLaw.xmlGetNodeList('variable'):
                        if XMLGWScalar:
                            XMLGWScalarProp = XMLGWScalar.xmlGetNode('property')
                            if XMLGWScalarProp:
                                if XMLGWScalarProp.xmlGetNode('formula'):
                                    XMLGWScalarProp.xmlRemoveChild('formula')
                                XMLGWScalar.xmlRemoveChild('property')
                            XMLGWLaw.xmlRemoveChild('variable')

                            print("Warning: xml parameter file must be updated manually.")
                            print("Settings of diffusivity and delay for the "
                                  "defined scalars have been removed "
                                  "to ensure a partial compatibility.")

            XMLVelPrNode = XMLThermoPhysicalModelNode.xmlGetNode('velocity_pressure')
            for node in XMLVelPrNode.xmlGetNodeList('variable', 'name'):
                if node['name'] == 'pressure':
                    node['name'] = 'hydraulic_head'
                if node['label'] == 'Pressure':
                    node['label'] = 'HydraulicHead'

            XMLBoundaryNode = self.case.xmlGetNode('boundary_conditions')
            XMLGWNodes = XMLBoundaryNode.xmlGetNodeList('groundwater')
            for XMLGWNode in XMLGWNodes:
                for node in (XMLGWNode.xmlGetNodeList('neumann') + XMLGWNode.xmlGetNodeList('dirichlet')):
                    if node['name'] == 'pressure':
                        node['name'] = 'hydraulic_head'

        XMLPhysicalModelALE = XMLThermoPhysicalModelNode.xmlGetNode('ale_method')
        if XMLPhysicalModelALE:
            node = XMLPhysicalModelALE.xmlGetNode('external_coupling_post_synchronizationre')
        if node:
            node.xmlRemoveNode()

        node = self.case.xmlGetNode('variable', name='vec_potential')
        if node:
            node['dim'] = '3'

        for var in ['current_re', 'laplace_force', 'current_im']:
            node = self.case.xmlGetNode('property', name=var)
            if node:
                node['dim'] = '3'

        node = self.case.xmlGetNode('variable', name='vec_potential_01')
        if node:
            node.xmlRemoveNode()
        node = self.case.xmlGetNode('variable', name='vec_potential_02')
        if node:
            node.xmlRemoveNode()
        node = self.case.xmlGetNode('variable', name='vec_potential_03')
        if node:
            node.xmlRemoveNode()
        node = self.case.xmlGetNode('property', name='current_re_1')
        if node:
            node.xmlRemoveNode()
            nModels         = self.case.xmlGetNode('thermophysical_models')
            self.node_joule = nModels.xmlInitNode('joule_effect',      'model')
            self.setNewProperty(self.node_joule, 'current_re', dim = '3')
            self.setNewProperty(self.node_joule, 'laplace_force', dim = '3')
            self.setNewProperty(self.node_joule, 'electric_field', dim = '3')
        node = self.case.xmlGetNode('property', name='current_re_2')
        if node:
            node.xmlRemoveNode()
        node = self.case.xmlGetNode('property', name='current_re_3')
        if node:
            node.xmlRemoveNode()

        rd_list = None
        XMLBoundaryNode = self.case.xmlGetNode('boundary_conditions')
        if XMLBoundaryNode:
            rd_list = XMLBoundaryNode.xmlGetNodeList('radiative_data')
        if rd_list:
            for rd in rd_list:
                node = rd.xmlGetNode('output_zone')
                if node:
                    node.xmlRemoveNode()

        node_joule = XMLThermoPhysicalModelNode.xmlGetNode('joule_effect', 'model')
        if node_joule:
            if node_joule['model'] != 'off':
                node_gas = XMLThermoPhysicalModelNode.xmlGetNode('gas_combustion',
                                                                 'model')
                if node_gas:
                    node = node_gas.xmlGetNode('data_file')
                    if node:
                        node.xmlRemoveNode()

        # Update some postprocessing outputs

        npr = XMLThermoPhysicalModelNode.xmlGetNode('thermal_scalar')
        if npr:
            if npr['model'] != 'off':
                node = npr.xmlGetNode('property', name="input_thermal_flux")
                node2 = npr.xmlGetNode('property', name="thermal_flux")
                if node2 and node:
                    node.xmlRemoveNode()
                elif node:
                    node['name'] = "thermal_flux"
            if npr['model']:
                ThermalScalarModel(self.case).setThermalModelOutputs(npr['model'])

        npr = XMLThermoPhysicalModelNode.xmlGetNode('radiative_transfer')
        if npr:
            for name in ("wall_temp", "qrad_y", "qrad_z"):
                node = npr.xmlGetNode('property', name=name)
                if node:
                    node.xmlRemoveNode()
            node = npr.xmlGetNode('property', name="qrad_x")
            if node:
                node['name'] = "qrad"
                node['label'] = "Qrad"

        npr = XMLThermoPhysicalModelNode.xmlGetNode('velocity_pressure')
        if npr:
            node = npr.xmlGetNode('property', name="all_variables")
            if node:
                node.xmlRemoveNode()


        # renames

        self.__renameSingle('thermophysical_models', 'heads_losses', 'head_losses')
        self.__renameSingle('solution_domain', 'extrude_meshes', 'extrusion')

        # renumber boundary and volume zones if required

        m = LocalizationModel('BoundaryZone', self.case)
        m.renumberZones()

        m = LocalizationModel('VolumicZone', self.case)
        m.renumberZones()

        # missing info from some old files

        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                if not node['name']:
                    node['name'] = node['label']

        # fix name of thermal conductivity in radiative transfer node
        npr = XMLThermoPhysicalModelNode.xmlGetNode('radiative_transfer')
        if npr:
            node = npr.xmlGetNode('property', name="thermal_conductivity")
            if node:
                node['name'] = "wall_thermal_conductivity"


    def __backwardCompatibilityFrom_5_0(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        XMLThermoPhysicalModelNode = self.case.xmlGetNode('thermophysical_models')

        nch = XMLThermoPhysicalModelNode.xmlGetNode('solid_fuels')
        if nch:
            if nch['model'] == 'homogeneous_fuel_moisture_lagr':
                nch['model'] = 'homogeneous_fuel_moisture'

        # fix name of Reynolds stress tensor
        ntur = XMLThermoPhysicalModelNode.xmlGetNode('turbulence')
        if ntur:
            if ntur['model'] in  ('Rij-SSG', 'Rij-epsilon', 'Rij-EBRSM'):

                node = ntur.xmlGetNode('variable', name="r11")
                if node:
                    node['label'] = "Rij"
                    node['name']  = "rij"
                    node['dimension'] = "6"

                node = ntur.xmlGetNode('variable', name="r22")
                if node:
                    node.xmlRemoveNode()
                node = ntur.xmlGetNode('variable', name="r33")
                if node:
                    node.xmlRemoveNode()
                node = ntur.xmlGetNode('variable', name="r12")
                if node:
                    node.xmlRemoveNode()
                node = ntur.xmlGetNode('variable', name="r23")
                if node:
                    node.xmlRemoveNode()
                node = ntur.xmlGetNode('variable', name="r13")
                if node:
                    node.xmlRemoveNode()

            # Update Time Averages using Rij if needed !
            rij_lbls = ['r11', 'r22', 'r33', 'r12', 'r23', 'r13']
            for node in self.case.xmlGetNodeList('time_average'):
                if node:
                    for var in node.xmlGetChildNodeList('var_prop'):
                        varName = var['name']
                        if varName in rij_lbls:
                            var['name'] = 'rij'
                            var['component'] = rij_lbls.index(varName)

            # Update profiles using Rij
            for node in self.case.xmlGetNodeList('profile'):
                if node:
                    for var in node.xmlGetChildNodeList('var_prop'):
                        varName = var['name']
                        if varName in rij_lbls:
                            var['name'] = 'rij'
                            var['component'] = rij_lbls.index(varName)


            # Lagrangian model setup renames and cleanup
            XMLBoundaryNode = self.case.xmlGetNode('boundary_conditions')
            XMLInletNodes = XMLBoundaryNode.xmlGetNodeList('inlet')
            for XMLInletNode in XMLInletNodes:
                XMLParticeNodes = XMLInletNode.xmlGetNodeList('particles')
                for XMLParticleNode in XMLParticeNodes:
                    XMLClassNodes = XMLInletNode.xmlGetNodeList('class')
                    for XMLClassNode in XMLClassNodes:
                        node = XMLClassNode.xmlGetNode('velocity', choice='subroutine')
                        if node:
                            node['choice'] = 'components'
                        else:
                            node = XMLClassNode.xmlGetNode('velocity', choice='components')
                            if node:
                                for d in (('u', 'x'), ('v', 'y'), ('w', 'z')):
                                    onode = node.xmlGetNode('velocity_' + d[0])
                                    if onode:
                                        nnode = node.xmlInitNode('velocity_' + d[1])
                                        nnode.xmlChildsCopy(onode)
                                        onode.xmlRemoveNode()
                        if node:
                            node['choice'] = 'components'
                        for attr in ('statistical_weight', 'temperature'):
                            node = XMLClassNode.xmlGetNode(attr, choice='subroutine')
                            if node:
                                node['choice'] = 'prescribed'
                        node = XMLClassNode.xmlGetNode('diameter')
                        if node:
                            node.xmlDelAttribute('choice')
                        for attr in ('coal_composition', 'statitical_weight'):
                            node = XMLClassNode.xmlGetNode(attr)
                            if node:
                                node.xmlRemoveNode()

            XMLLagrangianNode = self.case.xmlGetNode('lagrangian')
            if XMLLagrangianNode:
                for attr in ('continuous_injection', 'particles_max_number'):
                    node = XMLLagrangianNode.xmlGetNode(attr)
                    if node:
                        node.xmlRemoveNode()
                stats_node = XMLLagrangianNode.xmlGetNode('statistics')
                if stats_node:
                    vol_node = stats_node.xmlGetNode('volume')
                    if vol_node:
                        for attr in ('Part_velocity_X', 'Part_velocity_Y',
                                     'Part_velocity_Z'):
                            node = vol_node.xmlGetNode('property', name=attr)
                            if node:
                                node.xmlRemoveNode()
                output_node = XMLLagrangianNode.xmlGetNode('output')
                if output_node:
                    for attr in ('trajectory', 'particles', 'number_of_particles',
                                 'postprocessing_format', 'postprocessing_options',
                                 'postprocessing_frequency'):
                        node = output_node.xmlGetNode(attr)
                        if node:
                            node.xmlRemoveNode()
        return

    def __backwardCompatibilityFrom_5_1(self):
        """
        Change XML in order to ensure backward compatibility.
        """

    def __backwardCompatibilityFrom_5_2(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        XMLAnaControl = self.case.xmlGetNode('analysis_control')
        node_steady = XMLAnaControl.xmlGetNode('steady_management')
        if node_steady:
            if node_steady['status'] == "on":
                node_time = XMLAnaControl.xmlGetNode('time_parameters')
                idtvar = 2
                n = self.case.xmlGetNode('numerical_parameters')
                if n:
                    n = n.xmlGetNode('velocity_pressure_algo', 'choice')
                    if n:
                        if n['choice'] == 'simple':
                            idtvar = -1
                if idtvar == -1:
                    iterations = node_steady.xmlGetString('iterations')
                    if iterations != None:
                        node_time.xmlSetData('iterations', iterations)
                    relax_c = node_steady.xmlGetString('relaxation_coefficient')
                    if relax_c != None:
                        node_time.xmlSetData('relaxation_coefficient', relax_c)
                node_time.xmlSetData('time_passing', idtvar)
            node_steady.xmlRemoveNode()

        XMLThermoPhysicalModelNode = self.case.xmlGetNode('thermophysical_models')
        npr = XMLThermoPhysicalModelNode.xmlGetNode('radiative_transfer')
        if npr:
            node = npr.xmlGetNode('property', name="qrad")
            if node:
                node['name'] = "radiative_flux"
            node = npr.xmlGetNode('property', name="srad")
            if node:
                node['name'] = "rad_st"
            node = npr.xmlGetNode('property', name="flux_convectif")
            if node:
                node['name'] = "rad_convective_flux"
            node = npr.xmlGetNode('property', name="flux_incident")
            if node:
                node['name'] = "rad_incident_flux"
            node = npr.xmlGetNode('property', name="flux_net")
            if node:
                node['name'] = "rad_net_flux"
            node = npr.xmlGetNode('property', name="coeff_ech_conv")
            if node:
                node['name'] = "rad_exchange_coefficient"
            node = npr.xmlGetNode('property', name="absorption")
            if node:
                node['name'] = "rad_absorption"
            node = npr.xmlGetNode('property', name="absorption_coefficient")
            if node:
                node['name'] = "rad_absorption_coeff"
            node = npr.xmlGetNode('property', name="emission")
            if node:
                node['name'] = "rad_emission"

    def __backwardCompatibilityFrom_5_3(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        XMLAnaControl = self.case.xmlGetNode('analysis_control')
        node_time = XMLAnaControl.xmlGetNode('time_parameters')
        if node_time:
            node = node_time.xmlGetNode('zero_time_step')
            if node:
                node.xmlRemoveNode()

        # Place ALE boundary definitions in the correct type of node
        # instead of a wall node

        XMLBoundaryNode = self.case.xmlInitNode('boundary_conditions')
        for node in XMLBoundaryNode.xmlGetNodeList('wall'):
            nn = node.xmlGetNode("ale")
            if nn:
                label = node['label']
                for nature in ('inlet', 'outlet', 'symmetry',
                               'free_inlet_outlet', 'groundwater'):
                    for nc in XMLBoundaryNode.xmlGetNodeList(nature):
                        if nc['label'] == label:
                            nc.xmlChildsCopy(node)
                            node.xmlRemoveNode()


    def __backwardCompatibilityCurrentVersion(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        # Reference values
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        nodeRefValues = XMLThermoPhysicalNode.xmlInitNode('reference_values')
        nodeTurb = XMLThermoPhysicalNode.xmlInitNode('turbulence', 'model')

        for ref_str in ['velocity', 'length']:
            value = nodeRefValues.xmlGetDouble(ref_str)
            if value:
                nodeTurb.xmlSetData('reference_' + ref_str, value)
                nodeRefValues.xmlRemoveChild(ref_str)


#-------------------------------------------------------------------------------
# XMLinit test case
#-------------------------------------------------------------------------------

class XMLinitTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from code_saturne.Base import XMLengine
        Toolbox.GuiParam.lang = 'en'
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
        doc = \
        '<Code_Saturne_GUI case="" study="" version="1.0">'\
        '<solution_domain/>'\
        '<thermophysical_models>'\
                '<velocity_pressure>'\
                        '<variable label="Pressure" name="pressure"/>'\
                        '<variable label="Velocity" name="velocity"/>'\
                        '<property label="total_pressure" name="total_pressure"/>'\
                '</velocity_pressure>'\
                '<turbulence model="k-epsilon">'\
                        '<variable label="TurbEner" name="turb_k"/>'\
                        '<variable label="Dissip" name="turb_eps"/>'\
                        '<property label="TurbVisc" name="turbulent_viscosity"/>'\
                        '<initialization choice="reference_velocity">'\
                                '<reference_velocity>1.0</reference_velocity>'\
                        '</initialization>'\
                '</turbulence>'\
                '<initialization>'\
                        '<zone name="1">0</zone>'\
                '</initialization>'\
                '<thermal_scalar model="off"/>'\
                '<gas_combustion model="off"/>'\
                '<solid_fuels model="off"/>'\
                '<joule_effect model="off"/>'\
                '<radiative_transfer model="off"/>'\
        '</thermophysical_models>'\
        '<numerical_parameters/>'\
        '<physical_properties>'\
                '<fluid_properties>'\
                        '<property choice="constant" label="Density" name="density">'\
                                '<initial_value>1.17862</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="LamVisc" name="molecular_viscosity">'\
                                '<initial_value>1.83e-05</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="SpecHeat" name="specific_heat">'\
                                '<initial_value>1017.24</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="ThermalCond" name="thermal_conductivity">'\
                                '<initial_value>0.02495</initial_value>'\
                        '</property>'\
                '</fluid_properties>'\
        '</physical_properties>'\
        '<additional_scalars/>'\
        '<boundary_conditions/>'\
        '<analysis_control>'\
                '<time_parameters>'\
                        '<time_step_ref>0.1</time_step_ref>'\
                        '<property label="CourantNb" name="courant_number">'\
                        '<property label="FourierNb" name="fourier_number">'\
                '</time_parameters>'\
        '</analysis_control>'\
        '<calculation_management/>'\
        '</Code_Saturne_GUI>'

        XMLinit(self.case).initialize()

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
