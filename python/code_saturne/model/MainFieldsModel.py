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

import sys, unittest, copy

from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.Common import LABEL_LENGTH_MAX
from code_saturne.model.NeptuneFieldModel import NeptuneField
from code_saturne.model.GlobalNumericalParametersModel import GlobalNumericalParametersModel

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

from code_saturne.model.EosWrapper import eosWrapper

class PredefinedFlowsModel:
    """
    This class manages the Field objects for a predefined flow in the XML file
    """

    fieldsCouple = ["None",
                    "free_surface",
                    "boiling_flow",
                    "droplet_flow",
                    "particles_flow",
                    "multiregime"]

    fieldsCoupleProperties = {}
    fieldsCoupleProperties[fieldsCouple[0]] = ("", "")
    fieldsCoupleProperties[fieldsCouple[1]] = (("continuous", "continuous"),
                                               ("liquid", "gas"))
    fieldsCoupleProperties[fieldsCouple[2]] = (("continuous", "dispersed"),
                                               ("liquid", "gas"))
    fieldsCoupleProperties[fieldsCouple[3]] = (("continuous", "dispersed"),
                                               ("gas", "liquid"))
    fieldsCoupleProperties[fieldsCouple[4]] = (("continuous", "dispersed"),
                                               ("gas", "solid"))
    fieldsCoupleProperties[fieldsCouple[5]] = (("continuous", "continuous"),
                                               ("liquid", "gas"))

#-------------------------------------------------------------------------------
# Model for main fields
#-------------------------------------------------------------------------------

class MainFieldsModel(Variables, Model):
    """
    This class manages the Field objects in the XML file
    """

    def __init__(self, case):
        """
        Constructor.
        """
        # XML file parameters
        self.case = case

        self.XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
        self.__XMLNodefields = self.XMLNodethermo.xmlInitNode('fields')
        self.XMLNodeVariable = self.XMLNodethermo.xmlInitNode('variables')
        self.XMLNodeproperty = self.XMLNodethermo.xmlInitNode('properties')

        self.XMLClosure      = self.case.xmlGetNode('closure_modeling')

        self.XMLturbulence   = self.XMLClosure.xmlInitNode('turbulence')
        self.XMLforces       = self.XMLClosure.xmlInitNode('interfacial_forces')
        self.node_anal       = self.case.xmlInitNode('analysis_control')
        self.node_average    = self.node_anal.xmlInitNode('time_averages')
        self.node_profile    = self.node_anal.xmlInitNode('profiles')

        pressure_node = self.XMLNodethermo.xmlGetNode('variable',
                                                      name='pressure')
        if pressure_node is None:
            pressure_node = self.XMLNodethermo.xmlGetNode('variable',
                                                          name='Pressure')
        if pressure_node is None:
            self.setNewVariableProperty("variable",
                                                        "",
                                                        self.XMLNodeVariable,
                                                        "none",
                                                        "pressure",
                                                        "Pressure",
                                                        post = True)
        porosity_node = self.XMLNodethermo.xmlGetNode('property',
                                                      name='porosity')
        if porosity_node is None:
            self.setNewVariableProperty("property",
                                                        "",
                                                        self.XMLNodeproperty,
                                                        "none",
                                                        "porosity",
                                                        "porosity")

        #EOS
        self.eos = eosWrapper()

        # Test of NeptuneField
        self.list_of_fields = self.load_fields_from_xml()

    def load_fields_from_xml(self):
        fields = []
        for field_id in self.getFieldIdList():
            field = NeptuneField(self.case, field_id)
            field.load_from_xml()
            fields.append(field)

        return fields

    def defaultValues(self):
        default = {}
        default['id']                         = ""
        default['label']                      = "defaultLabel"
        if self.getPhaseChangeTransferStatus() == "on":
            default['enthalpyResolutionStatus']   = "on"
            default['enthalpyResolutionModel']    = "total_enthalpy"
        else:
            default['enthalpyResolutionStatus']   = "off"
            default['enthalpyResolutionModel']    = "off"
        default['typeChoice']                 = "continuous"
        default['phase']                      = "liquid"
        default['carrierField']               = "off"
        default['compressibleStatus']         = "off"
        default['defaultPredefinedFlow'] = "None"
        default["phaseChangeTransfer"] = "off"

        return default


    def getFieldIdList(self, include_none=False):
        """
        Return the id field list
        """
        lst = []

        if include_none:
            lst.append('none')

# TODO : see if we can move this XML manipulation to NeptuneField Object
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            lst.append(node['field_id'])

        return lst

    def getFieldFromId(self, fieldId):
        for field in self.list_of_fields:
            if str(field.f_id) == str(fieldId):
                return field
        if fieldId == "none":
            return NeptuneField(self.case, "none")
        elif fieldId == "off":
            return NeptuneField(self.case, "off")
        elif fieldId == "all":
            return NeptuneField(self.case, "all")
        return None

    @Variables.undoGlobal
    def addField(self,fieldId=None):
        """
        add a new field
        """
        label = ""
        labNum= 0

        if fieldId not in self.getFieldIdList():
            label = "Field" + str(len(self.getFieldIdList())+1)
            labNum = len(self.getFieldIdList())+1
            if label in self.getFieldLabelsList():
               labelNumber = 1
               label = "Field" + str(labelNumber)
               labNum = labelNumber
               while label in self.getFieldLabelsList():
                  labelNumber += 1
                  label = "Field" + str(labelNumber)
                  labNum = labelNumber
            fieldId = len(self.getFieldIdList())+1

            type         = self.defaultValues()['typeChoice']
            phase        = self.defaultValues()['phase']
            hmodel       = self.defaultValues()['enthalpyResolutionModel']
            compressible = self.defaultValues()['compressibleStatus']
            carrierFieldId = self.defaultValues()['carrierField']

            # TODO replace these lines by call to setDefinedField
            new_field = NeptuneField(self.case, fieldId)
            new_field.label = label
            new_field.flow_type = type
            new_field.phase = phase
            new_field.carrier_id = carrierFieldId
            new_field.enthalpy_model = hmodel
            new_field.compressible = compressible

            new_field.createVariableProperties()
            self.list_of_fields.append(new_field)
            # Update default carrier field id
            NeptuneField.default_carrier_id = self.getFirstContinuousField()

        return new_field


    @Variables.undoGlobal
    def addDefinedField(self, fieldId, label, typeChoice, phase, carrierField="off",
                        hmodel="off", compressible="off"):
        """
        add field with initiaization attributes (for test purposes mostly)
        """

        self.isInList(phase, ['solid', 'liquid', 'gas'])
        self.isInList(typeChoice, ['dispersed', 'continuous'])

        # Check that the field does not already exist
        field = self.getFieldFromId(fieldId)
        if field != None:
            field.label = label
        else:
            fieldId = (self.addField()).f_id
        self.setDefinedField(fieldId, typeChoice, phase, carrierField, hmodel, compressible)

    def setDefinedField(self, fieldId, typeChoice, phase, carrierFieldId, energyModel, compressible):
        """
        Set properties of an already defined field
        """

        field = NeptuneField(self.case, fieldId)
        field.flow_type = typeChoice
        field.phase = phase
        field.carrier_id = carrierFieldId
        field.enthalpy_model = energyModel
        field.compressible = compressible
        is_multiregime = (self.getPredefinedFlow() == "multiregime")
        field.createVariableProperties(is_multiregime)
        NeptuneField.default_carrier_id = self.getFirstContinuousField()


    def getFieldLabelsList(self, include_none=False):
        """
        return list of label for field
        """
        lst = []

        if include_none:
            lst.append('none')

        for field in self.list_of_fields:
            lst.append(field.label)
        return lst


    def getContinuousFieldList(self):
        """ return list of continuous fields """
        return [f for f in self.list_of_fields if f.flow_type == "continuous"]

    def getDispersedFieldList(self):
        """ return list of dispersed fields """
        return [f for f in self.list_of_fields if f.flow_type == "dispersed"]

    def getGasPhaseList(self):
        """ return list of gas fields """
        return [f for f in self.list_of_fields if f.phase == "gas"]

    def getSolidPhaseList(self):
        """ return list of solid fields """
        return [f for f in self.list_of_fields if f.phase == "solid"]

    def getFirstContinuousField(self):
        """
        return id of first continuous field
        """
        f_id = 0
        if len(self.getContinuousFieldList()) > 0:
            f_id = self.getContinuousFieldList()[0].f_id
        return f_id


    def getFirstGasField(self) :
        """
        return id of first continuous field
        """
        fid = 0
        if len(self.getGasPhaseList()) > 0:
            fid = int(self.getGasPhaseList()[0].f_id)
        return fid


    @Variables.noUndo
    def hasEnthalpyResolvedField(self) :
        """ return list of fields with enthalpy resolution """
        return ([f for f in self.list_of_fields if f.enthalpy_model != "off"] != [])


    @Variables.undoLocal
    def setCriterion(self, fieldId, ftype):
        """
        Put type of field
        """
        #TODO this method is too complex and should be refactored, but it might imply a complete change of architecture
        #     for MainFieldsModel, InterfacialForcesModel, InterfacialEnthalpyModel, possibly TurbulenceNeptuneModel.
        #TODO move to NeptuneFieldModel if possible (it might be difficult to update
        # the carrier field properly without the knowledge of other fields)
        self.isFieldIdValid(fieldId)

        field = self.getFieldFromId(fieldId)
        oldtype = field.flow_type
        field.flow_type = ftype
        field_name = field.label

        # update carrier field  if current field is not continous anymore
        if ftype != "continuous":
            field.carrier_id = self.getFirstContinuousField()
            # Modify carrier field
            if oldtype != ftype:
                for disp_field in self.getDispersedFieldList():
                    if disp_field.carrier_id == str(fieldId):
                        disp_field.carrier_id = self.getFirstContinuousField()

        # if ftype is changed from continuous to dispersed and vice versa, delete old coupling information
        if oldtype != ftype and oldtype == "continuous":
            node = self.XMLturbulence.xmlGetNode("field", field_id=fieldId)
            if node:
                node.xmlRemoveNode()
            node_list = self.XMLforces.xmlGetChildNodeList("continuous_field_momentum_transfer", field_id_a=fieldId)
            node_list += self.XMLforces.xmlGetChildNodeList("continuous_field_momentum_transfer", field_id_b=fieldId)
            for node in node_list:
                if node:
                    node.xmlRemoveNode()
        if oldtype != ftype and oldtype == "dispersed":
            node = self.XMLturbulence.xmlGetNode("field", field_id=fieldId)
            if node:
                node.xmlRemoveNode()
            node_list = self.XMLforces.xmlGetChildNodeList("force",
                                                           field_id_b=fieldId)  # dispersed phase is always in second position
            for node in node_list:
                if node:
                    node.xmlRemoveNode()

        # if number of continuous field < 2 suppress coupling and momentum forces
        if len(self.getContinuousFieldList()) < 2:
            node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
            if node:
                node.xmlRemoveNode()
            node = self.XMLforces.xmlGetNode('continuous_field_momentum_transfer')
            if node:
                node.xmlRemoveNode()

        # remove closure law linked to field for coherence of XML
        if oldtype != ftype :
           for node in self.__nodesWithFieldIDAttribute():
               try :
                   if (node['field_id_a'] == str(fieldId) or node['field_id_b'] == str(fieldId)):
                       node.xmlRemoveNode()
               except :
                   pass

        # TODO mettre en coherence pour les aires interf., tout ce qui est closure law a faire aussi pour la nature.
        # Activated if dispersed or second continuous phase of GLIM
        if field.flow_type == "dispersed" or self.getPredefinedFlow() == "multiregime":
           self.setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "diameter", "diam_"+field_name)
           self.setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "drift_component", "drift_component_"+field_name, dim='3')
        else :
           self.removeVariableProperty("property", self.XMLNodeproperty, fieldId, "diameter")
           self.removeVariableProperty("property", self.XMLNodeproperty, fieldId, "drift_component")

        self.updateXML()

    @Variables.undoGlobal
    def deleteField(self, fieldId):
        """
        delete a field in XML and update
        """
        self.isFieldIdValid(fieldId)
        field_to_delete = self.getFieldFromId(fieldId)

        # Update for field Id
        # Note : this works only because the first field is always continuous !!!
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            nodec = node.xmlGetNode('carrier_field')
            if nodec != None :
                if str(nodec['field_id']) == str(fieldId) :
                    id = int(self.getFirstContinuousField())
                    if id > 0 :
                       nodec['field_id'] = id
                    else :
                       # only dispersed field -> change to continuous field
                       currentId = node['field_id']
                       self.setCriterion(currentId, "continuous")
                       self.getFieldFromId(currentId).carrier_id = self.defaultValues()["carrierField"]

        # Remove all previous occurences of fieldId
        # Warning : this needs to be done after carrier fields have been updated
        nodes_to_remove = self.case.xmlGetNodeWithAttrList("field_id", field_id=fieldId)
        nodes_to_remove += self.case.xmlGetNodeWithAttrList("field_id_a", field_id_a=fieldId)
        nodes_to_remove += self.case.xmlGetNodeWithAttrList("field_id_b", field_id_b=fieldId)
        if str(fieldId) == "1":
            nodes_to_remove += self.case.xmlGetNodeWithAttrList("field_id", field_id="none")
        for node in nodes_to_remove:
            node.xmlRemoveNode()

        self.list_of_fields = [fld for fld in self.list_of_fields if fld.f_id != str(fieldId)]
        # Renumber fields after deleted field (this is necessary because the solver works
        # with a consecutive numbering of the fields)
        for field in self.list_of_fields[int(fieldId)-1:]:
            field.f_id = int(field.f_id) - 1

        self.updateXML()
        NeptuneField.default_carrier_id = self.getFirstContinuousField()
        from code_saturne.model.SpeciesModel import SpeciesModel # WARNING : potential circular dependency here
        SpeciesModel(self.case).forceConsecutiveScalarIds()


    def updateXML(self):
        """
        method for update in case of suppress or change attribute
        """
        # suppress solid information
        if (len(self.getSolidPhaseList()) < 1) :
            node = self.XMLNodethermo.xmlGetNode('solid_compaction')
            if node:
                node.xmlRemoveNode()

        # suppress continuous-continuous information
        if (len(self.getContinuousFieldList()) < 2) :
            node = self.XMLforces.xmlGetNode('continuous_field_coupling')
            if node:
                node.xmlRemoveNode()
            node = self.XMLforces.xmlGetNode('continuous_field_momentum_transfer')
            if node:
                node.xmlRemoveNode()

        # suppress continuous-dispersed information
        if (len(self.getDispersedFieldList()) < 1) and self.getPredefinedFlow() != "multiregime":
            node = self.XMLClosure.xmlGetNode('interfacial_area_diameter')
            if node:
                node.xmlRemoveNode()


    def __nodesWithFieldIDAttribute(self) :
        """
        Return a list of nodes whose one of atributes is 'field_id'
        """
        root = self.case.root()
        list = []
        tagNodesWithFileIdAttribute = ('field', 'inlet', 'outlet', 'variable',
                                       'property', 'scalar',
                                       'initialization', 'var_prop',
                                       'wall', 'symetry', 'transfer', 'force',
                                       'closure_modeling', 'interfacial_area_diameter')
        for node in tagNodesWithFileIdAttribute :
            list1 = root.xmlGetNodeList(node)
            for n in list1:
                list.append(n)
        return list

    @Variables.undoGlobal
    def setPredefinedFlow(self, flow_choice, overwriteEnergy=True):
        """
        """

        self.isInList(flow_choice, PredefinedFlowsModel.fieldsCouple)

        node = self.XMLNodethermo.xmlInitChildNode('predefined_flow')
        node.xmlSetAttribute(choice=flow_choice)
        self.XMLNodeclosure = self.case.xmlGetNode('closure_modeling')
        self.XMLMassTrans = self.XMLNodeclosure.xmlInitNode('mass_transfer_model')

        # Variables : neptune_cfd.core.XMLvariables.Variables
        self.setNewVariableProperty("property", "constant", self.XMLNodeproperty, "none",
                                                    "surface_tension", "Surf_tens")
        energyModel = "total_enthalpy"
        if self.getPhaseChangeTransferStatus() == "off":
            energyModel = "off"
        if flow_choice == "particles_flow":
            energyModel = "specific_enthalpy"

        field_id_list = self.getFieldIdList()
        label_list = self.getFieldLabelsList()

        create_fields = False
        if len(field_id_list) < 2:
            field_id_list = ["1", "2"]
            label_list = ["Field1", "Field2"]
            create_fields = True

        if flow_choice != "None":

            for id in range(2):  # Always two fields in predefined regimes
                fieldId = field_id_list[id]
                typeChoice = PredefinedFlowsModel.fieldsCoupleProperties[flow_choice][0][id]
                phase = PredefinedFlowsModel.fieldsCoupleProperties[flow_choice][1][id]
                field = self.getFieldFromId(fieldId)
                if typeChoice == "dispersed":
                    carrierField = field_id_list[1 - id]  # other field
                else:
                    carrierField = "off"
                if not overwriteEnergy:
                    energyModel = field.enthalpy_model
                self.case.undoStop()
                if not (create_fields):
                    compressible = field.compressible
                    self.setDefinedField(fieldId, typeChoice, phase, carrierField, energyModel,
                                         compressible)
                else:  # Otherwise create a new field
                    compressible = "off"
                    self.addDefinedField(fieldId, label_list[id], typeChoice, phase, carrierField,
                                         energyModel,
                                         compressible)

            # Remove remnant fields from previous flow choice, starting from the last one
            for fieldId in field_id_list[-1:1:-1]:
                self.deleteField(fieldId)

            # modification du choix pour le predicteur de vitesse
            self.case.undoStop()
            if flow_choice == "boiling_flow":
                GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(
                    GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_bubble'])
            else:
                GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(
                    GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_std'])

            if flow_choice != "boiling_flow" and flow_choice != "free_surface" and flow_choice != "multiregime":
                self.XMLMassTrans.xmlRemoveChild('nucleate_boiling')

            if flow_choice == "particles_flow":
                self._deleteFieldsProperties()
            else:
                # Recreate fields in previous flow choice is "particles_flow"
                status = self.getPhaseChangeTransferStatus()
                self.setPhaseChangeTransferStatus(status)

        else :
            GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(
                GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_std'])
            self.XMLMassTrans.xmlRemoveChild('nucleate_boiling')

            # create two field flow continuous/continuous by default
            # Create the first field
            if create_fields:
                fieldId = "1"
                label = "Field1"
                typeChoice = "continuous"
                phase = "liquid"
                carrierField = "off"
                if not overwriteEnergy:
                    energyModel = self.getFieldFromId(fieldId).enthalpy_model
                compressible = "off"
                self.case.undoStop()
                self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, energyModel, compressible)
                # Create the second field
                fieldId = "2"
                label = "Field2"
                typeChoice = "continuous"
                phase = "gas"
                carrierField = "off"
                if not overwriteEnergy:
                    energyModel = self.getFieldFromId(fieldId).enthalpy_model
                compressible = "off"
                self.case.undoStop()
                self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, energyModel, compressible)

        self.case.undoStart()

        return node

    def _setFieldMaterial(self, fieldId, material):
        from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
        self.case.undoStop()
        ThermodynamicsModel(self.case).setMaterials(fieldId, material)
        fls = self.eos.getFluidMethods(material)
        if "Cathare2" in fls:
            ThermodynamicsModel(self.case).setMethod(fieldId, "Cathare")
        elif "Cathare" in fls:
            ThermodynamicsModel(self.case).setMethod(fieldId, "Cathare")
        else:
            ThermodynamicsModel(self.case).setMethod(fieldId, fls[0])

    def _deleteFieldsProperties(self):
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationTemperature")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none",
                                                    "SaturationEnthalpyLiquid")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationEnthalpyGas")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Liquid")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Gas")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Tsat_d_P")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "LatentHeat")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_total_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_total_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_evaporation_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_quenching_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none",
                                                    "wall_liquid_convective_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none",
                                                    "wall_steam_convective_flux")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "boundary_temperature")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_temperature")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_oversaturation")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "unal_diameter")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none",
                                                    "wall_diameter_mesh_independancy")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_roughness")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none",
                                                    "wall_dispersed_phase_mass_source_term")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "boiling_criteria")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "exchange_coefficient")
        self.removeVariableProperty("property", self.XMLNodeproperty, "none", "uninfluenced_part")

    def _createSaturationProperties(self):
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "SaturationTemperature", "TsatK", post=True)
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "SaturationEnthalpyLiquid", "Hsat_Liquid")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "SaturationEnthalpyGas", "Hsat_Gas")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Hsat_d_P_Liquid",
                                                    "dHsat_dp_Liquid")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Hsat_d_P_Gas",
                                                    "dHsat_dp_Gas")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Tsat_d_P",
                                                    "dTsat_dp")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "LatentHeat", "Hlat")

    def _createWallFieldsProperties(self):
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_total_flux",
                                                    "wall_total_flux", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_liquid_total_flux", "wall_liquid_total_flux",
                                                    support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_evaporation_flux", "wall_evaporation_flux",
                                                    support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_quenching_flux",
                                                    "wall_quenching_flux", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_liquid_convective_flux", "wall_liquid_convective_flux",
                                                    support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_steam_convective_flux", "wall_steam_convective_flux",
                                                    support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "boundary_temperature", "wall_temperature", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_liquid_temperature", "wall_liquid_temperature",
                                                    support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_oversaturation",
                                                    "wall_oversaturation", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "unal_diameter",
                                                    "unal_diameter", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_diameter_mesh_independancy",
                                                    "wall_diameter_mesh_independancy", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_roughness",
                                                    "wall_roughness", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "wall_dispersed_phase_mass_source_term",
                                                    "wall_dispersed_phase_mass_source_term", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "boiling_criteria",
                                                    "boiling_criteria", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none",
                                                    "exchange_coefficient", "exchange_coefficient", support="boundary")
        self.setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "uninfluenced_part",
                                                    "uninfluenced_part", support="boundary")

    @Variables.noUndo
    def detectFlowType(self):
        """
        Detect the type of flow set by the user.
        This method is similar to but different from getPredefinedFlow.
        Similar: because it determines a flow type (stored and set by the user in the case of getPredefinedFlow).
        Different: because it is not associated with locked input options.
        This method is used to determine the default values in the flow modeling.
        For now this method is implemented only for two-phase flow.
        """
        if len(self.list_of_fields) != 2:
            return "unknown"
        # Discriminate between stratified / multiphase flow and dispersed flows
        field2 = self.list_of_fields[1]
        if field2.flow_type == "continuous":
            return "free_surface"
        else:
            if field2.phase == "solid":
                return "particles_flow"
            elif field2.phase == "liquid":
                return "droplet_flow"
            elif field2.phase == "gas":
                return "boiling_flow" # TODO rename into bubbly_flow (need to check backward compatibility issues)
        return "unknown"

    @Variables.noUndo
    def getPredefinedFlow(self):
        """
        Retrieve the predefined flow set by the user.
        """
        node = self.XMLNodethermo.xmlGetNode('predefined_flow')
        if node is None:
            node = self.setPredefinedFlow(self.defaultValues()['defaultPredefinedFlow'])

        return node['choice']

    @Variables.noUndo
    def setPhaseChangeTransferStatus(self, status):
        self.isOnOff(status)
        node = self.XMLNodethermo.xmlInitChildNode('phase_change_transfer')
        node.xmlSetAttribute(status=status)
        if status == "on":
            self._createSaturationProperties()
            self._createWallFieldsProperties()
        else:
            self._deleteFieldsProperties()

    @Variables.noUndo
    def getPhaseChangeTransferStatus(self):
        node = self.XMLNodethermo.xmlGetNode('phase_change_transfer')
        if node is None:
            self.setPhaseChangeTransferStatus(self.defaultValues()['phaseChangeTransfer'])
            node = self.XMLNodethermo.xmlGetNode('phase_change_transfer')
        return node["status"]

    @Variables.noUndo
    def getFieldId(self, label, include_none=False):
        """
        return field id  for a label
        """
        if label in ["none", "off", "all"]:
            return label
        self.isInList(label, self.getFieldLabelsList(include_none=include_none))
        # TODO : refactor if function is to be kept
        field = [fld for fld in self.list_of_fields if fld.label == label][0]
        return field.f_id


    def isFieldIdValid(self, fieldId, strict_check=True):
        """
        Check if fieldId corresponds to an existing field OR if it has
        an acceptable value ("none", "off", "all").
        If strict_check = True, fieldId has to match an existing field.
        """
        if strict_check:
            return self.isInList(str(fieldId), self.getFieldIdList())
        else:
            if fieldId in ["none", "all", "off"]:
                return True
            else:
                return self.isInList(str(fieldId), self.getFieldIdList())

    def getXMLNodefields(self):
        return self.__XMLNodefields


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class MainFieldsTestCase(ModelTest):
    """
    """
    def checkMainFieldsInstantiation(self):
        """Check whether the MainFieldsModel class could be instantiated"""
        model = None
        model = MainFieldsModel(self.case)
        assert model != None, 'Could not instantiate MainFieldsModel'


    def checkGetFieldIdList(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get FieldIdList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addField()
        assert mdl.getFieldIdList() == ['1', '2'],\
            'Could not get Field id list'


    def checkaddField(self):
        """Check whether the MainFieldsModel class could add field"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not add field'


    def checkaddDefinedField(self):
        """Check whether the MainFieldsModel class could add defined field"""
        mdl = MainFieldsModel(self.case)
        mdl.addDefinedField('1', 'field1', 'dispersed', 'gas', 'on', 'on', 'off')
        doc = '''<fields>
                         <field field_id="1" label="field1">
                                 <type choice="dispersed"/>
                                 <carrier_field field_id="on"/>
                                 <phase choice="gas"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not defined field'


    def checkGetFieldLabelsList(self):
        """Check whether the MainFieldsModel class could set and get FieldLabelsList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addField()
        assert mdl.getFieldLabelsList() == ['Field1', 'Field2'],\
            'Could not get Field label list'


    def checkGetContinuousFieldList(self):
        """Check whether the MainFieldsModel class could set and get ContinuousFieldList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off')
        assert mdl.getContinuousFieldList().f_id == ['1'],\
            'Could not get continuous field list'


    def checkGetDispersedFieldList(self):
        """Check whether the MainFieldsModel class could set and get DispersedFieldList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off')
        assert mdl.getDispersedFieldList().f_id == ['2'],\
            'Could not get dispersed field list'


    def checkGetGasPhaseList(self):
        """Check whether the MainFieldsModel class could set and get FieldIdList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off')
        assert mdl.getGasPhaseList()[0].f_id == '2',\
            'Could not get GasPhaseList'


    def checkGetSolidPhaseList(self):
        """Check whether the MainFieldsModel class could set and get SolidPhaseList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off')
        assert mdl.getSolidPhaseList()[0].f_id == '2',\
            'Could not get SolidPhaseList'


    def checkGetFirstContinuousField(self):
        """Check whether the MainFieldsModel class could set and get FirstContinuousField"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addField()
        assert mdl.getFirstContinuousField() == '1',\
            'Could not get FirstContinuousField'


    def checkGetFirstGasField(self):
        """Check whether the MainFieldsModel class could set and get FirstGasField"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off')
        assert mdl.getFirstGasField() == '2' ,\
            'Could not get FirstGasField'


    def checkGetandSetCriterion(self):
        """Check whether the MainFieldsModel class could set and get Criterion"""
        mdl = MainFieldsModel(self.case)
        field = mdl.addField()
        mdl.setCriterion('1','dispersed')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="dispersed"/>
                                 <carrier_field field_id="0"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set Criterion'
        assert field.flow_type == 'dispersed',\
            'Could not get Criterion'


    def checkGetandSetFieldNature(self):
        """Check whether the MainFieldsModel class could set and get FieldNature"""
        mdl = MainFieldsModel(self.case)
        field1 = mdl.addField()
        field1.phase = "gas"
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="gas"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set FieldNature'
        assert field1.phase == 'gas',\
            'Could not get FieldNature'

    def checkDeleteField(self):
        """Check whether the MainFieldsModel class could DeleteField"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.deleteField(0)
        doc = '''<fields/>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not delete field'


    def checkGetandSetPredefinedFlow(self):
        """Check whether the MainFieldsModel class could set and get PredefinedFlow"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.setPredefinedFlow('free_surface')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                                 <method choice="Thetis"/>
                                 <material choice="Water"/>
                                 <reference choice="WaterLiquid"/>
                         </field>
                         <field field_id="2" label="Field2">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="gas"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                                 <method choice="Thetis"/>
                                 <material choice="Water"/>
                                 <reference choice="WaterVapor"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set PredefinedFlow'
        assert mdl.getPredefinedFlow() == 'free_surface',\
            'Could not get PredefinedFlow'


    def checkGetFieldId(self):
        """Check whether the MainFieldsModel class could set and get field id"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        assert mdl.getFieldId('Field1') == '1',\
            'Could not get field id'


def suite():
    testSuite = unittest.makeSuite(MainFieldsTestCase, "check")
    return testSuite


def runTest():
    print("MainFieldsTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
