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

import sys, unittest, copy

from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.Common import LABEL_LENGTH_MAX
from code_saturne.model.ProfilesModel import ProfilesModel
from code_saturne.model.TimeAveragesModel import TimeAveragesModel
from code_saturne.model.GlobalNumericalParametersModel import GlobalNumericalParametersModel

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

EOS = 1
try:
   import eosAva
except:
   EOS = 0

class PredefinedFlowsModel:
    """
    This class manages the Field objects for a predefined flow in the XML file
    """

    fieldsCouple = ["None" ,
                    "free_surface",
                    "boiling_flow",
                    "droplet_flow",
                    "particles_flow"]

    fieldsCoupleProperties = {}
    fieldsCoupleProperties[fieldsCouple[0]] = ("","")
    fieldsCoupleProperties[fieldsCouple[1]] = (("continuous", "continuous"),
                                               ("liquid", "gas"))
    fieldsCoupleProperties[fieldsCouple[2]] = (("continuous", "dispersed"),
                                               ("liquid", "gas"))
    fieldsCoupleProperties[fieldsCouple[3]] = (("continuous", "dispersed"),
                                               ("gas", "liquid"))
    fieldsCoupleProperties[fieldsCouple[4]] = (("continuous", "dispersed"),
                                               ("gas", "solid"))


#-------------------------------------------------------------------------------
# Description of fields attribute
#-------------------------------------------------------------------------------

class FieldAttributesDescription:
    """
    """
    typeChoiceValues = ['continuous', 'dispersed', 'auto']
    phaseValues = ['liquid', 'gas', 'particle']


#-------------------------------------------------------------------------------
# Model for main fields
#-------------------------------------------------------------------------------

class MainFieldsModel(Variables, Model):
    """
    This class manages the Field objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
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
        if pressure_node == None:
            pressure_node = self.XMLNodethermo.xmlGetNode('variable',
                                                          name='Pressure')
        if pressure_node == None:
            Variables(self.case).setNewVariableProperty("variable",
                                                        "",
                                                        self.XMLNodeVariable,
                                                        "none",
                                                        "pressure",
                                                        "Pressure",
                                                        post = True)
        porosity_node = self.XMLNodethermo.xmlGetNode('property',
                                                      name='porosity')
        if porosity_node == None:
            Variables(self.case).setNewVariableProperty("property",
                                                        "",
                                                        self.XMLNodeproperty,
                                                        "none",
                                                        "porosity",
                                                        "porosity")


    def defaultValues(self):
        default = {}
        default['id']                         = ""
        default['label']                      = "defaultLabel"
        default['enthalpyResolutionStatus']   = "on"
        default['enthalpyResolutionModel']    = "total_enthalpy"
        default['typeChoice']                 = FieldAttributesDescription.typeChoiceValues[0]
        default['phase']                      = FieldAttributesDescription.phaseValues[0]
        default['carrierField']               = "off"
        default['compressibleStatus']         = "off"
        default['defaultPredefinedFlow']      = "None"

        return default


    def getFieldIdList(self):
        """
        Return the id field list
        """
        list = []

        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            list.append(node['field_id'])
        return list


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
            carrierfield = self.defaultValues()['carrierField']

            self.addDefinedField(fieldId,
                                 label,
                                 type,
                                 phase,
                                 carrierfield,
                                 hmodel,
                                 compressible,
                                 labNum)

        return fieldId


    @Variables.undoGlobal
    def addDefinedField(self, fieldId, label, typeChoice, phase, carrierField,
                        hmodel, compressible, labNum):
        """
        add field for predefined flow
        """

        self.isInList(phase,['solid','liquid','gas'])
        self.isInList(typeChoice,['dispersed','continuous'])

        # Check that the field does not allreay exist
        nf = self.__XMLNodefields.xmlGetNode('field', field_id=fieldId)
        if nf == None:
            self.__XMLNodefields.xmlInitChildNode('field', field_id = fieldId, label = label)
        else:
            nf['label'] = label

        self.setCriterion(fieldId, typeChoice)
        self.setFieldNature(fieldId, phase)
        self.setEnergyModel(fieldId, hmodel)
        self.setCompressibleStatus(fieldId, compressible)
        self.setCarrierField(fieldId, carrierField)
        self.iniVariableProperties(fieldId, labNum)


    @Variables.undoLocal
    def iniVariableProperties(self, fieldNumber, labNum):
        """
        add XML variable and properties
        """

        field_name = self.getFieldLabelsList()[int(fieldNumber)-1]

        Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldNumber, "volume_fraction", "vol_f_"+field_name, post = True)
        Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldNumber, "velocity", "U_"+field_name, dim='3', post = True)
        Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldNumber, "enthalpy", "enthalpy_"+field_name, post = True)

        # Physical properties are set by default to "constant" to avoid uninitialized states with the GUI
        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "density", "density_"+field_name)
        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "molecular_viscosity", "molecular_viscosity_"+field_name)
        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "specific_heat", "specific_heat_"+field_name)
        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "thermal_conductivity", "thermal_conductivity_"+field_name)

        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "mass_trans", "mass_trans_"+field_name)
        Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "wall_distance", "y_plus_"+field_name, support = "boundary")
        if self.getCompressibleStatus(fieldNumber) == "on":
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "drodp", "drodp_"+field_name)
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "drodh", "drodh_"+field_name)
        if self.getFieldNature(fieldNumber) == "solid":
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "emissivity", "emissivity_"+field_name)
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "elasticity", "elasticity_"+field_name)
        if self.getEnergyResolution(fieldNumber) == "on":
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "temperature", "temp_"+field_name, post = True)
        if self.getCriterion(fieldNumber) == "dispersed":
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "diameter", "diam_"+field_name)
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldNumber, "drift_component", "drift_component_"+field_name, dim='3')


    def getFieldLabelsList(self):
        """
        return list of label for field
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            list.append(node['label'])
        return list


    def getContinuousFieldList(self):
        """
        return list of id for field
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            childNode = node.xmlInitChildNode('type')
            if childNode['choice'] == "continuous" :
               list.append(node['field_id'])
        return list


    def getDispersedFieldList(self):
        """
        return list of id for field
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            childNode = node.xmlInitChildNode('type')
            if childNode['choice'] == "dispersed" :
               list.append(node['field_id'])
        return list


    def getGasPhaseList(self):
        """
        return list of id for field
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            childNode = node.xmlInitChildNode('phase')
            if childNode['choice'] == "gas" :
               list.append(node['field_id'])
        return list


    def getSolidFieldIdList(self):
        """
        return list of id for solid field
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            childNode = node.xmlInitChildNode('phase')
            if childNode['choice'] == "solid" :
               list.append(node['field_id'])
        return list


    def getFirstContinuousField(self) :
        """
        return id of first continuous field
        """
        id = 0
        if len(self.getContinuousFieldList()) > 0:
            id = self.getContinuousFieldList()[0]
        return id


    def getFirstGasField(self) :
        """
        return id of first continuous field
        """
        id = 0
        if len(self.getGasPhaseList()) > 0:
            id = self.getGasPhaseList()[0]
        return id


    @Variables.noUndo
    def getEnthalpyResolvedField(self) :
        """
        return list of id for field with enthalpy resolution
        """
        list = []
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            childNode = node.xmlInitChildNode('hresolution')
            if childNode['status'] == "on" :
               list.append(node['field_id'])
        return list


    @Variables.undoLocal
    def setLabel(self, fieldId, label):
        """
        Put label
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        old_label = label
        label_new = label[:LABEL_LENGTH_MAX]
        if label_new not in self.getFieldLabelsList():
            node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
            if node :
               old_label = node['label']
               node['label'] = label

        # Renaming of variables and properties after field label change
        for node in self.case.xmlGetNodeList('variable') \
                + self.case.xmlGetNodeList('property'):
            if node['field_id'] == str(fieldId):
                li = node['label'].rsplit(old_label, 1)
                node['label'] = label.join(li)


    @Variables.noUndo
    def getLabel(self, fieldId):
        """
        get label
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        label = ""
        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        if node:
            label = node['label']
        return label


    @Variables.undoLocal
    def setCriterion(self, fieldId, type):
        """
        Put type of field
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('type')
        oldtype = childNode['choice']
        childNode.xmlSetAttribute(choice = type)
        # update carrier field
        if type == "continuous" :
            self.setCarrierField(fieldId,self.defaultValues()['carrierField'])
        else :
            self.setCarrierField(fieldId,self.getFirstContinuousField())
            if oldtype != type :
               for id in self.getDispersedFieldList() :
                  if self.getCarrierField(id) == str(fieldId) :
                     self.setCarrierField(id,self.getFirstContinuousField())
        # if number of continuous field < 2 suppress coupling and momentum forces
        if len(self.getContinuousFieldList()) < 2 :
            node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
            if node :
                node.xmlRemoveNode()
            node = self.XMLforces.xmlGetNode('continuous_field_momentum_transfer')
            if node :
                node.xmlRemoveNode()

        # remove closure law linked to field for coherence of XML
        if oldtype != type :
           for node in self.__nodesWithFieldIDAttribute():
               try :
                   if (node['field_id_a'] == str(fieldId) or node['field_id_b'] == str(fieldId)):
                       node.xmlRemoveNode()
               except :
                   pass

        # TODO mettre en coherence pour les aires interf., tout ce qui est closure law a faire aussi pour la nature.
        # Activated if dispersed or second continuous phase of GLIM
        if self.getCriterion(fieldId) == "dispersed":
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "diameter", "diam_"+field_name)
           Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "drift_component", "drift_component_"+field_name, dim='3')
        else :
           Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "diameter")
           Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "drift_component")

        self.updateXML()


    @Variables.noUndo
    def getCriterion(self, fieldId):
        """
        get type of field
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodet = node.xmlGetNode('type')
        if nodet == None :
            type = self.defaultValues()['typeChoice']
            self.setCriterion(fieldId,type)
        type = node.xmlGetNode('type')['choice']
        return type


    @Variables.undoLocal
    def setFieldNature(self, fieldId, phase):
        """
        put nature of field
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(phase,('liquid','solid','gas'))

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('phase')
        if childNode != None:
            oldstatus = childNode['choice']
            if phase != oldstatus:
               if phase == "solid":
                  Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "emissivity", "emissivity_"+field_name)
                  Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "elasticity", "elasticity_"+field_name)
                  self.setCriterion(fieldId, "dispersed")
               else :
                  Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "emissivity")
                  Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "elasticity")

        childNode.xmlSetAttribute(choice = phase)
        self.updateXML()


    @Variables.noUndo
    def getFieldNature(self, fieldId):
        """
        get nature of field
        """
        #~ self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodep = node.xmlGetNode('phase')
        if nodep == None :
            phase = self.defaultValues()['phase']
            self.setFieldNature(fieldId,phase)
        phase = node.xmlGetNode('phase')['choice']
        return phase


    @Variables.undoLocal
    def setEnergyResolution(self, fieldId, status):
        """
        set status for energy resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isOnOff(status)

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('hresolution')
        if childNode != None:
            oldstatus = childNode['status']
            if status != oldstatus:
               if status == "on":
                  Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldId, "enthalpy", "enthalpy_"+field_name, post = True)
                  Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "temperature", "temp_"+field_name, post = True)
               else :
                  Variables(self.case).removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "enthalpy")
                  Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "temperature")
        childNode.xmlSetAttribute(status = status)


    @Variables.noUndo
    def getEnergyResolution(self, fieldId):
        """
        get status for energy resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodeh= node.xmlGetNode('hresolution')
        if nodeh == None :
            hres = self.defaultValues()['enthalpyResolutionStatus']
            self.setEnergyResolution(fieldId,hres)
        hres= node.xmlGetNode('hresolution')['status']
        return hres


    @Variables.undoLocal
    def setEnergyModel(self, fieldId, mdl):
        """
        set model for energy resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(mdl, ('off', 'total_enthalpy', 'specific_enthalpy'))

        if mdl != 'off':
            self.setEnergyResolution(fieldId, 'on')
        else:
            self.setEnergyResolution(fieldId, 'off')

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('hresolution')
        childNode.xmlSetAttribute(model = mdl)


    @Variables.noUndo
    def getEnergyModel(self, fieldId):
        """
        get model for energy resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodeh= node.xmlGetNode('hresolution')
        if nodeh == None :
            self.getEnergyResolution(fieldId)
            hmdl = self.defaultValues()['enthalpyResolutionModel']
            self.setEnergyModel(fieldId,hmdl)
        hres= node.xmlGetNode('hresolution')['model']
        return hres


    @Variables.undoLocal
    def setCompressibleStatus(self, fieldId, status):
        """
        set status for compressible resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isOnOff(status)

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('compressible')
        if childNode != None:
            oldstatus = childNode['status']
            if status != oldstatus:
               if status == "on":
                  Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "d_rho_d_P", "drho_dP_"+field_name)
                  Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, "d_rho_d_h", "drho_dh_"+field_name)
               else :
                  Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "d_rho_d_P")
                  Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, "d_rho_d_h")
        childNode.xmlSetAttribute(status = status)


    @Variables.noUndo
    def getCompressibleStatus(self, fieldId):
        """
        get status for compressible resolution
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodec = node.xmlGetNode('compressible')
        if nodec == None :
            compress = self.defaultValues()['compressibleStatus']
            self.setCompressibleStatus(fieldId,compress)
        compress = node.xmlGetNode('compressible')['status']
        return compress


    @Variables.undoLocal
    def setCarrierField(self, fieldId, carrierfield):
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('carrier_field')
        childNode.xmlSetAttribute(field_id = str(carrierfield))


    @Variables.noUndo
    def getCarrierField(self, fieldId):
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodec =  node.xmlGetNode('carrier_field')
        if nodec == None :
            carrier = self.defaultValues()['carrierField']
            self.setCarrierField(fieldId,carrier)
        carrier =  node.xmlGetNode('carrier_field')['field_id']
        return carrier


    @Variables.undoGlobal
    def deleteField(self, row):
        """
        delete a field in XML and update
        """
        fieldId = row + 1
        self.isInList(str(fieldId),self.getFieldIdList())

        #suppress profile
        for node in reversed(self.node_profile.xmlGetNodeList('profile')) :
            suppress = 0
            for child in node.xmlGetNodeList('var_prop') :
                if (child['field_id'] == str(fieldId)) :
                    suppress = 1
            if suppress == 1 :
                label = node['label']
                ProfilesModel(self.case).deleteProfile(label)

        #suppress average
        for node in reversed(self.node_average.xmlGetNodeList('time_average')) :
            suppress = 0
            for child in node.xmlGetNodeList('var_prop') :
                if (child['field_id'] == str(fieldId)) :
                    suppress = 1
                    break
            if suppress == 1 :
                label = node['label']
                TimeAveragesModel(self.case).deleteTimeAverage(label)

        for node in self.__nodesWithFieldIDAttribute():
            try :
               if (node['field_id'] == str(fieldId) or ( fieldId == 1 and node['field_id'] == "none")):
                    node.xmlRemoveNode()
               elif (node['field_id_a'] == str(fieldId) or node['field_id_b'] == str(fieldId)):
                    node.xmlRemoveNode()
            except :
               pass

        # Update 'field_id' attributes for other fields in XML file
        for node in self.__nodesWithFieldIDAttribute():
            try :
                if node['field_id'] != "none" and node['field_id'] != "off":
                    if int(node['field_id']) > fieldId :
                        node['field_id'] = str(int(node['field_id']) - 1)
                    else:
                       if node['field_id_a'] > str(fieldId):
                          node['field_id_a'] = str(int(node['field_id_a']) - 1)
                       if node['field_id_b'] > str(fieldId):
                          node['field_id_b'] = str(int(node['field_id_b']) - 1)
            except:
                pass

        # Update for field Id
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            nodec = node.xmlGetNode('carrier_field')
            if nodec != None :
                if str(nodec['field_id']) == str(fieldId) :
                    id = self.getFirstContinuousField()
                    if id > 0 :
                       nodec['field_id'] = id
                    else :
                       # only dispersed field -> change to continuous field
                       currentId = node['field_id']
                       self.setCriterion(currentId, "continuous")
                       self.setCarrierField(currentId, self.defaultValues()['carrierField'])

        self.updateXML()


    def updateXML(self):
        """
        method for update in case of suppress or change attribute
        """
        # suppress solid information
        if (len(self.getSolidFieldIdList()) < 1) :
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
        if (len(self.getDispersedFieldList()) < 1) :
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
    def setPredefinedFlow(self, value):    # <=> setPredefinedFlow
        """
        """

        self.isInList(value, PredefinedFlowsModel.fieldsCouple)

        node = self.XMLNodethermo.xmlInitChildNode('predefined_flow')
        node.xmlSetAttribute(choice = value)
        from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
        self.XMLNodeclosure = self.case.xmlGetNode('closure_modeling')
        self.XMLMassTrans   = self.XMLNodeclosure.xmlInitNode('mass_transfer_model')

        #Variables : neptune_cfd.core.XMLvariables.Variables
        Variables(self.case).setNewVariableProperty("property", "constant", self.XMLNodeproperty, "none", "surface_tension", "Surf_tens")

        if value != "None":  #avant: "None"

            # Create the first field
            fieldId      = "1"
            label        = "Field1"
            typeChoice   = PredefinedFlowsModel.fieldsCoupleProperties[value][0][0]
            phase        = PredefinedFlowsModel.fieldsCoupleProperties[value][1][0]
            if typeChoice == "dispersed" :
                carrierField = "2"
            else :
                carrierField = "off"
            hmodel       = "total_enthalpy"
            compressible = "off"
            self.case.undoStop()

            self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, hmodel, compressible, 1)

            if value != "particles_flow":
                # initialisation for thermodynamique material
                if EOS == 1 :
                    material = 'Water'
                    self.ava = eosAva.EosAvailable()
                    self.ava.setMethods(material)
                    self.case.undoStop()
                    ThermodynamicsModel(self.case).setMaterials(fieldId, material)

                    fls = self.ava.whichMethods()
                    if "Cathare" in fls:
                        ThermodynamicsModel(self.case).setMethod(fieldId, "Cathare")
                    else:
                        for fli in fls:
                            if fli != "Ovap" and fli != "Flica4" and fli != "StiffenedGas":
                                ThermodynamicsModel(self.case).setMethod(fieldId, fli)
                            break

                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_total_flux", "wall_total_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_liquid_total_flux", "wall_liquid_total_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_evaporation_flux", "wall_evaporation_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_quenching_flux", "wall_quenching_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_liquid_convective_flux", "wall_liquid_convective_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_steam_convective_flux", "wall_steam_convective_flux", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "boundary_temperature", "wall_temperature", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_liquid_temperature", "wall_liquid_temperature", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_oversaturation", "wall_oversaturation", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "unal_diameter", "unal_diameter", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_diameter_mesh_independancy", "wall_diameter_mesh_independancy", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_roughness", "wall_roughness", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "wall_dispersed_phase_mass_source_term", "wall_dispersed_phase_mass_source_term", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "boiling_criteria", "boiling_criteria", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "exchange_coefficient", "exchange_coefficient", support = "boundary")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "uninfluenced_part", "uninfluenced_part", support = "boundary")

            # Create the second field
            fieldId      = "2"
            label        = "Field2"
            typeChoice   = PredefinedFlowsModel.fieldsCoupleProperties[value][0][1]
            phase        = PredefinedFlowsModel.fieldsCoupleProperties[value][1][1]
            if typeChoice == "dispersed" :
               carrierField = "1"
            else :
               carrierField = "off"
            hmodel       = "total_enthalpy"
            compressible = "off"
            self.case.undoStop()
            self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, hmodel, compressible, 2)

            if value != "particles_flow":
                # initialisation for thermodynamique material
                if EOS == 1 :
                    self.case.undoStop()
                    ThermodynamicsModel(self.case).setMaterials(fieldId, material)

                    fls = self.ava.whichMethods()
                    if "Cathare" in fls:
                        ThermodynamicsModel(self.case).setMethod(fieldId, "Cathare")
                    else:
                        for fli in fls:
                            if fli != "Ovap" and fli != "Flica4" and fli != "StiffenedGas":
                                ThermodynamicsModel(self.case).setMethod(fieldId, fli)
                            break

                #TODO attention a modifier lors de la generalisation des modeles
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "SaturationTemperature", "TsatK", post = True)
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "SaturationEnthalpyLiquid", "Hsat_Liquid")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "SaturationEnthalpyGas", "Hsat_Gas")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Hsat_d_P_Liquid", "dHsat_dp_Liquid")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Hsat_d_P_Gas", "dHsat_dp_Gas")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "d_Tsat_d_P", "dTsat_dp")
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, "none", "LatentHeat", "Hlat")
            else :
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationTemperature")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationEnthalpyLiquid")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationEnthalpyGas")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Liquid")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Gas")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Tsat_d_P")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "LatentHeat")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_total_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_total_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_evaporation_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_quenching_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_convective_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_steam_convective_flux")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "boundary_temperature")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_temperature")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_oversaturation")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "unal_diameter")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_diameter_mesh_independancy")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_roughness")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_dispersed_phase_mass_source_term")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "boiling_criteria")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "exchange_coefficient")
                Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "uninfluenced_part")

            # modification du choix pour le predicteur de vitesse
            self.case.undoStop()
            if (value == "boiling_flow") :
                GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_bubble'])
            else :
                GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_std'])

            if value != "boiling_flow" and value != "free_surface":
                self.XMLMassTrans.xmlRemoveChild('nucleate_boiling')

        else :
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationTemperature")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationEnthalpyLiquid")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "SaturationEnthalpyGas")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Liquid")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Hsat_d_P_Gas")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "d_Tsat_d_P")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "LatentHeat")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_total_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_total_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_evaporation_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_quenching_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_convective_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_steam_convective_flux")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "boundary_temperature")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_liquid_temperature")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_oversaturation")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "unal_diameter")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_diameter_mesh_independancy")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_roughness")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "wall_dispersed_phase_mass_source_term")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "boiling_criteria")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "exchange_coefficient")
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, "none", "uninfluenced_part")

            GlobalNumericalParametersModel(self.case).setVelocityPredictorAlgo(GlobalNumericalParametersModel(self.case).defaultValues()['velocity_predictor_algorithm_std'])
            self.XMLMassTrans.xmlRemoveChild('nucleate_boiling')

            # create two field flow continuous/continuous by default
            # Create the first field
            fieldId      = "1"
            label        = "Field1"
            typeChoice   = "continuous"
            phase        = "liquid"
            carrierField = "off"
            hmodel       = "total_enthalpy"
            compressible = "off"
            self.case.undoStop()
            self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, hmodel, compressible, 1)
            # Create the second field
            fieldId      = "2"
            label        = "Field2"
            typeChoice   = "continuous"
            phase        = "gas"
            carrierField = "off"
            hmodel       = "total_enthalpy"
            compressible = "off"
            self.case.undoStop()
            self.addDefinedField(fieldId, label, typeChoice, phase, carrierField, hmodel, compressible, 2)
        self.case.undoStart()

        return node

    @Variables.noUndo
    def getPredefinedFlow(self) :
        """
        """
        node = self.XMLNodethermo.xmlGetNode('predefined_flow')
        if node == None:
            node = self.setPredefinedFlow(self.defaultValues()['defaultPredefinedFlow'])

        return node['choice']


    @Variables.noUndo
    def getFieldId(self, label):
        """
        return field id  for a label
        """
        self.isInList(label,self.getFieldLabelsList())

        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            if node['label'] == label :
               return node['field_id']


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
        mdl.addDefinedField('1', 'field1', 'dispersed', 'gas', 'on', 'on', 'off', 1)
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
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        assert mdl.getContinuousFieldList() == ['1'],\
            'Could not get continuous field list'


    def checkGetDispersedFieldList(self):
        """Check whether the MainFieldsModel class could set and get DispersedFieldList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        assert mdl.getDispersedFieldList() == ['2'],\
            'Could not get dispersed field list'


    def checkGetGasPhaseList(self):
        """Check whether the MainFieldsModel class could set and get FieldIdList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        assert mdl.getGasPhaseList() == ['2'],\
            'Could not get GasPhaseList'


    def checkGetSolidFieldIdList(self):
        """Check whether the MainFieldsModel class could set and get SolidFieldIdList"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        assert mdl.getSolidFieldIdList() == ['2'],\
            'Could not get SolidFieldIdList'


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
        mdl.addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        assert mdl.getFirstGasField() == '2' ,\
            'Could not get FirstGasField'


    def checkGetandSetLabel(self):
        """Check whether the MainFieldsModel class could set and get Label"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.setLabel('1','example_label')
        doc = '''<fields>
                         <field field_id="1" label="example_label">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set Label'
        assert mdl.getLabel('1') == 'example_label',\
            'Could not get Label'


    def checkGetandSetCriterion(self):
        """Check whether the MainFieldsModel class could set and get Criterion"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
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
        assert mdl.getCriterion('1') == 'dispersed',\
            'Could not get Criterion'


    def checkGetandSetFieldNature(self):
        """Check whether the MainFieldsModel class could set and get FieldNature"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.setFieldNature('1','gas')
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
        assert mdl.getFieldNature('1') == 'gas',\
            'Could not get FieldNature'


    def checkGetandSetEnergyResolution(self):
        """Check whether the MainFieldsModel class could set and get EnergyResolution"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.setEnergyResolution('1','off')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="off"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set EnergyResolution'
        assert mdl.getEnergyResolution('1') == 'off',\
            'Could not get EnergyResolution'


    def checkGetandSetCompressibleStatus(self):
        """Check whether the MainFieldsModel class could set and get CompressibleStatus"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.setCompressibleStatus('1','on')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="on"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set CompressibleStatus'
        assert mdl.getCompressibleStatus('1') == 'on',\
            'Could not get CompressibleStatus'


    def checkGetandSetCarrierField(self):
        """Check whether the MainFieldsModel class could set and get CarrierField"""
        mdl = MainFieldsModel(self.case)
        mdl.addField()
        mdl.addField()
        mdl.setFieldNature('2','gas')
        mdl.setCarrierField('2','1')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                         <field field_id="2" label="Field2">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="1"/>
                                 <phase choice="gas"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefields() == self.xmlNodeFromString(doc),\
            'Could not set CarrierField'
        assert mdl.getCarrierField('2') == '1',\
            'Could not get CarrierField'


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
