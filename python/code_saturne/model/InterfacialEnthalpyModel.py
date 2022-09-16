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

# -------------------------------------------------------------------------------

import unittest
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import *

# -------------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------------

class InterfacialEnthalpyModel(Variables, Model):
    """
    This class manages the wall tranfer model objects in the XML file
    """

    def __init__(self, case):
        """
        Constructor.
        """
        #
        # XML file parameters
        self.mainFieldsModel = MainFieldsModel(case)
        self.case = case
        self.XMLClosure            = self.case.xmlGetNode('closure_modeling')
        self.XMLInterfEnthalpyNode = self.XMLClosure.xmlInitNode('interfacial_enthalpy')

        self.__availableLiquidBubble     = ['relaxation_time',
                                            'no_source_term',
                                            'bulk',
                                            'flashing',
                                            'bubble_model_for_liquid']

        self.__availableVaporBubble      = ['relaxation_time',
                                            'no_source_term',
                                            'relaxation_time_subcooled',
                                            'bubble_model_for_vapour']

        self.__availableLiquidDroplet    = ['relaxation_time',
                                            'no_source_term',
                                            'droplet_model_for_liquid']

        self.__availableVaporDroplet     = ['relaxation_time',
                                            'no_source_term',
                                            'bulk',
                                            'droplet_model_for_vapour']

        self.__availableLiquidContinuous = ['relaxation_time',
                                            'no_source_term',
                                            'wall_law_type_model']

        self.__availableVaporContinuous  = ['relaxation_time',
                                            'no_source_term',
                                            'sublayer_LI3C']

        # Init freeCouples for enthalpy : criterion checking !
        self.__liquidVaporCouples = []

        for field_a in self.mainFieldsModel.getContinuousFieldList() :
            if field_a.enthalpy_model != 'off' :
                for field_b in self.mainFieldsModel.getContinuousFieldList():
                    if field_a.phase == field_b.phase:
                        continue
                    if field_b.enthalpy_model != 'off' and field_b.f_id > field_a.f_id:
                        self.__liquidVaporCouples.append((field_a.f_id, field_b.f_id))
                for field_b in self.mainFieldsModel.getDispersedFieldList():
                    if field_a.phase == field_b.phase:
                        continue
                    if field_b.enthalpy_model != 'off' and field_b.phase != "solid":
                        self.__liquidVaporCouples.append((field_a.f_id, field_b.f_id))

    def getLiquidVaporCouples(self):
        """
        return list of free couples
        """
        return self.__liquidVaporCouples

    def getAvailableModels(self, fieldId):
        """
        Return available models for fieldId depending on the nature of field a and b
        """
        field_id_list = self.getEnthalpyCoupleFieldId()
        if field_id_list is None:
            return []
        fieldaId, fieldbId = field_id_list
        field_a = self.mainFieldsModel.getFieldFromId(fieldaId)
        field_b = self.mainFieldsModel.getFieldFromId(fieldbId)
        if field_a.phase == "liquid":
            if field_b.flow_type == "continuous":
                if fieldId == fieldaId:
                    return self.__availableLiquidContinuous
                else:
                    return self.__availableVaporContinuous
            else:
                if fieldId == fieldaId:
                    return self.__availableLiquidBubble
                else:
                    return self.__availableVaporBubble
        else:
            if field_b.flow_type == "continuous":
                if fieldId == fieldaId:
                    return self.__availableVaporContinuous
                else:
                    return self.__availableLiquidContinuous
            else:
                if fieldId == fieldaId:
                    return self.__availableVaporDroplet
                else:
                    return self.__availableLiquidDroplet

    def defaultValues(self):
        default = {}

        # Liquid/vapour
        default['continuousliquid'] = "relaxation_time"
        default['dispersedliquid']  = "relaxation_time"
        default['continuousgas']    = "relaxation_time"
        default['dispersedgas']     = "relaxation_time"
        default['solid']            = "none"
        default['relaxation_time']  = 0.01
        default['ponderation_coef'] = 'alp1'
        default['pool_boiling']     = "off"

        return default

    @Variables.undoGlobal
    def addLiquidVaporEnthalpyTransfer(self, field_id_a, field_id_b):
        """
        add a new interfacial enthalpy couple model
        """
        node = self.XMLInterfEnthalpyNode.xmlInitChildNode('enthalpy', field_id_a=field_id_a, field_id_b=field_id_b)

        for field_id in [field_id_a, field_id_b]:
            model = ""
            field = self.mainFieldsModel.getFieldFromId(field_id)
            if field.phase == "gas":
                if field.flow_type == "continuous":
                    model = self.defaultValues()['continuousgas']
                else:
                    model = self.defaultValues()['dispersedgas']
            elif field.phase == "liquid":
                if field.flow_type == "continuous":
                    model = self.defaultValues()['continuousliquid']
                else:
                    model = self.defaultValues()['dispersedliquid']
            else:
                raise ValueError("Field nature is neither gas nor liquid. No liquid-vapor enthalpy transfer added.")

            relaxation_time = self.defaultValues()['relaxation_time']
            ponderation_coef = self.defaultValues()['ponderation_coef']

            model_node = node.xmlInitChildNode('enthalpy_model', field_id=field_id)
            model_node['model'] = model
            if model == 'relaxation_time':
                ponderation_node = node.xmlInitChildNode('ponderation', field_id=field_id)
                ponderation_node['ponderation'] = ponderation_coef
                ponderation_node['relaxation'] = relaxation_time
        return

    @Variables.undoGlobal
    def deleteLiquidVaporEnthalpyTransfer(self):
        """
        suppress enthalpy interfacial couple
        """
        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        if node is not None:
            node.xmlRemoveNode()

        node_pool = self.XMLInterfEnthalpyNode.xmlGetNode('pool_boiling_model')
        if node_pool is not None:
            node_pool.xmlRemoveNode()

    def getEnthalpyCoupleFieldId(self):
        XMLNodesEnthal = self.XMLInterfEnthalpyNode.xmlGetNodeList('enthalpy')
        if len(XMLNodesEnthal) > 1:
            raise Exception("XML file contains more than one enthalpy transfer couple. "
                            "Correct your XML data before relaunching the GUI.")
        elif len(XMLNodesEnthal) == 1:
            fieldaId = XMLNodesEnthal[0]['field_id_a']
            fieldbId = XMLNodesEnthal[0]['field_id_b']
            return fieldaId, fieldbId
        else:
            return None

    def setEnthalpyCoupleFieldId(self, field_id_a=None, field_id_b=None):
        node = self.XMLInterfEnthalpyNode.xmlInitNode('enthalpy')
        old_id_a, old_id_b = node["field_id_a"], node["field_id_b"]
        children_a = [node.xmlGetNode(child, field_id=old_id_a) for child in ["enthalpy_model", "ponderation"]]
        children_b = [node.xmlGetNode(child, field_id=old_id_b) for child in ["enthalpy_model", "ponderation"]]
        if field_id_a is not None:
            node["field_id_a"] = field_id_a
            for child_node in children_a:
                if child_node is not None:
                    child_node["field_id"] = field_id_a
        if field_id_b is not None:
            node["field_id_b"] = field_id_b
            for child_node in children_b:
                if child_node is not None:
                    child_node["field_id"] = field_id_b
        return

    def __updatemodel(self, fieldaId, fieldbId):
        """
        update model after fieldIdChange
        """
        modela = ""
        modelb = ""
        field_a = self.mainFieldsModel.getFieldFromId(fieldaId)
        field_b = self.mainFieldsModel.getFieldFromId(fieldbId)

        if field_a.phase == "liquid":
            if field_b.flow_type == "continuous":
               modela = self.defaultValues()['continuousliquid']
               modelb = self.defaultValues()['continuousgas']
            else:
                modela = self.defaultValues()['continuousliquid']
                modelb = self.defaultValues()['dispersedgas']
        else:
            if field_b.flow_type == "continuous":
                modela = self.defaultValues()['continuousgas']
                modelb = self.defaultValues()['continuousliquid']
            else:
                modela = self.defaultValues()['continuousgas']
                modelb = self.defaultValues()['dispersedliquid']

        self.setFieldModel(fieldaId, modela)
        self.setFieldModel(fieldbId, modelb)

    @Variables.undoGlobal
    def setFieldModel(self, fieldId, choice):
        """
        """
        if self.getEnthalpyCoupleFieldId() is not None:
            fieldaId, fieldbId = self.getEnthalpyCoupleFieldId()
        self.isInList(str(fieldId), (fieldaId, fieldbId))
        self.isInList(choice, self.getAvailableModels(fieldId))

        node = self.XMLInterfEnthalpyNode.xmlInitNode('enthalpy')
        childNode = node.xmlInitNode('enthalpy_model', field_id=fieldId)
        oldchoice = childNode['model']

        childNode['model'] = choice

        if choice == 'relaxation_time':
            if oldchoice != choice:
                relaxation_time = self.defaultValues()['relaxation_time']
                ponderation_coef = self.defaultValues()['ponderation_coef']

                node.xmlInitChildNode('ponderation', field_id=fieldId, ponderation=ponderation_coef,
                                      relaxation = relaxation_time)
        else :
            childNode =  node.xmlGetNode('ponderation', field_id = fieldId)
            if childNode != None :
                childNode.xmlRemoveNode()

    @Variables.noUndo
    def getFieldModel(self, fieldId):
        """
        """
        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        if node:
            return node.xmlGetNode('enthalpy_model', field_id=fieldId)["model"]
        return None

    @Variables.undoLocal
    def setRelaxationTime(self, fieldId, value):
        """
        set relaxation time for fieldId
        """
        self.isFloat(value)
        if self.getEnthalpyCoupleFieldId() is not None:
            fieldaId, fieldbId = self.getEnthalpyCoupleFieldId()
        self.isInList(str(fieldId), (fieldaId, fieldbId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        childNode = node.xmlGetNode('ponderation', field_id=fieldId)

        childNode['relaxation'] = value

    @Variables.noUndo
    def getRelaxationTime(self, fieldId):
        """
        return relaxation time for fieldId
        """

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        childNode = node.xmlGetNode('ponderation', field_id=fieldId)

        return childNode['relaxation']

    @Variables.undoLocal
    def setPonderationCoef(self, fieldId, choice):
        """
        """
        self.isInList(choice, ('alp1', 'alp2', 'alp1_alp2'))
        if self.getEnthalpyCoupleFieldId() is not None:
            fieldaId, fieldbId = self.getEnthalpyCoupleFieldId()
        self.isInList(str(fieldId), (fieldaId, fieldbId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        childNode = node.xmlGetNode('ponderation', field_id=fieldId)

        childNode['ponderation'] = choice

    @Variables.noUndo
    def getPonderationCoef(self, fieldId):
        """
        """

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy')
        childNode = node.xmlGetNode('ponderation', field_id=fieldId)

        return childNode['ponderation']

    @Variables.undoLocal
    def setSolidEnergyTransfer(self, model):
        """
        set model for solid interfacial transfer
        """
        self.isInList(model, ('none', 'gas_particule'))
        n = self.XMLInterfEnthalpyNode.xmlInitNode('solid_enthalpy_transfer')
        n['model'] = model

    @Variables.noUndo
    def getSolidEnergyTransfer(self):
        """
        return model for solid interfacial transfer
        """
        value = self.defaultValues()['solid']
        n = self.XMLInterfEnthalpyNode.xmlGetNode('solid_enthalpy_transfer')
        if n != None:
            value = n['model']
        else:
            self.setSolidEnergyTransfer(value)

        return value

    @Variables.undoLocal
    def setSolidEnergyTransferStatus(self, fieldaId, fieldbId, status):
        """
        set model for solid interfacial transfer
        TODO check whether this method is useful (remove otherwise)
        """
        n = self.XMLInterfEnthalpyNode.xmlInitNode('enthalpy', field_id_a=fieldaId, field_id_b=fieldbId)
        for fieldId in [fieldaId, fieldbId]:
            model_node = n.xmlInitNode("enthalpy_model", field_id=fieldId)
            previous_model = model_node['model']
            if status == "off" and previous_model == "gas_particule":
                model_node['model'] = "none"
            if status == "on":
                model_node["model"] = "gas_particule"

    @Variables.undoLocal
    def setPoolBoiling(self, state):
        """
        Activate or deactivate pool boiling model.
        """
        n = self.XMLInterfEnthalpyNode.xmlInitNode('pool_boiling_model')
        if state == True:
            n['state'] = 'on'
        else:
            n['state'] = 'off'


    @Variables.noUndo
    def getPoolBoiling(self):
        value = self.defaultValues()['pool_boiling']
        n = self.XMLInterfEnthalpyNode.xmlGetNode('pool_boiling_model')
        if n != None:
            value = n['state']
        else:
            self.setPoolBoiling(value)

        return value

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class InterfacialEnthalpyTestCase(ModelTest):
    """
    """
    def checkInterfacialEnthalpyInstantiation(self):
        """Check whether the InterfacialEnthlalpyModel class could be instantiated"""
        model = None
        model = InterfacialEnthalpyModel(self.case)
        assert model != None, 'Could not instantiate InterfacialEnthalpyModel'

    def checkGetLiquidVaporCouples(self):
        """Check whether the InterfacialEnthalpyModel class could get FreeCouples"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getLiquidVaporCouples() == [('1', '2'), ('1', '3'), ('2', '3')], \
            'Could not get FreeCouples'

    def checkaddLiquidVaporEnthalpyTransfer(self):
        """Check whether the InterfacialEnthalpyModel class could addEnthalpyCouple"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not addEnthalpyCouple'


    def checkSetFielda(self):
        """Check whether the InterfacialEnthalpyModel class could SetFielda"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'continuous', 'liquid', 'on', 'on', 'off', 3)
        MainFieldsModel(self.case).addDefinedField('4', 'field4', 'continuous', 'gas', 'on', 'on', 'off', 4)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('3', '4')
        mdl.setEnthalpyCoupleFieldId(field_id_a='3')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="3" field_id_b="4">
                                 <enthalpy_model field_id="3" model="relaxation_time"/>
                                 <ponderation field_id="3" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="4" model="relaxation_time"/>
                                 <ponderation field_id="4" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not SetFielda'


    def checkSetFieldb(self):
        """Check whether the InterfacialEnthalpyModel class could SetFieldb"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'continuous', 'liquid', 'on', 'on', 'off', 3)
        MainFieldsModel(self.case).addDefinedField('4', 'field4', 'continuous', 'gas', 'on', 'on', 'off', 4)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '3')
        mdl.setEnthalpyCoupleFieldId(field_id_b='3')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="3">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="3" model="relaxation_time"/>
                                 <ponderation field_id="3" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not SetFieldb'


    def checkgetandsetFieldModel(self):
        """Check whether the InterfacialEnthalpyModel class could get and set FieldModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        mdl.setFieldModel('1', 'wall_law_type_model')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="wall_law_type_model"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not SetFieldModel'
        assert mdl.getFieldModel('1') == 'wall_law_type_model', \
            'Could not getFieldModel'


    def checkGetandSetRelaxationTime(self):
        """Check whether the InterfacialEnthalpyModel class could set and get RelaxationTime"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        mdl.setRelaxationTime('1', 123.25)
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="123.25"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not set RelaxationTime'
        assert mdl.getRelaxationTime('1') == '123.25', \
            'Could not get RelaxationTime'


    def checkGetandSetPonderationCoef(self):
        """Check whether the InterfacialEnthalpyModel class could set and get PonderationCoef"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        mdl.setPonderationCoef('1', 'alp2')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp2" relaxation="0.01"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not set PonderationCoef'
        assert mdl.getPonderationCoef('1') == 'alp2', \
            'Could not get PonderationCoef'

    def checkdeleteLiquidVaporEnthalpyTransfer(self):
        """Check whether the InterfacialEnthalpyModel class could deleteEnthalpyCouple"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        mdl.deleteLiquidVaporEnthalpyTransfer()
        doc = '''<interfacial_enthalpy/>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc), \
            'Could not deleteEnthalpyCouple'

    def checkSetSolidEnergyTransferStatus(self):
        """Check whether the InterfacialEnthalpyModel class could set SolidEnergyTransfer status"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addLiquidVaporEnthalpyTransfer('1', '2')
        mdl.setSolidEnergyTransfer('gas_particule')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not set SolidEnergyTransfer'
        assert mdl.getSolidEnergyTransfer() == 'gas_particule',\
            'Could not get SolidEnergyTransfer'


def suite():
    testSuite = unittest.makeSuite(InterfacialEnthalpyTestCase, "check")
    return testSuite


def runTest():
    print("InterfacialEnthalpyTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
