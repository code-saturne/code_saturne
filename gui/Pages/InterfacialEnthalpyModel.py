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
from code_saturne.Base.XMLvariables import Model
from code_saturne.Base.XMLengine import *
from code_saturne.Base.XMLmodel import *
from code_saturne.Pages.MainFieldsModel import *

#-------------------------------------------------------------------------------
# Constructor
#-------------------------------------------------------------------------------


class InterfacialEnthalpyModel(MainFieldsModel, Variables, Model):
    """
    This class manages the wall tranfer model objects in the XML file
    """


    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
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
        self.__freeCouples = []

        for fieldaId in self.getContinuousFieldList() :
            if self.getEnergyResolution(fieldaId) == 'on' :
                for fieldbId in self.getContinuousFieldList() :
                    if self.getEnergyResolution(fieldbId) == 'on' and fieldbId > fieldaId :
                        self.__freeCouples.append((fieldaId, fieldbId))
                for fieldbId in self.getDispersedFieldList() :
                    if self.getEnergyResolution(fieldbId) == 'on' and self.getFieldNature(fieldbId) != "solid" :
                        self.__freeCouples.append((fieldaId, fieldbId))

        XMLNodesEnthal = self.XMLInterfEnthalpyNode.xmlGetNodeList('enthalpy', 'field_id_a', 'field_id_b')
        for node in XMLNodesEnthal :
            # Get fields
            fieldaId = node['field_id_a']
            fieldbId = node['field_id_b']
            #
            # Update free couples
            try :
                self.__freeCouples.remove((fieldaId, fieldbId))
            except:
                pass


    def getFreeCouples(self) :
         """
         return list of free couples
         """
         return self.__freeCouples


    def getFieldIdaList (self, fieldaId) :
         """
         return list of free couples
         """
         list = []
         list.append(fieldaId)
         for fielda, fieldb in self.__freeCouples :
             if fielda not in list :
                 list.append(fielda)
         return list


    def getFieldIdbList (self, fieldaId) :
         """
         return list of free couples
         """
         list = []
         for fielda, fieldb in self.__freeCouples :
             if str(fielda) == str(fieldaId) :
                 list.append(fieldb)
         return list


    def getAvailableLiquidBubble(self) :
        """
        return list of available enthalpy model for liquid phase in case of bubble
        """
        return self.__availableLiquidBubble


    def getAvailableVaporBubble(self) :
        """
        return list of available enthalpy model for gas phase in case of bubble
        """
        return self.__availableVaporBubble


    def getAvailableLiquidDroplet(self) :
        """
        return list of available enthalpy model for liquid phase in case of droplets
        """
        return self.__availableLiquidDroplet


    def getAvailableVaporDroplet(self) :
        """
        return list of available enthalpy model for gas phase in case of droplets
        """
        return self.__availableVaporDroplet


    def getAvailableLiquidContinuous(self) :
        """
        return list of available enthalpy model for liquid phase in case of LIM
        """
        return self.__availableLiquidContinuous


    def getAvailableVaporContinuous(self) :
        """
        return list of available enthalpy model for gas phase in case of LIM
        """
        return self.__availableVaporContinuous


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


    def getEnthalpyCoupleList(self) :
        """
        return list of field couple for interfacial enthalpy
        """
        list = []
        for node in self.XMLInterfEnthalpyNode.xmlGetNodeList('enthalpy') :
            couple=[node['field_id_a'], node['field_id_b']]
            list.append(couple)
        return list


    @Variables.undoGlobal
    def addEnthalpyCouple(self) :
        """
        add a new interfacial enthalpy couple model
        """
        couple = []
        if len(self.getFreeCouples()) > 0 :
            couple = self.getFreeCouples()[0]
            node = self.XMLInterfEnthalpyNode.xmlInitChildNode('enthalpy', field_id_a = couple[0], field_id_b = couple[1])

            model = ""
            if self.getFieldNature(couple[0]) == "gas" :
                if self.getCriterion(couple[0]) == "continuous" :
                    model = self.defaultValues()['continuousgas']
                else :
                    model = self.defaultValues()['dispersedgas']
            else :
                if self.getCriterion(couple[0]) == "continuous" :
                    model = self.defaultValues()['continuousliquid']
                else :
                    model = self.defaultValues()['dispersedliquid']

            relaxation_time = self.defaultValues()['relaxation_time']
            ponderation_coef = self.defaultValues()['ponderation_coef']

            node.xmlInitChildNode('enthalpy_model', field_id = couple[0], model = model)
            if model == 'relaxation_time' :
                node.xmlInitChildNode('ponderation', field_id = couple[0], ponderation = ponderation_coef,
                                      relaxation = relaxation_time)

            model = ""
            if self.getFieldNature(couple[1]) == "gas" :
                if self.getCriterion == "continuous" :
                    model = self.defaultValues()['continuousgas']
                else :
                    model = self.defaultValues()['dispersedgas']
            else :
                if self.getCriterion == "continuous" :
                    model = self.defaultValues()['continuousliquid']
                else :
                    model = self.defaultValues()['dispersedliquid']

            relaxation_time = self.defaultValues()['relaxation_time']
            ponderation_coef = self.defaultValues()['ponderation_coef']

            node.xmlInitChildNode('enthalpy_model', field_id = couple[1], model = model)
            if model == 'relaxation_time' :
                node.xmlInitChildNode('ponderation', field_id = couple[1], ponderation = ponderation_coef,
                                      relaxation = relaxation_time)

            try :
                self.__freeCouples.remove((couple[0], couple[1]))
            except:
                pass
        return couple


    @Variables.undoGlobal
    def setFielda(self, oldfieldaId, fieldbId, fieldaId):
        """
        put field id for field a
        """
        self.isInList([oldfieldaId, fieldbId], self.getEnthalpyCoupleList())

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = oldfieldaId, field_id_b = fieldbId)
        node['field_id_a'] = fieldaId

        if oldfieldaId != fieldaId :
            childNode = node.xmlGetNode('enthalpy_model', field_id = oldfieldaId)
            if childNode != None :
                childNode['field_id'] = fieldaId
            childNode = node.xmlGetNode('ponderation', field_id = oldfieldaId)
            if childNode != None :
                childNode['field_id'] = fieldaId

            if (fieldaId, fieldbId) not in self.__freeCouples :
                new_fieldb = self.getFieldIdbList(fieldaId)[0]
                node['field_id_b'] = new_fieldb
                childNode = node.xmlGetNode('enthalpy_model', field_id = fieldbId)
                if childNode != None :
                    childNode['field_id'] = new_fieldb
                childNode = node.xmlGetNode('ponderation', field_id = fieldbId)
                if childNode != None :
                    childNode['field_id'] = new_fieldb
                if self.getFieldNature(fieldbId) != self.getFieldNature(new_fieldb) :
                    self.__updatemodel(fieldaId, new_fieldb)
                self.__freeCouples.remove((fieldaId, new_fieldb))
            else :
                self.__freeCouples.remove((fieldaId, fieldbId))
                self.__updatemodel(fieldaId, fieldbId)

            self.__freeCouples.append((oldfieldaId, fieldbId))


    @Variables.undoGlobal
    def setFieldb(self, fieldaId, oldfieldbId, fieldbId):
        """
        put field id for field b
        """
        self.isInList([fieldaId, oldfieldbId], self.getEnthalpyCoupleList())

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = oldfieldbId)
        node['field_id_b'] = fieldbId

        if oldfieldbId != fieldbId :
            self.__freeCouples.remove((fieldaId, fieldbId))
            self.__freeCouples.append((fieldaId, oldfieldbId))
            childNode = node.xmlGetNode('enthalpy_model', field_id = oldfieldbId)
            if childNode != None :
                childNode['field_id'] = fieldbId
            childNode = node.xmlGetNode('ponderation', field_id = oldfieldbId)
            if childNode != None :
                childNode['field_id'] = fieldbId
            self.__updatemodel(fieldaId, fieldbId)


    def __updatemodel(self, fieldaId, fieldbId) :
        """
        update model after fieldIdChange
        """
        modela = ""
        modelb = ""

        if self.getFieldNature(fieldaId) == "liquid" :
            if self.getCriterion(fieldbId) == "continuous" :
               modela = self.defaultValues()['continuousliquid']
               modelb = self.defaultValues()['continuousgas']
            else :
               modela = self.defaultValues()['continuousliquid']
               modelb = self.defaultValues()['dispersedgas']
        else :
            if self.getCriterion(fieldbId) == "continuous" :
               modela = self.defaultValues()['continuousgas']
               modelb = self.defaultValues()['continuousliquid']
            else :
               modela = self.defaultValues()['continuousgas']
               modelb = self.defaultValues()['dispersedliquid']

        self.setFieldModel(fieldaId, fieldbId, fieldaId, modela)
        self.setFieldModel(fieldaId, fieldbId, fieldbId, modelb)


    @Variables.undoGlobal
    def setFieldModel(self, fieldaId, fieldbId, fieldId, choice) :
        """
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))

        if self.getFieldNature(fieldaId) == "liquid" :
            if self.getCriterion(fieldbId) == "continuous" :
                if fieldId == fieldaId :
                   self.isInList(choice, self.getAvailableLiquidContinuous())
                else :
                   self.isInList(choice, self.getAvailableVaporContinuous())
            else :
                if fieldId == fieldaId :
                   self.isInList(choice, self.getAvailableLiquidBubble())
                else :
                   self.isInList(choice, self.getAvailableVaporBubble())
        else :
            if self.getCriterion(fieldbId) == "continuous" :
                if fieldId == fieldaId :
                   self.isInList(choice, self.getAvailableVaporContinuous())
                else :
                   self.isInList(choice, self.getAvailableLiquidContinuous())
            else :
                if fieldId == fieldaId :
                   self.isInList(choice, self.getAvailableVaporDroplet())
                else :
                   self.isInList(choice, self.getAvailableLiquidDroplet())

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('enthalpy_model', field_id = fieldId)
        oldchoice = childNode['model']

        childNode['model'] = choice

        if choice == 'relaxation_time' :
            if oldchoice != choice :
                relaxation_time = self.defaultValues()['relaxation_time']
                ponderation_coef = self.defaultValues()['ponderation_coef']

                node.xmlInitChildNode('ponderation', field_id = fieldId, ponderation = ponderation_coef,
                                      relaxation = relaxation_time)
        else :
            childNode =  node.xmlGetNode('ponderation', field_id = fieldId)
            if childNode != None :
                childNode.xmlRemoveNode()


    @Variables.noUndo
    def getFieldModel(self, fieldaId, fieldbId, fieldId) :
        """
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('enthalpy_model', field_id = fieldId)

        return childNode['model']


    @Variables.undoLocal
    def setRelaxationTime(self, fieldaId, fieldbId, fieldId, value) :
        """
        set relaxation time for fieldId
        """
        self.isFloat(value)
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('ponderation', field_id = fieldId)

        childNode['relaxation'] = value


    @Variables.noUndo
    def getRelaxationTime(self, fieldaId, fieldbId, fieldId) :
        """
        return relaxation time for fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('ponderation', field_id = fieldId)

        return childNode['relaxation']


    @Variables.undoLocal
    def setPonderationCoef(self, fieldaId, fieldbId, fieldId, choice) :
        """
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))
        self.isInList(choice, ('alp1','alp2','alp1_alp2'))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('ponderation', field_id = fieldId)

        childNode['ponderation'] = choice


    @Variables.noUndo
    def getPonderationCoef(self, fieldaId, fieldbId, fieldId) :
        """
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        self.isInList(str(fieldId),(fieldaId, fieldId))

        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('ponderation', field_id = fieldId)

        return childNode['ponderation']


    @Variables.undoGlobal
    def deleteEnthalpyCouple(self, fieldaId, fieldbId) :
        """
        suppress enthalpy interfacial couple
        """
        self.isInList([fieldaId, fieldbId], self.getEnthalpyCoupleList())
        node = self.XMLInterfEnthalpyNode.xmlGetNode('enthalpy', field_id_a = fieldaId, field_id_b = fieldbId)
        node.xmlRemoveNode()

        node_pool = self.XMLInterfEnthalpyNode.xmlGetNode('pool_boiling_model')
        node_pool.xmlRemoveNode()

        # update free couples
        self.__freeCouples.append((fieldaId, fieldbId))


    @Variables.undoLocal
    def setSolidEnergyTransfer(self, model) :
        """
        set model for solid interfacial transfer
        """
        self.isInList(model,('none','gas_particule'))
        n = self.XMLInterfEnthalpyNode.xmlInitNode('solid_enthalpy_transfer')
        n['model'] = model


    @Variables.noUndo
    def getSolidEnergyTransfer(self) :
        """
        return model for solid interfacial transfer
        """
        value = self.defaultValues()['solid']
        n = self.XMLInterfEnthalpyNode.xmlGetNode('solid_enthalpy_transfer')
        if n != None :
            value = n['model']
        else :
            self.setSolidEnergyTransfer(value)

        return value

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


    def checkGetFreeCouples(self):
        """Check whether the InterfacialEnthalpyModel class could get FreeCouples"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getFreeCouples() == [('1', '2'), ('1', '3'), ('2', '3')],\
            'Could not get FreeCouples'


    def checkGetFieldIdaList(self):
        """Check whether the InterfacialEnthalpyModel class could get FieldIdaList"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getFieldIdaList('1') == ['1', '2'],\
            'Could not get FieldIdaList'
        assert mdl.getFieldIdaList('2') == ['2', '1'],\
            'Could not get FieldIdaList'
        assert mdl.getFieldIdaList('3') == ['3', '1', '2'],\
            'Could not get FieldIdaList'


    def checkGetFieldIdbList(self):
        """Check whether the InterfacialEnthalpyModel class could get FieldIdbList"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getFieldIdbList('1') ==['2', '3'] ,\
            'Could not get FieldIdbList'
        assert mdl.getFieldIdbList('2') == ['3'] ,\
            'Could not get FieldIdbList'
        assert mdl.getFieldIdbList('3') == [],\
            'Could not get FieldIdbList'


    def checkGetAvailableLiquidBubble(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableLiquidBubble list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableLiquidBubble() == ['relaxation_time', 'no_source_term', 'user_function', 'bulk', 'flashing', 'bubble_model_for_liquid'],\
            'Could not get AvailableLiquidBubbleList'


    def checkGetAvailableVaporBubble(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableVaporBubble list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableVaporBubble() == ['relaxation_time', 'no_source_term', 'user_function', 'relaxation_time_subcooled', 'bubble_model_for_vapour'],\
            'Could not get AvailableVaporBubbleList'


    def checkGetAvailableLiquidDroplet(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableLiquidDroplet list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableLiquidDroplet() == ['relaxation_time', 'no_source_term', 'user_function'],\
            'Could not get AvailableLiquidDropletList'


    def checkGetAvailableVaporDroplet(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableVaporDroplet list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableVaporDroplet() == ['relaxation_time', 'no_source_term', 'user_function', 'bulk'],\
            'Could not get AvailableVaporDropletList'


    def checkGetAvailableLiquidContinuous(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableLiquidContinuous list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableLiquidContinuous() == ['relaxation_time', 'no_source_term', 'user_function', 'wall_law_type_model'],\
            'Could not get AvailableLiquidContinuousList'


    def checkGetAvailableVaporContinuous(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableVaporContinuous list"""
        mdl = InterfacialEnthalpyModel(self.case)
        assert mdl.getAvailableVaporContinuous() == ['relaxation_time', 'no_source_term', 'user_function', 'sublayer_LI3C'],\
            'Could not get AvailableVaporContinuousList'


    def checkGetEnthalpyCoupleList(self):
        """Check whether the InterfacialEnthalpyModel class could get the EnthalpyCoupleList"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addField()
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.addEnthalpyCouple()
        assert mdl.getEnthalpyCoupleList() == [['1', '2'], ['1', '3']],\
            'Could not get EnthalpyCoupleList'


    def checkaddEnthalpyCouple(self):
        """Check whether the InterfacialEnthalpyModel class could addEnthalpyCouple"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
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
        mdl.addEnthalpyCouple()
        mdl.setFielda('1','2','3')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="3" field_id_b="4">
                                 <enthalpy_model field_id="3" model="relaxation_time"/>
                                 <ponderation field_id="3" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="4" model="relaxation_time"/>
                                 <ponderation field_id="4" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not SetFielda'


    def checkSetFieldb(self):
        """Check whether the InterfacialEnthalpyModel class could SetFieldb"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'continuous', 'liquid', 'on', 'on', 'off', 3)
        MainFieldsModel(self.case).addDefinedField('4', 'field4', 'continuous', 'gas', 'on', 'on', 'off', 4)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.setFieldb('1','2','3')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="3">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="3" model="relaxation_time"/>
                                 <ponderation field_id="3" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not SetFieldb'


    def checkgetandsetFieldModel(self):
        """Check whether the InterfacialEnthalpyModel class could get and set FieldModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.setFieldModel('1','2','1','wall_law_type_model')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="wall_law_type_model"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not SetFieldModel'
        assert mdl.getFieldModel('1','2','1') == 'wall_law_type_model',\
            'Could not getFieldModel'


    def checkGetandSetRelaxationTime(self):
        """Check whether the InterfacialEnthalpyModel class could set and get RelaxationTime"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.setRelaxationTime('1','2','1',123.25)
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="123.25"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not set RelaxationTime'
        assert mdl.getRelaxationTime('1','2','1') == '123.25',\
            'Could not get RelaxationTime'


    def checkGetandSetPonderationCoef(self):
        """Check whether the InterfacialEnthalpyModel class could set and get PonderationCoef"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.setPonderationCoef('1','2','1','alp2')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp2" relaxation="0.01"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                 </interfacial_enthalpy>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not set PonderationCoef'
        assert mdl.getPonderationCoef('1','2','1') == 'alp2',\
            'Could not get PonderationCoef'


    def checkdeleteEnthalpyCouple(self):
        """Check whether the InterfacialEnthalpyModel class could deleteEnthalpyCouple"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.deleteEnthalpyCouple('1','2')
        doc = '''<interfacial_enthalpy/>'''
        assert mdl.XMLInterfEnthalpyNode == self.xmlNodeFromString(doc),\
            'Could not deleteEnthalpyCouple'


    def checkGetandSetSolidEnergyTransfer(self):
        """Check whether the InterfacialEnthalpyModel class could set and get SolidEnergyTransfer"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialEnthalpyModel(self.case)
        mdl.addEnthalpyCouple()
        mdl.setSolidEnergyTransfer('gas_particule')
        doc = '''<interfacial_enthalpy>
                         <enthalpy field_id_a="1" field_id_b="2">
                                 <enthalpy_model field_id="1" model="relaxation_time"/>
                                 <ponderation field_id="1" ponderation="alp1" relaxation="0.01"/>
                                 <enthalpy_model field_id="2" model="relaxation_time"/>
                                 <ponderation field_id="2" ponderation="alp1" relaxation="0.01"/>
                         </enthalpy>
                         <solid_enthalpy_transfer model="gas_particule"/>
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
