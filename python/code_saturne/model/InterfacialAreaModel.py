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

import sys, unittest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import *
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel


# TODO : try to include this model in "InterfaceForcesModel" directly ?
class InterfacialAreaModel(Variables, Model):

    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constructor.
        """
        #
        # XML file parameters
        self.mainFieldsModel = MainFieldsModel(case)
        self.case              = case
        self.XMLclosure        = self.case.xmlGetNode('closure_modeling')
        self.XMLAreaDiam       = self.XMLclosure.xmlInitNode('interfacial_area_diameter')
        self.XMLNodethermo     = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeVariable   = self.XMLNodethermo.xmlInitNode('variables')
        self.__AreaModel       = ['constant', 'interfacial_area_transport']
        self.__GasSourceTerm   = ['no_coalescence_no_fragmentation','wei_yao','kamp_colin','ruyer_seiler']
        self.__SolidSourceTerm = ['no_coalescence_no_fragmentation']


    def defaultValues(self):
        default = {}

        default['areamodel'] = "constant"
        default['diameter'] = 0.001
        default['sourceterm'] = "no_coalescence_no_fragmentation"
        default['mindiam'] = 1.e-6
        default['maxdiam'] = 0.1

        flow_type = self.mainFieldsModel.getPredefinedFlow()
        if flow_type == "None":
            flow_type = self.mainFieldsModel.detectFlowType()

        if flow_type in ["boiling_flow", "multiregime_flow"]:
            default['areamodel'] = "interfacial_area_transport"
            default['sourceterm'] = "ruyer_seiler"
        elif flow_type == "droplet_flow":
            default['areamodel'] = "interfacial_area_transport"
            default['sourceterm'] = "no_coalescence_no_fragmentation"

        return default

    def setDefaultParameters(self, field_id):
        """
        Reset all parameters to default when activating a predefined flow (non user)
        """
        predefined_flow = self.mainFieldsModel.getPredefinedFlow()
        if predefined_flow == "None":
            return

        default = self.defaultValues()
        self.setAreaModel(field_id, default['areamodel'])
        self.setInitialDiameter(field_id, default['diameter'])
        if default['areamodel'] == "constant":
            self.setSourceTerm(field_id, default['sourceterm'])
            self.setMinDiameter(field_id, default['mindiam'])
            self.setMaxDiameter(field_id, default['maxdiam'])

    def getAreaModelList(self) :
        """
        """
        return self.__AreaModel


    def getSourceTermList(self, fieldId) :
        """
        """
        lst = []
        field = self.mainFieldsModel.getFieldFromId(fieldId)
        if field.phase == "gas" :
            lst = self.__GasSourceTerm
        else :
            lst = self.__SolidSourceTerm
        return lst


    def getVariableAIList(self) :
        """
        """
        list = []
        for field in self.mainFieldsModel.getDispersedFieldList() :
            if self.getAreaModel(field.f_id) != "constant" :
                list.append(field)
        return list


    @Variables.undoGlobal
    def setAreaModel(self, fieldId, model) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isInList(model, self.getAreaModelList())

        oldmodel = None
        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        if node != None :
            oldmodel = node['model']

        node = self.XMLAreaDiam.xmlInitNode('field', field_id = fieldId)
        node.xmlSetAttribute(model = model)

        # Remove interfacial area parameters if needed
        if model == "constant" :
            for elt in ['source_term', 'solmeth', 'min_diameter', 'max_diameter']:
                noden = node.xmlGetChildNode(elt)
                if noden:
                    noden.xmlRemoveNode()

        if oldmodel != model :
            if model == "constant" :
                self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "Xd")
                self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "X2")
            else:
                field = self.mainFieldsModel.getFieldFromId(fieldId)
                field_name = field.label
                self.setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldId, "Xd", "Xd_"+field_name)


    @Variables.noUndo
    def getAreaModel(self, fieldId) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        cte_field_1 = self.mainFieldsModel.getContinuousFieldList()[0].f_id

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        if node is None :
            model = self.defaultValues()['areamodel']
            if int(fieldId) == int(cte_field_1):
                model = 'constant'
            self.setAreaModel(fieldId, model)
            node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        model = node['model']

        return model


    @Variables.undoGlobal
    def setSourceTerm(self, fieldId, model) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isInList(model, self.getSourceTermList(fieldId))

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)

        oldmodel = None
        childNode = node.xmlGetNode('source_term')
        if childNode != None :
            oldmodel = childNode['model']

        childNode = node.xmlInitChildNode('source_term')
        childNode.xmlSetAttribute(model = model)

        field = self.mainFieldsModel.getFieldFromId(fieldId)
        field_name = field.label
        if oldmodel != model :
            if model == 'kamp_colin' :
                self.setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldId, "X2", "X2_"+field_name)
            else:
                self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "X2")


    @Variables.noUndo
    def getSourceTerm(self, fieldId) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('source_term')
        if noden is None :
            model = ""
            model = self.defaultValues()['sourceterm']
            self.setSourceTerm(fieldId, model)
        model = node.xmlGetNode('source_term')['model']

        return model


    @Variables.undoLocal
    def setInitialDiameter(self, fieldId, value) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isPositiveFloat(value)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        node.xmlInitNode('diameter')
        node.xmlSetData('diameter', value)


    @Variables.noUndo
    def getInitialDiameter(self, fieldId) :
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        value = node.xmlGetDouble('diameter')
        if value is None :
            self.setInitialDiameter(fieldId, self.defaultValues()['diameter'])
            value = node.xmlGetDouble('diameter')
        return value


    @Variables.undoLocal
    def setMinDiameter(self, fieldId, value):
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isPositiveFloat(value)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
        node.xmlInitNode('min_diameter')
        node.xmlSetData('min_diameter', value)

    @Variables.noUndo
    def getMinDiameter(self, fieldId):
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
        value = node.xmlGetDouble('min_diameter')
        if value is None:
            self.setMinDiameter(fieldId, self.defaultValues()['mindiam'])
            node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
            value = node.xmlGetDouble('min_diameter')
        return value

    @Variables.undoLocal
    def setMaxDiameter(self, fieldId, value):
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isPositiveFloat(value)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
        node.xmlInitNode('max_diameter')
        node.xmlSetData('max_diameter', value)

    @Variables.noUndo
    def getMaxDiameter(self, fieldId):
        """
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
        value = node.xmlGetDouble('max_diameter')
        if value is None:
            self.setMaxDiameter(fieldId, self.defaultValues()['maxdiam'])
            node = self.XMLAreaDiam.xmlGetNode('field', field_id=fieldId)
            value = node.xmlGetDouble('max_diameter')
        return value

    def remove(self):
        self.XMLAreaDiam.xmlRemoveChildren()
        for fieldId in self.mainFieldsModel.getFieldIdList():
            self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "Xd")
            self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "X2")


# -------------------------------------------------------------------------------
# DefineUsersScalars test case
# -------------------------------------------------------------------------------
class InterfacialAreaTestCase(ModelTest):
    """
    """

    def checkInterfacialAreaInstantiation(self):
        """Check whether the InterfacialAreaModel class could be instantiated"""
        model = None
        model = InterfacialAreaModel(self.case)
        assert model != None, 'Could not instantiate InterfacialAreaModel'


    def checkGetAreaModelList(self):
        """Check whether the InterfacialAreaModel class could get the AreaModel list"""
        mdl = InterfacialAreaModel(self.case)
        assert mdl.getAreaModelList() == ['constant', 'interfacial_area_transport'],\
            'Could not get AreaModelList'


    def checkGetSourceTermList(self):
        """Check whether the InterfacialAreaModel class could get the AreaModel list"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field1', 'continuous', 'solid', 'on', 'on', 'off', 2)
        mdl = InterfacialAreaModel(self.case)
        assert mdl.getSourceTermList('1') == ['no_coalescence_no_fragmentation', 'wei_yao', 'kamp_colin', 'ruyer_seiler'],\
            'Could not get SourceTermList'
        assert mdl.getSourceTermList('2') == ['wei_yao', 'kamp_colin', 'ruyer_seiler'],\
            'Could not get SourceTermList'


    def checkGetandSetAreaModel(self):
        """Check whether the InterfacialAreaModel class could set and get AreaModel"""
        MainFieldsModel(self.case).addField()
        mdl = InterfacialAreaModel(self.case)
        mdl.setAreaModel('1','interfacial_area_transport')
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport"/>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc),\
            'Could not set AreaModel'
        assert mdl.getAreaModel('1') == 'interfacial_area_transport',\
            'Could not get AreaModel'


    def checkGetandSetSourceTerm(self):
        """Check whether the InterfacialAreaModel class could set and get SourceTerm"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialAreaModel(self.case)
        mdl.setAreaModel('1','interfacial_area_transport')
        mdl.setSourceTerm('1','wei_yao')
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport">
                                 <source_term model="wei_yao"/>
                         </field>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc),\
            'Could not set SourceTerm'
        assert mdl.getSourceTerm('1') == 'wei_yao',\
            'Could not get SourceTerm'


    def checkGetandSetInitialDiameter(self):
        """Check whether the InterfacialAreaModel class could set and get InitialDiameter"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialAreaModel(self.case)
        mdl.setAreaModel('1','interfacial_area_transport')
        mdl.setInitialDiameter('1',0.6)
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport">
                                 <diameter>
                                         0.6
                                 </diameter>
                         </field>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc),\
            'Could not set InitialDiameter'
        assert mdl.getInitialDiameter('1') == 0.6,\
            'Could not get InitialDiameter'


    def checkGetandSetMinMaxDiameter(self):
        """Check whether the InterfacialAreaModel class could set and get MinMaxDiameter"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialAreaModel(self.case)
        mdl.setAreaModel('1', 'interfacial_area_transport')
        mdl.setMaxDiameter("1", 4.5)
        mdl.setMinDiameter("1", 8.6)
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport"/>
                         <max_diameter>
                                 4.5
                         </max_diameter>
                         <min_diameter>
                                 8.6
                         </min_diameter>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc), \
            'Could not set MinMaxDiameter'
        assert mdl.getMaxDiameter("1") == 4.5, \
            'Could not get MaxDiameter'
        assert mdl.getMinDiameter("1") == 8.6, \
            'Could not get MinDiameter'


def suite():
    testSuite = unittest.makeSuite(InterfacialAreaTestCase, "check")
    return testSuite


def runTest():
    print("InterfacialAreaTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
