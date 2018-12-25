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
from code_saturne.Pages.InterfacialForcesModel import InterfacialForcesModel


class InterfacialAreaModel(MainFieldsModel, Variables, Model):

    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case              = case
        self.XMLclosure        = self.case.xmlGetNode('closure_modeling')
        self.XMLAreaDiam       = self.XMLclosure.xmlInitNode('interfacial_area_diameter')
        self.XMLNodethermo     = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeVariable   = self.XMLNodethermo.xmlInitNode('variables')
        self.__AreaModel       = ['constant', 'interfacial_area_transport']
        self.__GasSourceTerm   = ['no_coalescence_no_fragmentation','wei_yao','kamp_colin','ruyer_seiler']
        self.__SolidSourceTerm = ['no_coalescence_no_fragmentation']
        self.__Solution        = ['uncoupled','coupled']


    def defaultValues(self):
        default = {}

        default['diameter']        = 0.001
        default['mindiam']         = 1.e-6
        default['maxdiam']         = 0.1
        default['coupling']        = "uncoupled"
        default['areamodel']       = "constant"
        default['sourcetermgas']   = "no_coalescence_no_fragmentation"
        default['sourcetermsolid'] = "no_coalescence_no_fragmentation"
        if InterfacialForcesModel(self.case).getBubblesForLIMStatus() == 'on':
            default['areamodel']     = "interfacial_area_transport"
            default['sourcetermgas'] = "ruyer_seiler"

        return default


    def getAreaModelList(self) :
        """
        """
        return self.__AreaModel


    def getSourceTermList(self, fieldId) :
        """
        """
        list = []
        if self.getFieldNature(fieldId) == "gas" :
            list = self.__GasSourceTerm
        else :
            list = self.__SolidSourceTerm
        return list


    def getSolutionMethodList(self) :
        """
        """
        return self.__Solution


    def getVariableAIList(self) :
        """
        """
        list = []
        for field in self.getDispersedFieldList() :
            if self.getAreaModel(field) != "constant" :
                list.append(field)
        return list


    @Variables.undoGlobal
    def setAreaModel(self, fieldId, model) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, self.getAreaModelList())

        oldmodel = None
        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        if node != None :
            oldmodel = node['model']

        node = self.XMLAreaDiam.xmlInitNode('field', field_id = fieldId)
        node.xmlSetAttribute(model = model)

        if model == "constant" :
            noden = node.xmlGetNode('source_term')
            if noden != None :
                noden.xmlRemoveNode()
            noden = node.xmlGetNode('solmeth')
            if noden != None :
                noden.xmlRemoveNode()
            if len(self.getVariableAIList() ) == 0 :
                noden = self.XMLAreaDiam.xmlGetNode('min_diameter')
                if noden != None :
                    noden.xmlRemoveNode()
                noden = self.XMLAreaDiam.xmlGetNode('max_diameter')
                if noden != None :
                    noden.xmlRemoveNode()

        if oldmodel != model :
            if model == "constant" :
                Variables(self.case).removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "Xd")
                Variables(self.case).removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "X2")
            else:
                Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldId, "Xd", "Xd_"+str(fieldId))


    @Variables.noUndo
    def getAreaModel(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        if node == None :
            model = self.defaultValues()['areamodel']
            self.setAreaModel(fieldId, model)
            node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        model = node['model']

        return model


    @Variables.undoGlobal
    def setSourceTerm(self, fieldId, model) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, self.getSourceTermList(fieldId))

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)

        oldmodel = None
        childNode = node.xmlGetNode('source_term')
        if childNode != None :
            oldmodel = childNode['model']

        childNode = node.xmlInitChildNode('source_term')
        childNode.xmlSetAttribute(model = model)

        if oldmodel != model :
            if model == 'kamp_colin' :
                Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeVariable, fieldId, "X2", "X2_"+str(fieldId))
            else:
                Variables(self.case).removeVariableProperty("variable", self.XMLNodeVariable, fieldId, "X2")


    @Variables.noUndo
    def getSourceTerm(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('source_term')
        if noden == None :
            model = ""
            if self.getFieldNature(fieldId) == "gas" :
               model = self.defaultValues()['sourcetermgas']
            else :
               model = self.defaultValues()['sourcetermsolid']
            self.setSourceTerm(fieldId, model)
        model = node.xmlGetNode('source_term')['model']

        return model


    @Variables.undoLocal
    def setSolutionMethod(self, fieldId, model) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, self.getSolutionMethodList())

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('solmeth')
        childNode.xmlSetAttribute(model = model)


    @Variables.noUndo
    def getSolutionMethod(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('solmeth')
        if noden == None :
            model = self.defaultValues()['coupling']
            self.setSolutionMethod(fieldId, model)
        model = node.xmlGetNode('solmeth')['model']

        return model


    @Variables.undoLocal
    def setInitialDiameter(self, fieldId, value) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isPositiveFloat(value)

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        node.xmlInitNode('diameter')
        node.xmlSetData('diameter', value)


    @Variables.noUndo
    def getInitialDiameter(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLAreaDiam.xmlGetNode('field', field_id = fieldId)
        value = node.xmlGetDouble('diameter')
        if value is None :
            self.setInitialDiameter(fieldId, self.defaultValues()['diameter'])
            value = node.xmlGetDouble('diameter')
        return value


    @Variables.undoLocal
    def setMinDiameter(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.XMLAreaDiam.xmlInitNode('min_diameter')
        self.XMLAreaDiam.xmlSetData('min_diameter', value)


    @Variables.noUndo
    def getMinDiameter(self) :
        """
        """
        value = self.XMLAreaDiam.xmlGetDouble('min_diameter')
        if value is None :
            self.setMinDiameter(self.defaultValues()['mindiam'])
            value = self.XMLAreaDiam.xmlGetDouble('min_diameter')
        return value


    @Variables.undoLocal
    def setMaxDiameter(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.XMLAreaDiam.xmlInitNode('max_diameter')
        self.XMLAreaDiam.xmlSetData('max_diameter', value)


    @Variables.noUndo
    def getMaxDiameter(self) :
        """
        """
        value = self.XMLAreaDiam.xmlGetDouble('max_diameter')
        if value is None :
            self.setMaxDiameter(self.defaultValues()['maxdiam'])
            value = self.XMLAreaDiam.xmlGetDouble('max_diameter')
        return value

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
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


    def checkGetSolutionMethodList(self):
        """Check whether the InterfacialAreaModel class could get the SolutionMethod list"""
        mdl = InterfacialAreaModel(self.case)
        assert mdl.getSolutionMethodList() == ['uncoupled', 'coupled'],\
            'Could not get SolutionMethodList'


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


    def checkGetandSetSolutionMethod(self):
        """Check whether the InterfacialAreaModel class could set and get SolutionMethod"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = InterfacialAreaModel(self.case)
        mdl.setAreaModel('1','interfacial_area_transport')
        mdl.setSolutionMethod('1','coupled')
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport">
                                 <solmeth model="coupled"/>
                         </field>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc),\
            'Could not set SolutionMethod'
        assert mdl.getSolutionMethod('1') == 'coupled',\
            'Could not set SolutionMethod'


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
        mdl.setAreaModel('1','interfacial_area_transport')
        mdl.setMaxDiameter(4.5)
        mdl.setMinDiameter(8.6)
        doc = '''<interfacial_area_diameter>
                         <field field_id="1" model="interfacial_area_transport"/>
                         <max_diameter>
                                 4.5
                         </max_diameter>
                         <min_diameter>
                                 8.6
                         </min_diameter>
                 </interfacial_area_diameter>'''
        assert mdl.XMLAreaDiam == self.xmlNodeFromString(doc),\
            'Could not set MinMaxDiameter'
        assert mdl.getMaxDiameter() == 4.5,\
            'Could not get MaxDiameter'
        assert mdl.getMinDiameter() == 8.6,\
            'Could not get MinDiameter'


def suite():
    testSuite = unittest.makeSuite(InterfacialAreaTestCase, "check")
    return testSuite


def runTest():
    print("InterfacialAreaTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
