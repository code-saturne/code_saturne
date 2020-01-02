# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import *
from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModel


class SolidModel(TurbulenceModel):

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
        TurbulenceModel.__init__(self, case)
        self.case = case
        self.XMLThermo = self.case.xmlGetNode('thermophysical_models')
        self.__XMLNodefields = self.XMLNodethermo.xmlInitNode('fields')


    def defaultValues(self):
        default = TurbulenceModel.defaultValues(self)

        default['compaction']    = 0.64
        default['friction']      = "none"
        default['granular']      = "none"
        default['kinetic']       = "none"
        default['polydispersed'] = "off"
        return default


    @Variables.undoLocal
    def setFrictionModel(self, fieldId, model) :
        """
        Public method.
        Put model for friction.
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, ('none', 'pressure', 'fluxes'))

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('friction')
        childNode.xmlSetAttribute(model = model)


    @Variables.noUndo
    def getFrictionModel(self, fieldId) :
        """
        Public method.
        Return model for friction.
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('friction')
        if noden == None :
            model = self.defaultValues()['friction']
            self.setFrictionModel(fieldId, model)
        model = node.xmlGetNode('friction')['model']

        return model


    @Variables.undoLocal
    def setGranularModel(self, fieldId, model) :
        """
        Public method.
        Put model for granular.
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, ('none', 'pressure', 'fluxes'))

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('granular')
        childNode.xmlSetAttribute(model = model)


    @Variables.noUndo
    def getGranularModel(self, fieldId) :
        """
        Public method.
        Return model for granular.
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('granular')
        if noden == None :
            model = self.defaultValues()['granular']
            self.setGranularModel(fieldId, model)
        model = node.xmlGetNode('granular')['model']

        return model


    @Variables.undoLocal
    def setKineticModel(self, fieldId, model) :
        """
        Public method.
        Put model for kinetic.
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isInList(model, ('none', 'uncorrelate_collision', 'correlate_collision'))

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('kinetic')
        childNode.xmlSetAttribute(model = model)


    @Variables.noUndo
    def getKineticModel(self, fieldId) :
        """
        Public method.
        Return model for kinetic.
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        noden = node.xmlGetNode('kinetic')
        if noden == None :
            model = self.defaultValues()['kinetic']
            self.setKineticModel(fieldId, model)
        model = node.xmlGetNode('kinetic')['model']

        return model


    @Variables.undoLocal
    def setCompaction(self, value) :
        """
        Public method.
        Return values for compaction.
        """
        self.isPositiveFloat(value)

        self.XMLThermo.xmlSetData('solid_compaction', value)


    @Variables.noUndo
    def getCompaction(self) :
        """
        Public method.
        put values for compaction.
        """
        value = self.XMLThermo.xmlGetDouble('solid_compaction')
        if value == None :
           value = self.defaultValues()['compaction']
           self.setCompaction(value)
        return value


    @Variables.undoLocal
    def setCouplingStatus(self, fieldId, status):
        """
        put polydispersed coupling status
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isOnOff(status)

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        node.xmlSetData('polydispersed', status)


    @Variables.noUndo
    def getCouplingStatus(self, fieldId):
        """
        get polydispersed coupling status
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)

        ChildNode = node.xmlGetChildNode('polydispersed')
        if ChildNode == None:
           status = self.defaultValues()['polydispersed']
           self.setCouplingStatus(fieldId, status)
        status = node.xmlGetString('polydispersed')
        return status


    def getXMLNodefields(self):
        return self.__XMLNodefields

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class SolidTestCase(ModelTest):
    """
    """
    def checkSolidInstantiation(self):
        """Check whether the SolidModel class could be instantiated"""
        model = None
        model = SolidModel(self.case)
        assert model != None, 'Could not instantiate SolidModel'


    def checkGetandSetFrictionModel(self):
        """Check whether the SolidModel class could set and get FrictionModel"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = SolidModel(self.case)
        mdl.setFrictionModel('2','pressure')
        doc = '''<field field_id="2" label="field2">
                         <type choice="dispersed"/>
                         <carrier_field field_id="on"/>
                         <phase choice="solid"/>
                         <hresolution status="on"/>
                         <compressible status="off"/>
                         <friction model="pressure"/>
                 </field>'''
        assert mdl.getXMLNodefields().xmlGetNode('field', field_id = '2') == self.xmlNodeFromString(doc),\
            'Could not set FrictionModel'
        assert mdl.getFrictionModel('2') == 'pressure',\
            'Could not get FrictionModel'


    def checkGetandSetGranularModel(self):
        """Check whether the SolidModel class could set and get GranularModel"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = SolidModel(self.case)
        mdl.setGranularModel('2','pressure')
        doc = '''<field field_id="2" label="field2">
                         <type choice="dispersed"/>
                         <carrier_field field_id="on"/>
                         <phase choice="solid"/>
                         <hresolution status="on"/>
                         <compressible status="off"/>
                         <granular model="pressure"/>
                 </field>'''
        assert mdl.getXMLNodefields().xmlGetNode('field', field_id = '2') == self.xmlNodeFromString(doc),\
            'Could not set GranularModel'
        assert mdl.getGranularModel('2') == 'pressure',\
            'Could not get GranularModel'


    def checkGetandSetKineticModel(self):
        """Check whether the SolidModel class could set and get KineticModel"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = SolidModel(self.case)
        mdl.setKineticModel('2','correlate_collision')
        doc = '''<field field_id="2" label="field2">
                         <type choice="dispersed"/>
                         <carrier_field field_id="on"/>
                         <phase choice="solid"/>
                         <hresolution status="on"/>
                         <compressible status="off"/>
                         <kinetic model="correlate_collision"/>
                 </field>'''
        assert mdl.getXMLNodefields().xmlGetNode('field', field_id = '2') == self.xmlNodeFromString(doc),\
            'Could not set KineticModel'
        assert mdl.getKineticModel('2') == 'correlate_collision',\
            'Could not get KineticModel'


    def checkGetandSetCompaction(self):
        """Check whether the SolidModel class could set and get Compaction"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = SolidModel(self.case)
        mdl.setCompaction(0.6)
        doc = '''<solid_compaction>
                         0.6
                 </solid_compaction>'''
        assert mdl.XMLThermo.xmlGetNode('solid_compaction') == self.xmlNodeFromString(doc),\
            'Could not set Compaction'
        assert mdl.getCompaction() == 0.6,\
            'Could not get Compaction'


    def checkGetandSetCouplingStatus(self):
        """Check whether the SolidModel class could set and get CouplingStatus"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = SolidModel(self.case)
        mdl.setCouplingStatus('2','on')
        doc = '''<field field_id="2" label="field2">
                         <type choice="dispersed"/>
                         <carrier_field field_id="on"/>
                         <phase choice="solid"/>
                         <hresolution status="on"/>
                         <compressible status="off"/>
                         <polydispersed>
                                 on
                         </polydispersed>
                 </field>'''
        assert mdl.getXMLNodefields().xmlGetNode('field', field_id = '2') == self.xmlNodeFromString(doc),\
            'Could not set CouplingStatus'
        assert mdl.getCouplingStatus('2') == 'on',\
            'Could not get CouplingStatus'


def suite():
    testSuite = unittest.makeSuite(SolidTestCase, "check")
    return testSuite


def runTest():
    print("SolidTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
