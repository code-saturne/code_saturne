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
from code_saturne.model.Common import LABEL_LENGTH_MAX
from code_saturne.model.ProfilesModel import ProfilesModel
from code_saturne.model.TimeAveragesModel import TimeAveragesModel

class SpeciesModel(MainFieldsModel, Variables, Model):

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
        self.case            = case
        self.XMLUserScalar   = self.case.xmlGetNode('additional_scalars')
        self.XMLUser         = self.XMLUserScalar.xmlInitNode('users')
        self.XMLScalar       = self.XMLUserScalar.xmlInitNode('scalars')
        self.node_anal       = self.case.xmlInitNode('analysis_control')
        self.node_average    = self.node_anal.xmlInitNode('time_averages')
        self.node_profile    = self.node_anal.xmlInitNode('profiles')


    def defaultValues(self):
        default = {}

        default['carrierField']     = 1
        default['diffusion_coefficient']    = 1e-5
        default['schmidt']          = 0.9
        default['minValue']         = 0
        default['maxValue']         = 1
        default['timeDependStatus'] = "on"
        default['diffusionStatus']  = "on"
        default['coupledModel']     = "coupled"
        return default


    def getScalarLabelList(self):
        """
        return list of label for scalar
        """
        scalarList = []
        for node in self.XMLScalar.xmlGetNodeList('variable'):
            scalarList.append(node['label'])
        return scalarList


    def getScalarNameList(self):
        """
        return list of name for scalar (scalar + i)
        """
        scalarList = []
        for node in self.XMLScalar.xmlGetNodeList('variable'):
            scalarList.append(node['name'])
        return scalarList


    def getScalarByFieldId(self, FieldId):
        """
        Return the scalar name list for a fieldId
        """
        self.isInList(str(FieldId),self.getFieldIdList())
        list = []
        for node in self.XMLScalar.xmlGetNodeList('variable'):
            if self.getScalarFieldIdByName(node['name']) == str(FieldId):
                list.append(node['name'])
        return list


    @Variables.undoGlobal
    def addScalar(self, scalarId=None):
        """
        add a new scalar
        """
        name = "scalar_" + str(scalarId)
        label = "species" + str(scalarId)
        if label in self.getScalarLabelList():
           labelNumber = 1
           label = "species" + str(labelNumber)
           while label in self.getScalarLabelList():
              labelNumber += 1
              label = "species" + str(labelNumber)

        carrierfield = self.defaultValues()['carrierField']

        Variables(self.case).setNewVariableProperty("variable", "", self.XMLScalar, carrierfield, name, label)

        return label


    @Variables.undoGlobal
    def deleteScalar(self, scalarId):
        """
        delete a scalar
        """
        # Update scalar Id
        name = "scalar_" + str(scalarId)
        self.isInList(name,self.getScalarNameList())

        # delete scalar
        for node in self.XMLScalar.xmlGetNodeList('variable'):
            try:
               if node['name'] == name:
                  node.xmlRemoveNode()
            except:
               pass

        #suppress profile
        for node in reversed(self.node_profile.xmlGetNodeList('profile')):
            suppress = 0
            for child in node.xmlGetNodeList('var_prop'):
                if (child['name'] == name):
                    suppress = 1
            if suppress == 1:
                name = node['name']
                ProfilesModel(self.case).deleteProfile(name)

        #suppress average
        for node in reversed(self.node_average.xmlGetNodeList('time_average')):
            suppress = 0
            for child in node.xmlGetNodeList('var_prop'):
                if (child['name'] == name):
                    suppress = 1
                    break
            if suppress == 1:
                name = node['name']
                TimeAveragesModel(self.case).deleteTimeAverage(name)

        # update name for other scalar in XML file
        index = 1
        for node in self.XMLScalar.xmlGetNodeList('variable'):
            try:
               if index >= scalarId:
                  oldname = "scalar_" + str(index + 1)
                  name = "scalar_" + str(index)
                  node['name'] = name
                  for n in self.case.xmlGetNodeList('var_prop', name=oldname):
                      n['name'] = name
            except:
               pass
            index += 1


    @Variables.undoLocal
    def setScalarLabel(self, scalarId, label):
        """
        Put label
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name,self.getScalarNameList())
        label_new = label[:LABEL_LENGTH_MAX]
        old_label = ''
        if label_new not in self.getScalarLabelList():
            node = self.XMLScalar.xmlGetNode('variable', name = name)
            if node:
                old_label = node['label']
                node['label'] = label_new
        # udpate label in formula if needed
        if old_label != '':
            for no in self.case.xmlGetNodeList('formula'):
                txt = no.xmlGetTextNode()
                if txt != None:
                    f = txt.replace(old_label, label_new)
                    no.xmlSetTextNode(f)


    @Variables.noUndo
    def getScalarLabel(self, scalarId):
        """
        Get label
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name,self.getScalarNameList())
        label = ""
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            label = node['label']
        return label


    @Variables.noUndo
    def getScalarLabelByName(self, name):
        """
        Get label
        """
        label = ""
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            label = node['label']
        return label


    @Variables.undoLocal
    def setScalarFieldId(self, scalarId, carrierfield):
        """
        put carrier field id for scalar
        """
        self.isInList(str(carrierfield),self.getFieldIdList())

        name = "scalar_" + str(scalarId)
        self.isInList(name,self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            node['field_id'] = carrierfield


    @Variables.noUndo
    def getScalarFieldId(self, scalarId):
        """
        get carrier field id for scalar
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        carrierField = self.defaultValues()['carrierField']
        if node:
            carrierField = node['field_id']
        return carrierField


    @Variables.noUndo
    def getScalarFieldIdByName(self, name):
        """
        get carrier field id for scalar
        """
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        carrierField = self.defaultValues()['carrierField']
        if node:
            carrierField = node['field_id']
        return carrierField


    @Variables.undoLocal
    def setDiffusionCoef(self, scalarId, value):
        """
        put diffusion coefficient for scalar
        """
        self.isFloat(value)

        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            node.xmlSetData('diffusion_coefficient', value)


    @Variables.noUndo
    def getDiffusionCoef(self, scalarId):
        """
        get diffusion coefficient for scalar
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('diffusion_coefficient')
        if ChildNode == None:
            value = self.defaultValues()['diffusion_coefficient']
            self.setDiffusionCoef(scalarId, value)
        value = node.xmlGetDouble('diffusion_coefficient')
        return value


    @Variables.undoLocal
    def setSchmidt(self, scalarId, value):
        """
        put schmidt for scalar
        """
        self.isFloat(value)

        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            node.xmlSetData('schmidt', value)


    @Variables.noUndo
    def getSchmidt(self, scalarId):
        """
        get schmidt for scalar
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('schmidt')
        if ChildNode == None:
            value = self.defaultValues()['schmidt']
            self.setSchmidt(scalarId, value)
        value = node.xmlGetDouble('schmidt')
        return value


    @Variables.undoLocal
    def setMinValue(self, scalarId, value):
        """
        put min value for scalar
        """
        self.isFloat(value)

        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            node.xmlSetData('minValue', value)


    @Variables.noUndo
    def getMinValue(self, scalarId):
        """
        get min value for scalar
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('minValue')
        if ChildNode == None:
            value = self.defaultValues()['minValue']
            self.setMinValue(scalarId, value)
        value = node.xmlGetDouble('minValue')
        return value


    @Variables.undoLocal
    def setMaxValue(self, scalarId, value):
        """
        put max value for scalar
        """
        self.isFloat(value)

        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        if node:
            node.xmlSetData('maxValue', value)


    @Variables.noUndo
    def getMaxValue(self, scalarId):
        """
        get max value for scalar
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('maxValue')
        if ChildNode == None:
            value = self.defaultValues()['maxValue']
            self.setMaxValue(scalarId, value)
        value = node.xmlGetDouble('maxValue')
        return value


    @Variables.undoLocal
    def setTimeDependStatus(self, scalarId, status):
        """
        put time depend status
        """
        self.isOnOff(status)
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        node.xmlSetData('timeDependStatus', status)


    @Variables.noUndo
    def getTimeDependStatus(self, scalarId):
        """
        get thickness status
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('timeDependStatus')
        if ChildNode == None:
           status = self.defaultValues()['timeDependStatus']
           self.setTimeDependStatus(scalarId, status)
        status = node.xmlGetString('timeDependStatus')
        return status


    @Variables.undoLocal
    def setDiffusionStatus(self, scalarId, status):
        """
        put diffusion resolution status
        """
        self.isOnOff(status)
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        node.xmlSetData('diffusionStatus', status)


    @Variables.noUndo
    def getDiffusionStatus(self, scalarId):
        """
        get diffusion resolution status
        """
        name = "scalar_" + str(scalarId)
        self.isInList(name, self.getScalarNameList())
        node = self.XMLScalar.xmlGetNode('variable', name = name)
        ChildNode = node.xmlGetChildNode('diffusionStatus')
        if ChildNode == None:
           status = self.defaultValues()['diffusionStatus']
           self.setDiffusionStatus(scalarId, status)
        status = node.xmlGetString('diffusionStatus')
        return status


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class SpeciesTestCase(ModelTest):
    """
    """
    def checkSpeciesInstantiation(self):
        """Check whether the SpeciesModel class could be instantiated"""
        model = None
        model = SpeciesModel(self.case)
        assert model != None, 'Could not instantiate SpeciesModel'


def suite():
    testSuite = unittest.makeSuite(SpeciesTestCase, "check")
    return testSuite


def runTest():
    print("SpeciesTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

