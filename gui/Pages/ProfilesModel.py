# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module defines the 1D profile management page.

This module defines the following classes:
- ProfilesModel
- ProfilesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import XMLmodel, ModelTest
from Base.XMLvariables import Model
from Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class ProfilesModel(Model):
    """
    1D profile management.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_anal      = self.case.xmlInitNode('analysis_control')
        self.node_prof      = self.node_anal.xmlInitNode('profiles')
        self.node_model     = self.case.xmlInitNode('thermophysical_models')
        self.node_model_vp  = self.node_model.xmlInitNode('velocity_pressure')
        self.node_var_vp    = self.node_model_vp.xmlGetNodeList('variable')
        self.node_pro_vp    = self.node_model_vp.xmlGetNodeList('property')

        self.suffixe = self.__defaultValues()['suffixe']
        self.__var_prop_list = self.getVariablesAndVolumeProperties()


    def __defaultValues(self):
        """
        Private method.
        Returns a dictionnary with default values.
        """
        value = {}
        value['nfreq'] = -1
        value['X1']    =  0.
        value['Y1']    =  0.
        value['Z1']    =  0.
        value['X2']    =  1.
        value['Y2']    =  1.
        value['Z2']    =  1.
        value['suffixe'] =  ""

        return value


    def __updateBatchScriptFile(self, param, profile):
        """
        Update the backup file if it's ready to run.
        """
        key = self.case['computer']
        if key:
            if not self.case['batchScript'][key]:
                return

            from BatchRunningModel import BatchRunningModel
            batch = BatchRunningModel(self.case)
            batch.initializeBatchScriptFile()

            if batch.dicoValues['USER_OUTPUT_FILES']:
                vlist = string.split(batch.dicoValues['USER_OUTPUT_FILES'])
            else:
                vlist = []

            if param == "delete":
                if profile in vlist:
                    vlist.remove(profile)
            elif param == "add":
                if profile not in vlist:
                    vlist.append(profile)

            batch.dicoValues['USER_OUTPUT_FILES'] = string.join(vlist, " ")
            batch.updateBatchScriptFile('USER_OUTPUT_FILES')
            del BatchRunningModel


    def getVariablesAndVolumeProperties(self):
        """
        Creates a dictionnary to connect name and label from
        variables and properties.
        """
        # FIXME: merge this method with the same one in TimeAveragesView
        self.dicoLabel2Name = {}
        model = XMLmodel(self.case)
        output = OutputVolumicVariablesModel(self.case)
        for nodeList in [self.node_var_vp,
                         self.node_pro_vp,
                         model.getTurbNodeList(),
                         output.getFluidProperty(),
                         output.getAdditionalScalarProperty(),
                         output.getTimeProperty(),
                         output.getListOfTimeMeans(),
                         output.getPuCoalScalProper(),
                         output.getMeteoScalProper(),
                         output.getThermalScalar(),
                         output.getAdditionalScalar()]:

            for node in nodeList:
                name = node['name']
                label = node['label']
                if not label:
                    raise ValueError, "Node has no label"

                if not (node['support'] and node['support'] == "boundary"):
                    if name != 'local_time_step':
                        self.dicoLabel2Name[label] = name

        return self.dicoLabel2Name.keys()


    def __setCoordinates(self, label, x1, y1, z1, x2, y2, z2):
        """
        Private method.
        Get coordinates for profile named I{label}
        """
        self.isInList(label, self.getProfilesLabelsList())
        for coord in (x1, y1, z1, x2, y2, z2):
            self.isFloat(coord)

        label_xml = label + self.suffixe
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        node.xmlSetData('x1', x1)
        node.xmlSetData('y1', y1)
        node.xmlSetData('z1', z1)
        node.xmlSetData('x2', x2)
        node.xmlSetData('y2', y2)
        node.xmlSetData('z2', z2)


    def __getCoordinates(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffixe
        node = self.node_prof.xmlGetNode('profile', label=label_xml)

        x1 = node.xmlGetDouble('x1')
        if x1 == None:
            x1 = self.__defaultValues()['X1']
            node.xmlSetData('x1', x1)

        y1 = node.xmlGetDouble('y1')
        if y1 == None:
            y1 = self.__defaultValues()['Y1']
            node.xmlSetData('y1', y1)

        z1 = node.xmlGetDouble('z1')
        if z1 == None:
            z1 = self.__defaultValues()['Z1']
            node.xmlSetData('z1', z1)

        x2 = node.xmlGetDouble('x2')
        if x2 == None:
            x2 = self.__defaultValues()['X2']
            node.xmlSetData('x2', x2)

        y2 = node.xmlGetDouble('y2')
        if y2 == None:
            y2 = self.__defaultValues()['Y2']
            node.xmlSetData('y2', y2)

        z2 = node.xmlGetDouble('z2')
        if z2 == None:
            z2 = self.__defaultValues()['Z2']
            node.xmlSetData('z2', z2)

        return x1, y1, z1, x2, y2, z2


    def getProfilesLabelsList(self):
        """
        Public method.
        Returns the profiles labels list.
        """
        list = []
        for node in self.node_prof.xmlGetNodeList('profile'):
            #label = node['label'][:-4]
            label = node['label']
            list.append(label)
        return list


    def setProfile(self, label, list, freq, x1, y1, z1, x2, y2, z2):
        """
        Public method.
        Sets data to create one profile named I{label}.
        """
        self.isNotInList(label, self.getProfilesLabelsList())
        self.isInt(freq)
        for coord in (x1, y1, z1, x2, y2, z2):
            self.isFloat(coord)

        label_xml = label + self.suffixe
        node = self.node_prof.xmlInitNode('profile', label=label_xml)
        for var in list:
            self.isInList(var, self.__var_prop_list)
            node.xmlAddChild('var_prop', name=self.dicoLabel2Name[var])
        node.xmlSetData('output_frequency', freq)
        self.__setCoordinates(label, x1, y1, z1, x2, y2, z2)
        self.__updateBatchScriptFile('add', label_xml)


    def replaceProfile(self, old_label, label, list, freq, x1, y1, z1, x2, y2, z2):
        """
        Public method.
        Replaces data from I{old_label} profile
        with label, frequency, and coordinates values.
        """
        self.isInList(old_label,self.getProfilesLabelsList())
        if label != old_label:
            self.isNotInList(label, self.getProfilesLabelsList())
        self.isInt(freq)
        for coord in (x1, y1, z1, x2, y2, z2):
            self.isFloat(coord)

        old_label_xml = old_label + self.suffixe
        label_xml = label + self.suffixe
        node = self.node_prof.xmlGetNode('profile', label=old_label_xml)
        if node:
            node.xmlRemoveChild('var_prop')
            node.xmlRemoveChild('output_frequency')
            for tag in ('x1', 'y1', 'z1', 'x2', 'y2', 'z2'):
                node.xmlRemoveChild(tag)
        for var in list:
            self.isInList(var, self.__var_prop_list)
            node.xmlAddChild('var_prop', name=self.dicoLabel2Name[var])
        node['label'] = label_xml
        node.xmlSetData('output_frequency', freq)
        self.__setCoordinates(label, x1, y1, z1, x2, y2, z2)

        if old_label != label:
            self.__updateBatchScriptFile('delete', old_label_xml)
            self.__updateBatchScriptFile('add', label_xml)


    def deleteProfile(self, label):
        """
        Public method.
        Deletes profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffixe
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        if node:
            node.xmlRemoveNode()
            self.__updateBatchScriptFile('delete', label_xml)


    def getProfileData(self, label):
        """
        Public method. Only for the GUI.
        Get profile named label and return list of variables or properties,
        frequency and coordinates.
        """
        self.isInList(label, self.getProfilesLabelsList())
        list = []
        label_xml = label + self.suffixe
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        freq = node.xmlGetInt('output_frequency')
        x1, y1, z1, x2, y2, z2 = self.__getCoordinates(label)
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == var['name']:
                    list.append(name)
        label_xml = node['label']
        #label = label_xml[:-4]
        label = label_xml

        return label, list, freq, x1, y1, z1, x2, y2, z2


#-------------------------------------------------------------------------------
# TimeAveragesModel test case
#-------------------------------------------------------------------------------


class ProfilesTestCase(ModelTest):
    """
    Unittest.
    """
    def checkProfilesModelInstantiation(self):
        """Check whether the ProfilesModel class could be instantiated"""
        model = None
        model = ProfilesModel(self.case)
        assert model != None, 'Could not instantiate ProfilesModel'


    def checkSetProfile(self):
        """Check whether the ProfilesModel class could set one profile"""
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        doc = '''<profiles>
                    <profile label="prof1">
                            <var_prop name="velocity_U"/>
                            <var_prop name="velocity_V"/>
                            <var_prop name="velocity_W"/>
                            <output_frequency>20</output_frequency>
                            <x1>0.1</x1>
                            <y1>0.2</y1>
                            <z1>0.3</z1>
                            <x2>2.1</x2>
                            <y2>2.2</y2>
                            <z2>2.3</z2>
                    </profile>
                 </profiles>'''

        assert mdl.node_prof == self.xmlNodeFromString(doc),\
            'Could not set profiles in ProfilesModel'


    def checkReplaceProfile(self):
        """Check whether the ProfilesModel class could replace profiles"""
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.replaceProfile('prof1', 'premier', ['VelocitU'], 30, 0.1, 0.2, 0.3, 2.0, 2.0, 2.0)
        doc = '''<profiles>
                    <profile label="premier">
                            <var_prop name="velocity_U"/>
                            <output_frequency>30</output_frequency>
                            <x1>0.1</x1>
                            <y1>0.2</y1>
                            <z1>0.3</z1>
                            <x2>2.0</x2>
                            <y2>2.0</y2>
                            <z2>2.0</z2>
                    </profile>
                 </profiles>'''

        assert mdl.node_prof == self.xmlNodeFromString(doc),\
            'Could not replace profiles in ProfilesModel'


    def checkDeleteProfile(self):
        """Check whether the ProfilesModel class could delete profiles"""
        mdl = ProfilesModel(self.case)

        mdl.setProfile('prof1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        mdl.deleteProfile('prof1')
        doc = '''<profiles>
                    <profile label="prof2">
                            <var_prop name="velocity_U"/>
                            <output_frequency>20</output_frequency>
                            <x1>0.9</x1>
                            <y1>0.8</y1>
                            <z1>0.7</z1>
                            <x2>9.1</x2>
                            <y2>9.2</y2>
                            <z2>9.3</z2>
                    </profile>
                 </profiles>'''
        assert mdl.node_prof == self.xmlNodeFromString(doc),\
            'Could not delete profiles in ProfilesModel'


    def checkGetVarProfile(self):
        """Check whether the ProfilesModel class could get one profile"""
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        prof = ('prof2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)

        assert mdl.getProfileData('prof2') == prof,\
            'Could not get values for profile named label in ProfilesModel'


def suite():
    testSuite = unittest.makeSuite(ProfilesTestCase, "check")
    return testSuite


def runTest():
    print __file__
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
