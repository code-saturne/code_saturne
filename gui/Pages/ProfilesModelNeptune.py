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
from code_saturne.Base.XMLvariablesNeptune import Variables
from code_saturne.Base.XMLmodelNeptune import ModelTest
from code_saturne.Base.XMLengine import *
from code_saturne.Base.Common import LABEL_LENGTH_MAX

#-------------------------------------------------------------------------------
# Constructor
#-------------------------------------------------------------------------------

class ProfilesModel(Variables, Model):
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

        self.__var_prop_list = self.getVariablesAndVolumeProperties()


    def __defaultValues(self):
        """
        Private method.
        Returns a dictionnary with default values.
        """
        value = {}
        value['nfreq']     = -1
        value['formula']   =  "x = 0;\ny = s;\nz = 0;\n"
        value['points']    =  200
        value['label']     =  "profile"
        value['choice']    =  "frequency"
        value['frequency'] =  1
        value['format']    =  "CSV"

        return value


    @Variables.noUndo
    def getVariablesAndVolumeProperties(self):
        """
        Creates a dictionnary to connect name and label from
        variables and properties.
        """
        self.dicoLabel2Name = {}
        # get all variable nodes include mean and constant
        nodeList = self.getVariablesPropertiesList('yes', 'yes', 'yes')

        for node in nodeList:

            name = node['name']
            label = node['label']
            fieldId = node['field_id']
            # for time average
            if name == None:
                name = label
            if fieldId == None:
                fieldId = 'none'
            if not label:
                raise ValueError("Node has no label")

            dim = node['dimension']
            if int(dim) > 1:
                for ii in range(int(dim)):
                    label1 = label + "[" + str(ii) + "]"
                    if not (node['support'] and node['support'] == "boundary"):
                        self.dicoLabel2Name[label1] = (name, fieldId, str(ii))
            else:
                if not (node['support'] and node['support'] == "boundary"):
                    if name != 'local_time_step':
                        self.dicoLabel2Name[label] = (name, fieldId, str(0))

        return list(self.dicoLabel2Name.keys())


    @Variables.undoGlobal
    def addProfile(self):
        """
        Public method.
        Add a new profile and return a default label
        """
        label = self.__defaultValues()['label']
        def_label = label + str(len(self.getProfilesLabelsList()) + 1)

        # define default label
        if def_label in self.getProfilesLabelsList():
            i = 2
            while def_label in self.getProfilesLabelsList():
                def_label = label + str(len(self.getProfilesLabelsList()) + i)
                i = i + 1

        node = self.node_prof.xmlInitNode('profile', label = def_label)

        return def_label


    @Variables.undoLocal
    def setFormula(self, label, str):
        """
        Public method.
        Get coordinates for profile named I{label}
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlSetData('formula', str)


    def getFormula(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        return node.xmlGetString('formula')


    @Variables.undoLocal
    def setNbPoint(self, label, NbPoint):
        """
        Public method.
        Get coordinates for profile named I{label}
        """
        self.isInt(NbPoint)
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlSetData('points', NbPoint)


    @Variables.noUndo
    def __getNbPoint(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        return node.xmlGetInt('points')


    @Variables.noUndo
    def getProfilesLabelsList(self):
        """
        Public method.
        Returns the profiles labels list.
        """
        lst = []
        for node in self.node_prof.xmlGetNodeList('profile'):
            label = node['label']
            lst.append(label)
        return lst


    @Variables.undoLocal
    def setOutputType(self, label, choice):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        self.isInList(choice, ["end", "frequency", "time_value"])
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlSetData('output_type', choice)


    @Variables.undoLocal
    def setOutputFrequency(self, label, freq):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlSetData('output_frequency', freq)


    @Variables.noUndo
    def getOutputFrequency(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        return node.xmlGetDouble('output_frequency')


    @Variables.undoLocal
    def setLabel(self, old_label, label):
        """
        Public method.
        """
        self.isInList(old_label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = old_label)
        node['label'] = label


    @Variables.undoLocal
    def setFormat(self, label, fmt):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlRemoveChild('format')
        node.xmlAddChild('format', name=fmt)


    @Variables.undoLocal
    def setVariable(self, label, lst):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        node.xmlRemoveChild('var_prop')
        for var in lst:
            self.isInList(var, self.__var_prop_list)
            (name, field, comp) = self.dicoLabel2Name[var]
            node.xmlAddChild('var_prop', field_id=field, name=name, component=comp)


    @Variables.noUndo
    def getVariable(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)

        lst = []
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['field_id'], var['component']) :
                    lst.append(name)
        return lst


    @Variables.undoLocal
    def deleteProfile(self, label):
        """
        Public method.
        Deletes profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label=label)
        if node:
            node.xmlRemoveNode()


    @Variables.noUndo
    def getProfileData(self, label):
        """
        Public method. Only for the GUI.
        Get profile named label and return list of variables or properties,
        frequency and coordinates.
        """
        self.isInList(label, self.getProfilesLabelsList())
        lst = []
        node = self.node_prof.xmlGetNode('profile', label = label)
        choice = node.xmlGetString('output_type')

        if not choice:
            choice = self.__defaultValues()['choice']
            self.setOutputType(label, choice)

        if choice == "time_value":
            freq = node.xmlGetDouble('output_frequency')
        else:
            freq = node.xmlGetInt('output_frequency')
        if not freq:
            freq = self.__defaultValues()['frequency']
            self.setOutputFrequency(label, freq)

        f_node = node.xmlGetChildNode('format')
        if not f_node:
            fmt = self.__defaultValues()['format']
            self.setFormat(label, fmt)
        else:
            fmt = f_node['name']

        formula = self.getFormula(label)
        if not formula:
            formula = self.__defaultValues()['formula']
            self.setFormula(label, formula)

        NbPoint = self.__getNbPoint(label)
        if not NbPoint:
            NbPoint = self.__defaultValues()['points']
            self.setNbPoint(label, NbPoint)

        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['field_id'], var['component']) :
                    lst.append(name)

        return label, fmt, lst, choice, freq, formula, NbPoint


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
        """Check whether the ProfilesModel class could set Profiles"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        doc = '''<profiles>
                         <profile label="prof1.dat" title="title">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <var_prop field_id="1" name="VelocityZ"/>
                                 <output_frequency>20</output_frequency>
                                <formula>
                                        x = 0;
y = 0;
z = 0;

                                </formula>
                         </profile>
                 </profiles>'''
        assert mdl.node_prof == self.xmlNodeFromString(doc),\
            'Could not set Profiles'


    def checkGetVariablesAndVolumeProperties(self):
        """Check whether the ProfilesModel class could getVariables and volume properties"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        assert mdl.getVariablesAndVolumeProperties() == ['Th_cond1', 'Lam_vis1', 'alpha1', 'V1', 'U1', 'enthalpy1', 'Pressure', 'Temp1', 'mass_trans1', 'W1', 'Sp_heat1', 'density1'],\
            'Could not get variables and volume properties'


    def checkGetProfilesLabelsList(self):
        """Check whether the ProfilesModel class could be get the profiles labels list"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        assert mdl.getProfilesLabelsList() == ['prof1.dat'],\
            'Could not get the profiles labels list'


    def checkReplaceProfile(self):
        """Check whether the ProfilesModel class could replace profile"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.replaceProfile('prof1.dat','profil1.dat', 'titre', ['U1'], 10, 0.01, 0.02, 0.03, 2.01, 2.02, 2.03)
        doc = '''<profiles>
                         <profile label="profil1.dat" title="titre">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <output_frequency>10</output_frequency>
                                 <x1>0.01</x1>
                                 <y1>0.02</y1>
                                 <z1>0.03</z1>
                                 <x2>2.01</x2>
                                 <y2>2.02</y2>
                                 <z2>2.03</z2>
                         </profile>
                 </profiles>'''

        assert mdl.node_prof == self.xmlNodeFromString(doc),\
            'Could not replace profile'


    def checkDeleteProfile(self):
        """Check whether the ProfilesModel class could delete profile"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2.dat', 'title2', ['U1', 'V1', 'W1'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        mdl.deleteProfile('prof1.dat')
        doc = '''<profiles>
                         <profile label="prof2.dat" title="title2">
                                 <var_prop field_id="1" name="VelocityX"/>
                                 <var_prop field_id="1" name="VelocityY"/>
                                 <var_prop field_id="1" name="VelocityZ"/>
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
            'Could not delete profile'


    def checkGetProfileData(self):
        """Check whether the ProfilesModel class could get profile data"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        assert mdl.getProfileData('prof1.dat') == ('prof1.dat', 'title', ['U1', 'V1', 'W1'], 20, 0.10000000000000001, 0.20000000000000001, 0.29999999999999999, 2.1000000000000001, 2.2000000000000002, 2.2999999999999998),\
            'Could not get the profile data'


def suite():
    testSuite = unittest.makeSuite(ProfilesTestCase, "check")
    return testSuite


def runTest():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
