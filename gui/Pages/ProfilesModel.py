# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

from code_saturne.Base.Common import *
import Base.Toolbox as Tool
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel

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

        self.suffix = self.__defaultValues()['suffix']
        self.__var_prop_list = self.getVariablesAndVolumeProperties()


    def __defaultValues(self):
        """
        Private method.
        Returns a dictionnary with default values.
        """
        value = {}
        value['nfreq']  = -1
        value['formula']=  " "
        value['points']=  200
        value['suffix'] = ""
        value['title']  = ""

        return value


    @Variables.noUndo
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
                         output.getListOfTimeAverage(),
                         output.getPuCoalScalProper(),
                         output.getGasCombScalProper(),
                         output.getMeteoScalProper(),
                         output.getElecScalProper(),
                         output.getThermalScalar(),
                         output.getAdditionalScalar()]:

            for node in nodeList:
                name = node['name']
                label = node['label']
                if not label:
                    raise ValueError("Node has no label")

                dim = node['dimension']
                if dim and int(dim) > 1:
                    for ii in range(int(dim)):
                        label1 = label + "[" + str(ii) + "]"
                        if not (node['support'] and node['support'] == "boundary"):
                            self.dicoLabel2Name[label1] = (name, str(ii))
                else:
                    if not (node['support'] and node['support'] == "boundary"):
                        if name != 'local_time_step':
                            self.dicoLabel2Name[label] = (name, str(0))

        return list(self.dicoLabel2Name.keys())


    def __setFormula(self, label, str):
        """
        Private method.
        Get coordinates for profile named I{label}
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        node.xmlSetData('formula', str)


    def __getFormula(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        return node.xmlGetString('formula')


    def __setNbPoint(self, label, NbPoint):
        """
        Private method.
        Get coordinates for profile named I{label}
        """
        self.isInt(NbPoint)
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        node.xmlSetData('points', NbPoint)


    def __getNbPoint(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
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


    @Variables.undoGlobal
    def setProfile(self, label, title, format, lst, choice, freq, formula, NbPoint):
        """
        Public method.
        Sets data to create one profile named I{label}.
        """
        self.isNotInList(label, self.getProfilesLabelsList())
        self.isInt(NbPoint)

        label_xml = label + self.suffix
        node = self.node_prof.xmlInitNode('profile', label=label_xml)
        node.xmlAddChild('format', name=format)
        for var in lst:
            self.isInList(var, self.__var_prop_list)
            (name, comp) = self.dicoLabel2Name[var]
            node.xmlAddChild('var_prop', name=name, component=comp)
        node.xmlSetData('output_type', choice)
        node.xmlSetData('output_frequency', freq)
        node['title'] = title
        self.__setFormula(label, formula)
        self.__setNbPoint(label, NbPoint)


    @Variables.undoGlobal
    def replaceProfile(self, old_label, label, title, format, lst, choice, freq, formula, NbPoint):
        """
        Public method.
        Replaces data from I{old_label} profile
        with label, frequency, and coordinates values.
        """
        self.isInList(old_label,self.getProfilesLabelsList())
        if label != old_label:
            self.isNotInList(label, self.getProfilesLabelsList())
        self.isInt(NbPoint)

        old_label_xml = old_label + self.suffix
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=old_label_xml)
        if node:
            node['title'] = ""
            for tag in ('format', 'var_prop', 'output_frequency', 'formula','points'):
                node.xmlRemoveChild(tag)
            node.xmlAddChild('format', name=format)
            for var in lst:
                self.isInList(var, self.__var_prop_list)
                (name, comp) = self.dicoLabel2Name[var]
                node.xmlAddChild('var_prop', name=name, component=comp)
            node['label'] = label_xml
            node.xmlSetData('output_type', choice)
            node.xmlSetData('output_frequency', freq)
            node['title'] = title
            self.__setFormula(label, formula)
            self.__setNbPoint(label, NbPoint)


    @Variables.undoLocal
    def deleteProfile(self, label):
        """
        Public method.
        Deletes profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
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
        label_xml = label + self.suffix
        node = self.node_prof.xmlGetNode('profile', label=label_xml)
        choice = node.xmlGetString('output_type')
        if choice == "time_value":
            freq = node.xmlGetDouble('output_frequency')
        else:
            freq = node.xmlGetInt('output_frequency')
        title = node['title']
        format = 'DAT'
        f_node = node.xmlGetChildNode('format')
        if f_node:
            format = f_node['name']
        formula = self.__getFormula(label)
        NbPoint = self.__getNbPoint(label)
        for var in node.xmlGetChildNodeList('var_prop'):
            for name in self.__var_prop_list:
                if self.dicoLabel2Name[name] == (var['name'], var['component']) :
                    lst.append(name)
        label_xml = node['label']
        label = label_xml

        return label, title, format, lst, choice, freq, formula, NbPoint


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
        mdl.setProfile('prof1.dat', 'title', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        doc = '''<profiles>
                    <profile label="prof1.dat" title="title">
                            <var_prop name="velocity" component="0"/>
                            <var_prop name="velocity" component="1"/>
                            <var_prop name="velocity" component="2"/>
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
        mdl.setProfile('prof1', 'title', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.replaceProfile('prof1', 'premier', 'title_bis', ['VelocitU'], 30, 0.1, 0.2, 0.3, 2.0, 2.0, 2.0)
        doc = '''<profiles>
                    <profile label="premier" title="title_bis">
                            <var_prop name="velocity" component="0"/>
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

        mdl.setProfile('prof1.dat', 'title1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2.dat', 'title2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        mdl.deleteProfile('prof1')
        doc = '''<profiles>
                    <profile label="prof2.dat" title="title2">
                            <var_prop name="velocity" component="0"/>
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
        mdl.setProfile('prof1.dat', 'title1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2.dat', 'title2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        prof = ('prof2.dat', 'title2', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)

        assert mdl.getProfileData('prof2') == prof,\
            'Could not get values for profile named label in ProfilesModel'


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
