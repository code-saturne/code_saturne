
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

"""
This module defines the 1D profile management page.

This module defines the following classes:
- ProfilesModel
- ProfilesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.OutputVolumicVariablesModel import OutputVolumicVariablesModel

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

        self.__var_prop_list = self.getVariablesAndVolumeProperties()
        self.__var_prop_list_compat = self.dicoLabel2Name_compat.keys()


    def __defaultValues(self):
        """
        Private method.
        Returns a dictionnary with default values.
        """
        value = {}
        value['nfreq']         = -1
        value['formula']       =  "x = 0;\ny = s;\nz = 0;\n"
        value['points']        =  200
        value['label']         =  "profile"
        value['choice']        =  "frequency"
        value['frequency']     =  1
        value['format']        =  "CSV"
        value['snap_mode']     =  "snap_to_center"
        value['interpolation'] =  "no"

        return value


    @Variables.noUndo
    def getVariablesAndVolumeProperties(self):
        """
        Creates a dictionnary to connect name and label from
        variables and properties.
        """
        mdl = OutputVolumicVariablesModel(self.case)

        self.dicoLabel2Name = mdl.getVolumeFieldsLabel2Name(time_averages=True)
        self.dicoLabel2Name_compat = mdl.getVolumeFieldsLabel2Name(time_averages=True,
                                                                   get_components=True)


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


    def getNbPoint(self, label):
        """
        Private method.
        Gets coordinates for profile named I{label}.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        return node.xmlGetInt('points')


    @Variables.undoLocal
    def setSnapMode(self, label, mode):
        """
        """
        node = self.node_prof.xmlGetNode('profile', label=label)
        if mode == self.__defaultValues()['snap_mode']:
            node.xmlRemoveChild('snap_mode')
        else:
            node.xmlSetData('snap_mode', mode)


    @Variables.noUndo
    def getSnapMode(self, label):
        """
        """
        node = self.node_prof.xmlGetNode('profile', label = label)
        mode = node.xmlGetString('snap_mode')
        if not mode:
            mode = self.__defaultValues()['snap_mode']
        return mode


    @Variables.undoLocal
    def setProfileInterpolation(self, label, state):
        """
        """
        node = self.node_prof.xmlGetNode('profile', label = label)
        if state == self.__defaultValues()["interpolation"]:
            node.xmlRemoveChild('interpolation')
        else:
            node.xmlSetData('interpolation', state)


    @Variables.noUndo
    def getProfileInterpolation(self, label):
        """
        """
        node = self.node_prof.xmlGetNode('profile', label = label)
        state = node.xmlGetString('interpolation')
        if not state:
            state = self.__defaultValues()["interpolation"]

        return state


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
        OutputVolumicVariablesModel(self.case).setVariablesAtNode(node, lst, time_averages=True, get_components=True)

    @Variables.noUndo
    def getVariable(self, label):
        """
        Public method.
        """
        self.isInList(label, self.getProfilesLabelsList())
        node = self.node_prof.xmlGetNode('profile', label = label)
        return OutputVolumicVariablesModel(self.case).getVariablesAtNode(node, time_averages=True, get_components=False)


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

        NbPoint = self.getNbPoint(label)
        if not NbPoint:
            NbPoint = self.__defaultValues()['points']
            self.setNbPoint(label, NbPoint)

        lst = OutputVolumicVariablesModel(self.case).getVariablesAtNode(node, time_averages=True, get_components=True)

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
        """Check whether the ProfilesModel class could set one profile"""
        mdl = ProfilesModel(self.case)
        mdl.setProfile('prof1.dat', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        doc = '''<profiles>
                    <profile label="prof1.dat">
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
        mdl.setProfile('prof1', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.replaceProfile('prof1', 'premier', ['VelocitU'], 30, 0.1, 0.2, 0.3, 2.0, 2.0, 2.0)
        doc = '''<profiles>
                    <profile label="premier">
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

        mdl.setProfile('prof1.dat', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2.dat', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        mdl.deleteProfile('prof1')
        doc = '''<profiles>
                    <profile label="prof2.dat">
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
        mdl.setProfile('prof1.dat', ['VelocitU', 'VelocitV', 'VelocitW'], 20, 0.1, 0.2, 0.3, 2.1, 2.2, 2.3)
        mdl.setProfile('prof2.dat', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)
        prof = ('prof2.dat', ['VelocitU'], 20, 0.9, 0.8, 0.7, 9.1, 9.2, 9.3)

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
