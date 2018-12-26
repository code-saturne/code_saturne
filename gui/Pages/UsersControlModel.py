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
from code_saturne.Base.Common import LABEL_LENGTH_MAX

class UsersControlModel(Variables, Model):

    """
    This class manages the users objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        # XML file parameters (TODO: rename matching tree nodes)
        self.case     = case
        XMLUserScalar = self.case.xmlInitNode('additional_scalars')
        self.XMLUser  = XMLUserScalar.xmlInitNode('users')
        self.f_type  = "property"


    def defaultValues(self):
        default = {}

        default['support'] = "cells"
        default['dim'] = "1"
        return default


    def getUsersList(self, location=None):
        """
        Return the user array list
        """
        lst = []
        for node in self.XMLUser.xmlGetNodeList(self.f_type):
            if location:
                if node['support'] == location:
                    lst.append(node)
            else:
                lst.append(node)

        return lst


    @Variables.undoGlobal
    def addUser(self, name):
        """
        add a new user field
        """
        if name not in self.getUserNamesList():
            self.XMLUser.xmlInitNode(self.f_type, name=name, type="user", label=name)

        n = self.XMLUser.xmlGetNode(self.f_type, name=name)
        n['support']   = self.defaultValues()['support']
        n['dimension'] = self.defaultValues()['dim']

        return name


    @Variables.noUndo
    def deleteUser(self, name):
        """
        delete a user field
        """
        # delete user
        for node in self.XMLUser.xmlGetNodeList(self.f_type):
            try :
               if node['name'] == name:
                  node.xmlRemoveNode()
            except :
               pass

        # suppress in profiles and averages
        node_profile = None
        node_moments = None
        node_ac = self.case.xmlGetNode('analysis_control')
        if node_ac:
            node_profile = node_ac.xmlGetNode('profiles')
        if node_profile:
            for node in reversed(node_profile.xmlGetNodeList('profile')):
                suppress = False
                for child in node.xmlGetNodeList('var_prop'):
                    if (child['name'] == name):
                        child.xmlRemoveNode()

        # suppress average
        if node_ac:
            node_m = node_ac.xmlGetNode('time_averages')
        if node_m:
            for node in reversed(node_m.xmlGetNodeList('time_average')):
                suppress = False
                for child in node.xmlGetNodeList('var_prop'):
                    if (child['name'] == name):
                        suppress = True
                    break
                if suppress:
                    node.xmlRemoveNode()


    @Variables.noUndo
    def getUserNamesList(self, location=None):
        """
        Get label list
        """
        lst = []
        for node in self.XMLUser.xmlGetNodeList(self.f_type):
            if location:
                if node['support'] == location:
                    lst.append(node['name'])
            else:
                lst.append(node['name'])
        return lst


    @Variables.undoLocal
    def setUsersName(self, old_name, new_name):
        """
        Put label
        """
        if new_name not in self.getUserNamesList():
            node = self.XMLUser.xmlGetNode(self.f_type, name=old_name)
            if node:
               node['name'] = new_name
               node['label'] = new_name


    @Variables.noUndo
    def getUsersLocation(self, name):
        """
        Get location
        """
        node = self.XMLUser.xmlGetNode(self.f_type, name=name)
        if node:
            location = node['support']
        if not location:
            location = self.defaultValues()['support']
        return location


    @Variables.undoLocal
    def setUsersLocation(self, name, support):
        """
        Set location
        """
        self.isInList(support, ('cells', 'internal', 'boundary', 'vertices'))
        node = self.XMLUser.xmlGetNode(self.f_type, name=name)
        if node :
            node['support'] = support


    @Variables.noUndo
    def getUsersDim(self, name):
        """
        Get array dimension
        """
        dim = 1
        node = self.XMLUser.xmlGetNode(self.f_type, name = name)
        if node:
            dim = node['dimension']
        return dim


    @Variables.undoLocal
    def setUsersDim(self, name, dim):
        """
        Set array dimension
        """
        self.isInList(dim, ("1", "2", "3", "6", "9"))
        node = self.XMLUser.xmlGetNode(self.f_type, name=name)
        if node:
            node['dimension'] = dim


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class AdditionalUserTestCase(ModelTest):
    """
    """
    def checkAdditionalUserInstantiation(self):
        """Check whether the AdditionalUserModel class could be instantiated"""
        model = None
        model = AdditionalUserModel(self.case)
        assert model != None, 'Could not instantiate AdditionalUserModel'


    def checkGetUserList(self):
        """Check whether the AdditionalUserModel class could get the user list"""
        mdl = AdditionalUserModel(self.case)
        mdl.addUser(1)
        assert mdl.getUsersList() == ['User_1'],\
            'Could not get UserList'


def suite():
    testSuite = unittest.makeSuite(AdditionalUserTestCase, "check")
    return testSuite


def runTest():
    print("AdditionalUserTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

