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
from code_saturne.Base.XMLmodelNeptune import *
from code_saturne.Base.Common import LABEL_LENGTH_MAX

class UsersControlModel(Variables, Model):

    """
    This class manages the users objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case          = case
        self.XMLUserScalar = self.case.xmlGetNode('additional_scalars')
        self.XMLUser       = self.XMLUserScalar.xmlInitNode('users')


    def defaultValues(self):
        default = {}

        default['support'] = "cells"
        default['dim'] = "1"
        return default


    def getUsersList(self):
        """
        Return the user array list
        """
        llst = {}
        for node in self.XMLUser.xmlGetNodeList('variable'):
            idx = node['name'].split('_')[1]
            tmp = idx.zfill(4)
            llst[tmp] = node['name']

        keys = list(llst.keys())
        keys.sort()

        lst = []
        for key in keys:
            lst.append(llst[key])
        return lst


    @Variables.undoGlobal
    def addUser(self, userId):
        """
        add a new user
        """
        self.isInt(userId)
        name  = "User_" + str(userId)
        label = "User_" + str(userId)
        # TODO control les names et labels pour noms differents

        Variables(self.case).setNewVariableProperty("variable",
                                                    "",
                                                    self.XMLUser,
                                                    "none",
                                                    name,
                                                    label,
                                                    dim = self.defaultValues()['dim'])

        n = self.XMLUser.xmlGetNode('variable', name = name)
        n['support']   = self.defaultValues()['support']
        n['dimension'] = self.defaultValues()['dim']

        return name



    @Variables.noUndo
    def deleteUser(self, userId):
        """
        delete a user
        """
        # delete user
        self.isInt(userId)
        name  = "User_" + str(userId)
        for node in self.XMLUser.xmlGetNodeList('variable'):
            try :
               if node['name'] == name :
                  node.xmlRemoveNode()
            except :
               pass

        #suppress profile
        for node in reversed(self.node_profile.xmlGetNodeList('profile')):
            suppress = 0
            for child in node.xmlGetNodeList('var_prop'):
                if (child['name'] == name):
                    suppress = 1
            if suppress == 1:
                label = node['label']
                ProfilesModel(self.case).deleteProfile(label)

        #suppress average
        for node in reversed(self.node_average.xmlGetNodeList('time_average')):
            suppress = 0
            for child in node.xmlGetNodeList('var_prop'):
                if (child['name'] == name):
                    suppress = 1
                    break
            if suppress == 1:
                label = node['label']
                TimeAveragesModel(self.case).deleteTimeAverage(label)

        # update name for other scalar in XML file
        index = 1
        for node in self.XMLUser.xmlGetNodeList('variable'):
            try :
               if index >= userId :
                  oldname = "User_" + str(index + 1)
                  name = "User_" + str(index)
                  node['name'] = name
            except :
               pass
            index += 1


    @Variables.noUndo
    def getUserLabelList(self):
        """
        Get label list
        """
        lst = []
        for node in self.XMLUser.xmlGetNodeList('variable'):
            if node['label']:
                lst.append(node['label'])
        return lst


    @Variables.undoLocal
    def setUsersLabel(self, idx, label):
        """
        Put label
        """
        label_new = label[:LABEL_LENGTH_MAX]
        name  = "User_" + str(idx)
        if label_new not in self.getUserLabelList():
            node = self.XMLUser.xmlGetNode('variable', name = name)
            if node :
               node['label'] = label


    @Variables.noUndo
    def getUsersLabel(self, name):
        """
        Get label
        """
        label = ""
        node = self.XMLUser.xmlGetNode('variable', name = name)
        if node:
            label = node['label']
        return label


    @Variables.noUndo
    def getUsersSupport(self, name):
        """
        Get support
        """
        node = self.XMLUser.xmlGetNode('variable', name = name)
        if node:
            support = node['support']
        return support


    @Variables.undoLocal
    def setUsersSupport(self, idx, support):
        """
        Put support
        """
        self.isInList(support, ('cells', 'internal', 'boundary'))
        name  = "User_" + str(idx)
        node = self.XMLUser.xmlGetNode('variable', name = name)
        if node :
            node['support'] = support


    @Variables.noUndo
    def getUsersDim(self, name):
        """
        Get array dimension
        """
        dim = 1
        node = self.XMLUser.xmlGetNode('variable', name = name)
        if node:
            dim = node['dimension']
        return dim


    @Variables.undoLocal
    def setUsersDim(self, idx, dim):
        """
        Set array dimension
        """
        self.isInList(dim, ("1", "2", "3", "6", "9"))
        name = "User_" + str(idx)
        node = self.XMLUser.xmlGetNode('variable', name = name)
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

