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

"""
This module contains the following classes and function:
- NotebookModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from math import sqrt
import os, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model, Variables

#-------------------------------------------------------------------------------
# Body Force model class
#-------------------------------------------------------------------------------

class NotebookModel(Model):

    def __init__(self, case):

        """
        Constuctor.
        """
        self.case = case

        self.node_pp   = self.case.xmlGetNode('physical_properties')
        if not self.node_pp:
            self.node_pp = self.case.root().xmlInitNode('physical_properties')
        self.node_note = self.node_pp.xmlInitNode('notebook')


    def defaultNotebookValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['val']  = 0.0
        default['name'] = "var"

        return default


    def getVarList(self):
        """
        Return list of id
        """
        lst = []
        for nn in self.node_note.xmlGetNodeList("var"):
            lst.append(int(nn['id']))
        return lst


    def getVarNameList(self):
        """
        Return list of name
        """
        lst = []
        for nn in self.node_note.xmlGetNodeList("var"):
            lst.append(nn['name'])
        return lst


    def getNotebookList(self):
        """
        Return list of name
        """
        lst = []
        for nn in self.node_note.xmlGetNodeList("var"):
            lst.append((nn['name'], nn['value']))
        return lst


    @Variables.undoGlobal
    def addVariable(self):
        """
        Return list of case
        """
        idx = len(self.node_note.xmlGetNodeList("var"))
        name = self.defaultNotebookValues()['name'] + str(idx)
        node = self.node_note.xmlInitChildNode("var", id = idx, name = name)


    @Variables.undoGlobal
    def deleteVariable(self, idx):
        """
        Return list of case
        """
        n = self.node_note.xmlGetNode("var", id = idx)
        if n:
            n.xmlRemoveNode()

            for node in self.node_note.xmlGetNodeList('var'):
                try:
                    if int(node['id']) > idx:
                        node['id'] = str(int(node['id']) - 1)
                except:
                    pass


    @Variables.undoGlobal
    def ImportVariableFromFile(self, name):
        """
        read a file and add/update variables
        """
        base, ext = os.path.splitext(name)
        ficIn = open(name, 'r')
        for line in ficIn.readlines():
            if ext == '.txt':
                content = line.split()
            else:
                content = line.split(',')

            if len(content) == 2:
                var = content[0]
                value = content[1]

                if var in self.getVarNameList():
                    node = self.node_note.xmlGetNode("var", name = var)
                    node['value'] = value
                else:
                    self.addVariable()
                    idx = len(self.getVarList())
                    self.setVariableValue(idx-1, value)
                    self.setVariableName(idx-1, var)


    @Variables.noUndo
    def getVariableName(self, idx):
        """
        Return name of variable
        """
        name = self.defaultNotebookValues()['name'] + str(idx)
        node = self.node_note.xmlInitChildNode("var", id = idx)
        if node:
            name = node['name']
        else:
            self.setVariableName(idx, name)

        return name


    @Variables.undoGlobal
    def setVariableName(self, idx, name):
        """
        Return list of case
        """
        node = self.node_note.xmlInitChildNode("var", id = idx)
        node['name'] = name


    @Variables.noUndo
    def getVariableValue(self, idx):
        """
        Return name of variable
        """
        node = self.node_note.xmlInitChildNode("var", id = idx)
        value = node['value']
        if not value:
            value = self.defaultNotebookValues()['val']
            self.setVariableValue(idx, value)

        return value


    @Variables.undoGlobal
    def setVariableValue(self, idx, val):
        """
        Return list of case
        """
        node = self.node_note.xmlInitChildNode("var", id = idx)
        node['value'] = val


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
