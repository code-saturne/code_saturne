# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This page is devoted to the time tables management.

This module contains the following classes and function:
- TimeTablesModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Time Tables Model class
#-------------------------------------------------------------------------------

class TimeTablesModel(Model):
    """
    Class to handle time tables page.
    """
    def __init__(self, case):
        """
        Constructor
        """
        self.case = case

        self.node_pp = self.case.xmlGetNode('physical_properties')
        if not self.node_pp:
            self.node_pp = self.case.root().xmlInitNode('physical_properties')

        self.node_time_tables = self.node_pp.xmlInitNode('time_tables')


    def defaultValues(self, k):
        """
        Return a dictionnary which contains default values
        """
        _default = {'name':'table',
                    'delimiter':' ',
                    'file_name':'',
                    'skip_rows':0,
                    'cols2import':'all',
                    'col_ids':'',
                    'headers_def':'user',
                    'headers_list':'',
                    'headers_line':1}

        return _default.get(k, '')


    def getNumberOfTables(self):
        """
        Get number of defined time tables.
        """
        lst = self.node_time_tables.xmlGetNodeList('table')

        return len(lst)


    def getTableNamesList(self):
        """
        Return list of tables names.
        """
        lst = []

        for n in self.node_time_tables.xmlGetNodeList('table'):
            lst.append(n['name'])

        return lst


    @Variables.noUndo
    def getTableIdFromName(self, name):
        """
        """
        retval = -1

        for i, n in enumerate(self.getTableNamesList()):
            if n == name:
                retval = i
                break

        return retval


    @Variables.noUndo
    def getTableHeadersList(self, idx):
        """
        """
        hl = [elt.strip() for elt in self.getTableProperty(idx, 'headers_list').split(',')]

        return hl


    @Variables.undoGlobal
    def addTable(self):
        """
        Create a new table entry
        """
        idx = self.getNumberOfTables()

        name = "_".join([self.defaultValues("name"), str(idx)])

        _fname = self.defaultValues("file_name")
        _d     = self.defaultValues("delimiter")
        node = self.node_time_tables.xmlInitChildNode('table',
                                                      id = idx,
                                                      name = name,
                                                      file_name = _fname,
                                                      delimiter = _d)

        self.setTableProperty(idx, 'headers_list', '')


    @Variables.undoGlobal
    def deleteTable(self, idx):
        """
        Delete a table
        """
        n = self.__getTableNode(idx)
        if n:
            n.xmlRemoveNode()

            for node in self.node_time_tables.xmlGetNodeList('table'):
                try:
                    if int(node['id']) > idx:
                            node['id'] = str(int(node['id']) - 1)
                except:
                    pass


    @Variables.undoGlobal
    def setTableName(self, idx, name):
        """
        Set table name
        """
        node = self.__getTableNode(idx)
        node['name'] = name


    @Variables.noUndo
    def getTableName(self, idx):
        """
        Get name of table
        """
        node = self.__getTableNode(idx)
        if node:
            name = node['name']
        else:
            name = "_".join([self.defaultValues('name'), idx])
            self.setTableName(idx, name)

        return name


    @Variables.undoGlobal
    def setTableDelimiter(self, idx, delimiter):
        """
        Set table delimiter
        """
        node = self.__getTableNode(idx)
        node['delimiter'] = delimiter


    @Variables.noUndo
    def getTableDelimiter(self, idx):
        """
        Get table delimiter.
        """
        node = self.__getTableNode(idx)
        if node:
            d = node['delimiter']
        else:
            d = self.defaultValues('delimiter')
            self.setTableDelimiter(idx, d)

        return d


    @Variables.undoGlobal
    def setTableFileName(self, idx, file_name):
        """
        Set table file name.
        """
        node = self.__getTableNode(idx)
        node['file_name'] = file_name


    @Variables.noUndo
    def getTableFileName(self, idx):
        """
        Get table file name
        """
        node = self.__getTableNode(idx)
        if node:
            _fname = node['file_name']
        else:
            _fname = self.defaultValues('file_name')
            self.setTableFileName(idx, _fname)

        return _fname


    @Variables.undoGlobal
    def setTableProperty(self, idx, ppty, val=None):
        """
        Set Value of a given property.
        """
        node = self.__getTableNode(idx)

        # Default values for properties are not stored in xml file
        # to avoid increasing files too much
        node.xmlRemoveChild(ppty)

        if val not in (None, self.defaultValues(ppty)):
            node.xmlSetData(ppty, val)


    @Variables.noUndo
    def getTableProperty(self, idx, ppty):
        """
        Get value of a given property
        """
        node = self.__getTableNode(idx)

        # Default values for properties are not stored in xml file
        # to avoid increasing files too much
        val = node.xmlGetString(ppty)
        if val in (None, ''):
            val = self.defaultValues(ppty)

        return val


    def __getTableNode(self, idx):
        """
        Get datamodel node based on table index
        """
        node = self.node_time_tables.xmlInitChildNode('table', id = idx)

        return node
#-------------------------------------------------------------------------------
