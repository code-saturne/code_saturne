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
This module defines the Cathare coupling model data management.

This module contains the following classes and function:
- CathareCouplingModel
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Cathare coupling model class
#-------------------------------------------------------------------------------

class CathareCouplingModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.__node_models = self.case.xmlGetNode('thermophysical_models')
        self.__node_nepcat = self.__node_models.xmlInitNode('cathare_coupling')

    def defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['cathare_elt']        = "TUBE1"
        default['cathare_first_cell'] = 1
        default['cathare_last_cell']  = 2
        default['neptune_bc']         = "off"
        default['neptune_1d_zone']    = "all[]"
        return default

    # ----------------------------------
    def __getNumberOfCathareCoupling(self):
        return len(self.__node_nepcat.xmlGetNodeList('coupling'))

    # ----------------------------------
    def getNumberOfFluids(self):

        nodef = self.__node_models.xmlGetNode('fields')

        np = len(nodef.xmlGetNodeList('field'))

        return np

    # ----------------------------------
    def getCathareActivationStatus(self):

        value = self.__node_nepcat.xmlGetInt('api_type')
        if value is None:
            value = 0

        return value

    # ----------------------------------
    def setCathareFile(self, value):

        self.__node_nepcat.xmlSetData('jdd_name', value)

    def getCathareFile(self):

        value = self.__node_nepcat.xmlGetString('jdd_name')
        if value is None or value == '':
            value = 'cathare_jdd.dat'

        return value

    # ----------------------------------
    def setCplName(self, value):

        self.__node_nepcat.xmlSetData('coupling_name', value)

    def getCplName(self):

        value = self.__node_nepcat.xmlGetString('coupling_name')
        if value is None or value == '':
            value = 'NEPCAT_CPL'

        return value

    # ----------------------------------
    def setCathareTime(self, value):

        self.__node_nepcat.xmlSetData('cathare_init_time', value)

    def getCathareTime(self):

        value = self.__node_nepcat.xmlGetDouble('cathare_init_time')
        if value is None:
            value = 0.0

        return value

    # ----------------------------------
    def setCathareInstanceName(self, name):

        self.__node_nepcat.xmlSetData('cathare_instance_name', name)


    def getCathareInstanceName(self):

        name = self.__node_nepcat.xmlGetString('cathare_instance_name')
        if name is None:
            name = 'CATHARE'

        return name

    # ----------------------------------
    def setNeptuneInstanceName(self, name):

        self.__node_nepcat.xmlSetData('neptune_instance_name', name)


    def getNeptuneInstanceName(self):

        name = self.__node_nepcat.xmlGetString('neptune_instance_name')
        if name is None:
            name = 'NEPTUNE'

        return name


    # ----------------------------------
    def setCplTime(self, value):

        self.__node_nepcat.xmlSetData('coupling_run_time', value)

    def getCplTime(self):

        value = self.__node_nepcat.xmlGetDouble('coupling_run_time')
        if value is None:
            value = 0.0

        return value

    @Variables.noUndo
    def getCathareCouplingList(self):
        """
        @return: list of Cathare coupling description.
        @rtype: C{List}
        """
        node_list = self.__node_nepcat.xmlGetNodeList('coupling')
        lst = []
        for index in range(len(node_list)):
            num = index + 1

            cathare_elt        = self.getCathareEltName(num)
            cathare_first_cell = self.getCathareFCell(num)
            cathare_last_cell  = self.getCathareLCell(num)
            neptune_bc         = self.getNeptuneBc(num)
            neptune_1d_zone    = self.getNeptune1dZone(num)

            lst.append([cathare_elt, cathare_first_cell, cathare_last_cell,
                        neptune_bc, neptune_1d_zone])

        return lst


    @Variables.undoGlobal
    def addCathareCoupling(self, cathare_elt,
                           cathare_first_cell, cathare_last_cell,
                           neptune_bc, neptune_1d_zone):

        """
        Add a new definition of a Cathare coupling.

        @type cathare_elt: C{String}
        @param cathare_elt: Cathare coupled element
        @type cathare_first_cell: C{Int}
        @param cathare_first_cell: first coupled cell for cathare
        @type cathare_last_cell: C{Int}
        @param cathare_last_cell: last coupled cell for cathare
        @type neptune_bc : C{String}
        @param neptune_bc: neptune coupled boundary name
        @type neptune_1d_zone: C{String}
        @param neptune_1d_zone: selection criteria
        @return: new number of Cathare coupling
        @rtype: C{Int}
        """
        num = len(self.__node_nepcat.xmlGetNodeList('coupling'))
        node_new = self.__node_nepcat.xmlAddChild('coupling')

        num = num + 1
        self.setCathareEltName(num, cathare_elt)
        self.setCathareFCell(num, cathare_first_cell)
        self.setCathareLCell(num, cathare_last_cell)
        self.setNeptuneBc(num, neptune_bc)
        self.setNeptune1dZone(num, neptune_1d_zone)

        return num


    @Variables.undoLocal
    def deleteCathareCoupling(self, num):
        """
        Delete a definition of a Cathare coupling.

        @type num: C{Int}
        @param num: Cathare coupling number
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        node_list = self.__node_nepcat.xmlGetNodeList('coupling')
        node = node_list[num-1]
        node.xmlRemoveNode()

    #------------------------------------------------------------------
    # Helper function
    #------------------------------------------------------------------
    def __getStringData(self, index, name, setFunction):
        """
        Get string value from xml file.
        """
        self.isLowerOrEqual(index+1, self.__getNumberOfCathareCoupling())
        node = self.__node_nepcat.xmlGetNodeList('coupling')[index]
        value = node.xmlGetString(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getIntData(self, index, name, setFunction):
        """
        Get int value from xml file.
        """
        self.isLowerOrEqual(index+1, self.__getNumberOfCathareCoupling())
        node = self.__node_nepcat.xmlGetNodeList('coupling')[index]
        value = node.xmlGetInt(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getDefaultDataIfNone(self, index, value, name, setFunction):
        """
        Get default value if value is none.
        """
        if value is None or value == "":
            value = self.defaultValues()[name]
            setFunction(index+1, value)
        return value

    #------------------------------------------------------------------
    # API type: Stand alone NCFD, coupled NCFD or coupled CATHARE
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setApiType(self, value):

        self.__node_nepcat.xmlSetData('api_type', value)

    @Variables.noUndo
    def getApiType(self):

        return self.__node_nepcat.xmlGetInt('api_type')

    #------------------------------------------------------------------
    # Number of coupled phases
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setNphases(self, value):

        self.__node_nepcat.xmlSetData('nphases', value)

    @Variables.noUndo
    def getNphases(self):

        return self.__node_nepcat.xmlGetInt('nphases')

    #------------------------------------------------------------------
    # Cathare coupled element name
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setCathareEltName(self, num, value):
        """
        Set value of Cathare element name.

        @type num: C{Int}
        @param num: Cathare coupling number
        @type value: C{String}
        @param value: Cathare element name
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isStr(value)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('cathare_element', value)


    @Variables.noUndo
    def getCathareEltName(self, num):
        """
        Get value of Cathare element name.

        @type num: C{Int}
        @param num: Cathare coupling number
        @return: Cathare element name
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'cathare_element',
                                    self.setCathareEltName)

    #------------------------------------------------------------------
    # Coupled cells in the Cathare element
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setCathareFCell(self, num, value):
        """
        Set value of Cathare first coupled cell.

        @type num: C{Int}
        @param num: Cathare coupling number
        @type value: C{Int}
        @param value: Cathare first coupled cell
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isInt(value)
        self.isGreaterOrEqual(value, 1)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('first_cell', value)


    @Variables.noUndo
    def getCathareFCell(self, num):
        """
        Get value of Cathare first coupled cell.

        @type num: C{Int}
        @param num: Cathare coupling number
        @return: Cathare first coupled cell
        @rtype: C{Int}
        """
        return self.__getIntData(num-1,
                                 'first_cell',
                                 self.setCathareFCell)

    #------------------------------------------------------------------
    # Coupled cells in the Cathare element
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setCathareLCell(self, num, value):
        """
        Set value of Cathare last coupled cell.

        @type num: C{Int}
        @param num: Cathare coupling number
        @type value: C{Int}
        @param value: Cathare last coupled cell
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isInt(value)
        self.isGreaterOrEqual(value, 1)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('last_cell', value)


    @Variables.noUndo
    def getCathareLCell(self, num):
        """
        Get value of Cathare last coupled cell.

        @type num: C{Int}
        @param num: Cathare coupling number
        @return: Cathare last coupled cell
        @rtype: C{Int}
        """
        return self.__getIntData(num-1,
                                 'last_cell',
                                 self.setCathareLCell)

    #------------------------------------------------------------------
    # Coupled Boundary condition
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setNeptuneBc(self, num, value):
        """
        Set value of Neptune coupled boundary name.

        @type num: C{Int}
        @param num: Cathare coupling number
        @type value: C{String}
        @param value: Neptune coupled boundary name
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isStr(value)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('neptune_bc', value)


    @Variables.noUndo
    def getNeptuneBc(self, num):
        """
        Get value of Neptune coupled boundary name.

        @type num: C{Int}
        @param num: Cathare coupling number
        @return: Neptune coupled boundary name
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'neptune_bc',
                                    self.setNeptuneBc)

    @Variables.undoLocal
    def setNeptune2dZone(self, num, value):

        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isStr(value)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('selection_criteria_2d', value)


    @Variables.undoLocal
    def setNeptuneBcType(self, num, value):

        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isStr(value)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]

        bc_type = -1
        if value == 'inlet':
            bc_type = 0
        elif value == 'outlet':
            bc_type = 1
        elif value == 'wall':
            bc_type = 2

        node.xmlSetData('bc_cpl_type', bc_type)


    #------------------------------------------------------------------
    # Selection criteria
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setNeptune1dZone(self, num, value):
        """
        Set value of equivalent coupled 1D zone.

        @type num: C{Int}
        @param num: Cathare coupling number
        @type value: C{String}
        @param value: equivalent coupled 1D zone
        """
        self.isLowerOrEqual(num, self.__getNumberOfCathareCoupling())
        self.isStr(value)
        node = self.__node_nepcat.xmlGetNodeList('coupling')[num-1]
        node.xmlSetData('selection_criteria_3d', value)


    @Variables.noUndo
    def getNeptune1dZone(self, num):
        """
        Get value of equivalent coupled 1D zone.

        @type num: C{Int}
        @param num: Cathare coupling number
        @return: equivalent coupled 1D zone
        @rtype: C{String}
        """
        return self.__getStringData(num-1,
                                    'selection_criteria_3d',
                                    self.setNeptune1dZone)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
