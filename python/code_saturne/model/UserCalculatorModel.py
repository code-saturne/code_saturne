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

import sys, unittest
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.Common import LABEL_LENGTH_MAX
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.OutputVolumicVariablesModel import OutputVolumicVariablesModel

class UserCalculatorModel(Variables, Model):
    """
    This class manages the user defined calculator operators.
    """

    def __init__(self, case):
        """
        Constructor
        """
        self.case = case
        _node = self.case.root().xmlInitNode('user_functions')
        self.calculator = _node.xmlInitNode('calculator')

        self.nb = NotebookModel(self.case)

    def defaultValues(self):
        _d = {}
        _d['support'] = "cells"
        _d['dimension'] = "1"

        return _d


    def getNumberOfFunctions(self, location=None):
        """
        """
        return len(self.__get_nodes_list(location))


    def getFunctionsNamesList(self, location=None):
        """
        """
        lst = []
        for n in self.__get_nodes_list(location):
            lst.append(n['name'])

        return lst


    @Variables.undoGlobal
    def addFunction(self, name):
        """
        Add a new user defined function
        """
        # Add only if needed
        n = self.__get_or_add_function_node(name)

        return n


    @Variables.undoGlobal
    def deleteFunction(self, name):
        """
        """
        self.__remove_function_node(name)

        # TODO: Add removing of time averages if needed


    @Variables.undoLocal
    def setFunctionName(self, old_name, new_name):
        """
        """
        n = self.__get_function_node(old_name)
        if n is not None:
            for elt in ('name', 'label'):
                n[elt] = new_name


    @Variables.noUndo
    def getFunctionLocation(self, name):
        """
        """
        n = self.__get_function_node(name)

        location = n['support'] if n is not None else self.defaultValues()['support']

        return location


    @Variables.undoLocal
    def setFunctionLocation(self, name, support):
        """
        """
        self.isInList(support, ('cells', 'internal', 'boundary', 'vertices'))
        n = self.__get_function_node(name)
        if n:
            n['support'] = support


    @Variables.noUndo
    def getFunctionDim(self, name):
        """
        """
        n = self.__get_function_node(name)
        dim = n['dimension'] if n else 1

        return dim


    @Variables.undoLocal
    def setFunctionDim(self, name, dim):
        """
        """
        self.isInList(dim, ("1", "2", "3", "6", "9"))
        n = self.__get_function_node(name)
        if n:
            n['dimension'] = dim


    @Variables.noUndo
    def getFunctionFormula(self, name):
        """
        """
        retval = None
        n = self.__get_function_node(name)
        if n:
            retval = n.xmlGetChildString('formula')

        # Generate default formula
        if not retval:
            dim = int(self.getFunctionDim(name))
            if dim == 1:
                retval = "{} = 1.;".format(name)
            else:
                retval = "{}[0] = 1.;".format(name)
                for i in range(1, dim):
                    retval += "\n{}[{}] = 1.;".format(name,i)

            self.setFunctionFormula(name, retval)

        return retval


    @Variables.undoLocal
    def setFunctionFormula(self, name, exp):
        """
        """
        n = self.__get_function_node(name)
        n.xmlSetData('formula', exp)


    @Variables.noUndo
    def getFunctionFormulaComponents(self, name):
        """
        """
        exp = self.getFunctionFormula(name)

        dim = int(self.getFunctionDim(name))
        if dim == 1:
            req = [(name, name)]
        else:
            req = []
            for i in range(dim):
                req.append(("{}[{}]".format(name,i), name))

        # Predefined symbols
        symbols = []

        # -> coordinates
        _coord_obj = ""
        _s = self.getFunctionLocation(name)
        if _s == "cells":
            _coord_obj = "cell center"
        elif _s == "boundary":
            _coord_obj = "boundary face center"
        elif _s == "internal":
            _coord_obj = "internal face center"
        elif _s == "vertices":
            _coord_obj = "vertex"

        for elt in ("x", "y", "z"):
            symbols.append((elt, "{} coordinate".format(_coord_obj)))


        # Known fields (for equation->C translation) within the code
        known_fields = []
        om = OutputVolumicVariablesModel(self.case)
        _data = om.getVolumeFieldsLabel2Name(time_averages=False, get_components=True)
        for k in _data.keys():
            d = _data[k]
            # Known fields is only fields!
            if d[1] == "-1":
                known_fields.append((k, d[0], d[2]))
            elif d[1] is None:
                known_fields.append((k, d[0]))

            if d[1] != "-1":
                symbols.append((k, d[0]))


        for (nme, val) in self.nb.getNotebookList():
            symbols.append((nme, 'value (notebook) = {}'.format(str(val))))

        return exp, req, known_fields, symbols


    @Variables.noUndo
    def getPrintingStatus(self, f_name):
        """
        """
        n = self.__get_function_node(f_name)
        return self.__get_on_off_sub_node_status(n, 'listing_recording')


    @Variables.undoLocal
    def setPrintingStatus(self, f_name, status):
        """
        Set status for printing in listing.
        """
        n = self.__get_function_node(f_name)
        self.__update_on_off_sub_node(n, 'listing_printing', status)


    @Variables.noUndo
    def getPostStatus(self, f_name):
        """
        """
        n = self.__get_function_node(f_name)
        return self.__get_on_off_sub_node_status(n, 'postprocessing_recording')


    @Variables.undoLocal
    def setPostStatus(self, f_name, status):
        """
        Update Postprocessing status
        """
        n = self.__get_function_node(f_name)
        self.__update_on_off_sub_node(n, 'postprocessing_recording', status)


    @Variables.noUndo
    def getMonitorStatus(self, f_name):
        """
        """
        n = self.__get_function_node(f_name)
        return self.__get_on_off_sub_node_status(n, 'probes_recording')


    @Variables.undoLocal
    def setMonitorStatus(self, f_name, status):
        """
        """
        n = self.__get_function_node(f_name)
        self.__update_on_off_sub_node(n, 'probes_recording', status)


    # -------------------------------------------------------------------------
    # Private functions
    # -------------------------------------------------------------------------

    def __get_function_node(self, f_name):
        """
        Get calculator function xml node
        """
        n = self.calculator.xmlGetNode("function", name=f_name)

        return n


    def __add_function_node(self, f_name):
        """
        Add a new calculator function node
        """

        self.calculator.xmlInitNode("function", name=f_name, label=f_name)

        n = self.__get_function_node(f_name)
        for key in self.defaultValues().keys():
            n[key] = self.defaultValues()[key]

        return n


    def __get_or_add_function_node(self, f_name):
        """
        """

        n = self.__get_function_node(f_name)
        if n is None:
            n = self.__add_function_node(f_name)

        return n

    def __remove_function_node(self, f_name):
        """
        """
        n = self.__get_function_node(f_name)
        if n is not None:
            n.xmlRemoveNode()


    def __update_on_off_sub_node(self, node, sub_node_name, status):
        """
        """
        if node is None:
            return

        self.isOnOff(status)

        if status is "off":
            node.xmlInitChildNode(sub_node_name)['status'] = status
        else:
            if node.xmlGetChildNode(sub_node_name):
                node.xmlRemoveChild(sub_node_name)


    def __get_on_off_sub_node_status(self, node, sub_node_name):
        """
        """
        retval = "on"

        if node:
            sn = node.xmlGetChildNode(sub_node_name, 'status')
            if sn:
                retval = sn['status']

        return retval


    def __get_nodes_list(self, location=None):
        """
        """

        lst = []
        for n in self.calculator.xmlGetNodeList("function"):
            if location:
                if n['support'] == location:
                    lst.append(n)
            else:
                lst.append(n)


        return lst
