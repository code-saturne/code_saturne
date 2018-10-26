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
This module defines the Variables class which creates the <variable>,
<scalar> and <property> markups.

This module contains the following classes and function:
- Variables
- VariablesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Common import *
import Toolbox
from code_saturne.Base.XMLvariables import Variables as saturneVariables

#-------------------------------------------------------------------------------
# class Variables : creates <variable>, <scalar> and <property> markups.
#-------------------------------------------------------------------------------

class Variables(saturneVariables):
    """
    This class creates <variable>, <scalar> and <property> markups.
    Each new markup has a 'name' and a 'label' attribute.
    Each new markup has <listing_printing status='on'>,
    <postprocessing_recording status='on'> and several
    <probe_recording name="XX"> (if any <probe> exists) as child markups.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        saturneVariables.__init__(self, case)


    def setOutputControl(self, variable, post=False):
        """
        Update the output markups <probe_recording name="XX">,
        <postprocessing_recording status='on'> and
        <listing_printing status='on'> for the new 'variable' markup.
        """
        analysis_ctrl = self.case.xmlGetNode('analysis_control')
        node_output   = analysis_ctrl.xmlInitChildNode('output')
        if not (variable['support'] and variable['support'] == "boundary"):
            for node in node_output.xmlGetChildNodeList('probe', 'name'):
                num = node['name']
                variable.xmlInitChildNode('probe_recording', name=num)

        variable.xmlInitChildNode('listing_printing', status='on')
        if post:
            variable.xmlInitChildNode('postprocessing_recording', status='on')
        else:
            variable.xmlInitChildNode('postprocessing_recording', status='off')


    def setNewVariableProperty(self, type, choice, node, num, name, label, dim=None, support=None, post=False):
        id = str(num)
        label = label
        if not node.xmlGetNode(type, field_id=id,  name=name):
            if type != 'property' :
                if dim != None:
                    n = node.xmlInitNode(type, field_id=id, name=name, label=label, dimension=dim)
                else:
                    n = node.xmlInitNode(type, field_id=id, name=name, label=label, dimension='1')
            else :
                if support != None:
                    n = node.xmlInitNode(type, field_id=id, choice=choice, name=name, label=label, dimension='1', support = support)
                else:
                    n = node.xmlInitNode(type, field_id=id, choice=choice, name=name, label=label, dimension='1')
            self.setOutputControl(n, post=post)
            self.updateLabel(n)


    def removeVariableProperty(self, type, node, id, name):
        """
        """
        try :
            noder = node.xmlGetNode(type, field_id=id,  name=name)
            noder.xmlRemoveNode()
        except :
            pass


    def removeVariablesProperties(self, type, node, id):
        """
        """
        nodesList = node.xmlGetNodeList(type, field_id=str(id))
        for variableNode in nodesList:
            variableNode.xmlRemoveNode()


    def setNewTurbField(self, node, id, model, coupling):
        """
        Input a new node
        <field id="my_id" model="my_model", two_way_coupling="coupling">
        in the xmldoc.
        """
        if not node.xmlGetNode('field', field_id=id, model=model):
            v1 = node.xmlInitNode('field', field_id=str(id),
                                        model=model,
                                        two_way_coupling=coupling)
            self.updateLabel(v1)


    def getVariablesPropertiesList(self, average, constant, user) :
        """
        return list of variables, properties (and scalar)
        for Output field, profiles and averages
        if constant == yes we take account constant variables
        """
        self.XMLNodethermo  = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeclosure = self.case.xmlGetNode('closure_modeling')
        self.XMLNodeTurb    = self.XMLNodeclosure.xmlInitNode('turbulence')
        self.XMLNodeAna     = self.case.xmlGetNode('analysis_control')
        self.XMLNodeAverage = self.XMLNodeAna.xmlGetNode('time_averages')
        self.XMLUserScalar  = self.case.xmlGetNode('additional_scalars')
        self.XMLScalar      = self.XMLUserScalar.xmlInitNode('scalars')
        self.XMLUsers       = self.XMLUserScalar.xmlGetNode('users')
        list = []
        #TODO for variableType in ('variable', 'property', 'scalar') :
        for node in self.XMLNodethermo.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLNodeTurb.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLScalar.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLNodethermo.xmlGetNodeList('property') :
            choice = node['choice']
            if constant == 'yes' :
                list.append(node)
            elif choice != 'constant' :
                list.append(node)
        for node in self.XMLNodeTurb.xmlGetNodeList('property') :
            choice = node['choice']
            if constant == 'yes' :
                list.append(node)
            elif choice != 'constant' :
                list.append(node)
        if average == 'yes':
            # Warning average node is different
            for node in self.XMLNodeAverage.xmlGetNodeList('time_average'):
                list.append(node)
        if user == 'yes' and self.XMLUsers:
            for node in self.XMLUsers.xmlGetNodeList('variable'):
                list.append(node)

        return list

#-------------------------------------------------------------------------------
# End of XMLvariables
#-------------------------------------------------------------------------------
