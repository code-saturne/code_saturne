# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2021 EDF S.A.
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
This module contains the creation of the dependancy graph between CASES.

This module contains the following classes:
- node_case
- dependency_graph
"""

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
import logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
#log.setLevel(logging.DEBUG)
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# class node_case
#-------------------------------------------------------------------------------

class node_case:

    def __init__(self, name, n_procs, n_iter, estim_wtime, tags, depends):
        """ initializes a Code_Saturne CASE
        """
        # name should be in the folling form: STUDY/CASE/run_id
        # ex: "81_HEAD_LOSSES_COMPONENT/CASE_BEND/mesh12"
        self.name        = name
        self.n_procs     = n_procs
        self.tags        = tags
        # TODO: verify that these parameters are well handled by smgr in general
        self.n_iter      = n_iter
        # estimated wall time for the submission of the case on cluster
        self.estim_wtime = estim_wtime
        # depends should be in the same format as name
        self.depends     = depends
        # level is modified in add_node
        self.level       = None

    def __str__(self):
        res  = '\nCase ' + self.name
        res += ' on ' + str(self.n_procs) + ' procs'
        res += ' on ' + str(self.n_iter) + ' iterations'
        res += ' with an estimated wall time of ' + str(self.estim_wtime)
        res += ' with tags ' + str(self.tags)
        res += ' with level ' + str(self.level)

        return res

#-------------------------------------------------------------------------------
# class dependency_graph
#-------------------------------------------------------------------------------

class dependency_graph(object):

    def __init__(self):
        """ Initializes a dependency graph object to an empty dictionary
        """
        self.graph_dict = {}
        self.add_node(node_case('root', n_procs=0, n_iter=0, estim_wtime=0,
                                tags=None, depends=None))

    def add_dependency(self, dependency):
        """ Defines dependency between two node_cases as an edge in the graph
        """
        (node1, node2) = dependency
        if node1 in self.graph_dict:
            self.graph_dict[node1].append(node2)
        else:
            self.graph_dict[node1] = [node2]

    def add_node(self, node):
        """ Add a node_case in the graph if not already there.
            Add a dependency when depends parameters is defined
        """
        if node not in self.graph_dict:
            self.graph_dict[node] = []

            root_node = self.root_node()
            if node is not root_node:
                if node.depends:
                    for neighbor in self.graph_dict:
                        if neighbor.name == node.depends:
                            # cases with dependency are level > 1 and connected to the dependency
                            self.add_dependency((node, neighbor))
                            node.level = neighbor.level + 1
                            break
                    if node.level is None:
                        msg = "Problem in graph construction : dependency " \
                              + node.depends + " is not found.\n"
                        sys.exit(msg)
                else:
                    # cases with no dependency are level 1 and connected to the root node
                    self.add_dependency((node, root_node))
                    node.level = 1
            else:
                # root node is level 0 and have no dependency
                node.level = 0

    def nodes(self):
        """ returns the cases of the dependency graph """
        return list(self.graph_dict.keys())

    def dependencies(self):
        """ returns the dependencies between the cases of the graph """
        dependencies = []
        for node in self.graph_dict:
            for neighbor in self.graph_dict[node]:
                if (neighbor, node) not in dependencies:
                    dependencies.append((node, neighbor))
        return dependencies

    def extract_sub_graph(self, level, n_procs):
        """ extracts a sub_graph based on level and n_procs criteria"""
        sub_graph = dependency_graph()
        for node in self.graph_dict:
            if node.level == level and int(node.n_procs) == n_procs:
                sub_graph.add_node(node)
        return sub_graph

    def root_node(self):
        """ returns the root node of the graph """
        for node in self.graph_dict:
            if node.name == 'root':
                return node

    def __str__(self):
        res = "\nList of cases: "
        for node in self.nodes():
            res += str(node) + " "
        res += "\nList of dependencies: "
        for dependency in self.dependencies():
            (node1, node2) = dependency
            res += '\n' + str(node1.name) + ' depends on ' + str(node2.name)
        return res

#-------------------------------------------------------------------------------

