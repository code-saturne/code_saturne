# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the XML data model.

This module contains the following class:
- XMLinit
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest, re

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables

#-------------------------------------------------------------------------------
# class XMLinit
#-------------------------------------------------------------------------------

class XMLinit(Variables):
    """
    This class initializes the XML contents of the case.
    """
    def __init__(self, case):
        """
        """
        self.case = case


    def initialize(self, reinit_indices = True):
        """
        Verify that all Headings exist only once in the XMLDocument and
        create the missing heading.
        """
        msg = self.__initHeading()
        if msg:
            return msg

        self.__backwardCompatibility()

        if reinit_indices:
            self.__reinitIndices()

        return msg


    def __initHeading(self):
        """
        Create if necessary headings from the root element of the case.
        """
        msg = ""
        tagList = ('repository','destination')

        for tag in tagList:
            nodeList = self.case.root().xmlInitChildNodeList(tag)

            if len(nodeList) > 1:
                msg = "More than one occurence of the tag: \n\n" + tag \
                      + "\n\nThe application will finish."

        for tag in tagList:
            nodeList = self.case.xmlInitNodeList(tag)

            if len(nodeList) > 1:
                msg = "More than one occurence of the tag: \n\n" + tag \
                      + "\n\nThe application will finish."

        return msg


    def __reinitIndices(self):
        """
        Insert indices to make xml compatible with GUI
        and reinitialize all indices.
        """
        for nn in self.case.xmlGetNodeList('study'):
            idx = 0
            for node in nn.xmlGetNodeList("case"):
                if not node['id']:
                    node['id'] = idx
                idx = idx + 1

        # ensure id for subplot 0 to n
        for nn in self.case.xmlGetNodeList('study'):
            idx = 0
            dico = {}
            idlst = []

            for node in nn.xmlGetNodeList("subplot"):
                dico[node['id']] = idx
                idlst.append(node['id'])
                node['id'] = idx
                idx = idx + 1

            idxx = 0
            for node in nn.xmlGetNodeList("figure"):
                lst = node['idlist']
                new_lst = ''
                if lst:
                    for idl in lst.split(" "):
                        if idl != " ":
                            if idl in idlst:
                                if new_lst != "":
                                    new_lst = new_lst + " " + str(dico[idl])
                                else:
                                    new_lst = str(dico[idl])
                    node['idlist'] = new_lst
                    if not node['id']:
                        node['id'] = idxx
                    idxx = idxx + 1

            idxx = 0
            for node in nn.xmlGetNodeList("plot"):
                lst = node['spids']
                new_lst = ''
                if lst:
                    for idl in lst.split(" "):
                        if idl != " ":
                            if new_lst != "":
                                new_lst = new_lst + " " + str(dico[idl])
                            else:
                                new_lst = str(dico[idl])
                    node['spids'] = new_lst
                    if not node['id']:
                        node['id'] = idxx
                    idxx = idxx + 1


    def __backwardCompatibility(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        cur_vers = self.case['package'].version

        if self.case.root()["solver_version"]:
            his_r = self.case.root()["solver_version"]
            history = his_r.split(";")
            last_vers = self.__clean_version(history[len(history) - 1])
            if last_vers == cur_vers:
                self.__backwardCompatibilityCurrentVersion()
            else:
                self.__backwardCompatibilityOldVersion(last_vers)
                self.__backwardCompatibilityCurrentVersion()
            his = ""
            vp = ""
            for v in history:
                vc = self.__clean_version(v)
                if vc != vp:
                    his += vc + ";"
                    vp = vc
            if cur_vers != vp:
                his += cur_vers + ";"
            his = his[:-1]
            if his != his_r:
                self.case.root().xmlSetAttribute(solver_version = his)

        else:
            vers = cur_vers
            self.case.root().xmlSetAttribute(solver_version = vers)

            # apply all backwardCompatibilities as we don't know
            # when it was created
            self.__backwardCompatibilityOldVersion("-1")
            self.__backwardCompatibilityCurrentVersion()


    def __backwardCompatibilityOldVersion(self, from_vers):
        """
        Change XML in order to ensure backward compatibility for old version
        """
        if from_vers <= "-1.0":
            self.__backwardCompatibilityBefore_6_0()


    def __backwardCompatibilityBefore_6_0(self):
        """
        Change XML to ensure backward compatibility from before 6.0 to 6.0
        """

    def __backwardCompatibilityCurrentVersion(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        # rename some atrributes of plot markup
        for o_attr in ['xfois', 'yfois', 'fig']:
            for node in self.case.xmlGetNodeList('plot', o_attr):
                val = node.xmlGetAttribute(o_attr)
                node.xmlDelAttribute(o_attr)
                if o_attr == "xfois":
                    node.xmlSetAttribute(xscale = val)
                elif o_attr == "yfois":
                    node.xmlSetAttribute(yscale = val)
                elif o_attr == "fig":
                    node.xmlSetAttribute(spids = val)

#-------------------------------------------------------------------------------
# End of XMLinit
#-------------------------------------------------------------------------------
