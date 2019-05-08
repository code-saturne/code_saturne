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
This module contains the XML file initialization class for studymanager.

This module contains the following class:
- smgr_xml_init
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest, re

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLinitialize import BaseXmlInit

#-------------------------------------------------------------------------------
# class smgr_xml_init
#-------------------------------------------------------------------------------

class smgr_xml_init(BaseXmlInit):
    """
    This class initializes the content of a smgr xml parameter file.
    """

    def initialize(self, reinit_indices = True):
        """
        Verify that all Headings exist only once in the XMLDocument and
        create the missing heading.
        """
        msg = self.__initHeading()
        if msg:
            return msg

        self._backwardCompatibility()

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


    def _backwardCompatibilityOldVersion(self, from_vers):
        """
        Change XML in order to ensure backward compatibility for old version
        """
        if from_vers <= "-1.0":
            self.__backwardCompatibilityBefore_6_0()


    def __backwardCompatibilityBefore_6_0(self):
        """
        Change XML to ensure backward compatibility from before 6.0 to 6.0
        """

    def _backwardCompatibilityCurrentVersion(self):
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
