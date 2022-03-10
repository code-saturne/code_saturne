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

    def initialize(self):
        """
        Verify that all Headings exist only once in the XMLDocument and
        create the missing heading.
        """
        msg = self.__initHeading()
        if msg:
            return msg

        self._backwardCompatibility()

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

        # Remove ids previously added by GUI.

        for nn in self.case.xmlGetNodeList('study'):
            for node in nn.xmlGetNodeList("case"):
                try:
                    if node['id'] != None:
                        del(node['id'])
                except Exception:
                    pass


    def countPreproNodes(self):
        """
        Check if XML has old "prepro" type nodes.
        """
        n_prepro_nodes = [0, 0]

        for node in self.case.xmlGetNodeList('prepro', 'label'):
            if node['status'] == 'on':
                n_prepro_nodes[0] += 1
            else:
                n_prepro_nodes[1] += 1

        return n_prepro_nodes


    def convertPreproNodes(self):
        """
        Convert old XML old "prepro" type nodes to "kw_args".
        """
        for sn in self.case.xmlGetNodeList('study'):
            for cn in sn.xmlGetNodeList("case"):
                for pn in cn.xmlGetNodeList('prepro', 'label'):
                    s = ' --prepro-script=' + pn['label'] + ' '
                    s += pn['args']
                    if pn['status'] == 'off':
                        s = ' --status=off'
                    lst = cn.xmlGetNodeList('kw_args')
                    if lst:
                        for i, n in enumerate(lst):
                            s += ' ' + n['args']
                            n.xmlRemoveNode()
                    kn = cn.xmlInitChildNode('kw_args', 'args')
                    s_new = kn['args'] + s
                    kn['args'] = s_new.strip()
                    pn.xmlRemoveNode()


#-------------------------------------------------------------------------------
# End of XMLinit
#-------------------------------------------------------------------------------
