#!/usr/bin/env python

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

import sys
import os.path
from xml.dom import minidom

#-------------------------------------------------------------------------------
# Utility functions
#-------------------------------------------------------------------------------

def getChildNode(node, tag, required=False):
    """
    Return a child node matching a tag.
    """
    childNode = None
    for child in node.childNodes:
        if child.nodeType == node.ELEMENT_NODE:
            if child.nodeName == tag:
                if childNode == None:
                    childNode = child
                else:
                    errStr = "Multiple instance of " + tag + "nodes"
                    raise Exception(errStr)

    if childNode == None and required:
        errStr = tag + " node not found under " + node.tagName
        raise Exception(errStr)

    return childNode

#-------------------------------------------------------------------------------

def childNodeList(node, tag):
    """
    Return a list of child nodes matching a tag.
    """
    childNodeList = []
    for child in node.childNodes:
        if child.nodeType == node.ELEMENT_NODE:
            if child.nodeName == tag:
                childNodeList.append(child)

    return childNodeList

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Default command: return simple string based on local information

    n_args = len(sys.argv) - 1

    if n_args != 1:
        print("Usage: %(prog) <salome_root_directory>" % {'prog':sys.argv[0]})
        sys.exit(1)

    salome_root_dir = sys.argv[1]

    root = None

    path = os.path.join(salome_root_dir, '.config_appli_template.xml')
    if not os.path.isfile(path):
        raise Exception('XML file: ' + path + ' not found')

    doc = minidom.parse(path)

    root = doc.documentElement

    env_modules_node = getChildNode(root, 'env_modules')
    nodeList = childNodeList(env_modules_node, 'env_module')

    for node in nodeList:
        name = str(node.getAttribute('name'))
        print(name)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
