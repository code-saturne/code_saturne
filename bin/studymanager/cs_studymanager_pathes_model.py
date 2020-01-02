# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
This module defines the Page in which the user defines the path of the treated
case. The Page verify that the case configuration directories is appropriate.

This module contents the following classes:
- PathesModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class PathesModel(Model):
    """
    Class to manipulate xml file.
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case


    def getRepositoryPath(self):
        """
        Get repository directory name
        """
        node = self.case.xmlGetNode("studymanager")
        path = node.xmlGetString('repository')
        return path


    def setRepositoryPath(self, dirname):
        """
        Put repository directory name
        """
        node = self.case.xmlGetNode("studymanager")
        node.xmlSetData('repository', dirname)


    def getStudyNumber(self):
        """
        Get number of study in XML
        """
        return len(self.case.xmlGetNodeList("study"))


    def loadStudyAndCases(self, dirname):
        """
        automatic load
        """
        # add study and case in xml
        node = self.case.xmlGetNode("studymanager")
        for tt in os.listdir(dirname):
            directory = os.path.abspath(os.path.join(dirname, tt))
            if os.path.isdir(directory) and "POST" in os.listdir(directory) \
                                        and "MESH" in os.listdir(directory):
                study = node.xmlInitNode('study', label = tt, status = 'on')
                # add case in study
                idx = 0
                for fl in os.listdir(directory):
                    rep = os.path.abspath(os.path.join(directory, fl))
                    if os.path.isdir(rep) and "DATA" in os.listdir(rep) \
                                          and "SRC" in os.listdir(rep):
                        study.xmlInitNode('case', label = fl, id = idx)
                        idx = idx + 1


    def getDestinationPath(self):
        """
        Get repository directory name
        """
        node = self.case.xmlGetNode("studymanager")
        path = node.xmlGetString('destination')
        return path


    def setDestinationPath(self, dirname):
        """
        Put repository directory name
        """
        node = self.case.xmlGetNode("studymanager")
        node.xmlSetData('destination', dirname)
