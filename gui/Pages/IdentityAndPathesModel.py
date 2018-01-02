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
This module defines the Page in which the user defines the path of the treated
case. The Page verify that the case configuration directories is appropriate.

This module contents the following classes:
- IdentityAndPathesModel
- IdentityAndPathesModelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, string, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class IdentityAndPathesModel(Model):
    """
    Class to manipulate xml file.
    """
    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case


    def getXmlFileName(self):
        """
        Get xml file'name
        """
        return self.case['xmlfile']


    def getCasePath(self):
        """
        Get xml file'name
        """
        return self.case['case_path']


    def setId(self, ncase, nstudy):
        """
        Put case and study's names
        """
        self.case.root().xmlSetAttribute(case=ncase, study=nstudy)


    @Variables.noUndo
    def setCasePath(self, dircase):
        """
        Put path of case into xml file
        """
        self.case['case_path'] = dircase


    @Variables.noUndo
    def setRelevantSubdir(self, val, directory):
        """
        Put relevant_subdir value into xml file
        """
        self.isInList(val, ('yes', 'no'))

        dirList = ["DATA", "RESU", "SRC", "SCRIPTS",
                   'data_path','resu_path',
                   'user_src_path','scripts_path','mesh_path']
        self.isInList(directory, dirList)

        self.case['relevant_subdir'] = val


    @Variables.noUndo
    def setPath(self, pathi, tag):
        """
        Put relevant_subdir value into xml file
        """
        self.isInList(pathi, ('data_path','resu_path',
                              'user_src_path',
                              'scripts_path','mesh_path'))
        self.case[pathi] = tag

#-------------------------------------------------------------------------------
# IdentityAndPathesModel test Case
#-------------------------------------------------------------------------------

class IdentityAndPathesModelTestCase(ModelTest):
    """
    Unittest.
    """
    def checkIdentityAndPathesModelInstantiation(self):
        """Check whether the IdentityAndPathesModel class could be instantiated"""
        model = None
        model = IdentityAndPathesModel(self.case)
        assert model != None, 'Could not instantiate IdentityAndPathesModel'


    def checkGetXmlFileName(self):
        model = IdentityAndPathesModel(self.case)
        model.setId('castoto', 'etude')
        self.case['xmlfile'] = 'trucmuche.xml'
        assert model.getXmlFileName() == 'trucmuche.xml',\
            'Could not get name of xml in IdentityAndPathesModel'


    def checkSetCaseandStudy(self):
        model = IdentityAndPathesModel(self.case)
        model.setId('castoto', 'etude')

        assert self.case.root()['case'] == 'castoto',\
            'Could not get case name in IdentityAndPathesModel'
        assert self.case.root()['study'] == 'etude',\
            'Could not get study name in IdentityAndPathesModel'


    def checkSetRelevantSubdir(self):
        model = IdentityAndPathesModel(self.case)
        model.setId('castoto', 'etude')
        model.setRelevantSubdir('yes', 'mesh_path')
        assert self.case['relevant_subdir']  == 'yes',\
            'Could not set relevant subdir in IdentityAndPathesModel'


    def checkSetandGetCase(self):
        model = IdentityAndPathesModel(self.case)
        model.setId('castoto', 'etude')
        directory = os.getcwd()
        model.setCasePath(directory)
        assert model.getCasePath() == directory,\
            'Could not set or get case path in IdentityAndPathesModel'


    def checkSetPathI(self):
        model = IdentityAndPathesModel(self.case)
        model.setId('castoto', 'etude')
        dir = os.getcwd()
        dir_resu = dir + '/RESU'
        model.setPath('resu_path', dir_resu)
        assert self.case['resu_path'] == dir_resu,\
            'Could not set pathI in IdentityAndPathesModel'


def suite():
    testSuite = unittest.makeSuite(IdentityAndPathesModelTestCase, "check")
    return testSuite


def runTest():
    print("IdentityAndPathesModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
