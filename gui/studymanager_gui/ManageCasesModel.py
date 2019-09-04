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
This module manages studies and cases :

This module defines the following classes:
- ManageCasesModel
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

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class ManageCasesModel(Model):

    def __init__(self, case):
        """
        Constuctor.
        """
        self.case = case
        self.list_study = self.case.xmlGetNodeList("study")
        self.list_case = self.case.xmlGetNodeList("case")
        self.repo = self.case.xmlGetNode("studymanager").xmlGetString('repository')

        self.StudyList = []
        for node in self.list_study:
            self.StudyList.append(node['label'])


    def _defaultValues(self):
        """
        Return in a dictionary which contains default values
        """
        default = {}
        default['compute']        = "on"
        default['post']           = "on"
        default['status']         = "on"
        default['run_id']         = ""
        default['prepro_status']  = "off"
        default['post_status']    = "off"
        default['compare_status'] = "off"

        return default


    def getCaseList(self, name):
        """
        Get list of case name for a study
        """
        node = self.case.xmlGetNode('study', label = name)
        lst = []
        for nn in node.xmlGetNodeList("case"):
            lst.append(int(nn['id']))
        return lst


    def addCase(self, study, name):
        """
        Add case name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        idx = len(study_node.xmlGetNodeList("case"))
        node = study_node.xmlInitChildNode("case", id = idx, label = name)
        self.list_case.append(node)


    def addStudy(self, name):
        """
        Add study name from node with index
        """
        node = self.case.xmlGetNode("studymanager")

        study_node = node.xmlInitChildNode('study', label = name)
        self.StudyList.append(name)
        self.list_study = self.case.xmlGetNodeList("study")


    def loadCases(self, study):
        """
        load cases
        """
        node = self.case.xmlGetNode("studymanager")
        idx = 0
        directory = os.path.abspath(os.path.join(self.repo, study))
        study_node = node.xmlGetNode('study', label = study)
        for fl in os.listdir(directory):
            rep = os.path.abspath(os.path.join(directory, fl))
            if os.path.isdir(rep) and "DATA" in os.listdir(rep) \
                                  and "SRC" in os.listdir(rep):
                study_node.xmlInitNode('case', label = fl, id = idx)
                idx = idx + 1
        self.list_case = self.case.xmlGetNodeList("case")


    def deleteStudy(self, study):
        """
        delete study name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        study_node.xmlRemoveNode()

        self.list_study = self.case.xmlGetNodeList("study")
        self.list_case = self.case.xmlGetNodeList("case")
        self.StudyList = []
        for node in self.list_study:
            self.StudyList.append(node['label'])


    def deleteCase(self, study, idx):
        """
        delete case name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        node = study_node.xmlGetNode("case", id = idx)
        node.xmlRemoveNode()

        for node in study_node.xmlGetNodeList('case'):
            try:
                if int(node['id']) > idx:
                    node['id'] = str(int(node['id']) - 1)
            except:
                pass
        self.list_case = self.case.xmlGetNodeList("case")


    def duplicateCase(self, study, idx):
        """
        duplicate case name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        ii = len(study_node.xmlGetNodeList("case"))
        node = study_node.xmlGetNode("case", id = idx)
        name = node['label']
        new_node = study_node.xmlInitChildNode("case", id = ii, label = name)

        compute = node['compute']
        if compute:
            new_node['compute'] = compute

        post = node['post']
        if post:
            new_node['post'] = post

        status = node['status']
        if status:
            new_node['status'] = status

        run_id  = node['run_id']
        if run_id:
            new_node['run_id'] = run_id

        self.list_case = self.case.xmlGetNodeList("case")


    def getCaseName(self, study_name, idx):
        """
        Get case name from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        return node['label']


    def getStudyStatus(self, study_name):
        """
        Get study status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        status = study_node['status']
        if not status:
            status = self._defaultValues()['status']
            self.setStudyStatus(study_name, status)
        return status


    def setStudyStatus(self, study_name, status):
        """
        Put study status from node with index
        """
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        study_node['status'] = status


    def getComputeStatus(self, study_name, idx):
        """
        Get compute status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        status = node['compute']
        if not status:
            status = self._defaultValues()['compute']
            self.setComputeStatus(study_name, idx, status)
        return status


    def setComputeStatus(self, study_name, idx, status):
        """
        Put compute status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        node['compute'] = status


    def getPostStatus(self, study_name, idx):
        """
        Get post status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        status = node['post']
        if not status:
            status = self._defaultValues()['post']
            self.setPostStatus(study_name, idx, status)
        return status


    def setPostStatus(self, study_name, idx, status):
        """
        Put post status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        node['post'] = status


    def getStatus(self, study_name, idx):
        """
        Get general status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        status = node['status']
        if not status:
            status = self._defaultValues()['status']
            self.setStatus(study_name, idx, status)
        return status


    def setStatus(self, study_name, idx, status):
        """
        Put general status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        node['status'] = status


    def getRunId(self, study_name, idx):
        """
        Get run_id from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        runid = node['run_id']
        if not runid:
            runid = self._defaultValues()['run_id']
            self.setRunId(study_name, idx, runid)
        return runid


    def setRunId(self, study_name, idx, run_id):
        """
        Put run_id from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        node['run_id'] = run_id


    def getPreproScriptStatus(self, study_name, idx):
        """
        Get prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("prepro")
        if nn:
            status = nn['status']
            if not status:
                status = self._defaultValues()['prepro_status']
                self.setPreproScriptStatus(study_name, idx, status)
        else:
            status = "off"
        return status


    def setPreproScriptStatus(self, study_name, idx, status):
        """
        Put prepro script status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("prepro")
        nn['status'] = status


    def getPostScriptStatus(self, study_name, idx):
        """
        Get post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("script")
        if nn:
            status = nn['status']
            if not status:
                status = self._defaultValues()['post_status']
                self.setPostScriptStatus(study_name, idx, status)
        else:
            status = "off"
        return status


    def setPostScriptStatus(self, study_name, idx, status):
        """
        Put post script status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("script")
        nn['status'] = status


    def getStudyPostScriptStatus(self, study_name):
        """
        Get post script status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlGetNode("postpro")
        if nn:
            status = nn['status']
            if not status:
                status = self._defaultValues()['post_status']
                self.setStudyPostScriptStatus(study_name, status)
        else:
            status = "off"
        return status


    def setStudyPostScriptStatus(self, study_name, status):
        """
        Put post script status from node with index
        """
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlInitChildNode("postpro")
        nn['status'] = status


    def getCompareStatus(self, study_name, idx):
        """
        Get post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("compare")
        if nn:
            status = nn['status']
            if not status:
                status = self._defaultValues()['compare_status']
                self.setCompareStatus(study_name, idx, status)
        else:
            status = "off"
        return status


    def setCompareStatus(self, study_name, idx, status):
        """
        Put post script status from node with index
        """
        self.isInt(idx)
        self.isOnOff(status)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("compare")
        nn['status'] = status


    def getPreproScriptArgs(self, study_name, idx):
        """
        Get prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("prepro")
        args = ""
        if nn:
            args = nn['args']
        return args


    def setPreproScriptArgs(self, study_name, idx, args):
        """
        Put prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("prepro")
        nn['args'] = args


    def getPreproScriptName(self, study_name, idx):
        """
        Get prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("prepro")
        name = ""
        if nn:
            name = nn['label']
        return name


    def setPreproScriptName(self, study_name, idx, name):
        """
        Put prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("prepro")
        nn['label'] = name


    def getPostScriptArgs(self, study_name, idx):
        """
        Get post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("script")
        args = ""
        if nn:
            args = nn['args']
        return args


    def setPostScriptArgs(self, study_name, idx, args):
        """
        Put post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("post")
        nn['args'] = args


    def getPostScriptName(self, study_name, idx):
        """
        Get post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("script")
        name = ""
        if nn:
            name = nn['label']
        return name


    def setPostScriptName(self, study_name, idx, name):
        """
        Put post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("script")
        nn['label'] = name


    def getPostScriptInput(self, study_name, case_idx):
        """
        Get post script status from node with index
        """
        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = case_idx)
        nn = node.xmlGetNode("input")
        name = ""
        if nn:
            name = nn['file']
        return name


    def setPostScriptInput(self, study_name, case_idx, name):
        """
        Put post script status from node with index
        """
        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = case_idx)
        nn = node.xmlInitChildNode("input")
        nn['file'] = name


    def getStudyPostScriptArgs(self, study_name):
        """
        Get post script status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlGetNode("postpro")
        args = ""
        if nn:
            args = nn['args']
        return args


    def setStudyPostScriptArgs(self, study_name, args):
        """
        Put post script status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlInitChildNode("postpro")
        nn['args'] = args


    def getStudyPostScriptName(self, study_name):
        """
        Get post script status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlGetNode("postpro")
        name = ""
        if nn:
            name = nn['label']
        return name


    def setStudyPostScriptName(self, study_name, name):
        """
        Put post script status from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        nn = study_node.xmlInitChildNode("postpro")
        nn['label'] = name


    def getCompareArgs(self, study_name, idx):
        """
        Get prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlGetNode("compare")
        args = ""
        if nn:
            args = nn['args']
        return args


    def setCompareArgs(self, study_name, idx, args):
        """
        Put prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNode("case", id = idx)
        nn = node.xmlInitChildNode("compare")
        nn['args'] = args
