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
        default['tags']           = ""
        default['prepro_status']  = "off"
        default['post_status']    = "off"
        default['compare_status'] = "off"

        return default


    def __get_post_script_node__(self, study_name, case_idx, script_idx):
        """
        Get post script node for a given script index.
        This assumes we are sure the node ies present.
        """

        self.isInt(case_idx)
        self.isInt(script_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        nl = None
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
            nl = node.xmlGetChildNodeList("script")
        else:
            nl = study_node.xmlGetChildNodeList("postpro")
        return nl[script_idx]


    def getCaseList(self, name):
        """
        Get list of case name for a study
        """
        node = self.case.xmlGetNode('study', label = name)
        lst = []
        for id, nn in enumerate(node.xmlGetNodeList("case")):
            lst.append(id)
        return lst


    def addCase(self, study, name):
        """
        Add case name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        idx = len(study_node.xmlGetNodeList("case"))
        node = study_node.xmlInitChildNode("case", id_tmp = 'new', label = name)
        del(node['id_tmp'])
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        node.xmlRemoveNode()

        self.list_case = self.case.xmlGetNodeList("case")


    def duplicateCase(self, study, idx):
        """
        duplicate case name from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study)
        ii = len(study_node.xmlGetNodeList("case"))
        node = study_node.xmlGetNodeByIdx("case", idx)
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

        tags  = node['tags']
        if tags:
            new_node['tags'] = tags

        self.list_case = self.case.xmlGetNodeList("case")


    def getCaseName(self, study_name, idx):
        """
        Get case name from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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


    def getStudyTags(self, study_name):
        """
        Get study tags from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        tags = study_node['tags']
        if not tags:
            tags = self._defaultValues()['tags']
        return tags


    def setStudyTags(self, study_name, tags):
        """
        Put study tags from node with index
        """
        study_node = self.case.xmlGetNode('study', label = study_name)
        study_node['tags'] = tags
        if tags == "":
            del(study_node['tags'])


    def getComputeStatus(self, study_name, idx):
        """
        Get compute status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        node['compute'] = status


    def getPostStatus(self, study_name, idx):
        """
        Get post status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        node['post'] = status


    def getStatus(self, study_name, idx):
        """
        Get general status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        node['status'] = status


    def getRunId(self, study_name, idx):
        """
        Get run_id from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        node['run_id'] = run_id


    def getTags(self, study_name, idx):
        """
        Get tags from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        tags = node['tags']
        if not tags:
            tags = self._defaultValues()['tags']
        return tags


    def setTags(self, study_name, idx, tags):
        """
        Put tags from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        node['tags'] = tags
        if tags == "":
            del(node['tags'])


    def getCompareStatus(self, study_name, idx):
        """
        Get post script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        nn = node.xmlInitChildNode("compare")
        nn['status'] = status


    def getNotebookArgs(self, study_name, idx):
        """
        Get notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        args = ""
        lst = node.xmlGetNodeList("notebook")
        if lst:
            for n in lst:
                if n['args']:
                    args += " " + n['args']
            args = args.strip()
        return args


    def setNotebookArgs(self, study_name, idx, args):
        """
        Put notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        lst = node.xmlGetNodeList("notebook")
        if lst:
            for i, n in enumerate(lst):
                if i > 0:
                    n.xmlRemoveNode()
        nn = node.xmlInitChildNode("notebook")
        nn['args'] = args


    def getParametricArgs(self, study_name, idx):
        """
        Get notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        args = ""
        lst = node.xmlGetNodeList("parametric")
        if lst:
            for n in lst:
                if n['args']:
                    args += " " + n['args']
            args = args.strip()
        return args


    def setParametricArgs(self, study_name, idx, args):
        """
        Put notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        lst = node.xmlGetNodeList("parametric")
        if lst:
            for i, n in enumerate(lst):
                if i > 0:
                    n.xmlRemoveNode()
        nn = node.xmlInitChildNode("parametric")
        nn['args'] = args


    def getKwArgs(self, study_name, idx):
        """
        Get notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        args = ""
        lst = node.xmlGetNodeList("kw_args")
        if lst:
            for n in lst:
                if n['args']:
                    args += " " + n['args']
            args = args.strip()
        return args


    def setKwArgs(self, study_name, idx, args):
        """
        Put notebook arguments from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
        # Remove extra nodes if data spread over multiple nodes.
        lst = node.xmlGetNodeList("kw_args")
        if lst:
            for i, n in enumerate(lst):
                if i > 0:
                    n.xmlRemoveNode()
        nn = node.xmlInitChildNode("kw_args")
        nn['args'] = args


    def getPostScripts(self, study_name, case_idx):
        """
        Get list of postprocessing scripts for a given study or case/run_id.
        """

        ps_list = []
        if study_name is None or study_name == '':
            return ps_list

        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        nl = None
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
            nl = node.xmlGetChildNodeList("script")
        else:
            nl = study_node.xmlGetChildNodeList("postpro")
        if nl:
            for n in nl:
                script = n['label']
                args = n['args']
                status = n['status']
                if script == None:
                    continue
                if args == None:
                    args =''
                if status == None:
                    status = 'on'
                ps_list.append([script, args, status])

        return ps_list


    def setPostScriptStatus(self, study_name, case_idx, script_idx, status):
        """
        Set post script status from node with index
        """
        self.isOnOff(status)

        nn = self.__get_post_script_node__(study_name, case_idx, script_idx)
        nn['status'] = status


    def setPostScriptArgs(self, study_name, case_idx, script_idx, args):
        """
        Put post script status from node with index
        """

        nn = self.__get_post_script_node__(study_name, case_idx, script_idx)
        nn['args'] = args


    def addPostScript(self, study_name, case_idx, label):
        """
        Put post script status from node with index
        """

        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        nl = None
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
            nn = node.xmlInitNode("script", label=label, args='',
                                  status='on', idx=-1)
        else:
            nn = study_node.xmlInitNode("postpro", label=label, args='',
                                        status='on', idx=-1)
        del(nn['idx'])


    def removePostScript(self, study_name, case_idx, script_idx):
        """
        Remove post script
        """

        nn = self.__get_post_script_node__(study_name, case_idx, script_idx)
        nn.xmlRemoveNode()


    def getPostInput(self, study_name, case_idx):
        """
        Get post script status from node with index
        """
        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
        else:
            node = study_node
        inputs = []
        lst = node.xmlGetChildNodeList("input")
        if lst:
            for n in lst:
                inputs.append(n['file'])
        return inputs


    def addPostInput(self, study_name, case_idx, name):
        """
        Put post script status from node with index
        """
        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
        else:
            node = study_node
        lst = node.xmlGetChildNodeList("input")
        if lst:
            for n in lst:
                if name == n['file']:
                    return False
        nn = node.xmlInitNode("input", file=name)
        return True


    def removePostInput(self, study_name, case_idx, name):
        """
        Put post script status from node with index
        """
        self.isInt(case_idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        if case_idx > -1:
            node = study_node.xmlGetNodeByIdx("case", case_idx)
        else:
            node = study_node
        lst = node.xmlGetChildNodeList("input")
        if lst:
            for n in lst:
                if name == n['file']:
                    n.xmlRemoveNode()


    def getCompareArgs(self, study_name, idx):
        """
        Get prepro script status from node with index
        """
        self.isInt(idx)
        study_node = self.case.xmlGetNode('study', label = study_name)
        node = study_node.xmlGetNodeByIdx("case", idx)
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
        node = study_node.xmlGetNodeByIdx("case", idx)
        nn = node.xmlInitChildNode("compare")
        nn['args'] = args
