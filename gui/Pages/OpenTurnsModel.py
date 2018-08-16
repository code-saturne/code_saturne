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
This module is used to handle the OpenTurns Study capacity for Code_Saturne
- OpenTurnsModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

from __future__ import print_function

import sys, unittest
import os, os.path, shutil, sys, types, re

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model
from code_saturne.Pages.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Class OpenTurnsModel
#-------------------------------------------------------------------------------

class OpenTurnsModel(Model):
    """
    This class modifies the OpenTurns study properties
    """
    # ---------------------------------------
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        case_dir      = self.case['case_path']
        self.ref_case_name = os.path.split(case_dir)[-1]
        self.otstudy_path  = os.path.split(case_dir)[0]

        # path to OpenTurns study config file
        self.cfg_ot   = None
        self.cfg_path = None

        self.__get_cfg()

        # Test for uninitialized cfg file
        if not self.cfg_ot.has_section('host_parameters'):
            self.__set_default_cfg(self.cfg)

        #FIXME:
        self.resfile_name = "output.dat"

        # ------------------------------------
        # Read cfg info:

        # HOST
        self.host_name     = self.cfg_ot.get('host_parameters', 'host_name')
        self.batch_manager = self.cfg_ot.get('host_parameters', 'batch_manager')
        self.build_name    = self.cfg_ot.get('host_parameters', 'build_name')
        self.arch_path     = self.cfg_ot.get('host_parameters', 'arch_path')

        # STUDY
        self.otstudy_name = self.cfg_ot.get('study_parameters', 'study_name')
        # Check if the study name needs update
        tmp_study_name = os.path.split(self.otstudy_path)[-1]
        if self.otstudy_name != tmp_study_name:
            self.otstudy_name = tmp_study_name

        self.resfile_name = self.cfg_ot.get('study_parameters', 'results_file')

        # SUBMISSION
        self.wall_clock   = self.cfg_ot.get('batch_parameters', 'wall_clock')
        self.nprocs       = int(self.cfg_ot.get('batch_parameters', 'nprocs'))
        self.wckey        = self.cfg_ot.get('batch_parameters', 'wckey')
        if self.batch_manager == 'none':
            self.nnodes    = 'none'
            self.ntasks    = 'none'
            self.nthreads  = 'none'
            self.dist_wdir = 'none'
        else:
            self.nnodes    = self.cfg_ot.get('batch_parameters', 'number_of_nodes')
            self.ntasks    = self.cfg_ot.get('batch_parameters', 'tasks_per_node')
            self.nthreads  = self.cfg_ot.get('batch_parameters', 'threads_per_process')
            self.dist_wdir = self.cfg_ot.get('batch_parameters', 'distant_workdir')


    # ---------------------------------------
    def __set_cfg_file_path(self):
        """
        Set the path to the OpenTurns study cfg file
        """

        self.cfg_path = os.path.join(self.otstudy_path, "openturns_study.cfg")


    # ---------------------------------------
    def __get_cfg(self):
        """
        get the config file of the study
        """
        if not self.cfg_path:
            self.__set_cfg_file_path()

        config = configparser.ConfigParser()

        cfg_found = config.read(self.cfg_path)

        if not cfg_found:
            f = open(self.cfg_path, "w")
            f.close()

            self.__set_default_cfg(config)

            f = open(self.cfg_path, "w")
            config.write(f)
            f.close()

        self.cfg_ot = config

    # ---------------------------------------
    def __set_default_cfg(self, cfg):
        """
        set default options for the config file
        """

        cfg.add_section('host_parameters')
        cfg.set('host_parameters', 'host_name', 'localhost')
        cfg.set('host_parameters', 'batch_manager', 'none')
        cfg.set('host_parameters', 'build_name', 'default')
        cfg.set('host_paramaters', 'arch_path', 'default')

        cfg.add_section('study_parameters')
        case_dir = self.case['case_path']
        cfg.set('study_parameters', 'code_name', self.case['package'].name)

        tmp_study_name = os.path.split(self.otstudy_path)[-1]
        cfg.set('study_parameters', 'study_name', tmp_study_name)
        cfg.set('study_parameters', 'ref_case', 'REF_CASE')
        cfg.set('study_parameters', 'xmlfile', 'setup.xml')
        cfg.set('study_parameters', 'results_file', 'output.dat')


        cfg.add_section('batch_parameters')
        cfg.set('batch_parameters', 'wall_clock', '0-00:10:00')
        cfg.set('batch_parameters', 'wckey', 'SATURNE')
        cfg.set('batch_parameters', 'nprocs', 1)
        cfg.set('batch_parameters', 'number_of_nodes', 1)
        cfg.set('batch_parameters', 'tasks_per_node', 1)
        cfg.set('batch_parameters', 'threads_per_process', 1)
        cfg.set('batch_parameters', 'distant_workdir', 'DEFAULT')

    # ---------------------------------------
    def setHostName(self, host_name):
        """
        set host name
        """

        self.host_name = host_name

    def getHostName(self):
        """
        get host name
        """

        return self.host_name

    # ---------------------------------------
    def setBuildName(self, build_name):
        """
        set build name
        """

        self.build_name = build_name


    def getBuildName(self):
        """
        Return build name
        """

        return self.build_name


    # ---------------------------------------
    def setBatchManager(self, bmgr):
        """
        set batch manager
        """

        self.batch_manager = bmgr

    def getBatchManager(self):
        """
        returns the host batch manager
        """

        return self.batch_manager

    # ---------------------------------------
    def setOtStudyName(self, sname):
        """
        set OpenTurns Study name
        """

        self.otstudy_name = sname


    def getOtStudyName(self):

        return self.otstudy_name

    # ---------------------------------------
    def setDistWorkdir(self, dwdir):
        """
        Set the distant workdir path
        """

        self.dist_wdir = dwdir


    def getDistWorkdir(self):

        return self.dist_wdir

    # ---------------------------------------
    def setWallClockTime(self, d, h, m, s='00'):
        """
        set the wall clock time
        """

        self.wall_clock = str(d)+'-'+str(h)+':'+str(m)+':'+str(s)

    def getWallClockTime(self):
        """
        get wall clock time
        """

        wct = self.wall_clock.split('-')
        d = wct[0]

        h = wct[1].split(':')[0]
        m = wct[1].split(':')[1]
        s = wct[1].split(':')[2]

        return d, h, m, s

    # ---------------------------------------
    def setInputVars(self, vlist):
        """
        sets the input variables list
        """

        if type(vlist) is list:
            self.input_variables = vlist
        else:
            raise Exception('The input variables needed format is a list')


    def getInputVars(self):
        """
        returns the list of input variables
        """

        return self.input_variables


    # ---------------------------------------
    def setOutputVars(self, vlist):
        """
        sets the output variables list
        """

        if type(vlist) is list:
            self.output_variables = vlist
        else:
            raise Exception('The output variables needed format is a list')


    def getOutputVars(self):
        """
        returns the list of output variables
        """

        return self.output_variables

    # ---------------------------------------
    def setNprocs(self, nprocs):
        """
        Set the number of wanter processors
        """

        self.nprocs = nprocs


    def getNprocs(self):
        """
        Get the number of wanter processors
        """

        return self.nprocs

    # ---------------------------------------
    def setClusterParams(self, nnodes=None, ntasks=None, nthreads=None):
        """
        Set cluster paramaters: nnodes, ntasks and/or nthreads
        """

        if nnodes:
            self.nnodes = nnodes

        if ntasks:
            self.ntasks = ntasks

        if nthreads:
            self.nthreads = nthreads

        self.nprocs = nnodes * ntasks * nthreads

    def getClusterParams(self):
        """
        Get cluster paramaters: nnodes, ntasks and nthreads
        """

        return self.nnodes, self.ntasks, self.nthreads

    # ---------------------------------------
    def update_ot_variables(self):
        """
        Retrieve from the notebook the list of input and output OpenTurns
        variables
        """

        nb = NotebookModel(self.case)

        iv = []
        ov = []

        for v in nb.getVarList():
            otv = nb.getVariableOt(v)

            if otv[0:3] == "Yes":
                if otv[5:] == "Input":
                    iv.append(nb.getVariableName(v))
                elif otv[5:] == "Output":
                    ov.append(nb.getVariableName(v))

        self.input_variables  = iv
        self.output_variables = ov

    # ---------------------------------------
    def update_cfg_file(self):
        """
        Update the OpenTurns study cfg file
        """

        self.update_ot_variables()

        if self.cfg_ot:

            self.cfg_ot.set('host_parameters', 'host_name', self.host_name)
            self.cfg_ot.set('host_parameters', 'batch_manager', self.batch_manager)
            self.cfg_ot.set('host_parameters', 'build_name', self.build_name)
            self.cfg_ot.set('host_parameters', 'arch_path', self.arch_path)

            self.cfg_ot.set('study_parameters', 'study_name', self.otstudy_name)
            self.cfg_ot.set('study_parameters', 'ref_case', self.ref_case_name)
            self.cfg_ot.set('study_parameters', 'xmlfile',
                            os.path.split(self.case['xmlfile'])[-1])

            self.cfg_ot.set('batch_parameters', 'wall_clock', self.wall_clock)
            self.cfg_ot.set('batch_parameters', 'wckey', self.wckey)
            self.cfg_ot.set('batch_parameters', 'nprocs', self.nprocs)
            self.cfg_ot.set('batch_parameters', 'number_of_nodes', self.nnodes)
            self.cfg_ot.set('batch_parameters', 'tasks_per_node', self.ntasks)
            self.cfg_ot.set('batch_parameters', 'threads_per_process', self.nthreads)
            self.cfg_ot.set('batch_parameters', 'distant_workdir', self.dist_wdir)

            f = open(self.cfg_path, 'w')
            self.cfg_ot.write(f)
            f.close()

#-------------------------------------------------------------------------------
# End of OpenTurnsModel
#-------------------------------------------------------------------------------
