# ==============================================================================
# -*- coding: utf-8 -*-
#
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
# ==============================================================================
"""
"""

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, os.path
try:
    import ConfigParser as configparser # Python2
except Exception:
    import configparser  # Python3

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


# ==============================================================================
# OpenTurns Study Model class
class cfd_openturns_study:
    """
    Class used for an OpenTurns study within SALOME_CFD
    Public methods:
    - init
    - run
    - study2code
    - code2study
    """
# ------------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# ------------------------------------------------------------------------------

    # ---------------------------------------
    def __init__(self, study_path, study_cfg, vars_dico):
        '''
        Classs constructor
        '''

        # --------------------------
        # Init from user values
        self.study_path = study_path

        self.cfg_name = study_cfg
        cfg_path      = os.path.join(study_path, study_cfg)
        self.cfg_path = cfg_path

        print("-----------------------------------")
        print(cfg_path)
        print("-----------------------------------")

        self.cfg = configparser.ConfigParser()
        self.cfg.read(cfg_path)

        self.vars_dico = vars_dico

        self.cs_launcher = None

        self.ref_case  = self.cfg.get('study_parameters', 'ref_case')
        if self.cfg.get('host_parameters', 'host_name') == 'localhost':
            self.run_type = 'local'
        else:
            self.run_type = 'distant'

        # --------------------------
        self.__setPackage__()

        self.__setCaseId__()
    # ---------------------------------------

    # ---------------------------------------
    def load_launcher(self, case_dir):

        new_launcher = None

        if self.run_type == 'local':
            from CFDSTUDYOTURNS_LocalLauncher import cfd_openturns_local_launcher

            new_launcher = cfd_openturns_local_launcher(case_dir   = case_dir,
                                                        params_cfg = self.cfg,
                                                        package    = self.pkg)
        elif self.run_type == 'distant':
            from CFDSTUDYOTURNS_DistantLauncher import CFDSTUDY_DistantLauncher

            new_launcher = CFDSTUDY_DistantLauncher(case_dir   = case_dir,
                                                    params_cfg = self.cfg,
                                                    package    = self.pkg)

        self.cs_launcher = new_launcher
    # ---------------------------------------

    # ---------------------------------------
    def study2code(self):
        '''
        A function which takes into account the variables values given by
        OpenTurns and creating a Code_Saturne case thanks to it.
        '''
        from glob import  glob

        case_dir = os.path.join(self.study_path, self.case_id)

        if os.path.isdir(case_dir):
            # Check if the case directory allready exists.
            # if it exists, check if a case has allready been run
            case_resu_dir = os.path.join(case_dir, "RESU")

            dirs = filter(os.path.isdir, glob(case_resu_dir + "/*"))

            self.__setCsCase__()
            self.case['case_path'] = case_dir
            self.load_launcher(case_dir)

            if self.run_type == 'distant' and len(dirs) > 0:
                dirs.sort(key=lambda x: os.path.getmtime(x))
                resu_dir = dirs[-1]
                self.cs_launcher.run_id = os.path.split(resu_dir)[-1]
                self.cs_launcher.reload_job(resu_dir)
        else:
            self.__createCase__()
            self.load_launcher(case_dir)

    # ---------------------------------------

    # ---------------------------------------
    def run(self):
        '''
        This method launches a Code_Saturne computation based on the
        values provided by OpenTurns
        '''

        self.cs_launcher.launch(force_submit=False)

        if self.run_type == 'distant':
            while True:
                if self.cs_launcher.need_restart():
                    old_run = self.cs_launcher.run_id

                    rst_dir = os.path.join(self.study_path,
                                           self.case_id,
                                           'RESU',
                                           old_run)
                    self.__setRestartPath__(rst_dir)
                    self.cs_launcher.__setRunId__()

                    self.cs_launcher.launch(force_submit=True)

                    if self.cs_launcher.job_failed:
                        break
                else:
                    break

            self.cs_launcher.sync_results()

    # ---------------------------------------

    # ---------------------------------------
    def code2study(self, n_values=1):
        """
        Opens the results file specified by the user.
        Returns a tuple of all the required values (needed by OpenTurns)
        """


        resfile = self.cfg.get('study_parameters', 'results_file')

        rspth = os.path.join(self.study_path,
                             self.case_id,
                             'RESU',
                             self.cs_launcher.run_id,
                             resfile)

        results = ()
        if os.path.exists(rspth):
            r = open(rspth, 'r').readlines()
            d = r[-1].split()
            for ed in d:
                results += (float(ed),)

        if len(results) == 0:
            results = (float('nan'),)*n_values


        return results;
    # ---------------------------------------

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------
# INTERNAL UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

    # ---------------------------------------
    def __setCaseId__(self):
        """
        Setting the case id based on the variables names and values
        """

        self.case_id = ''

        items = []
        for key in self.vars_dico.keys():
            items.append(str(key[:4])+'_'+str(self.vars_dico[key]))

        self.case_id = "_".join(items)
    # ---------------------------------------

    # ---------------------------------------
    def __setPackage__(self):
        """
        Setting the package according to the requested code by the user
        """

        pkg_name = self.cfg.get('study_parameters', 'code_name')

        if pkg_name == 'code_saturne':
            from code_saturne.cs_package import package
        elif pkg_name == 'neptune_cfd':
            from neptune_cfd.nc_package import package
        else:
            raise Exception("Uknown package: "+pkg_name)

        self.pkg = package()

    # ---------------------------------------

    # ---------------------------------------
    def __createCase__(self):
        """
        Creating the case for a given evaluation based on the user defined
        reference case
        """
        from code_saturne.model.NotebookModel import NotebookModel

        # Creating the case from the ref case
        from code_saturne.cs_script import master_script

        case_create_args = ['create',
                            '-c',
                            os.path.join(self.study_path, self.case_id),
                            '--copy-from',
                            os.path.join(self.study_path, self.ref_case) ]

        ms = master_script(case_create_args, self.pkg)
        retcode = ms.execute()

        self.__setCsCase__()
        self.__setOtVarsVals__()
    # ---------------------------------------

    # ---------------------------------------
    def __setCsCase__(self):
        """
        Initialize the cs_case structure according to the user provided
        data
        """

        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit

        paramfile = self.cfg.get('study_parameters', 'xmlfile')
        fp = os.path.join(self.study_path, self.case_id, 'DATA', paramfile)

        print("=======================================")
        print(fp)
        print("=======================================")

        self.case = Case(package=self.pkg, file_name=fp)
        self.case['xmlfile'] = fp
        self.case.xmlCleanAllBlank(self.case.xmlRootNode())
        if self.case['package'].name == 'code_saturne':
            from code_saturne.model.XMLinitialize import XMLinit
        else:
            from code_saturne.model.XMLinitializeNeptune import XMLinitNeptune as XMLinit
        XMLinit(self.case).initialize()


    # ---------------------------------------

    # ---------------------------------------
    def __setOtVarsVals__(self):
        """
        This method translates the values provided by OpenTURNS into the
        code xml input file
        """

        # Update xml file
        from code_saturne.model.NotebookModel import NotebookModel
        nb = NotebookModel(self.case)
        nb_ids = nb.getVarList()
        nb_names = nb.getVarNameList()
        for key in self.vars_dico.keys():
            try:
                ix = nb_names.index(key)
            except:
                ix = -1

            if ix != -1:
                idx = nb_ids[ix]
                nb.setVariableValue(idx=idx, val=str(self.vars_dico[key]))

        self.case.xmlSaveDocument()
    # ---------------------------------------

    # ---------------------------------------
    def __setRestartPath__(self, restart_path):
        """
        This method sets the restart path in the code xml input file
        """

        from code_saturne.model.StartRestartModel import StartRestartModel

        rs = StartRestartModel(self.case)

        rp = os.path.join(restart_path, 'checkpoint')
        rs.setRestartPath(rp)
        rs.setRestartWithAuxiliaryStatus('on')

        self.case.xmlSaveDocument()

    # ---------------------------------------

# ------------------------------------------------------------------------------

# ==============================================================================
