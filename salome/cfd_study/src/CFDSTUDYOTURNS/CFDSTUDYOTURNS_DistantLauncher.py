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

import salome
import os
from glob import glob
import launcher_proxy as jobmanager
"""
Distant launcher
================

Contains the functions needed for a distant launching of Code_Saturne on
another host (Computing cluster).
"""

# ==============================================================================
class CFDSTUDY_DistantLauncher:
    """
    Distant Code_Saturne launcher class. Provides all the necessary tools for
    distant runs.
    """

    # --------------------------------------------------------------------------
    def __init__(self, case_dir, params_cfg, package, run_prefix=None):
        """
        Class constructor
        """

        # Initialize the salome services
        salome.salome_init()

        self.case_dir        = case_dir
        self.case_name       = os.path.split(self.case_dir)[1]

        self.study_dir       = os.path.split(self.case_dir)[0]
        self.study_name      = os.path.split(self.study_dir)[1]

        self.pkg             = package

        self.cfg             = params_cfg
        self.host            = self.cfg.get('host_parameters', 'host_name')
        self.host_arch_path  = self.cfg.get('host_parameters', 'arch_path')
        self.host_build_name = self.cfg.get('host_parameters', 'build_name')
        self.host_bin_path   = os.path.join(self.host_arch_path,
                                            self.host_build_name,
                                            'bin')

        self.dist_wdir       = self.cfg.get('batch_parameters', 'distant_workdir')
        self.package_name    = self.cfg.get('study_parameters', 'code_name')
        self.paramfile       = self.cfg.get('study_parameters', 'xmlfile')
        self.results_file    = "cs_uncertain_output.dat"

        self.run_prefix = None
        if run_prefix:
            self.run_prefix = run_prefix + "_"

        self.run_id = None
        tmp = os.getcwd()
        os.chdir(self.case_dir)
        self.__setRunId__()
        os.chdir(tmp)

        self.job_params = self.__getJobParameters__()

        self.job_id = -1
        self.job_failed = False
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def launch(self, force_submit=False):
        """
        Launch and wait for job results
        """

        self.launch_job(force_submit=force_submit)
        self.wait_on_job()
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def launch_job(self, force_submit=False):
        """
        Launch the distant job
        """

        if self.job_id == -1 or force_submit:

            self.__createPrepareScript__()
            os.chmod(self.job_params.pre_command, 493)
            self.__createRuncaseScript__()
            self.job = jobmanager.Job.launch(self.job_params)
            self.job_id = self.job.job_id
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def wait_on_job(self):
        """
        Wait for job termination if needed
        """

        self.job.wait()
        exit_code = self.job.verify()

        if exit_code != 0:
            self.job_failed = True

    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def sync_results(self):
        """
        Sync the results files
        """

        self.job.getResults()
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def need_restart(self):
        """
        Check if the run requires a restart (insufficient time)
        """

        test_name = os.path.join(self.job_params.result_directory,
                                 "run_status.exceeded_time_limit")

        not_finished = os.path.isfile(test_name)

        return not_finished
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def relaunch_job(self):
        """
        Relaunch the job if restart is needed
        """

        self.__prepareRestart__()
        self.job.relaunch()

        return

    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def reload_job(self, resu_dir):
        """
        Reload the job into jobmanager
        """

        self.job = jobmanager.Job.reloadJob(resu_dir)
        self.job_id = self.job.job_id
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def __getJobParameters__(self):
        """
        Sets and returns the JobManager single job parameters.
        """

        # ---------------------------------
        resManager = salome.lcc.getResourcesManager()

        job_params = salome.JobParameters()
        # ---------------------------------

        # ---------------------------------
        job_params.resource_required         = salome.ResourceParameters()
        job_params.resource_required.name    = self.host
        job_params.resource_required.nb_proc = int(self.cfg.get('batch_parameters', 'nprocs'))
        job_params.resource_required.type    = 'rsync'

        # Jobmanager wall clock format is hh:mm !
        wall_clock = self.cfg.get('batch_parameters','wall_clock')
        days, hms = wall_clock.split("-")
        wch, wcm, wcs = hms.split(':')

        jp_wc = "%d:%s:%s" % (int(days)*24+int(wch), wcm, wcs)
        job_params.maximum_duration = jp_wc

        job_params.wckey            = self.cfg.get('batch_parameters', 'wckey')
        job_params.job_name         = "CS_OT"
        # ---------------------------------

        # ---------------------------------
        job_params.work_directory   = os.path.join(self.dist_wdir,
                                                   self.study_name,
                                                   self.case_name)

        job_params.local_directory = os.path.join(self.case_dir,
                                                  'RESU',
                                                  self.run_id)

        job_params.result_directory = job_params.local_directory
        # ---------------------------------

        # ---------------------------------
        job_params.in_files = [ os.path.join(self.case_dir, 'DATA'),
                                os.path.join(self.case_dir, 'SRC'),
                                os.path.join(self.case_dir, 'SCRIPTS'),
                                os.path.join(self.case_dir, 'RESU') ]
        # ---------------------------------

        # ---------------------------------
        job_params.out_files = []
        for f in (self.results_file, 'run_status.exceeded_time_limit'):
            df = os.path.join(job_params.work_directory,
                              'RESU',
                              self.run_id,
                              f)

            job_params.out_files.append(df)
        # ---------------------------------

        # ---------------------------------
        job_params.job_type = 'command'
        job_params.pre_command = os.path.join(self.case_dir, 'prepare_cs_case.sh')

        job_params.job_file = os.path.join(self.case_dir, 'run_cs_case.sh')
        # ---------------------------------

        return job_params
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def __setRunId__(self):
        """
        Set the current run id. Needed for the distant launching scripts
        """

        if self.run_prefix:
            idx = 0

            # Check if case exists if prefix is forced
            search_dir = os.path.join(self.case_dir,
                                      'RESU',
                                      self.run_prefix) + '*'
            dirlist = glob(search_dir)
            dirlist.sort()
            if len(dirlist) > 0:
                idx = int(dirlist[-1].split("_")[-1])

            self.run_id = self.run_prefix + '{:0>4}'.format(idx + 1)

        else:

            from code_saturne.cs_run import run as get_run_id

            id_args = ['run', '--suggest-id']
            self.run_id = get_run_id(id_args, self.pkg)[1]

    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def __createPrepareScript__(self):
        """
        Creating the prepro script for the cluster, which will run on the front
        nodes.
        """

        run_args = [self.package_name, 'run', '--stage',
                    '--id ', self.run_id,
                    '-p', self.paramfile]
        f = open('prepare_cs_case.sh', 'wt')

        f.write('cd SCRIPTS\n\n')
        export_line = 'export PATH=' + self.host_bin_path + ':$PATH\n\n'

        f.write('# Ensure the correct command is found:\n')
        f.write(export_line)
        f.write('# Run command:\n\n')

        run_line = ' '.join(run_args)
        f.write(run_line + '\n\n')
        f.write('cd ..\n')

        f.close()
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def __createRuncaseScript__(self):
        """
        Creating the main run script.
        """
        run_args = [self.package_name, 'run', '-p', self.paramfile,
                    '--initialize', '--finalize', '--id', self.run_id]

        f = open('run_cs_case.sh', 'wt')
        f.write('cd SCRIPTS\n\n')

        export_line = 'export PATH=' + self.host_bin_path + ':$PATH\n\n'
        f.write('# Ensure the correct command is found:\n')
        f.write(export_line)
        f.write('# Run command:\n')

        run_line = ' '.join(run_args)
        f.write(run_line)

        f.close()
    # --------------------------------------------------------------------------

# ==============================================================================
