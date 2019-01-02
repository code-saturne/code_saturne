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

import os

"""
Local launcher
==============

Contains the functions needed for a local launching of Code_Saturne
during an OpenTURNS study.
"""

# ==============================================================================
class cfd_openturns_local_launcher:
    """
    Local Code_Saturne launcher class. Provides all the necessary tools for
    the local study.
    """

    # --------------------------------------------------------------------------
    def __init__(self, case_dir, params_cfg, package, run_prefix=None):
        """
        Class constructor
        """

        self.case_dir  = case_dir
        self.case_name = os.path.split(self.case_dir)[-1]
        self.pkg       = package

        self.cfg       = params_cfg

        self.package_name = self.cfg.get('study_parameters', 'code_name')
        self.paramfile    = self.cfg.get('study_parameters', 'xmlfile')

        self.nprocs       = self.cfg.get('batch_parameters', 'nprocs')
#        self.study_path   = os.path.split(case_dir)[0]

        self.run_prefix = None
        if run_prefix:
            self.run_prefix = run_prefix + "_"

        self.run_id = None
        self.__setRunId__()
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    def launch(self, force_submit=False):

        script_dir = os.path.join(self.case_dir, 'SCRIPTS')
        os.chdir(script_dir)

        from code_saturne.cs_script import master_script

        run_args = ['run',
                    '-p',
                    self.paramfile,
                    '-n',
                    str(self.nprocs)]

        ms = master_script(run_args, self.pkg)

        retcode = ms.execute()
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
            from code_saturne.cs_script import master_script
            id_args = ['run', '--suggest-id']
            ms = master_script(id_args, self.pkg)
            retcode = ms.execute()

            self.run_id = retcode[0]

    # --------------------------------------------------------------------------

# ==============================================================================
