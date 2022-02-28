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
This module modifies the run_conf object.
- BatchRunningModel
"""
#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys, types

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from code_saturne.base import cs_batch
from code_saturne.base import cs_run_conf

#-------------------------------------------------------------------------------
# Class BatchRunningModel
#-------------------------------------------------------------------------------

class BatchRunningModel(object):
    """
    This class modifies the job info (in run.cfg)
    """
    def __init__(self, path=None, pkg=None, import_legacy=False):
        """
        Constructor.
        """

        self.pkg = pkg
        if self.pkg is None:
            from code_saturne.base.cs_package import package
            pkg = package()

        self.run_conf = None
        self.path = path

        # Configuration-based information

        i_c = cs_run_conf.get_install_config_info(self.pkg)

        self.resource_name = cs_run_conf.get_resource_name(i_c)
        self.compute_builds = i_c['compute_builds']

        self.batch = cs_batch.batch(self.pkg)

        # Convert from legacy runcase if not updated yet
        # (delaying application of file changes to save).
        # in this case, the run_conf object is pre-loaded so that
        # the "save" operation can apply the (conversion) changes
        # even when the configuration file has not been loaded.

        self.runcase_path = None
        if self.path:
            if import_legacy and not os.path.isfile(self.path):
                dirname = os.path.dirname(self.path)
                if os.path.basename(dirname) == 'DATA':
                    dirname = os.path.join(os.path.basename(dirname), 'SCRIPTS')
                runcase_path = os.path.join(dirname, 'runcase')
                if os.path.isfile(runcase_path):
                    self.runcase_path = runcase_path

        if self.runcase_path:
            from code_saturne.base import cs_runcase
            runcase = cs_runcase.runcase(runcase_path,
                                         package=self.pkg)
            sections = runcase.run_conf_sections(resource_name=self.resource_name,
                                                 batch_template=i_c['batch'])

            self.run_conf = cs_run_conf.run_conf(self.path,
                                                 package=self.pkg,
                                                 create_if_missing=True)
            for sn in sections:
                if not sn in self.run_conf.sections:
                    run_conf.sections[sn] = {}
                for kw in sections[sn]:
                    self.run_conf.sections[sn][kw] = sections[sn][kw]

    #---------------------------------------------------------------------------

    def __is_changed__(self):
        """
        Check if values have been changed relative to the initial load.
        """

        if not self.run_conf:
            return False

        changed = False

        if self.job_header_lines_ini != self.job_header_lines:
            changed = True

        if not changed:
            t = ('1', 't', 'true', 'y', 'yes')
            f = ('0', 'f', 'false', 'n', 'no')
            for k in self.run_dict:
                s1 = str(self.run_dict_ini[k]).lower()
                s2 = str(self.run_dict[k]).lower()
                if s1 in f and s2 in f:
                    continue
                elif s1 in t and s2 in t:
                    continue
                elif s1 != s2:
                    changed = True
                    break

        if not changed:
            for k in self.job_dict:
                if str(self.job_dict_ini[k]) != str(self.job_dict[k]):
                    changed = True
                    break

        return changed

    #---------------------------------------------------------------------------

    def __update__(self):
        """
        Update the associated dictionnaries and run_conf object.
        """

        if not self.have_mpi:
            self.job_dict['n_procs'] = None
        if not self.have_openmp:
            self.job_dict['n_threads'] = None

        if self.job_dict['debug_args'] == '':
            self.job_dict['debug_args'] = None

        if not self.run_conf:
            return

        # Run information

        for k in self.run_dict:
            # Special case for "initialize": if False, we mean we do not specify
            # this stage (since the others are not specified either)
            if k == 'initialize':
                if self.run_dict[k] == False:
                    self.run_conf.set('run', 'initialize', None)
                else:
                    self.run_conf.set('run', k, self.run_dict[k])
            else:
                self.run_conf.set('run', k, self.run_dict[k])

        # Resource-related information

        if self.job_header_lines != None:
            self.batch.update_lines(self.job_header_lines)
            self.run_conf.set(self.resource_name, 'job_header',
                              os.linesep.join(self.job_header_lines))

        for k in self.job_dict:
            if self.job_dict[k] is None:
                self.run_conf.set(self.resource_name, k, None)
            else:
                self.run_conf.set(self.resource_name, k, str(self.job_dict[k]))

    #---------------------------------------------------------------------------

    def load(self):
        """
        Load or the associated run_conf object if not already done.
        """

        if self.run_conf:
            return

        # Load or build run configuration

        self.run_conf = cs_run_conf.run_conf(self.path,
                                             package=self.pkg,
                                             create_if_missing=True)

        self.run_dict = {}
        self.job_dict = {}

        # Generic job running information (subset of possible "run" info)

        self.run_dict['id'] = self.run_conf.get('run', 'id')
        self.run_dict['compute_build'] = self.run_conf.get('run', 'compute_build')

        self.run_dict['initialize'] = self.run_conf.get_bool('run', 'initialize')
        self.run_dict['compute'] = None
        self.run_dict['finalize'] = None

        # Resource-specific info (subset of resource-based info, and batch)

        self.job_dict['n_procs'] = self.run_conf.get_int(self.resource_name,
                                                         'n_procs')
        self.job_dict['n_threads'] = self.run_conf.get_int(self.resource_name,
                                                           'n_threads')
        self.job_dict['debug_args'] = self.run_conf.get(self.resource_name,
                                                        'debug_args')

        self.job_header_lines = None

        if self.batch.rm_type:
            job_header = self.run_conf.get(self.resource_name, 'job_header')
            if not job_header:
                self.run_conf.rebuild_resource()
                job_header = self.run_conf.get(self.resource_name, 'job_header')
            if job_header != None:
                self.job_header_lines = job_header.split(os.linesep)

        if self.job_header_lines != None:
            self.batch.parse_lines(self.job_header_lines)

        # Save initial values (to determine which are changed)

        self.setup_ini = self.run_dict.get('setup', 'param')

        self.run_dict_ini = {}
        for k in self.run_dict:
            self.run_dict_ini[k] = self.run_dict[k]

        self.job_dict_ini = {}
        for k in self.job_dict:
            self.job_dict_ini[k] = self.job_dict[k]

        self.job_header_lines_ini = None
        if self.job_header_lines:
            self.job_header_lines_ini = list(self.job_header_lines)

        # Fix incompatible options if needed
        # (add entries to dictionnaries so __is_changed__ will detect update)

        if self.run_conf.get_bool('run', 'stage') == False:
            if not self.run_dict_ini['id']:
                self.run_dict_ini['stage'] = False
                self.run_dict['stage'] = None

        # Query info related to compute build

        self.updateComputeBuildInfo(self.run_dict['compute_build'])

    #---------------------------------------------------------------------------

    def updateComputeBuildInfo(self, compute_build=None):
        """
        Update relative to compute build
        """

        pkg_compute = None

        if self.compute_builds != None:
            if not compute_build:
                if len(self.compute_builds) > 0:
                    compute_build = self.compute_builds[0]
            if compute_build in self.compute_builds:
                pkg_compute = self.pkg.get_alternate_version(compute_build)

        if pkg_compute is None:
            pkg_compute = self.pkg

        config_features = pkg_compute.config.features
        if config_features['mpi'] == 'yes':
            self.have_mpi = True
        else:
            self.have_mpi = False
        if config_features['openmp'] == 'yes':
            self.have_openmp = True
        else:
            self.have_openmp = False

        if not self.have_mpi:
            self.job_dict['n_procs'] = None
        if not self.have_openmp:
            self.job_dict['n_threads'] = None

    #---------------------------------------------------------------------------

    def save(self, path=None, param=None, force=True):
        """
        Update the associated run_conf object and save file
        """

        if not self.run_conf:
            if  param and path:
                self.path = path
                if os.path.isfile(path):
                    self.load()
            if not self.run_conf:
                return
        else:
            try:
                self.__update__()
            except Exception:
                print("Error: can not update %s\n" % self.path)
                print("Probably not in a standard case directory structure.")
                return

        if param != None:
            param = os.path.basename(param)
        if param == 'setup.xml':
            param = None

        self.run_conf.set('setup', 'param', param)

        save_as = False

        if path != None:
            if self.path != None:
                if self.path != path:
                    save_as = True
                else:
                    save_as = False
        else:
            path = self.path

        if save_as:
            self.run_conf.save(path, new=True)

        else:

            if force == False:
                force = self.__is_changed__()

            if force:
                self.run_conf.save(path)

                if self.runcase_path != None:
                    os.remove(self.runcase_path)
                    try:
                        os.remove(os.path.dirname(self.runcase_path))
                    except Exception:
                        pass

        self.path = path

#-------------------------------------------------------------------------------
# End of BatchRunningModel
#-------------------------------------------------------------------------------
