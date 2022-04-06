# -*- coding: utf-8 -*-
# Package information, generated from cs_package.py.in by make.

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

import configparser
import os
import sys

try:
    from code_saturne.base import cs_config
except Exception:
    import cs_config

r_config_file_path = None
r_config_dict = {}

# Package information
#--------------------

class package:

    def __init__(self,
                 scriptdir=None,
                 alternate_config=False,
                 config_file=r_config_file_path,
                 install_prefix=None,
                 install_mode=False,
                 name='code_saturne'):

        # Global variables

        global r_config_file_path
        global r_config_dict

        # Load configration from file if needed and not already done.

        config_dict = r_config_dict

        if config_file and config_file != r_config_file_path:
            config_parser = configparser.ConfigParser()
            config_parser.read(config_file)
            config_dict = dict(config_parser)
            if alternate_config == False:  # Cache result in case of new call
                r_config_file_path = config_file
                r_config_dict = config_dict

        base_package = 'code_saturne'

        # System configuration (compilers, pre-requisites, ...)

        self.config = cs_config.config(config_dict)

        # Package information
        # -------------------

        self.name = name

        d = config_dict.get('package', {})

        for k in ('pkgversion', 'string', 'bugreport', 'url',
                  'version', 'version_full', 'version_short',
                  'revision'):
            fallback = '<' + str(k) + ':undefined>'
            s = d.get(k, fallback)
            setattr(self, k, s)

        self.code_name = "Code_Saturne"

        self.preprocessor = "cs_preprocess" + self.config.exeext
        self.solver = self.config.solver_modules[name]['solver']
        self.check_syntax = "cs_check_syntax" + self.config.exeext
        self.io_dump = "cs_io_dump" + self.config.exeext
        self.runcase = "runcase" + self.config.shext
        self.runsolver = "run_solver" + self.config.shext
        self.configfile = "code_saturne" + self.config.cfgext
        self.scratchdir = 'tmp_Saturne'

        # Installation directories
        # ------------------------

        if install_prefix is None:
            install_prefix = os.getenv('CS_ROOT_DIR')

        python_version = "%d.%d" % (sys.version_info.major, sys.version_info.minor)
        pythondir_rel = os.path.join('lib', 'python' + python_version,
                                     'site-packages')

        # Otherwise, assume that the standard tree structure is used
        if install_prefix is None:
            install_prefix = os.path.dirname(os.path.realpath(__file__))
            i = install_prefix.find(pythondir_rel)
            if i > -1:
                install_prefix = os.path.dirname(install_prefix[:i+1])

        exec_prefix = install_prefix
        pythondir = os.path.join(exec_prefix, pythondir_rel)
        datarootdir = os.path.join(install_prefix, 'share')

        self.dirs = {'prefix': install_prefix,
                     'exec_prefix': exec_prefix,
                     'bindir': os.path.join(exec_prefix, 'bin'),
                     'includedir': os.path.join(exec_prefix, 'include'),
                     'pkgincludedir': os.path.join(install_prefix, 'include',
                                                   base_package),
                     'libdir': os.path.join(exec_prefix, 'lib'),
                     'libexecdir': os.path.join(exec_prefix, 'libexec'),
                     'pkglibexecdir': os.path.join(exec_prefix, 'libexec',
                                                   base_package),
                     'pythondir': pythondir,
                     'pkgpythondir': os.path.join(pythondir, base_package),
                     'localedir': os.path.join(datarootdir, 'locale'),
                     'datarootdir': datarootdir,
                     'datadir': datarootdir,
                     'pkgdatadir': os.path.join(datarootdir, base_package),
                     'docdir': os.path.join(datarootdir, 'doc', base_package),
                     'sysconfdir': os.path.join(install_prefix, 'etc')}

        # If the build is not relocatable, read paths from configuration

        if self.config.features['relocatable'] == "no" or install_mode:
            d = config_dict.get('install', {})
            for k in d.keys():
                self.dirs[k] = d[k]

        # Adjust docdir for additional modules.

        # We should try to find a cleaner/more consistant solution,
        # perhaps storing all documents in the main module's
        # directory (code_saturne).

        docdir_1 = self.dirs['docdir']
        docdir_1 = os.path.join(os.path.split(docdir_1)[0], name)
        self.dirs['docdir'] = os.path.join(os.path.split(docdir_1)[0], name)


    def get_dir(self, installdir):

        return self.dirs[installdir]

    def get_preprocessor(self):

        if self.config.features['frontend'] == "no":
            raise Exception("This " + self.name + " build does not " + \
                            "include the front-end and Preprocessor.")
        return os.path.join(self.get_dir("pkglibexecdir"),
                            self.preprocessor)

    def get_io_dump(self):

        return os.path.join(self.get_dir("pkglibexecdir"),
                            self.io_dump)

    def get_solver(self):

        return os.path.join(self.get_dir("pkglibexecdir"),
                            self.solver)

    def get_check_syntax(self):

        return os.path.join(self.get_dir("pkglibexecdir"),
                            self.check_syntax)

    def get_io_dump(self):

        return os.path.join(self.get_dir("pkglibexecdir"),
                            self.io_dump)

    def get_global_configfile(self):

        # Windows:         C:\ProgramData
        # Linux and co:    /etc

        if sys.platform.startswith("win"):
            configdir = os.path.join(os.getenv("PROGRAMDATA"),
                                     self.code_name, self.version_short)
        else:
            configdir = self.get_dir("sysconfdir")

        return [os.path.join(configdir, self.configfile)]

    def get_user_configfile(self):

        # Windows:         C:\Users\{user}\AppData\Roaming
        # Linux and co:    /home/{user}   (or similar)

        if sys.platform.startswith("win"):
            configdir = os.path.join(os.getenv("APPDATA"),
                                     self.code_name, self.version_short)
            return [os.path.join(configdir, self.configfile)]
        else:
            configdir = os.path.expanduser("~")
            return [os.path.join(configdir, "." + self.configfile)]

    def get_configfiles(self):

        u_cfg = self.get_user_configfile()
        g_cfg = self.get_global_configfile()

        return g_cfg + u_cfg

    def get_batchdir(self):

        return os.path.join(self.get_dir("pkgdatadir"),
                            'batch')

    def get_pkgdatadir_script(self, script):

        return os.path.join(self.get_dir("pkgdatadir"),
                            script)

    def get_alternate_version(self, version):
        """
        Return alternate version package object
        """

        if not version:
            return self

        pkg = None

        # Determine path (by absolute or local name)
        pythondir = os.path.normpath(version)
        prefix = os.path.normpath(self.get_dir("exec_prefix"))
        if not os.path.isabs(version):
            version = os.path.join(os.path.split(prefix)[0], version)

        prefix = os.path.normpath(version)
        config_file = os.path.join(prefix, 'lib', 'code_saturne_build.cfg')

        # load alternate package
        pkg = package(alternate_config=True,
                      name=self.name,
                      config_file=config_file,
                      install_prefix=prefix)

        return pkg

    def get_cross_compile(self):
        """
        Return cross-conpilation info
        """
        return self.config.features['build_os']

#-------------------------------------------------------------------------------
