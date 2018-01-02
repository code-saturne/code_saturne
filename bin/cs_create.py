#!/usr/bin/env python
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
This module describes the script used to create a study/case for Code_Saturne.

This module defines the following functions:
- usage
- process_cmd_line
- set_executable
- unset_executable
and the following classes:
- study
- class
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from __future__ import print_function

import os, sys, shutil, stat
import types, string, re, fnmatch
from optparse import OptionParser
try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

import cs_exec_environment
import cs_runcase

#-------------------------------------------------------------------------------
# Process the passed command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog create [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-s", "--study", dest="study_name", type="string",
                      metavar="<study>",
                      help="create a new study")

    parser.add_option("-c", "--case", dest="case_names", type="string",
                      metavar="<case>", action="append",
                      help="create a new case")

    parser.add_option("--copy-from", dest="copy", type="string",
                      metavar="<case>",
                      help="create a case from another one")

    parser.add_option("--import-only", dest="import_only",
                      action="store_true",
                      help="rebuild scripts of existing case")

    parser.add_option("--noref", dest="use_ref",
                      action="store_false",
                      help="don't copy references")

    parser.add_option("-q", "--quiet",
                      action="store_const", const=0, dest="verbose",
                      help="do not output any information")

    parser.add_option("-v", "--verbose",
                      action="store_const", const=2, dest="verbose",
                      help="dump study creation parameters")

    parser.add_option("--syrthes", dest="syr_case_names", type="string",
                      metavar="<syr_cases>", action="append",
                      help="create new SYRTHES case(s).")

    parser.add_option("--aster", dest="ast_case_name", type="string",
                      metavar="<ast_case>",
                      help="create a new Code_Aster case.")

    parser.set_defaults(use_ref=True)
    parser.set_defaults(study_name=os.path.basename(os.getcwd()))
    parser.set_defaults(case_names=[])
    parser.set_defaults(copy=None)
    parser.set_defaults(verbose=1)
    parser.set_defaults(import_only=False)
    parser.set_defaults(n_sat=1)
    parser.set_defaults(syr_case_names=[])
    parser.set_defaults(ast_case_name=None)

    (options, args) = parser.parse_args(argv)

    if options.case_names == []:
        if len(args) > 0:
            options.case_names = args
        else:
            if not options.import_only:
                options.case_names = ["CASE1"]
            else:
                options.case_names = [""]

    return Study(pkg,
                 options.study_name,
                 options.case_names,
                 options.syr_case_names,
                 options.ast_case_name,
                 options.copy,
                 options.import_only,
                 options.use_ref,
                 options.verbose)

#-------------------------------------------------------------------------------
# Assign executable mode (chmod +x) to a file
#-------------------------------------------------------------------------------

def set_executable(filename):
    """
    Give executable permission to a given file.
    Equivalent to `chmod +x` shell function.
    """

    st   = os.stat(filename)
    mode = st[stat.ST_MODE] | stat.S_IXUSR
    if mode & stat.S_IRGRP:
        mode = mode | stat.S_IXGRP
    if mode & stat.S_IROTH:
        mode = mode | stat.S_IXOTH
    os.chmod(filename, mode)

    return

#-------------------------------------------------------------------------------
# Remove executable mode (chmod -x) from a file or files inside a directory
#-------------------------------------------------------------------------------

def unset_executable(path):
    """
    Remove executable permission from a given file or files inside a directory.
    Equivalent to `chmod -x` shell function.
    """

    if os.path.isfile(path):
        try:
            st   = os.stat(path)
            mode = st[stat.ST_MODE] | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            mode = mode ^ (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            os.chmod(path, mode)
        except Exception:
            pass

    elif os.path.isdir(path):
        l = os.listdir(path)
        for p in l:
            unset_executable(p)

    return

#-------------------------------------------------------------------------------
# Build lines necessary to import SYRTHES packages
#-------------------------------------------------------------------------------

def syrthes_path_line(pkg):
    """
    Build lines necessary to import SYRTHES packages
    """

    line = None

    config = configparser.ConfigParser()
    config.read(pkg.get_configfiles())

    if config.has_option('install', 'syrthes'):
        syr_datapath = os.path.join(config.get('install', 'syrthes'),
                                    os.path.join('share', 'syrthes'))
        line = 'sys.path.insert(1, \'' + syr_datapath + '\')\n'

    return line

#-------------------------------------------------------------------------------
# Definition of a class for a study
#-------------------------------------------------------------------------------

class Study:

    def __init__(self, package, name, cases, syr_case_names, ast_case_name,
                 copy, import_only, use_ref, verbose):
        """
        Initialize the structure for a study.
        """

        # Package specific information

        self.package = package

        self.name = name
        self.copy = copy
        if self.copy is not None:
            self.copy = os.path.abspath(self.copy)
        self.import_only = import_only
        self.use_ref = use_ref
        self.verbose = verbose

        self.cases = []
        for c in cases:
            self.cases.append(c)
        self.n_sat = len(cases)

        self.syr_case_names = []
        for c in syr_case_names:
            self.syr_case_names.append(c)

        self.ast_case_name = ast_case_name

        if self.import_only:
            self.use_ref = False
            self.copy = None


    def create(self):
        """
        Create a study.
        """

        if self.name != os.path.basename(os.getcwd()):

            if not self.import_only:
                if self.verbose > 0:
                    sys.stdout.write("  o Creating study '%s'...\n" % self.name)
                os.mkdir(self.name)
                os.chdir(self.name)
                os.mkdir('MESH')
                os.mkdir('POST')

            else:
                if self.verbose > 0:
                    sys.stdout.write("  o Importing study '%s'...\n" % self.name)
                os.chdir(self.name)

        # Creating Code_Saturne cases
        repbase = os.getcwd()
        for c in self.cases:
            os.chdir(repbase)
            self.create_case(c)

        # Creating SYRTHES cases
        if len(self.syr_case_names) > 0:
            self.create_syrthes_cases(repbase)

        # Creating Code_Aster case
        if self.ast_case_name is not None:
            config = configparser.ConfigParser()
            config.read(self.package.get_configfiles())
            if config.has_option('install', 'aster'):
                self.create_aster_case(repbase)

            else:
                sys.stderr.write("Cannot locate Code_Aster installation.")
                sys.exit(1)

        # Creating coupling structure
        if len(self.cases) + len(self.syr_case_names) > 1 or self.ast_case_name:
            self.create_coupling(repbase)


    def create_syrthes_cases(self, repbase):
        """
        Create and initialize SYRTHES case directories.
        """

        os_env_save = os.environ

        cs_exec_environment.source_syrthes_env(self.package)
        import syrthes

        for s in self.syr_case_names:
            os.chdir(repbase)
            if not self.import_only:
                retval = syrthes.create_syrcase(s)
            else:
                retval = 0
            if retval > 0:
                sys.stderr.write("Cannot create SYRTHES case: '%s'\n" % s)
                sys.exit(1)

        os_environ = os_env_save


    def create_aster_case(self, repbase):
        """
        Create and initialize Code_Aster case directory.
        """

        if self.verbose > 0:
            sys.stdout.write("  o Creating Code_Aster case  '%s'...\n" %
                             self.ast_case_name)

        c = os.path.join(repbase, self.ast_case_name)

        if not self.import_only:
            os.mkdir(c)

        datadir = self.package.get_dir("pkgdatadir")
        try:
            shutil.copy(os.path.join(datadir, 'salome', 'fsi.export'),
                        os.path.join(repbase, self.ast_case_name))
        except:
            sys.stderr.write("Cannot copy fsi.export file: " + \
                             os.path.join(datadir, 'salome', 'fsi.export') + ".\n")
            sys.exit(1)


    def create_coupling(self, repbase):
        """
        Create structure to enable code coupling.
        """

        if self.verbose > 0:
            sys.stdout.write("  o Creating coupling features ...\n")

        e_pkg = re.compile('PACKAGE')
        e_dom = re.compile('DOMAIN')

        solver_name = self.package.code_name

        header = \
"""# -*- coding: utf-8 -*-

#===============================================================================
# User variable settings to specify a coupling computation environnement.

# A coupling case is defined by a dictionnary, containing the following:

# Solver type ('Code_Saturne', 'SYRTHES', 'NEPTUNE_CFD' or 'Code_Aster')
# Domain directory name
# Run parameter setting file
# Number of processors (or None for automatic setting)
# Optional command line parameters. If not useful = None
#===============================================================================
"""

        dict_str = \
"""
# Define coupled domains

domains = [
"""
        sep = ""

        for c in self.cases:

            c = os.path.normpath(c)
            base_c = os.path.basename(c)

            template = sep + \
"""
    {'solver': 'PACKAGE',
     'domain': 'DOMAIN',
     'script': 'runcase',
     'n_procs_weight': None,
     'n_procs_min': 1,
     'n_procs_max': None}
"""
            template = re.sub(e_pkg, solver_name, template)
            template = re.sub(e_dom, base_c, template)

            dict_str += template

            if sep == "":
                sep = \
"""
    ,"""


        for c in self.syr_case_names:

            c = os.path.normpath(c)
            base_c = os.path.basename(c)

            template = \
"""
    ,
    {'solver': 'SYRTHES',
     'domain': 'DOMAIN',
     'script': 'syrthes_data.syd',
     'n_procs_weight': None,
     'n_procs_min': 1,
     'n_procs_max': None,
     'opt' : ''}               # Additional SYRTHES options
                               # (ex.: postprocessing with '-v ens' or '-v med')
"""
            template = re.sub(e_dom, base_c, template)
            dict_str += template

        if self.ast_case_name is not None:

            template = \
"""
    ,
    {'solver': 'Code_Aster',
     'domain': 'DOMAIN',
     'script': 'fsi.export'}
"""
            template = re.sub(e_dom, self.ast_case_name, template)
            dict_str += template

            template = \
"""
    ,
    {'coupler': 'FSI_coupler',
     'max_time_steps': 10,
     'n_sub_iterations': 1,
     'time_step': 0.0001,
     'start_time': 0.0,
     'epsilon': 0.00001}
"""
            dict_str += template

        # Now finish building dictionnary string

        dict_str += \
"""
    ]

"""
        # Result directory for coupling execution

        resu = os.path.join(repbase, 'RESU_COUPLING')
        if not self.import_only:
            os.mkdir(resu)

        coupling_base = 'coupling_parameters.py'
        coupling = os.path.join(repbase, coupling_base)

        fd = open(coupling, 'w')
        fd.write(header)

        if len(self.syr_case_names) > 0:
            syrthes_insert = syrthes_path_line(self.package)
            if syrthes_insert:
                fd.write("\n# Ensure the correct SYRTHES install is used.\n")
                fd.write(syrthes_insert)

        fd.write(dict_str)

        self.build_batch_file(distrep = repbase,
                              casename = 'coupling',
                              coupling = coupling_base)


    def create_case(self, casename):
        """
        Create a case for a Code_Saturne study.
        """

        casedirname = casename

        if self.verbose > 0:
            if not self.import_only:
                sys.stdout.write("  o Creating case  '%s'...\n" % casename)
            else:
                if not casename:
                    casedirname = "."
                    casename = os.path.basename(os.getcwd())
                sys.stdout.write("  o Importing case  '%s'...\n" % casename)

        datadir = self.package.get_dir("pkgdatadir")
        data_distpath  = os.path.join(datadir, 'data')

        if not self.import_only:
            try:
                os.mkdir(casedirname)
            except:
                sys.exit(1)

        os.chdir(casedirname)

        # Data directory

        data = 'DATA'

        if not self.import_only:
            os.mkdir(data)

        if self.use_ref:

            thch_distpath = os.path.join(data_distpath, 'thch')
            ref           = os.path.join(data, 'REFERENCE')
            os.mkdir(ref)
            for f in ['dp_C3P', 'dp_C3PSJ', 'dp_C4P', 'dp_ELE',
                      'dp_FCP.xml',
                      'dp_FUE', 'dp_transfo', 'meteo']:
                abs_f = os.path.join(thch_distpath, f)
                if os.path.isfile(abs_f):
                    shutil.copy(abs_f, ref)
                    unset_executable(ref)
            abs_f = os.path.join(datadir, 'cs_user_scripts.py')
            shutil.copy(abs_f, ref)
            unset_executable(ref)

        # Write a wrapper for GUI launching

        guiscript = os.path.join(data, self.package.guiname)

        fd = open(guiscript, 'w')
        cs_exec_environment.write_shell_shebang(fd)

        cs_exec_environment.write_script_comment(fd,
            'Ensure the correct command is found:\n')
        cs_exec_environment.write_prepend_path(fd, 'PATH',
                                               self.package.get_dir("bindir"))
        fd.write('\n')
        cs_exec_environment.write_script_comment(fd, 'Run command:\n')
        # On Linux systems, add a backslash to prevent aliases
        if sys.platform.startswith('linux'): fd.write('\\')
        fd.write(self.package.name + ' gui ' +
                 cs_exec_environment.get_script_positional_args() + '\n')

        fd.close()

        set_executable(guiscript)

        # User source files directory

        if not self.import_only:
            src = 'SRC'
            os.mkdir(src)

        if self.use_ref:

            user_distpath = os.path.join(datadir, 'user')
            user_examples_distpath = os.path.join(datadir, 'user_examples')

            user = os.path.join(src, 'REFERENCE')
            user_examples = os.path.join(src, 'EXAMPLES')
            shutil.copytree(user_distpath, user)
            shutil.copytree(user_examples_distpath, user_examples)

            add_datadirs = []
            if self.package.name == 'neptune_cfd' :
                add_datadirs.append(os.path.join(self.package.get_dir("datadir"),
                                                 self.package.name))

            for d in add_datadirs:
                user_distpath = os.path.join(d, 'user')
                user_examples_distpath = os.path.join(d, 'user_examples')

                if os.path.isdir(user_distpath):
                    s_files = os.listdir(user_distpath)
                    for f in s_files:
                        shutil.copy(os.path.join(user_distpath, f), user)

                if os.path.isdir(user_examples_distpath):
                    s_files = os.listdir(user_examples_distpath)
                    for f in s_files:
                        shutil.copy(os.path.join(user_examples_distpath, f), user_examples)

            unset_executable(user)
            unset_executable(user_examples)

        # Copy data and source files from another case

        if self.copy is not None:

            # Data files

            ref_data = os.path.join(self.copy, data)
            data_files = os.listdir(ref_data)

            for f in data_files:
                abs_f = os.path.join(ref_data, f)
                if os.path.isfile(abs_f) and \
                       f not in [self.package.guiname,
                                 'preprocessor_output']:
                    shutil.copy(os.path.join(ref_data, abs_f), data)
                    unset_executable(os.path.join(data, f))

            # Source files

            ref_src = os.path.join(self.copy, src)
            src_files = os.listdir(ref_src)

            c_files = fnmatch.filter(src_files, '*.c')
            cxx_files = fnmatch.filter(src_files, '*.cxx')
            cpp_files = fnmatch.filter(src_files, '*.cpp')
            h_files = fnmatch.filter(src_files, '*.h')
            hxx_files = fnmatch.filter(src_files, '*.hxx')
            hpp_files = fnmatch.filter(src_files, '*.hpp')
            f_files = fnmatch.filter(src_files, '*.[fF]90')

            for f in c_files + h_files + f_files + \
                    cxx_files + cpp_files + hxx_files + hpp_files:
                shutil.copy(os.path.join(ref_src, f), src)
                unset_executable(os.path.join(src, f))

        # Results directory (only one for all instances)

        resu = 'RESU'
        if not os.path.isdir(resu):
            os.mkdir(resu)

        # Script directory (only one for all instances)

        scripts = 'SCRIPTS'
        if not os.path.isdir(scripts):
            os.mkdir(scripts)

        self.build_batch_file(distrep = os.path.join(os.getcwd(), scripts),
                              casename = casename)


    def build_batch_file(self, distrep, casename, coupling=None):
        """
        Retrieve batch file for the current system
        Update batch file for the study
        """

        batch_file = os.path.join(distrep, 'runcase')
        if sys.platform.startswith('win'):
            batch_file = batch_file + '.bat'

        if self.copy is not None:
            ref_runcase_path = os.path.join(self.copy, 'SCRIPTS', 'runcase')
            if sys.platform.startswith('win'):
                ref_runcase_path += '.bat'
            try:
                shutil.copy(ref_runcase_path, batch_file)
            except Exception:
                pass

        # Add info from parent in case of copy

        runcase = cs_runcase.runcase(batch_file,
                                     package=self.package,
                                     rebuild=True,
                                     study_name=self.name,
                                     case_name=casename)

        if coupling:
            runcase.set_coupling(coupling)

        runcase.save()


    def dump(self):
        """
        Dump the structure of a study.
        """

        print()
        print("Name  of the study:", self.name)
        print("Names of the cases:", self.cases)
        if self.copy is not None:
            print("Copy from case:", self.copy)
        print("Copy references:", self.use_ref)
        if self.n_sat > 1:
            print("Number of instances:", self.n_sat)
        if self.syr_case_names != None:
            print("SYRTHES instances:")
            for c in self.syr_case_names:
                print("  " + c)
        if self.ast_case_name != None:
            print("Code_Aster instance:", self.ast_case_name)
        print()


#-------------------------------------------------------------------------------
# Creation of the study directory
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    welcome = """\
%(name)s %(vers)s study/case generation
"""

    study = process_cmd_line(argv, pkg)

    if study.verbose > 0:
        sys.stdout.write(welcome % {'name':pkg.name, 'vers': pkg.version})

    study.create()

    if study.verbose > 1:
        study.dump()


if __name__ == "__main__":
    main(sys.argv[1:], None)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
