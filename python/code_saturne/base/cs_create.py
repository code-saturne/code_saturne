#!/usr/bin/env python3
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
This module describes the script used to create a study/case for code_saturne.

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

import os, sys, shutil, stat
import types, string, re, fnmatch
from optparse import OptionParser
import configparser

from code_saturne.base import cs_exec_environment
from code_saturne.base import cs_run_conf
from code_saturne.base import cs_runcase

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

    parser.add_option("--noref", dest="use_ref",
                      action="store_false",
                      help="don't copy references")

    parser.add_option("--copy-ref", dest="use_ref",
                      action="store_true",
                      help="copy references and examples")

    parser.add_option("-q", "--quiet",
                      action="store_const", const=0, dest="verbose",
                      help="do not output any information")

    parser.add_option("-v", "--verbose",
                      action="store_const", const=2, dest="verbose",
                      help="dump study creation parameters")

    parser.add_option("--syrthes", dest="syr_case_names", type="string",
                      metavar="<syr_cases>", action="append",
                      help="create new SYRTHES case(s).")

    parser.add_option("--cathare", dest="cat_case_name", type="string",
                      metavar="<cat_case>",
                      help="create a new Cathare2 case.")

    parser.add_option("--python", dest="py_case_name", type="string",
                      metavar="<py_case>",
                      help="create a new Python script case.")

    parser.set_defaults(use_ref=False)
    parser.set_defaults(study_name=None)
    parser.set_defaults(case_names=[])
    parser.set_defaults(copy=None)
    parser.set_defaults(verbose=1)
    parser.set_defaults(n_sat=1)
    parser.set_defaults(syr_case_names=[])
    parser.set_defaults(cat_case_name=None)
    parser.set_defaults(py_case_name=None)

    (options, args) = parser.parse_args(argv)

    if options.case_names == []:
        if len(args) > 0:
            options.case_names = args
        else:
            options.case_names = ["CASE1"]

    return study(pkg,
                 options.study_name,
                 cases=options.case_names,
                 syr_case_names=options.syr_case_names,
                 cat_case_name=options.cat_case_name,
                 py_case_name=options.py_case_name,
                 copy=options.copy,
                 use_ref=options.use_ref,
                 verbose=options.verbose)

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
            unset_executable(os.path.join(path, p))

    return

#-------------------------------------------------------------------------------
# Create local launch script
#-------------------------------------------------------------------------------

def create_local_launcher(pkg, path=None):
    """
    Create a local launcher script.
    """

    # Use alternate wrapper if configured

    wrapper_postfix = None

    config = configparser.ConfigParser()
    config.read(pkg.get_configfiles())
    if config.has_option('install', 'wrapper_postfix'):
        wrapper_postfix = config.get('install', 'wrapper_postfix')

    # Generate script

    local_script = pkg.name + pkg.config.exeext
    if path:
        local_script = os.path.join(path, local_script)

    fd = open(local_script, 'w')
    cs_exec_environment.write_shell_shebang(fd)

    launcher_name = pkg.name

    if wrapper_postfix is None:
        cs_exec_environment.write_script_comment(fd,
                                                 'Ensure the correct command is found:\n')
        cs_exec_environment.write_prepend_path(fd, 'PATH',
                                               pkg.get_dir("bindir"))
        fd.write('\n')
    else:
        launcher_name += wrapper_postfix

    if sys.platform.startswith('win'):
        fd.write('\n')
        cs_exec_environment.write_script_comment(fd, 'Run command:\n')
        fd.write(launcher_name + ' gui ' +
                 cs_exec_environment.get_script_positional_args() + '\n')
    else:
        # On Linux or similar systems, add a backslash to prevent aliases
        # Start "code_saturne gui by default
        fd.write("""# Insert default command
cs_cmd=""
if test $# = 1; then
  if test -f $1; then
    cs_cmd=gui
  fi
elif test $# = 0; then
  cs_cmd=gui
fi

# Run command
"""
        )
        fd.write('\\')
        fd.write(launcher_name + ' $cs_cmd ' +
                 cs_exec_environment.get_script_positional_args() + '\n')

    fd.close()
    set_executable(local_script)

#-------------------------------------------------------------------------------
# Definition of a class for a study
#-------------------------------------------------------------------------------

class study:

    def __init__(self,
                 package,
                 study_name,
                 cases=[],
                 syr_case_names=[],
                 cat_case_name=None,
                 py_case_name=None,
                 copy=None,
                 use_ref=False,
                 verbose=0):
        """
        Initialize the structure for a study.
        """

        # Package specific information

        self.package = package

        self.study_name = study_name
        self.copy = copy
        if self.copy is not None:
            self.copy = os.path.abspath(self.copy)
        self.use_ref = use_ref
        self.verbose = verbose

        self.cases = []
        for c in cases:
            self.cases.append(c)
        self.n_sat = len(cases)

        self.syr_case_names = []
        for c in syr_case_names:
            self.syr_case_names.append(c)

        self.cat_case_name = cat_case_name

        self.py_case_name = py_case_name

    #---------------------------------------------------------------------------

    def create(self):
        """
        Create a study.
        """

        cur_dir = os.getcwd()

        if self.study_name:
            if self.verbose > 0:
                sys.stdout.write("  o Creating study '%s'...\n" % self.study_name)
            os.mkdir(self.study_name)
            os.chdir(self.study_name)
            os.mkdir('MESH')
            os.mkdir('POST')

        # Creating code_saturne cases
        repbase = os.getcwd()
        for c in self.cases:
            os.chdir(repbase)
            self.create_case(c)

        # Creating SYRTHES cases
        if len(self.syr_case_names) > 0:
            self.create_syrthes_cases(repbase)

        # Creating Cathare case
        if self.cat_case_name is not None:
            self.create_cathare_case(repbase)

        # Creating the Python script case
        if self.py_case_name is not None:
            self.create_python_case(repbase)

        # Creating coupling structures
        if len(self.cases) + len(self.syr_case_names) > 1 \
           or self.cat_case_name or self.py_case_name:
            self.create_coupling(repbase)
            create_local_launcher(self.package, repbase)

        if repbase != cur_dir:
            os.chdir(cur_dir)

    #---------------------------------------------------------------------------

    def create_syrthes_cases(self, repbase):
        """
        Create and initialize SYRTHES case directories.
        """

        os_env_save = os.environ

        cs_exec_environment.source_syrthes_env(self.package)
        import syrthes

        for s in self.syr_case_names:
            os.chdir(repbase)
            retval = syrthes.create_syrcase(s)
            if retval > 0:
                sys.stderr.write("Cannot create SYRTHES case: '%s'\n" % s)
                sys.exit(1)

        os_environ = os_env_save

    #---------------------------------------------------------------------------

    def create_cathare_case(self, repbase):
        """
        Create and initialize Cathare case directory for coupling.
        """

        os.chdir(repbase)
        self.create_case(self.cat_case_name)

    #---------------------------------------------------------------------------

    def create_python_case(self, repbase):
        """
        Create and initialize Python code case directory for coupling
        """

        os.chdir(repbase)

        if self.verbose > 0:
            sys.stdout.write("  o Creating Python code case  '%s'...\n" %
                             self.py_case_name)

        c = os.path.join(repbase, self.py_case_name)

        if not os.path.isdir(c):
            os.mkdir(c)

        os.chdir(c)

        for d in ["DATA", "SRC", "RESU"]:
            if not os.path.isdir(d):
                os.mkdir(d)

    #---------------------------------------------------------------------------

    def create_coupling(self, repbase):
        """
        Create structure to enable code coupling.
        """

        if self.verbose > 0:
            sys.stdout.write("  o Creating coupling features ...\n")

        e_pkg = re.compile('PACKAGE')
        e_dom = re.compile('DOMAIN')

        coupled_domains = []

        sep = ""

        for c in self.cases:

            c = os.path.normpath(c)
            base_c = os.path.basename(c)

            coupled_domains.append({'solver': 'Code_Saturne',
                                    'domain': base_c,
                                    'n_procs_weight': 'None'})

        for c in self.syr_case_names:

            c = os.path.normpath(c)
            base_c = os.path.basename(c)

            coupled_domains.append({'solver': 'SYRTHES',
                                    'domain': base_c,
                                    'param': 'syrthes_data.syd',
                                    'n_procs_weight': 'None',
                                    'opt': ''})
            # Last entry is for additional SYRTHES options
            # (ex.: postprocessing with '-v ens' or '-v med')

        if self.cat_case_name is not None:

            c = os.path.normpath(self.cat_case_name)
            base_c = os.path.basename(c)

            coupled_domains.append({'solver': 'CATHARE',
                                    'domain': base_c,
                                    'cathare_case_file':' jdd_case.dat',
                                    'neptune_cfd_domain': 'NEPTUNE'})


        if self.py_case_name is not None:

            c = os.path.normpath(self.py_case_name)
            base_c = os.path.basename(c)

            coupled_domains.append({'solver': 'PYTHON_CODE',
                                    'domain': self.py_case_name,
                                    'script': 'pycode.py',
                                    'command_line': '',
                                    'n_procs_weight': None,
                                    'n_procs_max': None})

        # Result directory for coupling execution

        resu = os.path.join(repbase, 'RESU_COUPLING')
        os.mkdir(resu)
        if not os.path.isdir(resu):
            os.mkdir(resu)

        self.__build_run_cfg__(distrep = repbase,
                               casename = 'coupling')

        self.__coupled_run_cfg__(distrep = repbase,
                                 coupled_domains = coupled_domains)

    #---------------------------------------------------------------------------

    def create_case(self, casename):
        """
        Create a case for a code_saturne study.
        """

        casedirname = casename

        datadir = self.package.get_dir("pkgdatadir")
        data_distpath  = os.path.join(datadir, 'data')

        if os.path.exists(casedirname):
            sys.stdout.write("  o Case  '%s' already exists\n" % casename)
            return

        if self.verbose > 0:
            sys.stdout.write("  o Creating case  '%s'...\n" % casename)

        os.mkdir(casedirname)
        os.chdir(casedirname)

        if self.copy is not None:
            if not (os.path.exists(os.path.join(self.copy, 'DATA', 'REFERENCE')) \
                    or os.path.exists(os.path.join(self.copy, 'SRC', 'REFERENCE'))):
                self.use_ref = False

        # Data directory

        data = 'DATA'

        os.mkdir(data)
        abs_setup_distpath = os.path.join(data_distpath, 'setup.xml')
        if os.path.isfile(abs_setup_distpath) and not self.copy:
            shutil.copy(abs_setup_distpath, data)
            unset_executable(data)

        if self.use_ref:
            thch_distpath = os.path.join(data_distpath, 'user')
            ref           = os.path.join(data, 'REFERENCE')
            shutil.copytree(thch_distpath, ref)
            unset_executable(ref)

        # Write a wrapper for code and launching

        create_local_launcher(self.package, data)

        # Generate run.cfg file or copy one

        run_conf = None
        run_conf_path = os.path.join(data, 'run.cfg')

        if not self.copy:
            run_conf = cs_run_conf.run_conf(run_conf_path,
                                            package=self.package,
                                            rebuild=True)

        # User source files directory

        src = 'SRC'

        if self.use_ref:
            user_distpath = os.path.join(datadir, 'user_sources')
            shutil.copytree(user_distpath, src)
            unset_executable(src)
        else:
            os.mkdir(src)

        # Copy data and source files from another case

        if self.copy is not None:

            # Data files

            ref_data = os.path.join(self.copy, data)
            data_files = os.listdir(ref_data)

            for f in data_files:
                abs_f = os.path.join(ref_data, f)
                if os.path.isfile(abs_f) and \
                       f not in ["SaturneGUI",
                                 "NeptuneGUI",
                                 self.package.name]:
                    shutil.copy(abs_f, data)
                    unset_executable(os.path.join(data, f))

            # Source files

            ref_src = os.path.join(self.copy, src)
            if os.path.exists(ref_src):
                src_files = os.listdir(ref_src)
            else:
                src_files = []

            for f in src_files:
                abs_f = os.path.join(ref_src, f)
                if os.path.isfile(abs_f):
                    shutil.copy(abs_f, src)
                    unset_executable(os.path.join(src, f))

            # If run.cfg was not present in initial case, generate it

            if not os.path.isfile(run_conf_path):

                run_conf = cs_run_conf.run_conf(run_conf_path,
                                                package=self.package,
                                                rebuild=True)

                # Runcase (for legacy structures)

                runcase_path = os.path.join(self.copy, 'SCRIPTS', 'runcase')

                if os.path.isfile(runcase_path):

                    i_c = cs_run_conf.get_install_config_info(self.package)
                    resource_name = cs_run_conf.get_resource_name(i_c)

                    runcase = cs_runcase.runcase(runcase_path,
                                                 package=self.package)
                    sections = runcase.run_conf_sections(resource_name=resource_name,
                                                         batch_template=i_c['batch'])
                    for sn in sections:
                        if not sn in run_conf.sections:
                            run_conf.sections[sn] = {}
                        for kw in sections[sn]:
                            run_conf.sections[sn][kw] = sections[sn][kw]

        # Now write run.cfg if not copied

        if run_conf != None:
            run_conf.save()

        # Results directory

        resu = 'RESU'
        if not os.path.isdir(resu):
            os.mkdir(resu)

    #---------------------------------------------------------------------------

    def __coupled_run_cfg__(self,
                          distrep,
                          coupled_domains=[]):
        """
        """
        if len(coupled_domains) == 0:
            return

        run_conf_path = os.path.join(distrep, 'run.cfg')
        run_conf = cs_run_conf.run_conf(run_conf_path,
                                        package=self.package,
                                        create_if_missing=True)

        domains = ""

        for i, dom in enumerate(coupled_domains):
            dom_name = dom["domain"]
            if i != 0:
                domains+=":"

            domains+=dom_name

            for key in dom.keys():
                run_conf.set(dom_name, key, dom[key])


        if domains != "":
            run_conf.set('setup', 'coupled_domains', domains)

        run_conf.save()

    #---------------------------------------------------------------------------

    def __build_run_cfg__(self,
                          distrep,
                          casename):
        """
        Retrieve batch file for the current system
        Update batch file for the study
        """

        run_conf_path = os.path.join(distrep, 'run.cfg')

        if self.copy is not None:
            ref_run_conf_path = os.path.join(self.copy, 'DATA', 'run.cfg')
            try:
                shutil.copy(ref_run_conf_path, run_conf_path)
            except Exception:
                pass

        # Add info from parent in case of copy

        run_conf = cs_run_conf.run_conf(run_conf_path,
                                        package=self.package,
                                        create_if_missing=True)

        run_conf.save()

    #---------------------------------------------------------------------------

    def dump(self):
        """
        Dump the structure of a study.
        """

        print()

        name = self.study_name
        if not name:
            name = os.path.basename(os.getcwd())

        print("Name  of the study:", name)
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
        if self.py_case_name != None:
            print("Python script instances:")
            for c in self.py_case_name:
                print("  " + c)
        print()

#-------------------------------------------------------------------------------
# Creation of the case of study directory tree
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
