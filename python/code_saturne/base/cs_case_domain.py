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

import configparser
import datetime
import fnmatch
import os
import os.path
import sys
import shutil
import stat

from code_saturne.base import cs_compile
from code_saturne.base import cs_xml_reader

from code_saturne.base.cs_exec_environment import run_command
from code_saturne.base.cs_exec_environment import enquote_arg, separate_args
from code_saturne.base.cs_exec_environment import get_ld_library_path_additions
from code_saturne.base.cs_exec_environment import source_syrthes_env

from code_saturne.base.cs_meg_to_c import meg_to_c_interpreter

#===============================================================================
# Utility functions
#===============================================================================

def any_to_str(arg):
    """Transform single values or lists to a whitespace-separated string"""

    s = ''

    if type(arg) == tuple or type(arg) == list:
        for e in arg:
            s += ' ' + str(e)
        return s[1:]

    else:
        return str(arg)

#-------------------------------------------------------------------------------

def make_clean_dir(path):
    """
    Create a directory, or remove files (not directories) in existing directory
    """

    if not os.path.isdir(path):
        os.mkdir(path)
    else:
        l = os.listdir(path)
        for f in l:
            os.remove(os.path.join(path, f))

#-------------------------------------------------------------------------------

class RunCaseError(Exception):
    """Base class for exception handling."""

    def __init__(self, *args):
        self.args = args

    def __str__(self):
        if len(self.args) == 1:
            return str(self.args[0])
        else:
            return str(self.args)

#   def __repr__(self):
#       return "%s(*%s)" % (self.__class__.__name__, repr(self.args))

#===============================================================================
# Classes
#===============================================================================

class base_domain:
    """
    Base class from which classes handling running case should inherit.
    """

    #---------------------------------------------------------------------------

    def __init__(self,
                 package,                    # main package
                 code_name = '<undefined>',  # code name
                 name = None,                # domain name
                 n_procs_weight = None,      # recommended number of processes
                 n_procs_min = 1,            # min. number of processes
                 n_procs_max = None):        # max. number of processes

        # Package specific information

        self.package = package

        self.code_name = code_name

        # User functions

        self.user = {}

        # Names, directories, and files in case structure

        self.case_dir = None
        self.case_root_dir = None

        self.name = name # used for multiple domains only

        self.data_dir = None
        self.result_dir = None
        self.src_dir = None

        # Notebook and parametric  definitions and additional user arguments

        self.notebook = None
        self.parametric_args = None
        self.kw_args = None

        # Working directory and executable

        self.exec_dir = None
        self.solver_path = None

        # Is the case already staged ?

        self.data_is_staged = False

        # Execution options

        self.n_procs = n_procs_weight
        if not n_procs_min:
            n_procs_min = 1
        self.n_procs_min = max(1, n_procs_min)
        self.n_procs_max = n_procs_max

        if self.n_procs is None:
            self.n_procs = 1
        self.n_procs = max(self.n_procs, self.n_procs_min)
        if self.n_procs_max != None:
            self.n_procs = min(self.n_procs, self.n_procs_max)

        # Error reporting
        self.error = ''

        self.error_long = ''

    #---------------------------------------------------------------------------

    def __input_path_abs_dir__(self, path):

        """
        Determine root directory based on given path.
        Paths including RESU should be based on dest_root_dir,
        those including MESH based on case_root_dir,
        others in case_dir.
        """

        prefix, base = os.path.split(path)
        while prefix != '':
            if base in ('RESU', 'RESU_COUPLING', 'MESH'):
                break
            prefix, base = os.path.split(prefix)

        if base in ('RESU', 'RESU_COUPLING'):
            if self.dest_root_dir:
                r_path = os.path.join(self.dest_root_dir, path)
                if os.path.exists(r_path):
                    return r_path
                c_path = os.path.join(self.case_root_dir, path)
                if os.path.exists(c_path):
                    return c_path
                # If path does not exist, assume it will be created
                # later, in chich case we use dest_root_dir
                return r_path
            else:
                return os.path.join(self.case_root_dir, path)
        elif base == 'MESH':
            return os.path.join(self.case_root_dir, path)

        return os.path.join(self.case_dir, path)

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir, staging_dir = None):

        # Names, directories, and files in case structure

        self.case_dir = case_dir
        self.case_root_dir = case_dir
        self.dest_root_dir = None

        if self.name != None:
            self.case_dir = os.path.join(self.case_dir, self.name)
            self.case_root_dir = os.path.join(case_dir, self.name)

        self.data_dir = os.path.join(self.case_dir, 'DATA')
        self.result_dir = os.path.join(self.case_dir, 'RESU')
        self.src_dir = os.path.join(self.case_dir, 'SRC')

        # If computation is already staged, avoid reading data in upstream
        # case directories, as it may have changed since staging.

        if staging_dir:
            self.case_dir = staging_dir
            if self.name != None:
                self.case_dir = os.path.join(self.case_dir, self.name)
            self.data_dir = self.case_dir
            self.src_dir = None

    #---------------------------------------------------------------------------

    def set_exec_dir(self, exec_dir):

        if os.path.isabs(exec_dir):
            self.exec_dir = exec_dir
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_dir)

        if self.name != None:
            self.exec_dir = os.path.join(self.exec_dir, self.name)

        if not os.path.isdir(self.exec_dir):
            os.makedirs(self.exec_dir)

    #---------------------------------------------------------------------------

    def set_result_dir(self, name, given_dir = None, dest_root_dir=None):
        """
        If suffix = true, add suffix to all names in result dir.
        Otherwise, create subdirectory
        """

        if given_dir is None:
            self.result_dir = os.path.join(self.case_dir, 'RESU', name)
        else:
            self.result_dir = given_dir

        if self.name != None:
            self.result_dir = os.path.join(self.result_dir, self.name)

        self.dest_root_dir = dest_root_dir

        if not os.path.isdir(self.result_dir):
            os.makedirs(self.result_dir)

    #---------------------------------------------------------------------------

    def copy_data(self):
        """
        Copy base data to the execution directory.
        """

        return

    #---------------------------------------------------------------------------

    def init_staged_data(self):
        """
        Initialize staged data in the execution directory.
        """

        self.data_is_staged = True

        return

    #---------------------------------------------------------------------------

    def prepare_data(self):
        """
        Prepare data in the execution directory prior to run
        """

    #---------------------------------------------------------------------------

    def copy_result(self, name, purge=False):
        """
        Copy a file or directory to the results directory,
        optionally removing it from the source.
        """

        # Determine absolute source and destination names

        if os.path.isabs(name):
            src = name
            dest = os.path.join(self.result_dir, os.path.basename(name))
        else:
            src = os.path.join(self.exec_dir, name)
            dest = os.path.join(self.result_dir, name)

        # If source and destination are identical, return
        if src == dest:
            return

        # Copy single file

        if os.path.isfile(src):
            shutil.copy2(src, dest)
            if purge:
                os.remove(src)

        # Copy single directory (possibly recursive)
        # Unlike os.path.copytree, the destination directory
        # may already exist.

        elif os.path.isdir(src):

            if not os.path.isdir(dest):
                os.mkdir(dest)
            l = os.listdir(src)
            for f in l:
                f_src = os.path.join(src, f)
                f_dest = os.path.join(dest, f)
                if os.path.isfile(f_src):
                    shutil.copy2(f_src, f_dest)
                elif os.path.isdir(f_src):
                    self.copy_result(f_src, f_dest)

            if purge:
                if os.path.islink(src):
                    os.remove(f)
                else:
                    shutil.rmtree(src)

    #---------------------------------------------------------------------------

    def purge_result(self, name):
        """
        Remove a file or directory from execution directory.
        """

        # Determine absolute name

        if os.path.isabs(name):
            f = name
        else:
            f = os.path.join(self.exec_dir, name)

        # Remove file or directory

        if os.path.isfile(f) or os.path.islink(f):
            os.remove(f)

        elif os.path.isdir(f):
            shutil.rmtree(f)

    #---------------------------------------------------------------------------

    def get_n_procs(self):
        """
        Returns an array (list) containing the current number of processes
        associated with a solver stage followed by the minimum and maximum
        number of processes.
        """

        return [self.n_procs, self.n_procs_min, self.n_procs_max]

    #---------------------------------------------------------------------------

    def set_n_procs(self, n_procs):
        """
        Assign a number of processes to a solver stage.
        """

        self.n_procs = n_procs

    #---------------------------------------------------------------------------

    def solver_command(self, need_abs_path=False):
        """
        Returns a tuple indicating the solver's working directory,
        executable path, and associated command-line arguments.
        """

        exec_path = self.solver_path
        if not os.path.isabs(exec_path) and need_abs_path:
            exec_path = os.path.join(self.exec_dir,
                                     os.path.basename(exec_path))

        return enquote_arg(self.exec_dir), enquote_arg(exec_path), ''

    #---------------------------------------------------------------------------

    def summary_info(self, s):
        """
        Output summary data into file s
        """

        if self.name:
            name = self.name
            exec_dir = os.path.join(self.exec_dir, name)
            result_dir = os.path.join(self.result_dir, name)
        else:
            name = os.path.basename(self.case_dir)
            exec_dir = self.exec_dir
            result_dir = self.result_dir

        s.write('  Case           : ' + name + '\n')
        s.write('    directory    : ' + self.case_dir + '\n')
        s.write('    results dir. : ' + self.result_dir + '\n')
        if exec_dir != result_dir:
            s.write('    exec. dir.   : ' + self.exec_dir + '\n')

#-------------------------------------------------------------------------------

class domain(base_domain):
    """Handle running case."""

    #---------------------------------------------------------------------------

    def __init__(self,
                 package,                     # main package
                 package_compute = None,      # package for compute environment
                 name = None,                 # domain name
                 n_procs_weight = None,       # recommended number of processes
                 n_procs_min = None,          # min. number of processes
                 n_procs_max = None,          # max. number of processes
                 logging_args = None,         # command-line options for logging
                 param = None,                # XML parameters file
                 prefix = None,               # installation prefix
                 adaptation = None):          # HOMARD adaptation script

        base_domain.__init__(self,
                             package,
                             'code_saturne',
                             name,
                             n_procs_weight,
                             n_procs_min,
                             n_procs_max)

        # Compute package if different from front-end

        if package_compute:
            self.package_compute = package_compute
        else:
            self.package_compute = self.package

        # Directories, and files in case structure

        self.restart_input = None
        self.restart_mesh_input = None
        self.mesh_input = None
        self.partition_input = None

        # Default executable

        self.solver_path = os.path.join(self.package_compute.get_dir("pkglibexecdir"),
                                        self.package.solver)

        # Preprocessor options

        self.mesh_dir = None
        self.meshes = None

        # Solver options

        self.preprocess_on_restart = False
        self.exec_solver = True

        if param:
            self.param = os.path.basename(param)
        else:
            self.param = None

        self.logging_args = logging_args
        self.solver_args = None

        # Additional data

        self.prefix = prefix

        self.compile_cflags = None
        self.compile_cxxflags = None
        self.compile_fcflags = None
        self.compile_nvccflags = None
        self.compile_libs = None

        # Adaptation using HOMARD

        self.adaptation = adaptation

        # MEG expression generator
        self.mci = None

    #---------------------------------------------------------------------------

    def __set_case_parameters__(self):

        update_xml = (self.data_is_staged == False)

        # We may now import user python script functions if present.

        self.user_locals = None

        user_scripts = os.path.join(self.exec_dir, 'cs_user_scripts.py')
        if os.path.isfile(user_scripts):
            try:
                exec(compile(open(user_scripts).read(), user_scripts, 'exec'),
                     locals(),
                     locals())
                self.user_locals = locals()
            except Exception:
                execfile(user_scripts, locals(), locals())
                self.user_locals = locals()

        # We may now parse the optional XML parameter file
        # now that its path may be built and checked.

        setup_path = os.path.join(self.exec_dir, "setup.xml")
        if os.path.isfile(setup_path):

            # Ensure XML file is up to date as a precaution,
            # and filter it in case of parametric or notebook arguments.

            case = self.__xml_case_initialize__(setup_path, apply_filters=True)
            case['xmlfile'] = setup_path

            P = cs_xml_reader.Parser(doc=case.doc)
            params = P.getParams()
            for k in list(params.keys()):
                self.__dict__[k] = params[k]

            case.xmlSaveDocument()

            if params['xml_root_name'] == 'NEPTUNE_CFD_GUI':
                solver_dir = self.package_compute.get_dir("pkglibexecdir")
                solver_name = "nc_solver" + self.package_compute.config.exeext
                self.solver_path = os.path.join(solver_dir, solver_name)
                self.code_name = 'neptune_cfd'

            self.param = "setup.xml"

        else:
            msg  = ('Remark:\n'
                    '  No setup.xml file was provided in the DATA folder.\n'
                    '  Default settings will be used.\n')
            print(msg, file = sys.stderr)

        # Now override or complete data from the XML file.

        if self.user_locals:
            m = 'define_domain_parameters'
            if m in self.user_locals.keys():
                eval(m + '(self)', globals(), self.user_locals)
                del self.user_locals[m]

        # Finally, ensure some fields are of the required types

        if type(self.meshes) != list:
            self.meshes = [self.meshes,]

    #---------------------------------------------------------------------------

    def __set_auto_restart__(self):
        """
        Select latest valid checkpoint directory for restart, based on name
        """

        self.restart_input = None

        from code_saturne.base.cs_exec_environment import get_command_output

        results_dir = os.path.abspath(os.path.join(self.result_dir, '..'))
        results = os.listdir(results_dir)
        results.sort(reverse=True)
        for r in results:
            m = os.path.join(results_dir, r, 'checkpoint', 'main.csc')
            if not os.path.isfile(m):
                m = os.path.join(results_dir, r, 'checkpoint', 'main')
            if os.path.isfile(m):
                try:
                    cmd = self.package.get_io_dump()
                    cmd += ' --location 0 ' + m
                    res = get_command_output(cmd)
                except Exception:
                    print('checkpoint of result: ' + r + ' does not seem usable')
                    continue
                self.restart_input = os.path.join(results_dir, r, 'checkpoint')
                break

        return

    #---------------------------------------------------------------------------

    def __xml_case_initialize__(self, path, apply_filters=False):
        """
        Build XML case object
        """

        from code_saturne.model.XMLengine import Case
        from code_saturne.model.SolutionDomainModel import getRunType
        case = Case(package=self.package, file_name=path)
        case['xmlfile'] = path
        case.xmlCleanAllBlank(case.xmlRootNode())
        preprocess_only = (getRunType(case) != 'standard')
        module_name = case.module_name()
        if module_name == 'code_saturne':
            from code_saturne.model.XMLinitialize import XMLinit
            XMLinit(case).initialize(preprocess_only)
        elif module_name == 'neptune_cfd':
            from code_saturne.model.XMLinitializeNeptune import XMLinitNeptune
            XMLinitNeptune(case).initialize(preprocess_only)

        if not apply_filters:
            return case

        # Apply changes defined through notebook or parametric options.

        if self.parametric_args:
            from code_saturne.base import cs_parametric_setup
            cs_parametric_setup.update_case_setup(case, self.parametric_args,
                                                  pkg=self.package)

        if self.notebook:
            from code_saturne.model.NotebookModel import NotebookModel
            notebookModel = NotebookModel(case)
            nbk_vars = notebookModel.getVarNameList()
            n_warnings = 0
            for k in self.notebook.keys():
                if k in nbk_vars:
                    vs = self.notebook[k]
                    v = None
                    try:
                        v = float(vs)
                    except Exception:
                        fmt = ("Warning: notebook variable '{0}'='{1}'"
                               " is not a real number.")
                        msg = fmt.format(k, vs)
                        print(msg, file = sys.stderr)
                        n_warnings += 1
                    if v != None:
                        notebookModel.setVariableValue(val=v, var=k)

                else:
                    fmt = ('Warning: {0} is not a known notebook variable.')
                    msg = fmt.format(k)
                    print(msg, file = sys.stderr)
                    n_warnings += 1
            if n_warnings > 0:
                print(file = sys.stderr)

        return case

    #---------------------------------------------------------------------------

    def for_domain_str(self):

        if self.name is None:
            return ''
        else:
            return 'for domain ' + str(self.name)

    #---------------------------------------------------------------------------

    def read_parameter_file(self, param):
        """
        Parse the optional XML parameter file.
        """

        if param is None:
            if os.path.isfile(os.path.join(self.data_dir, 'setup.xml')):
                param = 'setup.xml'

        if param != None:
            version_str = '2.0'
            P = cs_xml_reader.Parser(os.path.join(self.data_dir, param),
                                     version_str = version_str)
            params = P.getParams()
            for k in list(params.keys()):
                self.__dict__[k] = params[k]

            self.param = param

    #---------------------------------------------------------------------------

    def copy_data(self):
        """
        Copy base data to the execution directory
        """

        err_str = ""

        # Create the src folder if there are files to compile in source path

        src_files = []
        if self.src_dir != None:
            if os.path.exists(self.src_dir):
                src_files = cs_compile.files_to_compile(self.src_dir)
        if len(src_files) > 0:
            exec_src = os.path.join(self.exec_dir, 'src')

            make_clean_dir(exec_src)

            # Add header files to list so as not to forget to copy them
            dir_files = os.listdir(self.src_dir)
            src_files = src_files + (  fnmatch.filter(dir_files, '*.h')
                                     + fnmatch.filter(dir_files, '*.hxx')
                                     + fnmatch.filter(dir_files, '*.hpp'))

            # Copy source files to execution directory
            for f in src_files:
                src_file = os.path.join(self.src_dir, f)
                dest_file = os.path.join(exec_src, f)
                shutil.copy2(src_file, dest_file)

        # Copy data files

        dir_files = []
        if self.data_dir != self.exec_dir:
            dir_files = os.listdir(self.data_dir)

        if self.package.name in dir_files:
            dir_files.remove(self.package.name)

        for f in dir_files:
            src = os.path.join(self.data_dir, f)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(self.exec_dir, f))

        # Fixed parameter name

        setup_path = os.path.join(self.exec_dir, "setup.xml")

        if self.param != None:
            param_base = os.path.basename(self.param)
            if param_base != "setup.xml":
                src_path = os.path.join(self.exec_dir, param_base)
                if os.path.isfile(setup_path):
                    self.purge_result(setup_path) # in case of previous run here
                    fmt = ('Warning:\n'
                           '  Both {0} and {1} exist in\n'
                           '    {2}.\n'
                           '  {0} will be used for the computation.\n'
                           '  Be aware that to follow best practices '
                           'only one of the two should be present.\n\n')
                    msg = fmt.format(os.path.basename(self.param),
                                     os.path.basename(setup_path),
                                     self.data_dir)
                    print(msg, file = sys.stderr)
                try:
                    if os.path.islink(setup_path):
                        os.remove(setup_path)
                    os.symlink(self.param, setup_path)
                except Exception:
                    shutil.copy2(src_path, setup_path)

        # Now set parameters

        self.__set_case_parameters__()

    #---------------------------------------------------------------------------

    def init_staged_data(self):
        """
        Initialize staged data in the execution directory.
        """

        self.data_is_staged = True

        # Now set parameters

        self.__set_case_parameters__()

        # Ensure correct executable is used.

        exec_src = os.path.join(self.exec_dir, 'src')
        if os.path.isdir(exec_src):
            if len(cs_compile.files_to_compile(exec_src)) > 0:
                solver_name = os.path.basename(self.solver_path)
                self.solver_path = os.path.join('.', solver_name)

        # Additional precaution (in case source files were removed)

        solver = os.path.basename(self.solver_path)
        if os.path.isfile(os.path.join(self.exec_dir, solver)):
            self.solver_path = os.path.join('.', solver)

    #---------------------------------------------------------------------------

    def symlink(self, target, link=None, check_type=None):
        """
        Create a symbolic link to a file, or copy it if links are
        not possible
        """

        if target is None and link is None:
            return
        elif target is None:
            err_str = 'No target for link: ' + link
            raise RunCaseError(err_str)
        elif link is None:
            if self.exec_dir != None:
                link = os.path.join(self.exec_dir,
                                    os.path.basename(target))
            else:
                err_str = 'No path name given for link to: ' + target
                raise RunCaseError(err_str)

        if not os.path.exists(target):
            if check_type != 'allow_future':
                err_str = 'File: ' + target + ' does not exist.'
                raise RunCaseError(err_str)

        elif check_type == 'file':
            if not os.path.isfile(target):
                err_str = target + ' is not a regular file.'
                raise RunCaseError(err_str)

        elif check_type == 'dir':
            if not os.path.isdir(target):
                err_str = target + ' is not a directory.'
                raise RunCaseError(err_str)

        try:
            os.symlink(target, link)
        except AttributeError:
            if not os.path.isdir(target):
                shutil.copy2(target, link)
            else:
                if os.path.isdir(link):
                    shutil.rmtree(link)
                shutil.copytree(target, link)

    #---------------------------------------------------------------------------

    def needs_compile(self):
        """
        Compile and link user subroutines if necessary
        """
        # Check if there are files to compile in source path

        needs_comp = False

        exec_src = os.path.join(self.exec_dir, 'src')
        if os.path.isdir(exec_src):
            src_files = cs_compile.files_to_compile(exec_src)
            if self.exec_solver and len(src_files) > 0:
                needs_comp = True

        setup_path = os.path.join(self.exec_dir, "setup.xml")
        if os.path.isfile(setup_path):
            fp = os.path.join(self.exec_dir, self.param)
            case = self.__xml_case_initialize__(fp)

            # Do not call case.xmlSaveDocument() to avoid side effects in case
            # directory; is not required as meg_to_c_interpreter works from
            # case in memory

            module_name = case.module_name()
            self.mci = meg_to_c_interpreter(case,
                                            module_name=module_name,
                                            wdir = os.path.join(self.exec_dir, 'src'))

            if self.mci.has_meg_code():
                needs_comp = True
            else:
                self.mci = None

        return needs_comp

    #---------------------------------------------------------------------------

    def compile_and_link(self):
        """
        Compile and link user subroutines if necessary
        """

        src_files = []

        # Check if there are files to compile in source path
        # or if MEG functions need to be generated
        exec_src = os.path.join(self.exec_dir, 'src')
        if os.path.isdir(exec_src):
            src_files = cs_compile.files_to_compile(exec_src)

        if len(src_files) > 0 or self.mci != None:

            # Create the src folder if not done yet
            if not os.path.isdir(exec_src):
                os.mkdir(exec_src)

            if len(src_files) > 0:
                # Add header files to list so as not to forget to copy them
                dir_files = os.listdir(exec_src)
                src_files = src_files + (  fnmatch.filter(dir_files, '*.h')
                                         + fnmatch.filter(dir_files, '*.hxx')
                                         + fnmatch.filter(dir_files, '*.hpp'))

            if self.mci != None:
                mci_state = self.mci.save_all_functions()
                if mci_state['state'] == -1:
                    self.error = 'compile or link'
                    self.error_long = ' missing mathematical expressions:\n\n'
                    for i, eme in enumerate(mci_state['exps']):
                        self.error_long += " (%d/%d) %s is not provided for %s for zone %s\n" % (i+1, mci_state['nexps'], eme['func'], eme['var'], eme['zone'])

                    return
                elif mci_state['state'] == 2:
                    self.error = 'saving MEG generated sources'
                    self.error_long = ' Incorrect directory ?'

            log_name = os.path.join(self.exec_dir, 'compile.log')
            log = open(log_name, 'w')

            solver_name = os.path.basename(self.solver_path)

            retval = cs_compile.compile_and_link(self.package_compute,
                                                 solver_name,
                                                 exec_src,
                                                 self.exec_dir,
                                                 self.compile_cflags,
                                                 self.compile_cxxflags,
                                                 self.compile_fcflags,
                                                 self.compile_nvccflags,
                                                 self.compile_libs,
                                                 keep_going=True,
                                                 stdout=log,
                                                 stderr=log)

            log.close()

            if retval == 0:
                solver_dir = '.'
                self.solver_path = os.path.join(solver_dir, solver_name)
            else:
                # In case of error, copy source to results directory now,
                # as no calculation is possible, then raise exception
                for f in ['src', 'compile.log']:
                    self.copy_result(f)
                self.error = 'compile or link'

    #---------------------------------------------------------------------------

    def prepare_data(self):
        """
        Prepare data in the execution directory prior to run
        """

        err_str = ""

        if not self.exec_solver:
            return

        # Prepare some output directories
        # (Note: LUSTRE striping could be set here)

        try:
            os.mkdir(os.path.join(self.exec_dir, 'checkpoint'))
        except Exception:
            pass

        # Call user script if necessary

        if self.user_locals:
            m = 'domain_prepare_data_add'
            if m in self.user_locals.keys():
                eval(m + '(self)', globals(), self.user_locals)
                del self.user_locals[m]

        # Files or directories which migh be required but not available
        # yet (when the matching case may be staged but not run yet).

        upstream_pending = None

        # Restart files

        # Handle automatic case first

        ignore_checkpoint = False
        if self.solver_args:
            for a in ('--preprocess', '--quality', '-q'):
                if a in self.solver_args:
                    ignore_checkpoint = True

        if not ignore_checkpoint:
            if self.restart_input == '*':
                self.__set_auto_restart__()

            if self.restart_input != None:

                restart_input =  os.path.expanduser(self.restart_input)
                if not os.path.isabs(restart_input):
                    restart_input = self.__input_path_abs_dir__(restart_input)

                if not os.path.exists(restart_input):
                    err_str += restart_input + ' does not exist.\n\n'
                elif not os.path.isdir(restart_input):
                    err_str += restart_input + ' is not a directory.\n\n.'
                else:
                    self.symlink(restart_input,
                                 os.path.join(self.exec_dir, 'restart'))

                    print(' Restart from ' + self.restart_input)
                    print('')

                    r = os.path.join(restart_input, 'main.csc')
                    if not os.path.isfile(r):
                        if not os.path.isfile(os.path.join(self.restart_input, 'main')):
                            upstream_pending = [r]

            if self.restart_mesh_input != None and err_str == '':

                restart_mesh_input =  os.path.expanduser(self.restart_mesh_input)
                if not os.path.isabs(restart_mesh_input):
                    restart_mesh_input = self.__input_path_abs_dir__(restart_mesh_input)

                if not os.path.exists(restart_mesh_input):
                    if upstream_pending:
                        upstream_pending.append(restart_mesh_input)
                    else:
                        err_str += restart_mesh_input + ' does not exist.\n\n'

                elif not os.path.isfile(restart_mesh_input):
                    err_str += restart_mesh_input + ' is not a file.\n\n.'

                if err_str == '':
                    self.symlink(restart_mesh_input,
                                 os.path.join(self.exec_dir, 'restart_mesh_input'),
                                 check_type='allow_future')

                print(' Restart mesh ' + self.restart_mesh_input + '\n')

        # Mesh input file

        restart_input_mesh = None
        if self.restart_input:
            restart_input_mesh = os.path.join(self.exec_dir, 'restart', 'mesh_input.csm')
            if not os.path.exists(restart_input_mesh):
                restart_input_mesh = None

        if restart_input_mesh is None or self.preprocess_on_restart \
           or self.restart_mesh_input:
            if self.mesh_input:
                mesh_input = os.path.expanduser(self.mesh_input)
                if not os.path.isabs(mesh_input):
                    mesh_input = self.__input_path_abs_dir__(mesh_input)

                # Differentiate between a folder and file, since we now
                # have a file extension
                if os.path.isdir(mesh_input):
                    link_path = os.path.join(self.exec_dir, 'mesh_input')
                else:
                    link_path = os.path.join(self.exec_dir, 'mesh_input.csm')

                self.purge_result(link_path) # in case of previous run here
                self.symlink(mesh_input, link_path, check_type='allow_future')

                if not os.path.exists(mesh_input):
                    if not upstream_pending:
                        upstream_pending = []
                    upstream_pending.append(mesh_input)

        else:
            # use mesh from restart, with no symbolic link
            self.mesh_input = restart_input_mesh

        # Partition input files

        if self.partition_input != None and not err_str:

            partition_input = os.path.expanduser(self.partition_input)
            if not os.path.isabs(partition_input):
                partition_input = self.__input_path_abs_dir__(partition_input)

            if os.path.exists(partition_input):

                if not os.path.isdir(partition_input):
                    err_str += partition_input + ' is not a directory.\n\n'
                else:
                    self.symlink(partition_input,
                                 os.path.join(self.exec_dir, 'partition_input'))

        if upstream_pending:
            print(' Upstream computation might be staged but not run yet.')
            print(' Files required for computation stage but not yet present:')
            for r in upstream_pending:
                if os.path.islink(r):
                    print('   ' + (r))
                    print('   (-> ' + os.path.realpath(r) + ')')
                else:
                    print('   ' + os.path.normpath(r))
            print('')

        if len(err_str) > 0:
            self.error = 'data preparation'
            sys.stderr.write(err_str)
        else:
            try:
                link_path = os.path.join(self.exec_dir, "listing")
                os.symlink("run_solver.log", link_path)
            except Exception:
                pass

        self.data_is_staged = True

    #---------------------------------------------------------------------------

    def preprocess(self):
        """
        Runs the preprocessor in the execution directory
        """

        if self.mesh_input:
            return

        # Check if cartesian mesh is to be used
        if self.param != None:
            from code_saturne.model.XMLengine import Case
            from code_saturne.model.SolutionDomainModel import getMeshOriginType

            fp = os.path.join(self.exec_dir, self.param)
            case = Case(package=self.package, file_name=fp)
            case['xmlfile'] = fp
            case.xmlCleanAllBlank(case.xmlRootNode())

            if getMeshOriginType(case) == 'mesh_cartesian':
                return

        # If no mesh is provided return, since user can define mesh_input
        # using 'cs_user_mesh_input' user function.
        if len(self.meshes) == 1 and self.meshes[0] is None:
            return

        # Study directory
        study_dir = os.path.split(self.case_root_dir)[0]

        # User config file
        u_cfg = configparser.ConfigParser()
        u_cfg.read(self.package.get_user_configfile())

        # Global config file
        g_cfg = configparser.ConfigParser()
        g_cfg.read(self.package.get_global_configfile())

        # A mesh can be found in different mesh database directories
        # (case, study, user, global -- in this order)
        mesh_dirs = []
        if self.mesh_dir is not None:
            mesh_dir = os.path.expanduser(self.mesh_dir)
            if not os.path.isabs(mesh_dir):
                mesh_dir = os.path.join(self.case_root_dir, mesh_dir)
            mesh_dirs.append(mesh_dir)
        if os.path.isdir(os.path.join(study_dir, 'MESH')):
            mesh_dirs.append(os.path.join(study_dir, 'MESH'))
        if u_cfg.has_option('run', 'meshdir'):
            add_path = u_cfg.get('run', 'meshdir').split(':')
            for d in add_path:
                mesh_dirs.append(d)
        if g_cfg.has_option('run', 'meshdir'):
            add_path = g_cfg.get('run', 'meshdir').split(':')
            for d in add_path:
                mesh_dirs.append(d)

        # Switch to execution directory

        cur_dir = os.path.realpath(os.getcwd())
        if cur_dir != self.exec_dir:
            os.chdir(self.exec_dir)

        mesh_id = None

        if len(self.meshes) > 1:
            mesh_id = 0
            destdir = 'mesh_input'
            make_clean_dir(destdir)

        # Set environment

        ld_library_path_save = None
        add_lib_dirs = get_ld_library_path_additions(self.package)
        if add_lib_dirs:
            ld_library_path_save = os.getenv('LD_LIBRARY_PATH')
            ld_library_path = ""
            for d in add_lib_dirs:
                ld_library_path += d + ':'
            ld_library_path += ld_library_path_save
            os.environ['LD_LIBRARY_PATH'] = ld_library_path

        # Run once per mesh

        for m in self.meshes:

            # Get absolute mesh paths

            if m is None:
                err_str = 'Preprocessing stage required but no mesh is given'
                raise RunCaseError(err_str)

            if (type(m) == tuple):
                m0 = m[0]
            else:
                m0 = m

            m0 = os.path.expanduser(m0)

            mesh_path = m0
            if (not os.path.isabs(m0)) and len(mesh_dirs) > 0:
                for mesh_dir in mesh_dirs:
                    mesh_path = os.path.join(mesh_dir, m0)
                    if os.path.isfile(mesh_path):
                        break

            if not os.path.isfile(mesh_path):
                err_str = 'Mesh file ' + m0 + ' not found'
                if not (os.path.isabs(mesh_path) or mesh_dirs):
                    err_str += '(no mesh directory given)'
                raise RunCaseError(err_str)

            # Build command

            # Generate output mesh name
            if (mesh_id != None):
                mesh_id += 1
                _outputmesh = os.path.join('mesh_input',
                                           'mesh_%02d.csm' % (mesh_id))
            else:
                _outputmesh = 'mesh_input.csm'

            # Check if preprocessor is needed or not
            if mesh_path[-4:] == ".csm":
                # code_saturne mesh, no need to run preprocessor
                self.symlink(mesh_path,
                             os.path.join(self.exec_dir, _outputmesh))

                retcode = 0
            else:
                # run preprocessor if needed

                cmd = [self.package.get_preprocessor()]

                if (type(m) == tuple):
                    for opt in m[1:]:
                        cmd.append(opt)

                cmd = cmd + ['--out', _outputmesh]
                if (mesh_id != None):
                    cmd = cmd + ['--log', 'preprocessor_%02d.log' % (mesh_id)]
                    cmd = cmd + ['--case', 'preprocessor_%02d' % (mesh_id)]
                else:
                    cmd = cmd + ['--log']

                cmd.append(mesh_path)

                # Run command
                retcode = run_command(cmd, pkg=self.package)

            if retcode != 0:
                err_str = \
                    'Error running the preprocessor.\n' \
                    'Check the preprocessor.log file for details.\n\n'
                sys.stderr.write(err_str)

                self.exec_solver = False

                self.error = 'preprocess'

                break

        # Restore environment

        if ld_library_path_save:
            os.environ['LD_LIBRARY_PATH'] = ld_library_path_save

        # Revert to initial directory

        if cur_dir != self.exec_dir:
            os.chdir(cur_dir)

        return retcode

    #---------------------------------------------------------------------------

    def solver_command(self, need_abs_path=False):
        """
        Returns a tuple indicating the solver's working directory,
        executable path, and associated command-line arguments.
        """

        # Working directory and solver path from base class

        wd, exec_path, args = \
            base_domain.solver_command(self, need_abs_path)

        # Build kernel command-line arguments

        args = ''

        if self.logging_args != None:
            args += ' ' + self.logging_args

        if self.solver_args != None:
            args += ' ' + self.solver_args

        if self.package_compute.config.features['mpi'] == 'yes':
            if self.name != None:
                args += ' --mpi --app-name ' + enquote_arg(self.name)
            elif self.n_procs > 1:
                args += ' --mpi'

        return wd, exec_path, args

    #---------------------------------------------------------------------------

    def copy_results(self):
        """
        Retrieve solver results from the execution directory
        """

        # Call user script

        if self.user_locals:
            m = 'domain_copy_results_add'
            if m in self.user_locals.keys():
                eval(m + '(self)', globals(), self.user_locals)
                del self.user_locals[m]

        # Determine all files present in execution directory

        dir_files = os.listdir(self.exec_dir)

        # Only purge if we have no errors, to allow for future debugging.

        valid_dir = False
        purge = True

        if self.error != '':
            purge = False

        # Determine patterns from previous stages to ignore or possibly remove

        purge_list = []

        if not self.mesh_input and self.exec_solver:
            for f in ['mesh_input', 'mesh_input.csm']:
                if f in dir_files:
                    purge_list.append(f)

        for f in ['partition_input']:
            if f in dir_files:
                purge_list.append(f)

        # Determine files from this stage to ignore or to possibly remove

        solver_name = os.path.basename(self.solver_path)
        for f in [solver_name, self.package.runsolver]:
            if f in dir_files:
                purge_list.append(f)

        for f in purge_list:
            dir_files.remove(f)
            if purge:
                self.purge_result(f)

        if len(purge_list) > 0:
            valid_dir = True

        # Copy user sources, compile log, and xml file if present

        for f in ['src', 'compile.log', self.param]:
            if f in dir_files:
                valid_dir = True
                self.copy_result(f, purge)
                dir_files.remove(f)

        # Copy log files

        log_files = fnmatch.filter(dir_files, 'listing*')
        log_files.extend(fnmatch.filter(dir_files, '*.log'))
        log_files.extend(fnmatch.filter(dir_files, 'error*'))

        for f in log_files:
            self.copy_result(f, purge)
            dir_files.remove(f)

        if (len(log_files) > 0):
            valid_dir = True

        # Copy checkpoint files (in case of full disk, copying them
        # before other large data such as postprocessing output
        # increases chances of being able to continue).

        cpt = 'checkpoint'
        if cpt in dir_files:
            valid_dir = True
            self.copy_result(cpt, purge)
            dir_files.remove(cpt)

        # Now copy all other files

        if not valid_dir:
            return

        for f in dir_files:
            self.copy_result(f, purge)

    #---------------------------------------------------------------------------

    def summary_info(self, s):
        """
        Output summary data into file s
        """

        base_domain.summary_info(self, s)

        if not self.mesh_input:
            p = self.package.get_preprocessor()
            s.write('    preprocessor : ' + any_to_str(p) + '\n')

        if self.exec_solver:
            s.write('    solver       : ' + self.solver_path + '\n')

#-------------------------------------------------------------------------------

# SYRTHES coupling

class syrthes_domain(base_domain):

    def __init__(self,
                 package,
                 cmd_line = None,     # Command line to define optional syrthes4 behavior
                 name = None,
                 param = 'syrthes.data',
                 log_file = None,
                 n_procs_weight = None,
                 n_procs_min = 1,
                 n_procs_max = None,
                 n_procs_radiation = None,
                 verbose = True):

        base_domain.__init__(self,
                             package,
                             'SYRTHES',
                             name,
                             n_procs_weight,
                             n_procs_min,
                             n_procs_max)

        self.code_name = 'SYRTHES'

        self.n_procs_radiation = n_procs_radiation

        self.n_procs_ref = None
        self.n_procs_radiation = None

        # Additional parameters for code_saturne/SYRTHES coupling
        # Directories, and files in case structure

        self.cmd_line = cmd_line
        self.param = param

        self.logfile = log_file
        if self.logfile is None:
            self.logfile = 'syrthes.log'

        self.case_dir = None
        self.exec_dir = None
        self.data_dir = None
        self.src_dir = None
        self.result_dir = None
        self.echo_comm = None

        self.exec_solver = True

        # Always requires local buildn so relative path always known.

        self.solver_path = os.path.join('.', 'syrthes')

        # Generation of SYRTHES case deferred until we know how
        # many processors are really required

        self.syrthes_case = None

        # Determine environment which should be used for Syrthes
        # operations

        ld_library_path_save = os.getenv('LD_LIBRARY_PATH')

        source_syrthes_env(self.package, verbose)

        self.ld_library_path = os.getenv('LD_LIBRARY_PATH')

        self.__ld_library_path_restore__(ld_library_path_save)

    #---------------------------------------------------------------------------

    def __ld_library_path_restore__(self, ld_library_path_save):

        if ld_library_path_save:
            os.environ['LD_LIBRARY_PATH'] = ld_library_path_save
        else:
            try:
                os.environ.pop('LD_LIBRARY_PATH')
            except Exception:
                pass

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir, staging_dir = None):

        base_domain.set_case_dir(self, case_dir, staging_dir)

        # Names, directories, and files in case structure

        self.data_dir = self.case_dir
        self.src_dir = self.case_dir

    #---------------------------------------------------------------------------

    def set_exec_dir(self, exec_dir):

        if os.path.isabs(exec_dir):
            self.exec_dir = exec_dir
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_dir)

        self.exec_dir = os.path.join(self.exec_dir, self.name)

        if not os.path.isdir(self.exec_dir):
            os.mkdir(self.exec_dir)

    #---------------------------------------------------------------------------

    def solver_command(self, need_abs_path=False):
        """
        Returns a tuple indicating SYRTHES's working directory,
        executable path, and associated command-line arguments.
        """

        # Working directory and solver path from base class

        wd, exec_path, args = \
            base_domain.solver_command(self, need_abs_path)

        # Build Syrthes command-line arguments

        args = ''

        args += ' -d ' + enquote_arg(self.syrthes_case.data_file)
        args += ' -n ' + str(self.syrthes_case.n_procs)

        if self.syrthes_case.n_procs_ray > 0:
            args += ' -r ' + str(self.n_procs_ray)

        args += ' --name ' + enquote_arg(self.name)

        # Output to a logfile
        args += ' --log ' + enquote_arg(self.logfile)

        return wd, exec_path, args

    #---------------------------------------------------------------------------

    def __init_syrthes_case__(self):
        """
        Build or rebuild syrthes object
        """

        # Now try to read additional data set by Syrthes GUI

        part_tool_name = None

        f = open(os.path.join(self.case_dir, self.param), 'r')
        lines = f.readlines()
        del(f)
        for l in lines:
            if l[0] == '/':
                i = l.find('DOMAIN_POS=')
                if i > -1:
                    try:
                        id = int(l[i + 11:].strip())
                        if id == 0:
                            part_tool_name = 'scotch'
                        elif id == 1:
                            part_tool_name = 'metis'
                    except Exception:
                        pass
                i = l.find('LISTING=')
                if i > -1:
                    logfile = l[i + 8:].strip()
                    if logfile:
                        self.logfile = logfile
                i = l.find('NBPROC_COND=')
                if i > -1:
                     try:
                        self.n_procs_ref = int(l[i + 12:].strip())
                     except Exception:
                        pass
                i = l.find('NBPROC_RAD=')
                if i > -1:
                     try:
                        self.n_procs_radiation_ref = int(l[i + 11:].strip())
                     except Exception:
                        pass

        # Build command-line arguments

        args = '-d ' + os.path.join(self.case_root_dir, self.param)
        args += ' --name ' + self.name

        if self.n_procs != None and self.n_procs != 1:
            args += ' -n ' + str(self.n_procs)

        if self.n_procs_radiation != None:
            if self.n_procs_radiation > 0:
                args += ' -r ' + str(self.n_procs_radiation)

        if part_tool_name:
            args += ' -t ' + part_tool_name

        if self.data_dir != None:
            args += ' --data-dir ' + str(self.data_dir)

        if self.src_dir != None:
            args += ' --src-dir ' + str(self.src_dir)

        if self.exec_dir != None:
            args += ' --exec-dir ' + str(self.exec_dir)

        if self.cmd_line != None and len(self.cmd_line) > 0:
            args += ' ' + self.cmd_line

        # Define syrthes case structure

        import syrthes
        self.syrthes_case = syrthes.process_cmd_line(args.split())

        if self.syrthes_case.logfile is None:
            self.syrthes_case.set_logfile(self.logfile)
        else:
            self.logfile = self.syrthes_case.logfile

        # Read data file and store parameters

        self.syrthes_case.read_data_file()

    #---------------------------------------------------------------------------

    def prepare_data(self):
        """
        Fill SYRTHES domain structure
        Copy data to the execution directory
        Compile and link syrthes executable
        """

        # Build SYRTHES case; at this stage, the number of processes
        # may not be known, but this is not an issue here: for data
        # preparation, SYRTHES only uses the number of processes to
        # determine whether compilation should use MPI, and this is
        # always true anyways in the coupled case.

        ld_library_path_save = os.getenv('LD_LIBRARY_PATH')
        os.environ['LD_LIBRARY_PATH'] = self.ld_library_path

        self.__init_syrthes_case__()

        # Build exec_srcdir

        exec_srcdir = os.path.join(self.exec_dir, 'src')
        make_clean_dir(exec_srcdir)

        # Preparation of the execution directory and compile and link executable

        compile_logname = os.path.join(self.exec_dir, 'compile.log')

        retval = self.syrthes_case.prepare_run(exec_srcdir, compile_logname)

        self.copy_result(compile_logname)

        if retval != 0:
            err_str = '\n   Error during the SYRTHES preparation step\n'
            if retval == 1:
                err_str += '   Error during data copy\n'
            elif retval == 2:
                err_str += '   Error during syrthes compilation and link\n'
                # Copy source to results directory, as no run is possible
                for f in ['src', 'compile.log']:
                    self.copy_result(f)
            raise RunCaseError(err_str)

        sys.stdout.write('\n')

        # Restore environment

        self.__ld_library_path_restore__(ld_library_path_save)

    #---------------------------------------------------------------------------

    def preprocess(self):
        """
        Partition mesh for parallel run if required by user
        """

        ld_library_path_save = os.getenv('LD_LIBRARY_PATH')
        os.environ['LD_LIBRARY_PATH'] = self.ld_library_path

        # Rebuild Syrthes case now that number of processes is known
        self.__init_syrthes_case__()

        # Sumary of the parameters
        self.syrthes_case.dump()

        # Initialize output file if needed
        self.syrthes_case.logfile_init()

        # Pre-processing (including partitioning only if SYRTHES
        # computation is done in parallel)
        retval = self.syrthes_case.preprocessing()
        if retval != 0:
            err_str = '\n  Error during the SYRTHES preprocessing step\n'
            raise RunCaseError(err_str)

        self.__ld_library_path_restore__(ld_library_path_save)

    #---------------------------------------------------------------------------

    def copy_results(self):
        """
        Retrieve results from the execution directory
        """

        ld_library_path_save = os.getenv('LD_LIBRARY_PATH')
        os.environ['LD_LIBRARY_PATH'] = self.ld_library_path

        # Post-processing
        if self.syrthes_case.post_mode != None:
          retval = self.syrthes_case.postprocessing(mode = \
                   self.syrthes_case.post_mode)
        else:
          retval = 0

        if retval != 0:
            err_str = '\n   Error during SYRTHES postprocessing\n'
            raise RunCaseError(err_str)


        if self.exec_dir == self.result_dir:
            return

        retval = self.syrthes_case.save_results(save_dir = self.result_dir,
                                                horodat = False,
                                                overwrite = True)
        if retval != 0:
            err_str = '\n   Error saving SYRTHES results\n'
            raise RunCaseError(err_str)

        self.__ld_library_path_restore__(ld_library_path_save)

    #---------------------------------------------------------------------------

    def summary_info(self, s):
        """
        Output summary data into file s
        """

        base_domain.summary_info(self, s)

        s.write('    SYRTHES      : ' + self.solver_path + '\n')

#-------------------------------------------------------------------------------

# Cathare coupling

class cathare_domain(domain):
    """
    Class specific to coupling with CATHARE.
    """

    #---------------------------------------------------------------------------

    def __init__(self,
                 package,                     # main package
                 package_compute = None,      # package for compute environment
                 name = None,                 # domain name
                 n_procs_weight = None,       # recommended number of processes
                 n_procs_min = None,          # min. number of processes
                 n_procs_max = None,          # max. number of processes
                 logging_args = None,         # command-line options for logging
                 param = None,                # XML parameters file
                 prefix = None,               # installation prefix
                 cathare_case_file = None,    # Cathare case .dat file
                 neptune_cfd_dom   = None,    # neptune_cfd domain
                 adaptation = None):          # HOMARD adaptation script

        domain.__init__(self,
                        package,
                        package_compute,
                        name,
                        n_procs_weight,
                        n_procs_min,
                        n_procs_max,
                        logging_args,
                        param,
                        prefix,
                        adaptation)

        # If not provided, setup is automatically set to setup.xml
        if self.param is None:
            self.param = "setup.xml"

        self.cathare_case_file = cathare_case_file
        self.neptune_cfd_dom   = neptune_cfd_dom

    #---------------------------------------------------------------------------

    def read_parameter_file(self, param):

        # Getting the .xml from NEPTUNE_CFD, since all definitions are made
        # within its GUI and paramfile

        # First we check if code is launched using the run command,
        # hence in the main folder. If using submit, then  _spath
        # corresponds to the execution folder, so we remove the DATA
        _spath = os.path.split(self.case_dir)[0]
        nept_paramfile = os.path.join(_spath,
                                      self.neptune_cfd_dom,
                                      "DATA",
                                      param)

        if not os.path.exists(nept_paramfile):
            nept_paramfile = os.path.join(_spath,
                                          self.neptune_cfd_dom,
                                          param)
            if not os.path.exists(nept_paramfile):
                raise Exception("Could not find NEPTUNE_CFD data file.")

        ofile = open(nept_paramfile,'r').readlines()
        nfile = open(os.path.join(self.data_dir, param),'w')
        for line in ofile:
            if 'api_type' not in line:
                nfile.write(line)
            else:
                nfile.write('       <api_type>2</api_type>')

        nfile.close()

        if param != None:
            version_str = '2.0'
            P = cs_xml_reader.Parser(os.path.join(self.data_dir, param),
                                     version_str = version_str)
            params = P.getParams()
            for k in list(params.keys()):
                self.__dict__[k] = params[k]

            self.param = param

            solver_dir = self.package_compute.get_dir("pkglibexecdir")
            solver_name = "nc_solver" + self.package_compute.config.exeext
            self.solver_path = os.path.join(solver_dir, solver_name)

    #---------------------------------------------------------------------------

    def compile_cathare2_lib(self):
        """
        Copy the Cathare .dat file and the different .f files to the execution
        folder. Then generate the .so lib based on the ICoCo API of
        CATHARE2.
        """

        msg = "Generating CATHARE2 .so file\n" \
              "----------------------------\n\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

        import subprocess

        orig = os.getcwd()

        os.chdir(self.exec_dir)

        cathare_compilation_files = [self.cathare_case_file]

        dir_files = os.listdir(self.data_dir)
        for f in fnmatch.filter(dir_files, '*.f'):
            cathare_compilation_files.append(f)

        for f in cathare_compilation_files:
            src_file  = os.path.join(self.data_dir, f)
            dest_file = os.path.join(self.exec_dir, f)
            shutil.copy2(src_file, dest_file)

        # Generate the .so creation command (hence reducing dependencies)
        config = configparser.ConfigParser()
        config.read(self.package.get_configfiles())
        shell_cmd = 'export v25_3="%s"\n' % (config.get('install', 'cathare'))
        shell_cmd+= 'jdd_CATHARE="%s"\n' % (self.cathare_case_file)
        shell_cmd+= 'export LD_LIBRARY_PATH=${v25_3}/lib:$LD_LIBRARY_PATH\n'
        shell_cmd+= 'export LD_LIBRARY_PATH=${v25_3}/ICoCo/lib:$LD_LIBRARY_PATH\n'
        shell_cmd+= 'make -f ${v25_3}/ICoCo/Makefile_gad clean\n'
        shell_cmd+= '${v25_3}/unix-procedur/vers.unix\n'
        shell_cmd+= 'DATAFILE=${jdd_CATHARE} make -f ${v25_3}/ICoCo/Makefile_gad lib\n'

        # Shell
        user_shell = os.getenv('SHELL')
        if not user_shell:
            user_shell = '/bin/sh'

        # log
        log = open('cathare2_so_generation.log', 'w')

        p = subprocess.Popen(shell_cmd,
                             shell=True,
                             executable=user_shell,
                             stdout=log,
                             stderr=log,
                             universal_newlines=True)

        log.close()

        output, errors = p.communicate()

        if p.returncode != 0:
            self.error = 'compile cathare2 lib'
            self.error_long = 'Compilation of cathare2 .so library based on '
            self.error_long+= '%s file failed.' % (self.cathare_case_file)
            self.error_long+= ' Check "cathare2_so_generation.log' \


        os.chdir(orig)

    #---------------------------------------------------------------------------

    def compile_and_link(self):
        """
        Compile and link user subroutines if needed.
        """

        # Generate the CATHARE lib first
        self.compile_cathare2_lib()

        # Then compile NCFD source files if needed
        super(cathare_domain, self).compile_and_link()

    #---------------------------------------------------------------------------

    def solver_command(self, need_abs_path=False):
        """
        Returns a tuple indicating the script's working directory,
        executable path, and associated command-line arguments.
        """

        wd, exec_path, args \
            = super(cathare_domain, self).solver_command(need_abs_path)
        args += " --c2-wrapper"

        return wd, exec_path, args

#-------------------------------------------------------------------------------

# Coupling with a general Python-based code

class python_domain(base_domain):

    #---------------------------------------------------------------------------

    def __init__(self,
                 package,
                 cmd_line = None,
                 name = None,
                 script_name = 'partner_script.py',
                 log_file = None,
                 n_procs_weight = None,
                 n_procs_min = 1,
                 n_procs_max = None,
                 n_threads = 1):

        base_domain.__init__(self,
                             package,
                             'Python script',
                             name,
                             n_procs_weight,
                             n_procs_min,
                             n_procs_max)

        self.cmd_line = cmd_line
        self.logfile  = log_file
        if self.logfile is None:
            self.logfile = 'python.log'

        self.data_file = None

        self.solver_path = self.package.config.python
        self.exec_solver = True
        self.nthreads = n_threads

        self.script_name = script_name

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir, staging_dir = None):

        base_domain.set_case_dir(self, case_dir, staging_dir)

        self.data_dir = os.path.join(self.case_dir, "DATA")
        self.src_dir  = None
        self.result_dir = os.path.join(self.case_dir, "RESU")

    #---------------------------------------------------------------------------

    def set_exec_dir(self, exec_dir):

        if os.path.isabs(exec_dir):
            self.exec_dir = exec_dir
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_dir)

        self.exec_dir = os.path.join(self.exec_dir, self.name)

        if not os.path.isdir(self.exec_dir):
            os.mkdir(self.exec_dir)

    #---------------------------------------------------------------------------

    def copy_data(self):
        """
        Copy data to run directory
        """

        dir_files = os.listdir(self.data_dir)

        for f in dir_files:
            src = os.path.join(self.data_dir, f)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(self.exec_dir, f))

    #---------------------------------------------------------------------------

    def preprocess(self):
        """
        Preprocess dummy function: Does nothing for a standard python script
        """

        # Nothing to do

    #---------------------------------------------------------------------------

    def copy_results(self):
        """
        Copy results dummy function: Does nothing for a standard python script
        """
        # Nothing to do
        return

    #---------------------------------------------------------------------------

    def summary_info(self, s):
        """
        output summary data into file s
        """

        base_domain.summary_info(self, s)

    #---------------------------------------------------------------------------

    def solver_command(self, need_abs_path=False):
        """
        Returns a tuple indicating the script's working directory,
        executable path, and associated command-line arguments.
        """

        wd = enquote_arg(self.exec_dir)              # Working directory

        # Executable

        script_name = self.script_name
        if need_abs_path and not os.path.isabs(script_name):
            script_name = os.path.join(self.exec_dir,
                                       os.path.basename(script_name))

        exec_path = enquote_arg(self.solver_path) + ' ' + enquote_arg(script_name)

        # Build kernel command-line arguments

        args = ''

        if self.data_file:
            args += ' -d ' + enquote_arg(self.data_file)

        args += ' -n ' + str(self.n_procs)

        args += ' --name ' + enquote_arg(self.name)

        # Output to a logfile
        args += ' --log ' + enquote_arg(self.logfile)

        return wd, exec_path, args

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
