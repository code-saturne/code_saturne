#!/usr/bin/env python
#
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009-2010  EDF
#
#   Code_Saturne is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License,
#   or (at your option) any later version.
#
#   Code_Saturne is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import ConfigParser
import datetime
import fnmatch
import os
import os.path
import sys
import shutil
import stat

import cs_config
import cs_compile

from cs_exec_environment import run_command

#===============================================================================
# Constants
#===============================================================================

solver_base_name = 'cs_solver'

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
                 n_procs = None,          # recommended number of processes
                 n_procs_min = 1,         # min. number of processes
                 n_procs_max = None):     # max. number of processes

        # Names, directories, and files in case structure

        self.case_dir = None

        self.tag = None # used for multiple domains only

        self.data_dir = None
        self.result_dir = None
        self.src_dir = None

        self.result_suffix = None

        self.mesh_dir = None

        # Working directory and executable

        self.exec_dir = None
        self.solver_path = None

        # Execution and debugging options

        self.n_procs = n_procs
        self.n_procs_min = max(1, n_procs_min)
        self.n_procs_max = n_procs_max

        if self.n_procs == None:
            self.n_procs = 1
        self.n_procs = max(self.n_procs, self.n_procs_min)
        if self.n_procs_max != None:
            self.n_procs = min(self.n_procs, self.n_procs_max)

        self.valgrind = None

        # Error reporting
        self.error = ''

    #---------------------------------------------------------------------------

    def set_tag(self, tag):

        # Subdirectory or suffix in case of multiple domains
        # (must be called before set_case_dir)

        self.tag = str(tag)

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir):

        # Names, directories, and files in case structure

        study_dir = os.path.split(case_dir)[0]

        self.case_dir = case_dir

        self.data_dir = os.path.join(self.case_dir, 'DATA')
        self.result_dir = os.path.join(self.case_dir, 'RESU')
        self.src_dir = os.path.join(self.case_dir, 'SRC')

        if (self.tag != None):
            self.data_dir += '.' + self.tag
            self.src_dir += '.' + self.tag

    #---------------------------------------------------------------------------

    def set_exec_dir(self, exec_dir):

        if os.path.isabs(exec_dir):
            self.exec_dir = exec_dir
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_dir)

        if self.tag != None:
            self.exec_dir = os.path.join(self.exec_dir, self.tag)

        if not os.path.isdir(self.exec_dir):
            os.makedirs(self.exec_dir)

    #---------------------------------------------------------------------------

    def set_result_dir(self, name, by_suffix=True):
        """
        If suffix = true, add suffix to all names in result dir.
        Otherwise, create subdirectory
        """

        if by_suffix == True:
            self.result_dir = os.path.join(self.case_dir, 'RESU')
            self.result_suffix = name

        else:
            self.result_dir = os.path.join(self.case_dir, 'RESU', name)
            if (self.tag != None):
                self.result_dir = os.path.join(self.result_dir,
                                               self.tag)

        if not os.path.isdir(self.result_dir):
            os.makedirs(self.result_dir)

    #---------------------------------------------------------------------------

    def copy_data_file(self, name, copy_name=None, description=None):
        """
        Copy a data file to the execution directory.
        """
        if os.path.isabs(name):
            source = name
            if copy_name == None:
                dest = os.path.join(self.exec_dir, os.path.basename(name))
            elif os.path.isabs(copy_name):
                dest = copy_name
            else:
                dest = os.path.join(self.exec_dir, copy_name)
        else:
            source = os.path.join(self.data_dir, name)
            if copy_name == None:
                dest = os.path.join(self.exec_dir, name)
            elif os.path.isabs(copy_name):
                dest = copy_name
            else:
                dest = os.path.join(self.exec_dir, copy_name)

        if os.path.isfile(source):
            shutil.copy2(source, dest)
        else:
            if description != None:
                err_str = \
                    'The ' + description + ' file: ', name, '\n' \
                    'can not be accessed.'
            else:
                err_str = \
                    'File: ', name, '\n' \
                    'can not be accessed.'
            raise RunCaseError(err_str)

    #---------------------------------------------------------------------------

    def copy_result(self, name, destname=None):
        """
        Copy a file or directory to the results directory.
        """

        # Determine absolute source and destination names

        if os.path.isabs(name):
            src = name
        else:
            src = os.path.join(self.exec_dir, name)

        if destname == None:
            dest = os.path.basename(name)
        else:
            dest = destname

        if not os.path.isabs(dest):
            dest = os.path.join(self.result_dir, dest)
            if (self.result_suffix != None):
                if (self.tag != None):
                    dest += '.' + self.tag
                dest += '.' + self.result_suffix

        # Copy single file

        if os.path.isfile(src):
            shutil.copy2(src, dest)

        # Copy single directory (possibly recursive)
        # Unkike os.path.copytree, the destination directory
        # may already exist.

        elif os.path.isdir(src):
            if not os.path.isdir(dest):
                os.mkdir(dest)
            list = os.listdir(src)
            for f in list:
                f_src = os.path.join(src, f)
                f_dest = os.path.join(dest, f)
                if os.path.isfile(f_src):
                    shutil.copy2(f_src, f_dest)
                elif os.path.isdir(f_src):
                    self.copy_result(f_src, f_dest)

    #---------------------------------------------------------------------------

    def copy_results_to_dir(self, names, dirname):
        """
        Copy result files to a given directory (a subdirectory
        of the results directory if a relative path is given).
        """

        if os.path.isabs(dirname):
            dirpath = dirname
        else:
            dirpath = os.path.join(self.result_dir, dirname)
            if (self.result_suffix != None):
                if (self.tag != None):
                    dirpath += '.' + self.tag
                dirpath += '.' + self.result_suffix

        if not os.path.isdir(dirpath):
            os.mkdir(dirpath)

        for name in names:
            if os.path.isabs(name):
                src = name
                dest = os.path.join(dirpath, os.path.basename(name))
            else:
                src = os.path.join(self.exec_dir, name)
                dest = os.path.join(dirpath, name)
            if os.path.isfile(src):
                shutil.copy2(src, dest)
            elif os.path.isdir(src):
                self.copy_result(src, dest)

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

    def solver_args(self, **kw):
        """
        Returns a tuple indicating the solver's working directory,
        executable path, and associated command-line arguments.
        """

        return self.exec_dir, self.solver_path, ''

#-------------------------------------------------------------------------------

class domain(base_domain):
    """Handle running case."""

    #---------------------------------------------------------------------------

    def __init__(self,
                 n_procs = None,              # recommended number of processes
                 n_procs_min = None,          # min. number of processes
                 n_procs_max = None,          # max. number of processes
                 n_procs_partition = None,    # n. processes for partitioner
                 meshes = None,               # name or names of mesh files
                 mesh_dir = None,             # mesh database directory
                 reorient = False,            # reorient badly-oriented meshes
                 partition_list = None,       # list of partitions
                 partition_opts = None,       # partitioner options
                 mode_args = None,            # --quality or --benchmark ?
                 logging_args = None,         # command-line options for logging
                 param = None,                # XML parameters file
                 thermochemistry_data = None, # file name
                 meteo_data = None,           # meteo. profileFile name
                 user_input_files = None,     # file name or names
                 user_output_files = None,    # file or directory name or names
                 lib_add = None,              # linker command-line options
                 adaptation = None):          # HOMARD adaptation script

        base_domain.__init__(self, n_procs, n_procs_min, n_procs_max)

        # Names, directories, and files in case structure

        self.result_suffix = None

        self.restart_input_dir = None
        self.preprocess_output_in = None
        self.partition_output_in = None

        # Default executable

        self.solver_path = os.path.join(cs_config.dirs.bindir,
                                        solver_base_name)

        # Preprocessor options

        if mesh_dir is not None:
            self.mesh_dir = os.path.expanduser(mesh_dir)

        if type(meshes) == tuple or type(meshes) == list:
            self.meshes = meshes
        else:
            self.meshes = (meshes,)
        self.reorient = reorient

        # Partition options

        self.partition_n_procs = n_procs_partition
        self.partition_list = partition_list
        self.partition_opts = partition_opts

        # Solver options

        self.mode_args = mode_args

        self.logging_args = logging_args

        self.param = param

        # Additional data

        self.thermochemistry_data = thermochemistry_data
        self.meteo_data = meteo_data

        self.user_input_files = user_input_files
        self.user_output_files = user_output_files

        self.lib_add = lib_add

        # Adaptation using HOMARD
        self.adaptation = adaptation

        # Steps to execute
        self.exec_preprocess = True
        self.exec_partition = True
        self.exec_solver = True

    #---------------------------------------------------------------------------

    def for_domain_str(self):

        if self.tag == None:
            return ''
        else:
            return 'for domain ' + str(self.tag)

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir):

        # Names, directories, and files in case structure

        base_domain.set_case_dir(self, case_dir)

        self.restart_input_dir = os.path.join(self.data_dir, 'restart')
        self.preprocess_output_in = os.path.join(self.data_dir,
                                                 'preprocessor_output')
        self.partition_output_in = os.path.join(self.data_dir,
                                                'PARTITION_OUTPUT')

    #---------------------------------------------------------------------------

    def symlink(self, target, link=None, check_type=None):
        """
        Create a symbolic link to a file, or copy it if links are
        not possible
        """

        if target == None and link == None:
            return
        elif target == None:
            err_str = 'No target for link: ' + link
            raise RunCaseError(err_str)
        elif link == None:
            if self.exec_dir != None:
                link = os.path.join(self.exec_dir,
                                    os.path.basename(target))
            else:
                err_str = 'No path name given for link to: ' + target
                raise RunCaseError(err_str)

        if not os.path.exists(target):
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
            shutil.copy2(target, link)

    #---------------------------------------------------------------------------

    def needs_compile(self):
        """
        Compile and link user subroutines if necessary
        """
        # Check if there are files to compile in source path

        dir_files = os.listdir(self.src_dir)

        src_files = (fnmatch.filter(dir_files, '*.c')
                     + fnmatch.filter(dir_files, '*.[fF]90'))

        if len(src_files) > 0:
            return True
        else:
            return False

    #---------------------------------------------------------------------------

    def compile_and_link(self):
        """
        Compile and link user subroutines if necessary
        """
        # Check if there are files to compile in source path

        dir_files = os.listdir(self.src_dir)

        src_files = (fnmatch.filter(dir_files, '*.c')
                     + fnmatch.filter(dir_files, '*.[fF]90'))

        if len(src_files) > 0:

            # Add header files to list so as not to forget to copy them

            src_files = src_files + fnmatch.filter(dir_files, '*.h')

            # Copy source files to result directory

            copy_dir = os.path.join(self.result_dir, 'SRC')
            if (self.result_suffix != None):
                if (self.tag != None):
                    copy_dir += '.' + self.tag
                copy_dir += '.' + self.result_suffix

            os.makedirs(copy_dir)

            for f in src_files:
                src_file = os.path.join(self.src_dir, f)
                dest_file = os.path.join(copy_dir, f)
                shutil.copy2(src_file, dest_file)
                try:
                    oldmode = (os.stat(dest_file)).stmode
                    newmode = oldmode & (stat.IRUSR | stat.IRGRP | stat.IROTH)
                    os.chmod(dest_file, newmode)
                except Exception:
                    pass

            # Copy source files to execution directory

            if (self.exec_dir != self.result_dir):
                exec_src = os.path.join(self.exec_dir, 'src_saturne')
                os.mkdir(exec_src)
                for f in src_files:
                    src_file = os.path.join(self.src_dir, f)
                    dest_file = os.path.join(exec_src, f)
                    shutil.copy2(src_file, dest_file)

            log_name = os.path.join(self.exec_dir, 'compile.log')
            log = open(log_name, 'w')

            # Note: src_dir and copy_dir should contain identical source files,
            #       but we prefer to use the copied version, so that if the
            #       source is later modified, possible debug information in an
            #       executable file will reference the correct (saved) version.

            retval = cs_compile.compile_and_link(copy_dir,
                                                 self.exec_dir,
                                                 self.lib_add,
                                                 keep_going=True,
                                                 stdout=log,
                                                 stderr=log)

            log.close()
            self.copy_result(log_name)

            if retval == 0:
                self.solver_path = os.path.join(self.exec_dir,
                                                solver_base_name)
            else:
                raise RunCaseError('Compile or link error.')

    #---------------------------------------------------------------------------

    def check_model_consistency(self):
        """
        Check model user subroutine and xml options consistency
        """
        from cs_check_consistency import check_consistency
        return check_consistency(self.param, self.src_dir, self.n_procs)


    #---------------------------------------------------------------------------

    def copy_preprocessor_data(self):
        """
        Copy preprocessor data to execution directory
        """

        if self.exec_preprocess == False:
            return

        # Study directory
        study_dir = os.path.split(self.case_dir)[0]

        # User config file
        u_cfg = ConfigParser.ConfigParser()
        u_cfg.read(os.path.expanduser('~/.code_saturne.cfg'))

        # Global config file
        g_cfg = ConfigParser.ConfigParser()
        g_cfg.read(os.path.join(cs_config.dirs.sysconfdir, 'code_saturne.cfg'))

        # A mesh can be found in different mesh database directories
        # (case, study, user, global -- in this order)
        mesh_dirs = []
        if self.mesh_dir is not None:
            mesh_dirs.append(self.mesh_dir)
        if os.path.isdir(os.path.join(study_dir, 'MESH')):
            mesh_dirs.append(os.path.join(study_dir, 'MESH'))
        if u_cfg.has_option('run', 'meshdir'):
            mesh_dirs.append(u_cfg.get('run', 'meshdir'))
        if g_cfg.has_option('run', 'meshdir'):
            mesh_dirs.append(g_cfg.get('run', 'meshdir'))

        # Symlink the different meshes
        for mesh in self.meshes:

            if mesh is None:
                err_str = 'Preprocessing stage is asked but no mesh is given'
                raise RunCaseError(err_str)

            mesh = os.path.expanduser(mesh)

            if os.path.isabs(mesh):
                mesh_path = mesh
            elif len(mesh_dirs) > 0:
                for mesh_dir in mesh_dirs:
                    mesh_path = os.path.join(mesh_dir, mesh)
                    if os.path.isfile(mesh_path):
                        break
            else:
                err_str = 'No mesh directory given'
                raise RunCaseError(err_str)

            if not os.path.isfile(mesh_path):
                err_str = 'Mesh file ' + mesh + ' not found'
                raise RunCaseError(err_str)

            base_name = os.path.basename(mesh_path)

            link_path = os.path.join(self.exec_dir, base_name)
            self.symlink(mesh_path, link_path)

            # Special case for meshes in EnSight format: link to .geo file
            # necessary (retrieve name through .case file)
            base, ext = os.path.splitext(base_name)
            if ext == '.case':
                try:
                    f = open(mesh_path)
                    text = f.read(4096) # Should be largely sufficient
                    f.close()
                    m = re.search('^model:.*$', text, re.MULTILINE)
                    geo_name = (m.group()).split()[1]
                    mesh_path = os.path.join(self.mesh_dir, geo_name)
                    link_path = os.path.join(self.exec_dir, geo_name)
                    self.symlink(mesh_path, link_path)
                except Exception:
                    err_str = 'Model file name not found in ' + mesh_path
                    raise RunCaseError(err_str)

    #---------------------------------------------------------------------------

    def copy_preprocessor_output_data(self):
        """
        Copy preprocessor_output file to the execution directory,
        required both for the partitioner and the solver.
        """

        if self.exec_preprocess:
            return
        elif not (self.exec_partition or self.exec_solver):
            return

        if self.preprocess_output_in != None:
            self.symlink(self.preprocess_output_in,
                         os.path.join(self.exec_dir, 'preprocessor_output'),
                         'file')
        else:
            err_str = 'Error: no path name given for link to: ' + target
            raise RunCaseError(err_str)

    #---------------------------------------------------------------------------

    def copy_solver_data(self):
        """
        Copy solver data to the execution directory
        """

        if self.exec_solver == False:
            return

        if self.n_procs < 2:
            self.exec_partition = False
        elif self.exec_partition == False and self.partition_output_in != None:
            partition = os.path.join(self.partition_output_in,
                                     'domain_number_' + str(self.n_procs))
            if os.path.isfile(partition):
                self.symlink(partition)
            else:
                w_str = \
                    'Warning: no partitioning file is available\n' \
                    '         (no ' + partition + ').\n' \
                    '\n' \
                    '         Geometry-based partitioning will be used.\n'
                sys.stderr.write(w_str)

        # Parameters file

        if self.param != None:
            self.copy_data_file(self.param,
                                os.path.basename(self.param),
                                'parameters')

        # Restart files

        if self.restart_input_dir != None:

            if os.path.exists(self.restart_input_dir):

                if not os.path.isdir(self.restart_input_dir):
                    err_str = self.restart_input_dir + ' is not a directory.'
                    raise RunCaseError(err_str)
                else:
                    self.symlink(self.restart_input_dir,
                                 os.path.join(self.exec_dir, 'restart'))

        # Data for specific physics

        if self.thermochemistry_data != None:
            self.copy_data_file(self.thermochemistry_data,
                                'dp_thch',
                                'thermochemistry')
            if not os.path.isfile('JANAF'):
                self.copy_data_file(os.path.join(cs_config.dirs.pkgdatadir,
                                                 'data',
                                                 'thch',
                                                 'JANAF'),
                                    'JANAF')

        if self.meteo_data != None:
            self.copy_data_file(self.meteo_data,
                                'meteo',
                                'meteo profile')
            # Second copy so as to have correct name upon backup
            if self.meteo_data != 'meteo':
                self.copy_data_file(self.meteo_data)

        # Presence of user input files

        if self.user_input_files != None:
            for f in self.user_input_files:
                self.copy_data_file(f)

    #---------------------------------------------------------------------------

    def run_preprocessor(self):
        """
        Runs the preprocessor in the execution directory
        """

        if self.exec_preprocess == False:
            return

        # Build command

        cmd = os.path.join(cs_config.dirs.bindir, 'cs_preprocess')
        for m in self.meshes:
            cmd += ' --mesh ' + os.path.basename(m)

        if self.reorient:
            cmd += ' --reorient'

        cmd += ' --log'

        # Run command

        cur_dir = os.path.realpath(os.getcwd())
        if cur_dir != self.exec_dir:
            os.chdir(self.exec_dir)

        retcode = run_command(cmd)

        if retcode != 0:
            err_str = \
                'Error running the preprocessor.\n' \
                'Check the preprocessor.log file for details.\n\n'
            sys.stderr.write(err_str)

            self.exec_partition = False
            self.exec_solver = False

            self.error = 'preprocess'

        if cur_dir != self.exec_dir:
            os.chdir(cur_dir)

        return retcode

    #---------------------------------------------------------------------------

    def check_partitioner(self):
        """
        Tests if the partitioner is available and partitioning is defined.
        """

        partitioner = os.path.join(cs_config.dirs.bindir, 'cs_partition')
        if not os.path.isfile(partitioner):
            if self.n_procs > 1:
                w_str = \
                    'Warning: ' + partitioner + ' not found.\n\n' \
                    'The partitioner may not have been installed' \
                    '  (this is the case if neither METIS nor SCOTCH ' \
                    ' are available).\n\n' \
                    'Partitioning by a space-filling curve will be used.\n\n'
                sys.stderr.write(w_str)
            self.exec_partition = False
            self.partition_n_procs = None

        if self.partition_list == None and not self.exec_solver:
            err_str = \
                'Unable to run the partitioner:\n' \
                'The list of required partitionings is not set.\n' \
                'It should contain the number of processors for which a\n' \
                'partition is required, or a list of such numbers.\n'
            raise RunCaseError(err_str)

    #---------------------------------------------------------------------------

    def run_partitioner(self):
        """
        Runs the partitioner in the execution directory
        """

        self.check_partitioner()
        if self.exec_partition == False:
            return

        # Build command

        cmd = os.path.join(cs_config.dirs.bindir, 'cs_partition')

        if self.partition_opts != None:
            cmd += ' ' + self.partition_opts

        if self.partition_list != None:
            cmd += ' ' + any_to_str(self.partition_list)

        if self.exec_solver and self.n_procs != None:
            np = self.n_procs
            if self.partition_list == None:
                cmd += ' ' + str(np)
            elif np > 1 and not np in self.partition_list:
                cmd += ' ' + str(np)

        # Run command

        cur_dir = os.path.realpath(os.getcwd())
        if cur_dir != self.exec_dir:
            os.chdir(self.exec_dir)

        retcode = run_command(cmd)

        if retcode != 0:
            err_str = \
                'Error running the partitioner.\n' \
                'Check the partition.log file for details.\n\n'
            sys.stderr.write(err_str)

            self.exec_solver = False

            self.error = 'partition'

        if cur_dir != self.exec_dir:
            os.chdir(cur_dir)

        return retcode

    #---------------------------------------------------------------------------

    def partitioner_args(self):
        """
        Returns a tuple indicating the partitioner's working directory,
        executable path, and associated command-line arguments.
        """
        # Working directory and executable path

        wd = self.exec_dir
        exec_path = os.path.join(cs_config.dirs.bindir, 'cs_partition')

        # Build kernel command-line arguments

        args = ''

        if self.partition_n_procs > 1:
            args += ' --mpi'

        if self.partition_opts != None:
            args += ' ' + self.partition_opts

        if self.partition_list != None:
            args += ' ' + any_to_str(self.partition_list)

        if self.exec_solver and self.n_procs != None:
            np = self.n_procs
            if self.partition_list == None:
                args += ' ' + str(np)
            elif np > 1 and not np in self.partition_list:
                args += ' ' + str(np)

        # Adjust for Valgrind if used

        if self.valgrind != None:
            args = self.solver_path + ' ' + args
            exec_path = self.valgrind

        return wd, exec_path, args

    #---------------------------------------------------------------------------

    def solver_args(self, **kw):
        """
        Returns a tuple indicating the solver's working directory,
        executable path, and associated command-line arguments.
        """

        wd = self.exec_dir              # Working directory
        exec_path = self.solver_path    # Executable

        # Build kernel command-line arguments

        args = ''

        if self.param != None:
            args += ' --param ' + self.param

        if self.logging_args != None:
            args += ' ' + self.logging_args

        if self.mode_args != None:
            args += ' ' + self.mode_args

        if 'app_id' in kw:
            if self.tag == None:
                args += ' --mpi --app-name ' + os.path.basename(self.case_dir)
            else:
                args += ' --mpi --app-name ' + os.path.basename(self.case_dir) \
                    + '.' + self.tag
        elif self.n_procs > 1:
            args += ' --mpi'

        if 'syr_port' in kw:
            args += ' --syr-socket ' + str(kw['syr_port'])

        # Adjust for Valgrind if used

        if self.valgrind != None:
            args = self.solver_path + ' ' + args
            exec_path = self.valgrind + ' '

        return wd, exec_path, args

    #---------------------------------------------------------------------------

    def copy_preprocessor_results(self):
        """
        Retrieve preprocessor results from the execution directory
        """

        if self.exec_dir == self.result_dir and self.result_suffix == None:
            return

        dir_files = os.listdir(self.exec_dir)

        # Copy log file first

        f = os.path.join(self.exec_dir, 'preprocessor.log')
        if os.path.isfile(f):
            self.copy_result(f)

        # copy output if required

        if not self.exec_solver:
            f = os.path.join(self.exec_dir, 'preprocessor_output')
            if os.path.isfile(f):
                self.copy_result(f)

    #---------------------------------------------------------------------------

    def copy_partition_results(self):
        """
        Retrieve partition results from the execution directory
        """

        if self.exec_dir == self.result_dir and self.result_suffix == None:
            return

        dir_files = os.listdir(self.exec_dir)

        # Copy log file first

        f = os.path.join(self.exec_dir, 'partition.log')
        if os.path.isfile(f):
            self.copy_result(f)

        # copy output if required

        if not self.exec_solver:

            part_files = fnmatch.filter(dir_files, 'domain_number_*')
            self.copy_results_to_dir(part_files, 'PARTITION_OUTPUT')

    #---------------------------------------------------------------------------

    def copy_solver_results(self, date=None):
        """
        Retrieve solver results from the execution directory
        """

        if self.exec_dir == self.result_dir and self.result_suffix == None:
            return

        dir_files = os.listdir(self.exec_dir)

        # Copy log files first

        log_files = fnmatch.filter(dir_files, 'listing*')
        log_files.extend(fnmatch.filter(dir_files, 'error*'))

        for f in log_files:
            self.copy_result(f)

        # Copy checkpoint files second (in case of full disk,
        # increases chances of being able to continue).

        self.copy_result('checkpoint')

        # User files

        if self.user_output_files != None:
            user_files = []
            for f in self.user_output_files:
                if f in dir_files:
                    user_files.append(f)
            if len(user_files) > 0:
                self.copy_results_to_dir(user_files, 'RES_USER')

        # Parameter or similar data files

        if self.param != None:
            self.copy_result(self.param)

        if self.thermochemistry_data != None:
            self.copy_result(self.thermochemistry_data)

        if self.meteo_data != None:
            self.copy_result(self.meteo_data)

        # Plot history (including user history) files

        self.copy_result('monitoring')
        history_files = fnmatch.filter(dir_files, 'ush*')
        if len(history_files) > 0:
            self.copy_results_to_dir(history_files, 'monitoring')

        post_list = fnmatch.filter(dir_files, '*.ensight')
        post_list.extend(fnmatch.filter(dir_files, '*.med'))
        post_list.extend(fnmatch.filter(dir_files, '*.cgns'))

        # Retrieve postprocessor output.

        # We handle both files and directories here; *.ensight
        # files are always directories, and *.med or *.cgns files
        # are always regular files, but the user could define
        # an output directory with a .med extension, especially if
        # there are many separate *.med files.

        # Directories in the form dir.extension (for example chr.ensight)
        # are copied to DIR.SUFFIX, or DIR.extension (for example
        # CHR.SUFFIX or CHR.ensight) depending on the suffix behavior

        for p in post_list:

            p_abs = os.path.join(self.exec_dir, p)

            if os.path.isfile(p_abs):
                self.copy_result(p_abs)

            elif os.path.isdir(p_abs):
                p_base, p_ext = os.path.splitext(p)
                p_dest = os.path.join(self.result_dir, p_base.upper()) \
                    + p_ext.upper()
                if self.result_suffix != None:
                    if (self.tag != None):
                        p_dest += '.' + self.tag
                    p_dest += '.' + self.result_suffix

                self.copy_result(p_abs, p_dest)

        # Handle solver-specific postprocessing output

        radiative_files = fnmatch.filter(dir_files, 'bord*')
        if len(radiative_files) > 0:
            self.copy_results_to_dir(radiative_files, 'CHR')

        lagr_files = fnmatch.filter(dir_files, 'deplacement*')
        lagr_files.extend(fnmatch.filter(dir_files, 'trajectoire*'))
        lagr_files.extend(fnmatch.filter(dir_files, 'frontiere*'))
        lagr_files.extend(fnmatch.filter(dir_files, 'debug*'))
        if len(lagr_files) > 0:
            self.copy_results_to_dir(lagr_files, 'LAGR')

        # Matisse output files

        res_matisse = os.path.join(self.exec_dir, 'resuMatisse')
        if os.path.isfile(res_matisse):
            matisse = get_command_output('grep -i matisse '
                                         + os.path.join(self.data_dir,
                                                        self.param))
            if date != None and len(matisse) > 0:
                # Add the date to the first line of resuMatisse
                f = open(res_matisse, 'r')
                s = f.read()
                f.close()
                f = open(res_matisse, 'w')
                f.write('Date of the case                                       : ' + date + '\n')
                f.write(s)
                f.close()
            self.copy_result('resuMatisse')

#-------------------------------------------------------------------------------

class syrthes_domain(base_domain):

    def __init__(self,
                 syrthes_env = 'syrthes.env',  # SYRTHES environment file
                 echo_comm = None,             # coupling verbosity
                 coupling_mode = 'MPI',        # 'MPI' or 'sockets'
                 coupled_apps = None):         # coupled domain names
                                               # if n_domains > 1


        base_domain.__init__(self, 1, 1, 1)

        self.log_file = 'syrthes.log'

        # Names, directories, and files in case structure
        self.data_dir = None
        self.result_dir = None
        self.src_dir = None
        self.result_suffix = None
        self.syrthes_env = 'syrthes.env'
        self.echo_comm = None

        self.set_coupling_mode(coupling_mode)

        self.coupled_apps = coupled_apps

    #---------------------------------------------------------------------------

    def set_case_dir(self, case_dir):

        base_domain.set_case_dir(self, case_dir)

        # Names, directories, and files in case structure

        self.data_dir = os.path.join(case_dir, 'DATA_SYR')
        self.src_dir = os.path.join(case_dir, 'SRC_SYR')

        if (self.tag != None):
            self.data_dir += '.' + self.tag
            self.src_dir += '.' + self.tag

    #---------------------------------------------------------------------------

    def set_coupling_mode(self, coupling_mode):

        # Check that coupling mode is either 'MPI' or 'sockets'
        coupling_modes = ('MPI', 'sockets')
        if coupling_mode not in coupling_modes:
            err_str = \
                'SYRTHES3 coupling mode "' + coupling_mode + '" unknown.\n' \
                + 'Allowed modes: ' + str(coupling_modes) + '.\n'
            raise RunCaseError(err_str)

        # Coupling mode
        self.coupling_mode = coupling_mode

    #---------------------------------------------------------------------------

    def set_verbosity(self, verbosity):

        self.echo_comm = verbosity

    #---------------------------------------------------------------------------

    def compile_and_link(self):
        """
        Compile and link user subroutines if necessary
        """

        # Check if there are files to compile in source path

        dir_files = os.listdir(self.src_dir)

        src_files = (fnmatch.filter(dir_files, '*.c')
                     + fnmatch.filter(dir_files, '*.[fF]'))

        copy_dir = None

        if len(src_files) > 0:

            # Add header files to list so as not to forget to copy them

            src_files = src_files + fnmatch.filter(dir_files, '*.h')

            # Copy source files to result directory

            copy_dir = os.path.join(self.result_dir, 'SRC_SYRTHES')
            if (self.result_suffix != None):
                if (self.tag != None):
                    copy_dir += '.' + self.tag
                copy_dir += '.' + self.result_suffix

            os.makedirs(copy_dir)

            for f in src_files:
                src_file = os.path.join(self.src_dir, f)
                dest_file = os.path.join(copy_dir, f)
                shutil.copy2(src_file, dest_file)
                try:
                    oldmode = (os.stat(dest_file)).stmode
                    newmode = oldmode & (stat.IRUSR | stat.IRGRP | stat.IROTH)
                    os.chmod(dest_file, newmode)
                except Exception:
                    pass

            # Copy source files to execution directory

            if (self.exec_dir != self.result_dir):
                exec_src = os.path.join(self.exec_dir, 'src_syrthes')
                os.mkdir(exec_src)
                for f in src_files:
                    src_file = os.path.join(self.src_dir, f)
                    dest_file = os.path.join(exec_src, f)
                    shutil.copy2(src_file, dest_file)

        log_name = os.path.join(self.exec_dir, 'compile_syrthes.log')
        log = open(log_name, 'w')

        # Note: src_dir and copy_dir should contain identical source files,
        #       but we prefer to use the copied version, so that if the
        #       source is later modified, possible debug information in an
        #       executable file will reference the correct (saved) version.

        retval = cs_compile.compile_and_link_syrthes(copy_dir,
                                                     self.exec_dir,
                                                     stdout=log,
                                                     stderr=log)

        log.close()
        self.copy_result(log_name)

        if retval != 0:
            raise RunCaseError('Compile or link error.')

        self.solver_path = os.path.join(self.exec_dir, 'syrthes')

    #---------------------------------------------------------------------------

    def set_exec_dir(self, exec_dir):

        if os.path.isabs(exec_dir):
            self.exec_dir = exec_dir
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_dir)

        if self.tag != None:
            self.exec_dir = os.path.join(self.exec_dir, 'SYR_' + self.tag)

        if not os.path.isdir(self.exec_dir):
            os.mkdir(self.exec_dir)

    #---------------------------------------------------------------------------

    def set_result_dir(self, name, by_suffix=True):
        """
        If suffix = true, add suffix to all names in result dir.
        Otherwise, create subdirectory
        """

        if by_suffix == True:
            self.result_dir = os.path.join(self.case_dir, 'RESU', 'RESU_SYR')
            self.result_dir += '.' + name

        else:
            self.result_dir = os.path.join(self.case_dir,
                                           'RESU',
                                           name,
                                           'RESU_SYR')
        if (self.tag != None):
            self.result_dir += '.' + self.tag

        if not os.path.isdir(self.result_dir):
            os.makedirs(self.result_dir)

    #---------------------------------------------------------------------------

    def solver_args(self, **kw):
        """
        Returns a tuple indicating SYRTHES's working directory,
        executable path, and associated command-line arguments.
        """

        wd = self.exec_dir              # Working directory
        exec_path = self.solver_path    # Executable

        # Build kernel command-line arguments

        args = ''
        if self.echo_comm != None:
            args += ' --echo-comm ' + str(self.echo_comm)

        if self.coupling_mode == 'MPI':
            args += ' --app-name ' + os.path.basename(self.case_dir)
            args += ' --comm-mpi ' + any_to_str(self.coupled_apps)

        elif self.coupling_mode == 'sockets':
            if 'host_port' in kw:
                args += ' --comm-socket ' + any_to_str(kw['host_port'])

        # handled directly

        # Adjust for Valgrind if used

        if self.valgrind != None:
            args = self.solver_path + ' ' + args
            exec_path = self.valgrind

        return wd, exec_path, args

    #---------------------------------------------------------------------------

    def copy_data(self):
        """
        Copy data to the execution directory
        """

        if os.path.isabs(self.syrthes_env):
            syrthes_env = self.syrthes_env
        else:
            syrthes_env = os.path.join(self.data_dir, self.syrthes_env)

        cmd = os.path.join(cs_config.dirs.pkgdatadir, 'runcase_syrthes')
        cmd += ' -copy-data -syrthes-env=' + syrthes_env

        if run_command(cmd) != 0:
            raise RunCaseError

    #---------------------------------------------------------------------------

    def copy_results(self):
        """
        Retrieve results from the execution directory
        """

        if self.exec_dir == self.result_dir and self.result_suffix == None:
            return

        cmd = os.path.join(cs_config.dirs.pkgdatadir, 'runcase_syrthes') \
            + ' -copy-results -result-dir=' + self.result_dir

        if self.result_suffix != None:
            if self.tag != None:
                cmd += '.' + self.tag
            cmd += '.' + self.result_suffix

        if run_command(cmd) != 0:
            raise RunCaseError

#-------------------------------------------------------------------------------
