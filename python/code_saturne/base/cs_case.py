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
import os
import os.path
import platform
import sys
import stat

from code_saturne.base import cs_exec_environment, cs_run_conf

from code_saturne.base.cs_case_domain import *

homard_prefix = None

#===============================================================================
# Utility functions
#===============================================================================

def check_exec_dir_stamp(d):

    """
    Return path to execution directory if present in a given directory.
    Otherwise, return given directory path
    """

    retval = d

    p = os.path.join(d, 'run_status.exec_dir')

    if os.path.isfile(p):
        f = open(p, 'r')
        l = f.readlines()
        f.close()
        d = l[0].rstrip()
        if os.path.isdir(d):
            retval = d

    return retval

#-------------------------------------------------------------------------------
# Check if an execution directory is inside a standard case structure
#
# This may not be the case if a separate staging directory or destination
# directory has been specified.
#-------------------------------------------------------------------------------

def is_exec_dir_in_case(case_dir, exec_dir):
    """
    Check if an execution directory is inside a standard case structure.
    """

    in_case = False

    if case_dir and exec_dir:
        r_case_dir = os.path.normpath(case_dir)
        r_exec_dir = os.path.normpath(exec_dir)
        if not os.path.isabs(r_case_dir):
            r_case_dir = os.path.abspath(r_case_dir)
        if not os.path.isabs(r_exec_dir):
            r_exec_dir = os.path.abspath(r_exec_dir)
        if r_exec_dir.find(r_case_dir) == 0:
            exec_parent_dir = os.path.split(r_exec_dir)[0]
            if os.path.basename(exec_parent_dir) in ('RESU', 'RESU_COUPLING'):
                if os.path.split(exec_parent_dir)[0] == case_dir:
                    in_case = True

    return in_case

#-------------------------------------------------------------------------------

def get_case_dir(case=None, param=None, coupling=None, id=None):

    """
    Check and return case path based on partial or sub-path.
    """

    casedir = None
    staging_dir = None
    data = None
    src = None

    if coupling:

        # Multiple domain case

        if param:
            err_str = 'Error: coupling and param options are incompatible.\n'
            raise RunCaseError(err_str)

        coupling = os.path.realpath(coupling)
        if not os.path.isfile(coupling):
            err_str = 'Error: coupling parameters: ' + coupling + '\n' \
                      'not found or not a file.\n'
            raise RunCaseError(err_str)

        if id:
            cwd = os.path.split(coupling)[0]
            if os.path.basename(cwd) == str(id):
                d = os.path.split(cwd)[0]
                if os.path.basename(d) == 'RESU_COUPLING':
                    staging_dir = cwd

        if case:
            casedir = os.path.realpath(case)
        else:
            casedir = os.path.split(coupling)[0]
            if staging_dir:
                casedir = os.path.split(os.path.split(staging_dir)[0])[0]

    else:

        cwd = os.getcwd()

        # Single domain case

        if param:
            param = os.path.basename(param)
            if param != param:
                datadir = os.path.split(os.path.realpath(param))[0]
                (casedir, data) = os.path.split(datadir)
                if data != 'DATA': # inconsistent paramaters location.
                    casedir = None

        if id:
            if os.path.basename(cwd) == str(id):
                d = os.path.split(cwd)[0]
                if os.path.basename(d) == 'RESU':
                    staging_dir = cwd

        if case:
            if os.path.isabs(case):
                casedir = os.path.realpath(case)
            else:
                casedir = os.path.realpath(os.path.join(cwd, case))
            data = os.path.join(casedir, 'DATA')
            src = os.path.join(casedir, 'SRC')
        else:
            testdir = cwd
            while os.path.basename(testdir):
                data = os.path.join(testdir, 'DATA')
                src = os.path.join(testdir, 'SRC')
                cfg = os.path.join(data, 'run.cfg')
                if (os.path.isdir(data) and os.path.isdir(src)) \
                   or os.path.isfile(cfg):
                    casedir = testdir
                    break
                testdir = os.path.split(testdir)[0]

        if not (os.path.isdir(data)):
            casedir = None

    # Return casedir or None

    return casedir, staging_dir


#===============================================================================
# Main class
#===============================================================================

class case:
    """Base class from which classes handling running case should inherit."""

    #---------------------------------------------------------------------------

    def __init__(self,
                 package,                     # main package
                 package_compute = None,      # package for compute environment
                 case_dir = None,
                 dest_dir = None,
                 staging_dir = None,
                 domains = None,
                 syr_domains = None,
                 py_domains = None):

        # Package specific information

        self.package = package

        if package_compute:
            self.package_compute = package_compute
        else:
            self.package_compute = self.package

        # Determine caller script path (which may contain batch directives)

        self.parent_script = cs_exec_environment.get_parent_process_path()

        # Set environment modules if present

        cs_exec_environment.set_modules(self.package_compute)
        cs_exec_environment.source_rcfile(self.package_compute)

        # Re-clean os.environment as the above calls may have indirectly
        # restored some unwanted variables (such as lmod macros).

        cs_exec_environment.clean_os_environ_for_shell()

        # Ensure we have tuples or lists to simplify later tests

        if type(domains) == list:
            self.domains = domains
        elif type(domains) == tuple:
            self.domains = list(domains)
        else:
            self.domains = [domains]

        self.n_domains_native = self.domains  # code_saturne based domains

        if type(syr_domains) == tuple :
            self.domains += list(syr_domains)
        elif type(syr_domains) == list:
            self.domains += syr_domains
        elif syr_domains:
            self.domains += [syr_domains]

        n_domains = len(self.domains)

        if type(py_domains) == tuple:
            self.domains += list(py_domains)
        elif type(py_domains) == list:
            self.domains += py_domains
        elif py_domains:
            self.domains = [py_domains]

        self.n_domains_python = len(self.domains) - n_domains

        n_domains = len(self.domains)

        # Check names in case of multiple domains (coupled by MPI)

        if n_domains > 1:
            err_str = 'In case of multiple domains (i.e. code coupling), ' \
                + 'each domain must have a name.\n'
            for d in self.domains:
                if d.name is None:
                    raise RunCaseError(err_str)

        # Names, directories, and files in case structure:
        # associate case domains and set case directory

        self.case_dir = case_dir
        self.dest_dir = dest_dir
        self.__define_dest_dir__()

        self.result_dir = None

        # Determine parent (or coupling group) study directory

        if (os.path.isdir(os.path.join(case_dir, 'DATA')) and n_domains == 1):
            # Simple domain case from standard directory structure
            (self.study_dir, self.name) = os.path.split(case_dir)
            self.domains[0].name = None

        else:
            # Coupling or single domain run from coupling script
            self.study_dir = case_dir
            self.name = os.path.split(case_dir)[1]

        if staging_dir:
            staging_dir = check_exec_dir_stamp(staging_dir)

        n_nc_solver = 0

        for d in self.domains:
            d.set_case_dir(self.case_dir, staging_dir)
            try:
                solver_name = os.path.basename(d.solver_path)
                if solver_name == 'nc_solver':
                    n_nc_solver += 1
            except Exception:
                pass

        self.module_name = 'code_saturne'
        if n_nc_solver == self.n_domains_native:
            self.module_name = 'neptune_cfd'

        # Working directory

        self.exec_dir = None
        self.exec_prefix = None

        self.exec_solver = True

        n_exec_solver = 0
        for d in self.domains:
            if d.exec_solver:
                n_exec_solver += 1

        if n_exec_solver == 0:
            self.exec_solver = False
        elif n_exec_solver < n_domains:
            err_str = 'In case of multiple domains (i.e. code coupling), ' \
                + 'all or no domains must execute their solver.\n'
            raise RunCaseError(err_str)

        # Script hooks

        self.run_prologue = None
        self.run_epilogue = None
        self.compute_prologue = None
        self.compute_epilogue = None

        # Tool hooks

        self.debug_args = None
        self.tool_args = None
        self.mpi_tool_args = None

        # Date or other name

        self.run_id = None

        # Time limit

        self.time_limit = None

        # Error reporting
        self.error = ''
        self.error_long = ''

    #---------------------------------------------------------------------------

    def __define_dest_dir__(self):

        if self.dest_dir != None:
            if not os.path.isabs(self.dest_dir):
                dest_dir = os.path.join(os.path.split(self.case_dir)[0],
                                        self.dest_dir)
                self.dest_dir = os.path.abspath(dest_dir)

    #---------------------------------------------------------------------------

    def print_procs_distribution(self):

        """
        Print info on the processor distribution.
        """

        # Print process info

        name = self.module_name

        for d in self.domains:
            if len(self.domains) == 1:
                if d.n_procs > 1:
                    msg = ' Parallel ' + d.code_name + ' on ' \
                          + str(d.n_procs) + ' processes.\n'
                else:
                    msg = ' Single processor ' + d.code_name + ' simulation.\n'
                sys.stdout.write(msg)
            else:
                msg = ' ' + d.code_name + ' domain "' + d.name + '" on ' \
                    + str(d.n_procs) + ' process(es).\n'
                sys.stdout.write(msg)

        sys.stdout.write('\n')

    #---------------------------------------------------------------------------

    def distribute_procs(self, n_procs=None):

        """
        Distribute number of processors to domains,
        and returns the total number of associated processes.
        """

        # Initial count

        np_list = []

        for d in self.domains:
            np_list.append(d.get_n_procs())

        n_procs_tot = 0
        n_procs_min = 0
        for np in np_list:
            n_procs_tot += np[0]
            n_procs_min += np[1]

        # If no process count is given or everything fits:

        if n_procs is None or n_procs == n_procs_tot:
            return n_procs_tot

        # If process count is given and not sufficient, abort.

        elif n_procs < n_procs_min:
            err_str = ' Error:\n' \
                '   The current calculation scheme requires at least ' \
                + str(n_procs_min) + ' processes,\n' \
                + '   while the execution environment provides only ' \
                + str(n_procs) + '.\n' \
                + '   You may either allocate more processes or try to\n' \
                + '   oversubscribe by changing the number of processes\n' \
                + '   in the run.cfg configuration file.'
            raise RunCaseError(err_str)

        # Otherwise, rebalance process counts.

        n_passes = 0
        n_fixed_procs = 0
        n_fixed_apps = 0

        while (    n_procs_tot != n_procs and n_fixed_apps != len(np_list)
               and n_passes < 5):

            mult = (n_procs - n_fixed_procs) *1.0 \
                // (n_procs_tot - n_fixed_procs)

            n_procs_tot = 0

            for np in np_list:
                if np[2] is None:
                    np[2] = n_procs
                if np[2] != 0: # marked as fixed here
                    np[0] = int(np[0]*mult)
                    if np[0] < np[1]:
                        np[0] = np[1]
                        np[2] = 0 # used as marker here
                        n_fixed_procs += np[0]
                        n_fixed_apps += 1
                    elif np[0] > np[2]:
                        np[0] = np[2]
                        np[2] = 0 # used as marker here
                        n_fixed_procs += np[0]
                        n_fixed_apps += 1
                n_procs_tot += np[0]

            n_passes += 1

        # Minor adjustments may help in case of rounding

        n_passes = 0

        while (    n_procs_tot != n_procs and n_fixed_apps != len(np_list)
               and n_passes < 5):

            delta = int(round(  (n_procs_tot - n_procs)*1.0
                              // (len(np_list) - n_fixed_apps)))

            if delta == 0:
                if n_procs_tot < n_procs:
                    delta = 1
                else:
                    delta = -1

            for np in np_list:
                if np[2] != 0: # marked as fixed here
                    n_procs_prev = np[0]
                    np[0] -= delta
                    if np[0] < np[1]:
                        np[0] = np[1]
                        np[2] = 0 # used as marker here
                        n_fixed_procs += np[0]
                        n_fixed_apps += 1
                    elif np[0] > np[2]:
                        np[0] = np[2]
                        np[2] = 0 # used as marker here
                        n_fixed_procs += np[0]
                        n_fixed_apps += 1
                    n_procs_tot += (np[0] - n_procs_prev)
                if n_procs_tot == n_procs:
                    break

            n_passes += 1

        # We are now are ready to set return values

        app_id = 0

        for d in self.domains:
            d.set_n_procs(np_list[app_id][0])
            app_id += 1

        # Warn in case of over or under-subscription

        if n_procs_tot != n_procs:
            msg =  \
                'Warning:\n' \
                + '  total number of processes used (' + str(n_procs_tot) \
                + ') is different from the\n' \
                + '  number provided by the environment. (' + str(n_procs) \
                + ').\n\n'
            sys.stderr.write(msg)

        return n_procs_tot

    #---------------------------------------------------------------------------

    def define_exec_dir(self):
        """
        Define execution directory.
        """

        if self.exec_prefix != None:
            if self.case_dir != self.study_dir:
                study_name = os.path.split(self.study_dir)[1]
                exec_dir_name = study_name + '.' + self.name
            else:
                exec_dir_name = self.name
            exec_dir_name += '.' + self.run_id
            self.exec_dir = os.path.join(self.exec_prefix, exec_dir_name)
        else:
            if not self.result_dir:
                self.define_result_dir()
            self.exec_dir = self.result_dir

    #---------------------------------------------------------------------------

    def set_exec_dir(self, force=False):
        """
        Create execution directory.
        """

        # Define execution directory name

        self.define_exec_dir()

        # Check that directory does not exist.

        if os.path.isdir(self.exec_dir):
            if not force and self.result_dir != self.exec_dir:
                err_str = \
                    '\nWorking directory: ' + self.exec_dir \
                    + ' already exists.\n' \
                    + 'Calculation will not be run.\n'
                raise RunCaseError(err_str)

        else:
            os.makedirs(self.exec_dir)

        # Set execution directory

        for d in self.domains:
            d.set_exec_dir(self.exec_dir)

    #---------------------------------------------------------------------------

    def define_result_dir(self):

        if self.dest_dir != None:
            base_dir = os.path.join(self.dest_dir,
                                    os.path.basename(self.case_dir))
        else:
            base_dir = self.case_dir

        r = os.path.join(base_dir, 'RESU')

        # Coupled case
        if len(self.domains) > 1:
            r += '_COUPLING'

        if not os.path.isdir(r):
            os.makedirs(r)

        self.result_dir = os.path.join(r, self.run_id)

    #---------------------------------------------------------------------------

    def set_result_dir(self, force=True):

        self.define_result_dir()

        if os.path.isdir(self.result_dir):
            if not force:
                err_str = \
                    '\nResults directory: ' + self.result_dir \
                    + ' already exists.\n' \
                    + 'Calculation will not be run.\n'
                raise RunCaseError(err_str)

        else:
            os.mkdir(self.result_dir)

        dest_root_dir = None
        if self.dest_dir != None:
            dest_root_dir = os.path.join(self.dest_dir,
                                         os.path.basename(self.case_dir))

        for d in self.domains:
            d.set_result_dir(self.run_id,
                             self.result_dir,
                             dest_root_dir=dest_root_dir)

    #---------------------------------------------------------------------------

    def summary_init(self, exec_env):

        """
        Build summary start.
        """

        s_path = os.path.join(self.exec_dir, 'summary')
        s = open(s_path, 'w')

        r = exec_env.resources

        date = (datetime.datetime.now()).strftime("%A %B %d %H:%M:%S CEST %Y")
        t_uname = platform.uname()
        s_uname = ''
        for t in t_uname:
            s_uname = s_uname + t + ' '
        n_procs = r.n_procs
        if n_procs is None:
            n_procs = 1

        dhline = '========================================================\n'
        hline =  '  ----------------------------------------------------\n'

        s.write(dhline)
        s.write('Start time       : ' + date + '\n')
        s.write(dhline)
        cmd = sys.argv[0]
        for arg in sys.argv[1:]:
            cmd += ' ' + arg
        s.write('  Command        : ' + cmd + '\n')
        s.write(dhline)
        s.write('  Package path   : ' + self.package.get_dir('exec_prefix') + '\n')
        s.write(hline)
        if homard_prefix != None:
            s.write('  HOMARD          : ' + homard_prefix + '\n')
            s.write(hline)
        s.write('  MPI path       : ' + self.package_compute.config.libs['mpi'].bindir + '\n')
        if len(self.package_compute.config.libs['mpi'].variant) > 0:
            s.write('  MPI type       : ' + self.package_compute.config.libs['mpi'].variant + '\n')
        s.write(hline)
        s.write('  User           : ' + exec_env.user  + '\n')
        s.write(hline)
        s.write('  Machine        : ' + s_uname  + '\n')
        s.write('  N Procs        : ' + str(n_procs)  + '\n')
        if r.manager is None and r.hosts_list != None:
            s.write('  Hosts          :')
            for p in r.hosts_list:
                s.write(' ' + p)
            s.write('\n')
        s.write(hline)

        if len(self.domains) > 1:
            s.write('  Exec. dir.     : ' + self.exec_dir + '\n')
            s.write(hline)

        s.write('  Environment variables\n')
        ek = list(os.environ.keys())
        ek.sort()
        for k in ek:
            s.write('    ' + k + '=' + os.environ[k] + '\n')
        s.write(hline)

        for d in self.domains:
            d.summary_info(s)
            s.write(hline)

        s.close()

    #---------------------------------------------------------------------------

    def summary_finalize(self):

        """
        Output summary body.
        """

        date = (datetime.datetime.now()).strftime("%A %B %d %H:%M:%S CEST %Y")

        dhline = '========================================================\n'
        hline =  '  ----------------------------------------------------\n'

        s_path = os.path.join(self.exec_dir, 'summary')
        s = open(s_path, 'a')

        if self.error:
            s.write('  ' + self.error + ' failed\n')

        s.write(dhline)
        s.write('Finish time      : ' + date + '\n')
        s.write(dhline)

        s.close()

    #---------------------------------------------------------------------------

    def copy_log(self, name):
        """
        Retrieve single log file from the execution directory
        """

        if self.result_dir == self.exec_dir:
            return

        src = os.path.join(self.exec_dir, name)
        dest = os.path.join(self.result_dir, name)

        if os.path.isfile(src):
            shutil.copy2(src, dest)
            if self.error == '':
                os.remove(src)

    #---------------------------------------------------------------------------

    def copy_script(self):
        """
        Retrieve script (if present) from the execution directory
        """

        src = sys.argv[0]

        if os.path.basename(src) == self.package.name:
            src = self.parent_script

        if not src:
            src = self.package.name

        # Copy single file

        dest = os.path.join(self.result_dir, os.path.basename(src))
        if os.path.isfile(src) and os.path.abspath(src) != os.path.abspath(dest):
            shutil.copy2(src, dest)

    #---------------------------------------------------------------------------

    def copy_top_run_conf(self):
        """
        Retrieve top-level run.cfg (if present) from the top directory
        """

        n = "run.cfg"

        src = os.path.join(self.case_dir, n)
        dest = os.path.join(self.result_dir, n)

        if os.path.isfile(src) and src != dest:
            shutil.copy2(src, dest)

    #---------------------------------------------------------------------------

    def solver_script_path(self):
        """
        Determine name of solver script file.
        """
        return os.path.join(self.exec_dir, self.package.runsolver)

    #---------------------------------------------------------------------------

    def debug_wrapper_args(self):
        """
        Additional arguments to use debug wrapper
        """
        debug_args = ''
        python_exec = self.package.config.python
        if os.path.isfile(python_exec) or os.path.islink(python_exec):
            debug_args += python_exec + ' '
        else:
            debug_args += 'python '
        cs_pkg_dir = self.package.get_dir('pkgpythondir')
        if self.package.name != 'code_saturne':
            cs_pkg_dir = os.path.join(cs_pkg_dir, '../code_saturne')
        dbg_wrapper_path = os.path.join(cs_pkg_dir,
                                        'base', 'cs_debug_wrapper.py')
        debug_args += dbg_wrapper_path + ' '
        return debug_args

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_mpiexec(self, n_procs, mpi_env, tool_args=''):
        """
        Generate MPMD mpiexec command.
        """

        cmd = ''

        app_id = 0

        for d in self.domains:
            s_args = d.solver_command()
            if len(cmd) > 0:
                cmd += ' : '
            cmd += '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + tool_args + s_args[1] + s_args[2]
            app_id += 1

        return cmd

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile(self,
                                        n_procs,
                                        mpi_env,
                                        tool_args=''):
        """
        Generate MPMD mpiexec config file.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')
        e = open(e_path, 'w')

        app_id = 0

        for d in self.domains:
            s_args = d.solver_command()
            cmd = '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + tool_args + s_args[1] + s_args[2] + '\n'
            e.write(cmd)
            app_id += 1

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile_srun(self,
                                             n_procs,
                                             mpi_env,
                                             tool_args=''):
        """
        Generate MPMD mpiexec config file for SLURM's srun.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')
        e = open(e_path, 'w')

        e.write('# srun multiple program configuration file\n#\n')
        e.write('# > srun -n ' + str(n_procs) \
                + ' --multi-prog mpmd_configfile\n#\n')

        rank_id = 0

        for d in self.domains:
            s_args = d.solver_command(need_abs_path=True)

            cmd = '%d-%d\t' % (rank_id, rank_id + d.n_procs - 1) \
                   + tool_args + s_args[1] + s_args[2] \
                   + ' -wdir ' + os.path.basename(s_args[0]) + '\n'

            e.write(cmd)
            rank_id += d.n_procs

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile_ccc_mprun(self,
                                                  n_procs,
                                                  mpi_env,
                                                  tool_args=''):
        """
        Generate MPMD mpiexec config file for CCRT's ccc_mprun.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')

        e = open(e_path, 'w')

        app_id = 0

        for d in self.domains:
            s_args = d.solver_command()
            cmd = str(d.n_procs) +' ' \
                    + ' bash -c "' \
                    + ' cd ' + os.path.basename(s_args[0]) + ' &&' \
                    + ' ' + s_args[1] + s_args[2] + '"\n'

            e.write(cmd)
            app_id += 1

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_script(self, n_procs,
                                    mpi_env,
                                    tool_args='',
                                    use_mps=False):
        """
        Generate MPMD dispatch file.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_exec.sh')
        e = open(e_path, 'w')

        cs_exec_environment.write_shell_shebang(e)

        e.write('# Make sure to transmit possible additional '
                + 'arguments assigned by mpirun.\n'
                + '# we use $@ to forward arguments passed to this script'
                + ' to the executable files.\n\n')

        e.write('MPI_RANK=`'
                + self.package.get_pkgdatadir_script('runcase_mpi_rank')
                + ' $@`\n')

        if use_mps:
            cs_exec_environment.write_mps_start(e)

        app_id = 0
        nr = 0
        test_pf = 'if [ $MPI_RANK -lt '
        test_sf = ' ] ; then\n'

        for d in self.domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_command()
            e.write('  cd ' + s_args[0] + '\n')
            e.write('  ' + tool_args + s_args[1] + s_args[2] + ' $@\n')
            if app_id == 0:
                test_pf = 'el' + test_pf
            app_id += 1

        e.write('fi\n'
                + 'CS_RET=$?\n')

        if use_mps:
            cs_exec_environment.write_mps_stop(e)

        e.write('exit $CS_RET\n')

        e.close()

        oldmode = (os.stat(e_path)).st_mode
        newmode = oldmode | (stat.S_IXUSR)
        os.chmod(e_path, newmode)

        return e_path

    #---------------------------------------------------------------------------

    def solver_script_body(self, n_procs, mpi_env, s):
        """
        Write commands to solver script.
        """

        # Determine if an MPMD syntax (mpiexec variant) will be used

        mpiexec_mpmd = False
        if len(self.domains) > 1:
            if mpi_env.mpmd & cs_exec_environment.MPI_MPMD_mpiexec:
                mpiexec_mpmd = True
            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_configfile:
                mpiexec_mpmd = True

        use_srun = False

        # Start assembling command

        mpi_cmd = ''
        mpi_cmd_exe = ''
        mpi_cmd_args = ''

        if self.mpi_tool_args:
            mpi_cmd = self.mpi_tool_args + ' '

        if mpi_env.mpiexec != None:
            if os.path.basename(mpi_env.mpiexec)[:4] == 'srun':
                use_srun = True
            if (n_procs > 1 or use_srun):
                mpi_cmd += mpi_env.mpiexec

        if mpi_cmd:
            if mpi_env.mpiexec_opts != None:
                mpi_cmd += ' ' + mpi_env.mpiexec_opts
            if mpiexec_mpmd == False:
                if mpi_env.mpiexec_n != None and n_procs:
                    mpi_cmd += mpi_env.mpiexec_n + str(n_procs)
                if mpi_env.mpiexec_n_per_node != None:
                    mpi_cmd += mpi_env.mpiexec_n_per_node
                if mpi_env.mpiexec_separator != None:
                    mpi_cmd += ' ' + mpi_env.mpiexec_separator
                mpi_cmd += ' '
                if mpi_env.mpiexec_exe != None:
                    mpi_cmd += mpi_env.mpiexec_exe + ' '
                if mpi_env.mpiexec_args != None:
                    mpi_cmd_args = mpi_env.mpiexec_args + ' '
            else:
                mpi_cmd += ' '

        mpmd = mpi_env.mpmd
        if use_srun:
            mpmd = cs_exec_environment.MPI_MPMD_configfile

        # Additional (rank local) tool-related arguments

        tool_args = ''

        if self.tool_args:
            tool_args = self.tool_args + ' '
        if self.debug_args:
            tool_args += self.debug_wrapper_args() + self.debug_args + ' '

        # Check if we need MPS. If this is the case, force the
        # MPMD mode to script.
        use_mps = os.getenv('CS_USE_MPS')
        if use_mps:
            if use_mps.lower() not in ('0', 'false', 'no'):
                use_mps = True
                mpmd = cs_exec_environment.MPI_MPMD_script
        if use_mps != True:
            use_mps = False

        # Case with only one cs_solver instance possibly under MPI

        n_domains = len(self.domains)

        if n_domains == 1 and not use_mps:

            s_args = self.domains[0].solver_command()

            cs_exec_environment.write_script_comment(s, 'Run solver.\n')
            s.write(mpi_cmd + mpi_cmd_args + tool_args + s_args[1] + s_args[2])
            s.write(' ' + cs_exec_environment.get_script_positional_args() +
                    '\n')

        # General case

        else:

            if mpmd & cs_exec_environment.MPI_MPMD_mpiexec:

                if mpi_env.mpiexec_separator != None:
                    mpi_cmd += mpi_env.mpiexec_separator + ' '

                e_path = self.generate_solver_mpmd_mpiexec(n_procs,
                                                           mpi_env,
                                                           tool_args)

            elif mpmd & cs_exec_environment.MPI_MPMD_configfile:

                if mpi_env.mpiexec == 'srun':
                    e_path = self.generate_solver_mpmd_configfile_srun(n_procs,
                                                                       mpi_env,
                                                                       tool_args)
                    mpi_cmd += '--multi-prog ' + e_path

                elif mpi_env.mpiexec == 'ccc_mprun':
                    e_path = self.generate_solver_mpmd_configfile_ccc_mprun(n_procs,
                                                                            mpi_env,
                                                                            tool_args)
                    mpi_cmd += '-f ' + e_path
                else:
                    e_path = self.generate_solver_mpmd_configfile(n_procs,
                                                                  mpi_env,
                                                                  tool_args)
                    mpi_cmd += '-configfile ' + e_path

                e_path = ''

            elif mpmd & cs_exec_environment.MPI_MPMD_script:

                if mpi_env.mpiexec_separator != None:
                    mpi_cmd += mpi_env.mpiexec_separator + ' '

                e_path = self.generate_solver_mpmd_script(n_procs, mpi_env,
                                                          tool_args, use_mps)

            else:
                raise RunCaseError(' No allowed MPI MPMD mode defined.\n')

            s.write(mpi_cmd + e_path + '\n')

    #---------------------------------------------------------------------------

    def generate_solver_script(self, exec_env):
        """
        Generate localexec file.
        """

        # Initialize simple solver command script

        s_path = self.solver_script_path()

        s = open(s_path, 'w')

        cs_exec_environment.write_shell_shebang(s)

        # If n_procs not already given by environment, determine it

        n_procs = exec_env.resources.n_procs
        mpi_env = exec_env.mpi_env

        if n_procs is None:
            n_procs = 0
            for d in self.domains:
                n_procs += d.n_procs

        # Set PATH for Windows DLL search PATH

        if sys.platform.startswith('win'):
            cs_exec_environment.write_script_comment(s,
                'Ensure the correct command is found:\n')
            cs_exec_environment.write_prepend_path(s,
                                                   'PATH',
                                                   self.package.get_dir("bindir"))
            cs_exec_environment.write_prepend_path(s,
                                                   'PATH',
                                                   self.package.get_dir("libdir"))
            s.write('set CS_ROOT_DIR=' + self.package.get_dir("exec_prefix") + '\n')

        # Add MPI directories to PATH if in nonstandard path

        cs_exec_environment.write_script_comment(s, \
            'Export paths here if necessary or recommended.\n')
        mpi_bindir = self.package_compute.config.libs['mpi'].bindir
        if len(mpi_bindir) > 0:
            cs_exec_environment.write_prepend_path(s, 'PATH', mpi_bindir)
        mpi_libdir = self.package_compute.config.libs['mpi'].libdir
        if len(mpi_libdir) > 0:
            cs_exec_environment.write_prepend_path(s,
                                                   'LD_LIBRARY_PATH',
                                                   mpi_libdir)
        s.write('\n')

        # Handle cathare coupling
        wrote_cathare_path = False
        if self.domains:
            for d in self.domains:
                if hasattr(d, "cathare_case_file"):
                    if not wrote_cathare_path:
                        cs_exec_environment.write_script_comment(s, \
                            'Export paths necessary for CATHARE coupling.\n')
                        config = configparser.ConfigParser()
                        config.read(self.package.get_configfiles())
                        cathare_path = config.get('install', 'cathare')
                        for p in ['lib', 'ICoCo/lib']:
                            lp = os.path.join(cathare_path, p)
                            cs_exec_environment.write_prepend_path( \
                                    s, "LD_LIBRARY_PATH", lp)

                        wrote_cathare_path = True
                        s.write('\n')

        # Handle python coupling
        if self.n_domains_python > 0:
            cs_exec_environment.write_script_comment(s, \
                'Export paths necessary for python coupling.\n')
            pydir = self.package_compute.get_dir("pythondir")
            cs_exec_environment.write_prepend_path(s, "PYTHONPATH", pydir)
            # This export is necessary fr PyPLE to work correctly
            libdir = self.package_compute.get_dir("libdir")
            cs_exec_environment.write_prepend_path(s,
                                                   "LD_LIBRARY_PATH",
                                                   libdir)
            s.write('\n')

        # Handle environment modules if used

        if self.package_compute.config.env_modules != "no":
            cs_exec_environment.write_script_comment(s, \
               'Load environment if this script is run directly.\n')
            s.write('if test "$CS_ENVIRONMENT_SET" != "true" ; then\n')
            if self.package_compute.config.env_modules != "no":
                s.write('  module purge\n')
                for m in self.package_compute.config.env_modules.strip().split():
                    s.write('  module load ' + m + '\n')
            s.write('fi\n\n')

        # Add additional library search paths in case DT_RUNPATH
        # is used instead of DT_RPATH in ELF library

        add_lib_dirs = get_ld_library_path_additions(self.package)
        while add_lib_dirs:
            d = add_lib_dirs.pop()
            cs_exec_environment.write_prepend_path(s,
                                                   'LD_LIBRARY_PATH',
                                                   d)

        # Add paths for plugins or dynamic library dependencies

        lib_dirs, plugin_pythonpath_dirs, plugin_env_vars \
            = self.package_compute.config.get_run_environment_dependencies()
        for d in lib_dirs:
            cs_exec_environment.write_prepend_path(s,
                                                   'LD_LIBRARY_PATH',
                                                   d)
        for d in plugin_pythonpath_dirs:
            cs_exec_environment.write_prepend_path(s,
                                                   'PYTHONPATH',
                                                   d)

        # Add paths for local libraries
        if self.package_compute.config.features['relocatable'] == "yes":
            d = self.package_compute.get_dir('libdir')
            if d != '/usr/lib' and d != '/usr/local/lib':
                cs_exec_environment.write_prepend_path(s,
                                                       'LD_LIBRARY_PATH',
                                                       d)

        # Add additional environment variables

        for v in plugin_env_vars:
            cs_exec_environment.write_export_env(s, v, plugin_env_vars[v])

        s.write('\n')

        # Handle OpenMP if needed

        n_threads = exec_env.resources.n_threads
        if self.package_compute.config.features['openmp'] == 'yes' or n_threads:
            if not n_threads:
                n_threads = 1
            cs_exec_environment.write_export_env(s, 'OMP_NUM_THREADS',
                                                 str(n_threads))
            s.write('\n')

        # Handle rcfile if used

        rcfile = cs_exec_environment.get_rcfile(self.package_compute)

        if rcfile and rcfile != "no":
            cs_exec_environment.write_script_comment(s, \
               'Load environment if this script is run directly.\n')
            s.write('if test "$CS_ENVIRONMENT_SET" != "true" ; then\n')
            if rcfile:
                s.write('  source ' + rcfile + '\n')
            s.write('fi\n\n')

        # Boot MPI daemons if necessary

        if mpi_env.gen_hostsfile != None:
            cs_exec_environment.write_script_comment(s, 'Generate hostsfile.\n')
            s.write(mpi_env.gen_hostsfile + ' || exit $?\n\n')

        if n_procs > 1 and mpi_env.mpiboot != None:
            cs_exec_environment.write_script_comment(s, 'Boot MPI daemons.\n')
            s.write(mpi_env.mpiboot + ' || exit $?\n\n')

        # Ensure we are in the correct directory

        s.write('cd ' + cs_exec_environment.enquote_arg(self.exec_dir) + '\n\n')

        # Add user-defined prologue if defined

        if self.compute_prologue:
            s.write('# Compute prologue\n')
            s.write(self.compute_prologue)
            s.write('\n\n')

        # Generate script body

        self.solver_script_body(n_procs, mpi_env, s)

        # Obtain return value (or sum thereof)

        cs_exec_environment.write_export_env(s, 'CS_RET',
                                             cs_exec_environment.get_script_return_code())

        # Halt MPI daemons if necessary

        if n_procs > 1 and mpi_env.mpihalt != None:
            cs_exec_environment.write_script_comment(s, 'Halt MPI daemons.\n')
            s.write(mpi_env.mpihalt + '\n\n')

        if mpi_env.del_hostsfile != None:
            cs_exec_environment.write_script_comment(s, 'Remove hostsfile.\n')
            s.write(mpi_env.del_hostsfile + '\n\n')

        # Add user-defined epilogue if defined

        if self.compute_epilogue:
            s.write('\n# Compute epilogue\n')
            s.write(self.compute_epilogue)
            s.write('\n')

        if sys.platform.startswith('win'):
            s.write('\nexit %CS_RET%\n')
        else:
            s.write('\nexit $CS_RET\n')
        s.close()

        oldmode = (os.stat(s_path)).st_mode
        newmode = oldmode | (stat.S_IXUSR)
        os.chmod(s_path, newmode)

        return s_path

    #---------------------------------------------------------------------------

    def update_scripts_tmp(self, src, dest, caption=None):

        """
        Create a stamp file in the execution directory, rename it, or destroy it.
        """

        # Create a temporary file to determine status

        src_tmp_name = None
        dest_tmp_name = None

        base_path = os.path.join(self.exec_dir, 'run_status.')
        if not os.path.isdir(base_path):
            base_path = os.path.join(self.result_dir, 'run_status.')

        if src != None:
            if type(src) == tuple:
                for s in src:
                    p = base_path + s
                    if os.path.isfile(p):
                        src_tmp_name = p

            else:
                p = base_path + src
                if os.path.isfile(p):
                    src_tmp_name = p

        try:
            if dest != None:
                dest_tmp_name = base_path + dest
                if src_tmp_name is None:
                    scripts_tmp = open(dest_tmp_name, 'w')
                    if not caption:
                        caption = self.case_dir
                    scripts_tmp.write(caption + '\n')
                    scripts_tmp.close()
                else:
                    os.rename(src_tmp_name, dest_tmp_name)

            else:
                os.remove(src_tmp_name)

        except Exception:
            pass

    #---------------------------------------------------------------------------

    def add_exec_dir_stamp(self):

        """
        Create a file with the execution directory path in the results directory.
        """
        p = os.path.join(self.result_dir, 'run_status.exec_dir')

        f = open(p, 'w')
        f.write(self.exec_dir + '\n')
        f.close()

    #---------------------------------------------------------------------------

    def clear_exec_dir_stamp(self):

        """
        Remove file with the execution directory path from the results directory.
        """
        p = os.path.join(self.result_dir, 'run_status.exec_dir')

        try:
            os.remove(p)
        except Exception:
            pass

    #---------------------------------------------------------------------------

    def exceeded_time_limit(self):

        """
        Check if the allocated time limit was exceeded.
        """
        p = os.path.join(self.exec_dir, 'run_status.exceeded_time_limit')

        if os.path.isfile(os.path.join(p)):
            return True
        else:
            return False

    #---------------------------------------------------------------------------

    def set_run_id(self,
                   run_id = None):

        """
        Set run id for calculation.
        """

        if run_id != None:
            self.run_id = run_id
            self.define_result_dir()

        # When id not assigned, choose an id not already present

        if self.run_id is None:
            self.run_id = self.suggest_id()
            self.define_result_dir()
            j = 1
            run_id_base = self.run_id
            while os.path.isdir(self.result_dir):
                self.run_id = run_id_base + '_' + str(j)
                self.define_result_dir()

    #---------------------------------------------------------------------------

    def prepare_data(self,
                     force_id = False):

        """
        Prepare data for calculation.
        """

        # Before creating or generating file, create stage 'marker' file.

        self.update_scripts_tmp(None, 'preparing')

        # Copy script before changing to the working directory
        # (as after that, the relative path will not be up to date).

        self.copy_script()
        self.copy_top_run_conf()

        os.chdir(self.exec_dir)

        # Set (read and possibly modify) script parameters

        sys.stdout.write('Copying base setup data\n'
                         '-----------------------\n\n')

        for d in self.domains:
            d.copy_data()

        # Compile user subroutines if necessary
        # (for some domain types, such as for Syrthes, this may be done later,
        # during the general prepare_data stage).

        need_compile = False

        for d in self.domains:
            if not hasattr(d, 'needs_compile'):
                continue
            if d.needs_compile() == True:
                if need_compile == False: # Print banner on first pass
                    need_compile = True
                    msg = \
                        "Compiling and linking user-defined functions\n" \
                        "--------------------------------------------\n\n"
                    sys.stdout.write(msg)
                    sys.stdout.flush()
                d.compile_and_link()
                if len(d.error) > 0:
                    self.error = d.error
                    if len(d.error_long) > 0:
                        self.error_long = d.error_long

        if len(self.error) > 0:
            self.update_scripts_tmp('preparing', 'failed', self.error)
            return 1

        # Setup data
        #===========

        sys.stdout.write('Preparing calculation data\n'
                         '--------------------------\n\n')
        sys.stdout.flush()

        for d in self.domains:
            d.prepare_data()
            if len(d.error) > 0:
                self.error = d.error

        # Set run_id in run.cfg as a precaution

        run_conf_path = os.path.join(self.result_dir, "run.cfg")
        if os.path.isfile(run_conf_path):
            run_conf = cs_run_conf.run_conf(run_conf_path, package=self.package)
            run_conf.set('run', 'id', str(self.run_id))
            run_conf.save()

        # Rename temporary file to indicate new status

        if len(self.error) == 0:
            status = 'prepared'
        else:
            status = 'failed'

        self.update_scripts_tmp('preparing', status, self.error)

        # Standard or error exit

        if len(self.error) > 0:
            err_str = ' Error in ' + self.error + ' stage.\n\n'
            if len(self.error_long) > 0:
                err_str += self.error_long + '\n\n'
            sys.stderr.write(err_str)
            return 1
        else:
            return 0

    #---------------------------------------------------------------------------

    def init_prepared_data(self):

        """
        Initialize prepared data for calculation.
        """

        self.copy_script()

        os.chdir(self.exec_dir)

        for d in self.domains:
            d.init_staged_data()

    #---------------------------------------------------------------------------

    def preprocess(self,
                   n_procs = None,
                   n_threads = None,
                   mpiexec_options=None):

        """
        Preprocess data for and generate run script
        """

        os.chdir(self.exec_dir)

        # Adaptation should go in the future

        for d in self.domains:
            if hasattr(d, 'adaptation'):
                if d.adaptation:
                    adaptation(d.adaptation, saturne_script, self.case_dir)

        # Before creating or generating file, create stage 'marker' file.

        self.update_scripts_tmp('prepared', 'preprocessing')

        # Determine execution environment
        #================================
        # (priority for n_procs, in increasing order:
        # resource manager, method argument).

        n_procs_default = None
        if len(self.domains) == 1:
            d = self.domains[0]
            if hasattr(d, 'case_n_procs'):
                n_procs_default = int(d.case_n_procs)

        exec_env = cs_exec_environment.exec_environment(self.package_compute,
                                                        self.exec_dir,
                                                        n_procs,
                                                        n_procs_default,
                                                        n_threads)

        # Set user MPI options if required.

        if mpiexec_options != None:
            exec_env.mpi_env.mpiexec_opts = mpiexec_options

        # Compute number of processors.

        n_procs_tot = self.distribute_procs(exec_env.resources.n_procs)

        if n_procs_tot > 1:
            exec_env.resources.n_procs = n_procs_tot

        self.print_procs_distribution()

        # Preprocessing
        #==============

        sys.stdout.write('Preprocessing calculation\n'
                         '-------------------------\n\n')
        sys.stdout.flush()

        self.summary_init(exec_env)

        for d in self.domains:
            d.preprocess()
            if len(d.error) > 0:
                self.error = d.error

        if self.exec_solver:
            s_path = self.generate_solver_script(exec_env)

        # Rename temporary file to indicate new status

        if len(self.error) == 0:
            status = 'ready'
        else:
            status = 'failed'

        self.update_scripts_tmp('preprocessing', status, self.error)

        # Standard or error exit

        if len(self.error) > 0:
            if self.error == 'preprocess':
                error_stage = 'preprocessing'
            else:
                error_stage = self.error
            err_str = ' Error in ' + error_stage + ' stage.\n\n'
            sys.stderr.write(err_str)
            return 1
        else:
            return 0

    #---------------------------------------------------------------------------

    def run_solver(self):

        """
        Run solver proper.
        """

        if not self.exec_solver:
            return 0

        os.chdir(self.exec_dir)

        # Indicate status using temporary file for SALOME.

        self.update_scripts_tmp('ready', 'running')

        sys.stdout.write('Starting calculation\n'
                         '--------------------\n\n')
        sys.stdout.flush()

        # Maximum remaining time if defined through run configuration
        # or batch system.

        b = cs_exec_environment.batch_info()
        max_time = b.get_remaining_time()
        if self.time_limit:
            if max_time:
                if self.time_limit < max_time:
                    max_time = self.time_limit
            else:
                max_time = self.time_limit
        if max_time != None:
            os.putenv('CS_MAXTIME', str(max_time))

        # Tell the script it is being called through the main script
        # (implying environment modules are set and the environment loaded)

        rcfile = cs_exec_environment.get_rcfile(self.package_compute)

        if rcfile or self.package_compute.config.env_modules != "no":
            os.putenv('CS_ENVIRONMENT_SET', 'true')

        # Now run the calculation

        s_path = [self.solver_script_path()]

        retcode = cs_exec_environment.run_command(s_path)

        # Update error codes

        name = self.module_name

        if retcode != 0:
            self.error = 'solver'
            err_str = \
                ' solver script exited with status ' \
                + str(retcode) + '.\n\n'
            sys.stderr.write(err_str)

        if self.error == 'solver':

            if len(self.domains) == 1:
                err_str = \
                    'Error running the calculation.\n\n' \
                    'Check run_solver.log and error* files for details.\n\n'

            else:
                err_str = \
                    'Error running the coupled calculation.\n\n' \
                    'Either of the coupled codes may have failed.\n\n' \
                    'Check the following log files for details.\n\n'

            for d in self.domains:
                err_str += \
                    'Domain ' + str(d.name) + ' (' + d.code_name + '):\n'
                if d.code_name in ('code_saturne', 'neptune_cfd'):
                    err_str += '  run_solver.log, error*.\n\n'
                elif d.code_name == 'SYRTHES':
                    err_str += '  syrthes.log / listing\n\n'
                elif d.code_name == 'Python script':
                    err_str += '  python.log / listing\n\n'
                else:
                    err_str += '  available log files\n\n'

            sys.stderr.write(err_str)

            # Update error status for domains

            for d in self.domains:
                d.error = self.error

        # Indicate status using temporary file for SALOME.

        if retcode == 0:
            status = 'finished'
        else:
            status = 'failed'

        self.update_scripts_tmp(('running',), status, self.error)

        return retcode

    #---------------------------------------------------------------------------

    def save_results(self):

        """
        Save calcultation results from execution directory to result
        directory.
        """

        os.chdir(self.exec_dir)

        self.update_scripts_tmp(('ready', 'finished'), 'saving')

        # Now save results

        sys.stdout.write('Post-calculation operations\n'
                         '---------------------------\n\n')
        sys.stdout.flush()

        self.summary_finalize()
        self.copy_log('summary')

        n_domains = len(self.domains)
        if n_domains > 1 and self.error == '':
            dir_files = os.listdir(self.exec_dir)
            for f in [self.package.runsolver,
                      'mpmd_configfile', 'mpmd_exec.sh']:
                if f in dir_files:
                    try:
                        os.remove(f)
                    except Exception:
                        pass

        e_caption = None
        for d in self.domains:
            d.copy_results()
            if d.error:
                e_caption = d.error

        # Remove directories if empty

        try:
            os.removedirs(self.exec_dir)
            msg = ' Cleaned working directory:\n' \
                + '   ' +  str(self.exec_dir) + '\n'
            sys.stdout.write(msg)
        except Exception:
            pass

        if e_caption: # in case of error during copy/cleanup/postprocessing
            self.update_scripts_tmp('failed', 'failed', e_caption)

        self.update_scripts_tmp('saving', None)

    #---------------------------------------------------------------------------

    def run(self,
            n_procs = None,
            n_threads = None,
            time_limit = None,
            mpiexec_options=None,
            scratchdir = None,
            run_id = None,
            force_id = False,
            stages = None,
            notebook_args=None,
            parametric_args=None,
            kw_args=None):

        """
        Main script.
        """

        # Define run stages if not provided

        if not stages:
            stages = {'prepare_data':True,
                      'initialize':True,
                      'run_solver':True,
                      'save_results':True}

        # If preparation stage is not requested, it must have been done
        # previously, and the id must be forced.

        if not stages['prepare_data']:
            force_id = True

        # Run id and associated directories

        self.set_run_id(run_id)
        self.set_result_dir(force_id)

        # If preparation stage is missing, force it
        if stages['initialize'] and not stages['prepare_data']:
            self.define_exec_dir()
            if not os.path.isdir(self.exec_dir):
                stages['prepare_data'] = True

        # Define scratch directory
        # priority: argument, environment variable, preference setting.

        if scratchdir is None:
            scratchdir = os.getenv('CS_SCRATCHDIR')

        if scratchdir is None:

            # Read the possible config files

            config = configparser.ConfigParser()
            config.read(self.package.get_configfiles())

            # Determine default execution directory if not forced;
            # If the case is already in a sub-directory of the execution
            # directory determined in the configuration file, run in place.

            if config.has_option('run', 'scratchdir'):
                scratchdir = os.path.expanduser(config.get('run', 'scratchdir'))
                scratchdir = os.path.realpath(os.path.expandvars(scratchdir))
                if os.path.realpath(self.result_dir).find(scratchdir) == 0:
                    scratchdir = None

        if scratchdir != None:
            self.exec_prefix = os.path.join(scratchdir, self.package.scratchdir)

        # Define MPI execution options
        # priority: argument, environment variable, preference setting, defaults.

        if mpiexec_options is None:
            mpiexec_options = os.getenv('CS_MPIEXEC_OPTIONS')

        # Set working directory
        # (nonlocal filesystem, reachable by all the processes)

        self.set_exec_dir(force_id)

        # Optional time limit

        if time_limit != None:
            try:
                time_limit = int(time_limit)
                if time_limit >= 0:
                    self.time_limit = time_limit
            except Exception:
                pass

        # Greeting message.

        msg = \
            '\n' \
            + '                      ' + self.module_name + '\n' \
            + '                      ============\n' \
            + '\n' \
            + 'Version:   ' + self.package.version + '\n' \
            + 'Path:      ' + self.package.get_dir('exec_prefix') + '\n'
        if self.package_compute != self.package:
            msg += '  compute: ' + self.package_compute.get_dir('exec_prefix') + '\n\n'
        else:
            msg += '\n'
        msg += 'Result directory:\n' \
               + '  ' +  str(self.result_dir) + '\n\n'

        if self.exec_dir != self.result_dir:
            msg += 'Working directory (to be periodically cleaned):\n' \
                + '  ' +  str(self.exec_dir) + '\n\n'
            self.add_exec_dir_stamp()

        sys.stdout.write(msg)

        # Log user-defined notebook and keyword arguments
        # and pass them to domain structures

        if notebook_args:
            msg = 'Notebook variables:\n'
            keys = list(notebook_args.keys())
            keys.sort()
            for k in keys:
                msg += "  '" + k + "': '" + notebook_args[k] + "'\n"
            msg += '\n'
            sys.stdout.write(msg)

        if parametric_args:
            msg = 'cs_parametric filter arguments:\n'
            msg += '  ' + str(parametric_args) + '\n\n'
            sys.stdout.write(msg)

        if kw_args:
            msg = 'Additional user keyword arguments:\n'
            msg += '  ' + str(kw_args) + '\n\n'
            sys.stdout.write(msg)

        if notebook_args:
            for d in self.domains:
                if not hasattr(d, "notebook"):
                    continue
                if d.notebook is None:
                    d.notebook = notebook_args
                else:
                    err_str = 'domain: + ' + str(d.name) \
                              + ' notebook already defined\n'
                    sys.stderr.write(err_str)

        if parametric_args:
            for d in self.domains:
                if not hasattr(d, "parametric_args"):
                    continue
                if d.parametric_args is None:
                    d.parametric_args = parametric_args
                else:
                    err_str = 'domain: + ' + str(d.name) \
                              + ' parametric_args already defined\n'
                    sys.stderr.write(err_str)

        if kw_args:
            for d in self.domains:
                if d.kw_args is None:
                    d.kw_args = kw_args
                else:
                    err_str = 'domain: + ' + str(d.name) \
                              + ' kw_args already defined\n'
                    sys.stderr.write(err_str)

        sys.stdout.flush()

        # Remove possible status from previous run

        self.update_scripts_tmp(('failed,',
                                 'exceeded_time_limit'), None)

        # Now run

        if self.run_prologue:
            cwd = os.getcwd()
            os.chdir(self.exec_dir)
            retcode = cs_exec_environment.run_command(self.run_prologue)
            os.chdir(cwd)

        try:
            retcode = 0
            if stages['prepare_data']:
                retcode = self.prepare_data(force_id)
            else:
                self.init_prepared_data()

            if stages['initialize'] and  retcode == 0:
                retcode = self.preprocess(n_procs,
                                          n_threads,
                                          mpiexec_options)
            if stages['run_solver'] == True and retcode == 0:
                self.run_solver()

            if stages['save_results'] == True:
                self.save_results()
                self.clear_exec_dir_stamp()

        finally:
            if self.run_id != None:
                self.update_scripts_tmp(('preparing',
                                         'ready',
                                         'running',
                                         'finished'), None)

        # Standard or error exit

        retcode = 0

        if len(self.error) > 0:
            check_stage = {'preprocess':'preprocessing',
                           'solver':'calculation'}
            if self.error in check_stage:
                error_stage = check_stage[self.error]
            else:
                error_stage = self.error
            err_str = 'Run failed in ' + error_stage + ' stage.\n\n'
            if len(self.error_long) > 0:
                err_str += self.error_long + '\n\n'
            sys.stderr.write(err_str)
            retcode = 1

        if self.run_epilogue:
            cwd = os.getcwd()
            os.chdir(self.exec_dir)
            retcode = cs_exec_environment.run_command(self.run_epilogue)
            os.chdir(cwd)

        return retcode

    #---------------------------------------------------------------------------

    def suggest_id(self,
                   run_id_prefix = None,
                   run_id_suffix = None):

        """
        Suggest run id.
        """

        now = datetime.datetime.now()
        run_id_base = now.strftime('%Y%m%d-%H%M')

        if self.dest_dir != None:
            base_dir = os.path.join(self.dest_dir,
                                    os.path.basename(self.case_dir))
        else:
            base_dir = self.case_dir

        r = os.path.join(base_dir, 'RESU')

        if len(self.domains) > 1:
            r += '_COUPLING'
        elif len(self.domains) == 1 and not (run_id_prefix or run_id_suffix):
            try:
                d = self.domains[0]
                if not d.exec_solver:
                    run_id_prefix = 'import_'
                a = d.solver_args
                if a and a in ('--quality', '--preprocess', '--benchmark'):
                    run_id_prefix = a[2:] + '_'
            except Exception:
                pass

        # Make sure directory does not exist already

        j = 0

        while True:

            run_id = run_id_base
            if j > 0:
                run_id += '_' + str(j)

            if run_id_prefix:
                run_id = run_id_prefix + run_id
            if run_id_suffix:
                run_id = run_id + run_id_suffix

            result_dir = os.path.join(r, run_id)

            if os.path.isdir(result_dir):
                j += 1
            else:
                break

        return run_id

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
