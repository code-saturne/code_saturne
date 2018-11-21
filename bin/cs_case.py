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

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3
import datetime
import os
import os.path
import platform
import sys
import stat

import cs_exec_environment

from cs_case_domain import *

homard_prefix = None

#===============================================================================
# Utility functions
#===============================================================================

def adaptation(adaptation, saturne_script, case_dir):

    """
    Call HOMARD adaptation.
    """

    cmd = os.path.join(homard_prefix, 'saturne_homard')

    if adaptation == '-help':
        cmd.append('-help')
        cs_exec_environment.run_command(cmd)
        sys.exit(0)
    else:
        homard_options = '-v'
        cmd = cmd + ['-Saturne_Script', saturne_script, case_dir]
        cmd = cmd + ['-Pilotage_Adaptation', adaptation, homard_options]

    if cs_exec_environment.run_command(cmd) != 0:
        sys.exit(0)

#-------------------------------------------------------------------------------

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
                 staging_dir = None,
                 coupling_parameters = None,
                 domains = None,
                 syr_domains = None,
                 ast_domain = None,
                 fsi_coupler = None,
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

        # Ensure we have tuples or lists to simplify later tests

        if type(domains) == list:
            self.domains = domains
        elif type(domains) == tuple:
            self.domains = list(domains)
        else:
            self.domains = [domains]

        if syr_domains == None:
            self.syr_domains = []
        elif type(syr_domains) == tuple :
            self.syr_domains = list(syr_domains)
        elif type(syr_domains) == list:
            self.syr_domains = syr_domains
        else:
            self.syr_domains = [syr_domains]

        if ast_domain == None:
            self.ast_domains = []
        else:
            self.ast_domains = [ast_domain]

        self.fsi_coupler = fsi_coupler

        if py_domains == None:
            self.py_domains = []
        elif type(py_domains) == tuple:
            self.py_domains = list(py_domains)
        elif type(py_domains) == list:
            self.py_domains = py_domains
        else:
            self.py_domains = [py_domains]

        # Mark fluid domains as coupled with Code_Aster if coupling is present.
        # This should be improved by migrating tests to the executable,
        # to better handle cases were only some domains are coupled.
        # For now, we assume that the first CFD domain is the coupled one.

        if fsi_coupler:
            self.domains[0].fsi_aster = True

        # Check names in case of multiple domains (coupled by MPI)

        n_domains = len(self.domains) + len(self.syr_domains) + len(self.py_domains)

        if n_domains > 1:
            err_str = 'In case of multiple domains (i.e. code coupling), ' \
                + 'each domain must have a name.\n'
            for d in (self.domains + self.syr_domains + self.py_domains):
                if d.name == None:
                    raise RunCaseError(err_str)

        # Names, directories, and files in case structure:
        # test if we are using a master runcase or a single runcase,
        # associate case domains and set case directory

        self.case_dir = case_dir

        self.script_dir = os.path.join(self.case_dir, 'SCRIPTS')
        self.result_dir = None

        if (os.path.isdir(os.path.join(case_dir, 'DATA')) and n_domains == 1):
            # Simple domain case from standard directory structure
            (self.study_dir, self.name) = os.path.split(case_dir)
            if len(self.domains) == 1:
                self.domains[0].name = None
            elif len(self.syr_domains) == 1:
                self.syr_domains[0].name = None
            elif len(self.py_domains) == 1:
                self.py_domains[0].name = None

        else:
            # Coupling or single domain run from coupling script
            self.study_dir = case_dir
            self.name = os.path.split(case_dir)[1]
            self.script_dir = self.case_dir

        if staging_dir:
            staging_dir = check_exec_dir_stamp(staging_dir)

        for d in ( self.domains + self.syr_domains \
                 + self.ast_domains + self.py_domains ):
            d.set_case_dir(self.case_dir, staging_dir)

        # Working directory

        self.exec_dir = None
        self.exec_prefix = None

        self.exec_solver = True

        n_exec_solver = 0
        for d in ( self.domains + self.syr_domains \
                 + self.ast_domains + self.py_domains ):
            if d.exec_solver:
                n_exec_solver += 1

        if n_exec_solver == 0:
            self.exec_solver = False
        elif n_exec_solver <   len(self.domains) + len(self.syr_domains) \
                             + len(self.ast_domains) + len(self.py_domains):
            err_str = 'In case of multiple domains (i.e. code coupling), ' \
                + 'all or no domains must execute their solver.\n'
            raise RunCaseError(err_str)

        # Coupling parameters if present

        self.coupling_parameters = coupling_parameters

        # Date or other name

        self.run_id = None

        # Error reporting
        self.error = ''

    #---------------------------------------------------------------------------

    def print_procs_distribution(self):

        """
        Print info on the processor distribution.
        """

        # Print process info

        name = self.package.code_name

        if len(self.domains) == 1:
            if self.domains[0].n_procs > 1:
                msg = ' Parallel ' + name + ' on ' \
                    + str(self.domains[0].n_procs) + ' processes.\n'
            else:
                msg = ' Single processor ' + name + ' simulation.\n'
            sys.stdout.write(msg)
        else:
            for d in self.domains:
                msg = ' ' + name + ' domain ' + d.name + ' on ' \
                    + str(d.n_procs) + ' process(es).\n'
                sys.stdout.write(msg)

        if len(self.syr_domains) == 1:
            if self.syr_domains[0].n_procs > 1:
                msg = ' Parallel SYRTHES on ' \
                    + str(self.syr_domains[0].n_procs) + ' processes.\n'
            else:
                msg = ' Single processor SYRTHES simulation.\n'
            sys.stdout.write(msg)
        else:
            for d in self.syr_domains:
                msg = ' SYRTHES domain ' + d.name + ' on ' \
                    + str(d.n_procs) + ' processes.\n'
                sys.stdout.write(msg)

        if len(self.ast_domains) == 1:
            if self.ast_domains[0].n_procs > 1:
                msg = ' Parallel Code_Aster on ' \
                    + str(self.ast_domains[0].n_procs) + ' processes.\n'
            else:
                msg = ' Single processor Code_Aster simulation.\n'
            sys.stdout.write(msg)
        else:
            for d in self.ast_domains:
                msg = ' Code_Aster domain ' + d.name + ' on ' \
                    + str(d.n_procs) + ' processes.\n'
                sys.stdout.write(msg)

        if len(self.py_domains) == 1:
            if self.py_domains[0].n_procs > 1:
                msg = ' Parallel Python script on ' \
                    + str(self.py_domains[0].n_procs) + ' processes.\n'
            else:
                msg = ' Single processor Python script simulation.\n'
            sys.stdout.write(msg)
        else:
            for d in self.py_domains:
                msg = ' Python script domain ' + d.name + ' on ' \
                    + str(d.n_procs) + ' processes.\n'
                sys.stdout(msg)

        sys.stdout.write('\n')

    #---------------------------------------------------------------------------

    def distribute_procs(self, n_procs=None):

        """
        Distribute number of processors to domains,
        and returns the total number of associated processes.
        """

        # Initial count

        np_list = []

        for d in (self.domains + self.syr_domains + self.py_domains):
            np_list.append(d.get_n_procs())

        # Code_Aster is handled in a specific manner: processes are counted
        # against the total, but the coupler is not, and processes will not by
        # run by the same mpiexec instance (TODO: check/improve how this
        # interacts with resource managers, so that some nodes are not
        # oversubscribed while others are unused).

        for d in self.ast_domains:
            np_list.append(d.get_n_procs())

        n_procs_tot = 0
        n_procs_min = 0
        for np in np_list:
            n_procs_tot += np[0]
            n_procs_min += np[1]

        # If no process count is given or everything fits:

        if n_procs == None or n_procs == n_procs_tot:
            return n_procs_tot

        # If process count is given and not sufficient, abort.

        elif n_procs < n_procs_min:
            err_str = ' Error:\n' \
                '   The current calculation scheme requires at least ' \
                + str(n_procs_min) + 'processes,\n' \
                + '   while the execution environment provides only ' \
                + str(n_procs) + '.\n' \
                + '   You may either allocate more processes or try to\n' \
                + '   oversubscribe by forcing the number of processes\n' \
                + '   in the toplevel script.'
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
                if np[2] == None:
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

        for d in (self.domains + self.syr_domains + self.py_domains):
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
            r = os.path.join(self.case_dir, 'RESU')
            if len(self.domains) + len(self.syr_domains) \
             + len(self.ast_domains) + len(self.py_domains) > 1:
                r += '_COUPLING'
            self.exec_dir = os.path.join(r, self.run_id)

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

        for d in ( self.domains + self.syr_domains + self.ast_domains \
                 + self.py_domains):
            d.set_exec_dir(self.exec_dir)

    #---------------------------------------------------------------------------

    def define_result_dir(self):

        r = os.path.join(self.case_dir, 'RESU')

        # Coupled case
        if len(self.domains) + len(self.syr_domains) \
         + len(self.ast_domains) + len(self.py_domains) > 1:
            r += '_COUPLING'

        if os.path.isdir(r):
            self.result_dir = os.path.join(r, self.run_id)
        else:
            r = os.path.join(self.case_dir, 'RESU')
            err_str = \
                    '\nResults directory: ' + r + '\n' \
                    + 'does not exist.\n' \
                    + 'Calculation will not be run.\n'
            raise RunCaseError(err_str)

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

        for d in ( self.domains + self.syr_domains + self.ast_domains \
                 + self.py_domains):
            d.set_result_dir(self.run_id, self.result_dir)

    #---------------------------------------------------------------------------

    def summary_init(self, exec_env):

        """
        Build summary start.
        """

        s_path = os.path.join(self.exec_dir, 'summary')
        s = open(s_path, 'w')

        preprocessor = self.package.get_preprocessor()
        solver = os.path.join(self.exec_dir, self.package.solver)
        if not os.path.isfile(solver):
            solver = self.package.get_solver()

        r = exec_env.resources

        date = (datetime.datetime.now()).strftime("%A %B %d %H:%M:%S CEST %Y")
        t_uname = platform.uname()
        s_uname = ''
        for t in t_uname:
            s_uname = s_uname + t + ' '
        n_procs = r.n_procs
        if n_procs == None:
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
        if r.manager == None and r.hosts_list != None:
            s.write('  Hosts          :')
            for p in r.hosts_list:
                s.write(' ' + p)
            s.write('\n')
        s.write(hline)

        if len(self.domains) + len(self.syr_domains) + len(self.py_domains)  > 1:
            s.write('  Exec. dir.     : ' + self.exec_dir + '\n')
            s.write(hline)

        s.write('  Environment variables\n')
        ek = list(os.environ.keys())
        ek.sort()
        for k in ek:
            s.write('    ' + k + '=' + os.environ[k] + '\n')
        s.write(hline)

        for d in ( self.domains + self.syr_domains + self.ast_domains \
                 + self.py_domains):
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

        # Copy single file

        if src:
            dest = os.path.join(self.result_dir, os.path.basename(src))
            if os.path.isfile(src) and src != dest:
                shutil.copy2(src, dest)

    #---------------------------------------------------------------------------

    def solver_script_path(self):
        """
        Determine name of solver script file.
        """
        return os.path.join(self.exec_dir, self.package.runsolver)

    #---------------------------------------------------------------------------

    def generate_fsi_scheme(self):
        """
        "Creating YACS coupling scheme.
        """

        sys.stdout.write('Creating YACS coupling scheme.\n')

        cfd_d = self.domains[0]
        ast_d = self.ast_domains[0]

        ast_dir_name = os.path.basename(ast_d.exec_dir)
        ast_mesh_name = os.path.join(ast_dir_name, ast_d.mesh)
        ast_comm = os.path.join(ast_dir_name, ast_d.comm)

        c = self.fsi_coupler

        pkgdatadir = self.package.get_dir('pkgdatadir')
        s_path = os.path.join(pkgdatadir, 'salome', 'fsi_yacs_scheme.xml')
        s = open(s_path, 'r')
        template = s.read()
        s.close()

        import re

        e_re = re.compile('@AST_WORKINGDIR@')
        template = re.sub(e_re, ast_d.exec_dir, template)

        e_re = re.compile('@COCAS_WORKINGDIR@')
        template = re.sub(e_re, self.exec_dir, template)

        e_re = re.compile('@CFD_WORKINGDIR@')
        template = re.sub(e_re, cfd_d.exec_dir, template)

        e_re = re.compile('@AST_MAIL@')
        template = re.sub(e_re, ast_mesh_name, template)

        e_re = re.compile('@NBPDTM@')
        template = re.sub(e_re, str(c['max_time_steps']), template)

        e_re = re.compile('@NBSSIT@')
        template = re.sub(e_re, str(c['n_sub_iterations']), template)

        e_re = re.compile('@DTREF@')
        template = re.sub(e_re, str(c['time_step']), template)

        e_re = re.compile('@TTINIT@')
        template = re.sub(e_re, str(c['start_time']), template)

        e_re = re.compile('@EPSILO@')
        template = re.sub(e_re, str(c['epsilon']), template)

        e_re = re.compile('@COMM_FNAME@')
        template = re.sub(e_re, ast_comm, template)

        s_path = os.path.join(self.exec_dir, 'fsi_yacs_scheme.xml')
        s = open(s_path, 'w')
        s.write(template)
        s.close()

        # Now generate application

        if self.package.config.have_salome == "no":
            raise RunCaseError("SALOME is not available in this installation.\n")

        template = """\
%(salomeenv)s;
PYTHONPATH=%(pythondir)s/salome:%(pkgpythondir)s${PYTHONPATH:+:$PYTHONPATH};
export PYTHONPATH;
${KERNEL_ROOT_DIR}/bin/salome/appli_gen.py --prefix=%(applidir)s --config=%(configdir)s
"""

        cmd = template % {'salomeenv': self.package.config.salome_env,
                          'pythondir': self.package.get_dir('pythondir'),
                          'pkgpythondir': self.package.get_dir('pkgpythondir'),
                          'applidir': os.path.join(self.exec_dir, 'appli'),
                          'configdir': os.path.join(self.package.get_dir('sysconfdir'),
                                                    'salome',
                                                    'fsi_appli_config.xml')}

        retcode = cs_exec_environment.run_command(cmd,
                                                  stdout=None,
                                                  stderr=None)


    #---------------------------------------------------------------------------

    def generate_yacs_wrappers(self):
        """
        Generate wrappers for YACS.
        """

        apps = [('FSI_SATURNE.exe', os.path.join(self.exec_dir, 'cfd_by_yacs'))]

        i = 0
        for d in self.ast_domains:
            if len(self.ast_domains) == 1:
                apps.append(('FSI_ASTER.exe',
                             os.path.join(d.exec_dir, 'aster_by_yacs')))
            else:
                i = i+1
                apps.append(('FSI_ASTER' + str(i) + '.exe',
                             os.path.join(d.exec_dir, 'aster_by_yacs')))

        bin_dir = os.path.join(self.exec_dir, 'bin')
        if not os.path.isdir(bin_dir):
            os.mkdir(bin_dir)

        for a in apps:

            s_path = os.path.join(bin_dir, a[0])
            s = open(s_path, 'w')

            cs_exec_environment.write_shell_shebang(s)

            s_cmds = \
"""
# Get SALOME variables from launcher

export SALOME_CONTAINER=$1
export SALOME_CONTAINERNAME=$2
export SALOME_INSTANCE=$3

"""
            s.write(s_cmds)

            s.write(a[1] + '\n')

            s.close

            oldmode = (os.stat(s_path)).st_mode
            newmode = oldmode | (stat.S_IXUSR)
            os.chmod(s_path, newmode)

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_mpiexec(self, n_procs, mpi_env):
        """
        Generate MPMD mpiexec command.
        """

        cmd = ''

        app_id = 0

        for d in self.syr_domains:
            s_args = d.solver_command()
            if len(cmd) > 0:
                cmd += ' : '
            cmd += '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + s_args[1] + s_args[2]
            app_id += 1

        for d in self.domains:
            s_args = d.solver_command(app_id=app_id)
            if len(cmd) > 0:
                cmd += ' : '
            cmd += '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + s_args[1] + s_args[2]
            app_id += 1

        for d in self.py_domains:
            s_args = d.solver_command()
            if len(cmd) > 0:
                cmd += ' : '
            cmd += '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + s_args[1] + s_args[2]
            app_id += 1

        return cmd

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile(self, n_procs, mpi_env):
        """
        Generate MPMD mpiexec config file.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')
        e = open(e_path, 'w')

        app_id = 0

        for d in (self.syr_domains + self.domains + self.py_domains):
            s_args = d.solver_command()
            cmd = '-n ' + str(d.n_procs) \
                + ' -wdir ' + os.path.basename(s_args[0]) \
                + ' ' + s_args[1] + s_args[2] + '\n'
            e.write(cmd)
            app_id += 1

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile_srun(self, n_procs, mpi_env):
        """
        Generate MPMD mpiexec config file for SLURM's srun.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')
        e = open(e_path, 'w')

        e.write('# srun multiple program configuration file\n#\n')
        e.write('# > srun -n ' + str(n_procs) \
                + ' --multi-prog mpmd_configfile\n#\n')

        rank_id = 0

        for d in (self.syr_domains + self.domains + self.py_domains):
            s_args = d.solver_command()
            if s_args[1][0:2] == './':
                s_path = os.path.join(s_args[0], s_args[1])
                s_path = './' + os.path.relpath(s_path, self.exec_dir)
            else:
                s_path = s_args[1]
            cmd = '%d-%d\t' % (rank_id, rank_id + d.n_procs - 1) \
                   + s_path + s_args[2] \
                   + ' -wdir ' + os.path.basename(s_args[0]) + '\n'
            e.write(cmd)
            rank_id += d.n_procs

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_configfile_bgq(self, n_procs, mpi_env):
        """
        Generate MPMD mpiexec config file for BG/Q.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_configfile')
        e = open(e_path, 'w')

        rank_id = 0

        for d in (self.syr_domains + self.domains + self.py_domains):
            cmd = '#mpmdbegin %d-%d\n' % (rank_id, rank_id + d.n_procs - 1)
            e.write(cmd)
            s_args = d.solver_command()
            if s_args[1][0:2] == './':
                s_path = os.path.relpath(os.path.join(s_args[0], s_args[1]),
                                         self.exec_dir)
            else:
                s_path = s_args[1]
            cmd = '#mpmdcmd ' + s_path + s_args[2] \
                + ' -w ' + os.path.basename(s_args[0]) + '\n'
            e.write(cmd)
            e.write('#mpmdend\n')
            rank_id += d.n_procs

        e.write('#mapping ABCDET\n')

        # NOTE: adding the mapping line above seems to help in some cases where
        # the mapping file is not interpreted correctly, but seems not to be
        # always required (the documentation is minimalistic);
        # with driver V1R2M0_17, reading of the mapping file seems fragile
        # and subject to random failures unless this line is added.

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_script(self, n_procs, mpi_env):
        """
        Generate MPMD dispatch file.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_exec.sh')
        e = open(e_path, 'w')

        user_shell = cs_exec_environment.get_shell_type()

        e.write('#!' + user_shell + '\n\n')
        e.write('# Make sure to transmit possible additional '
                + 'arguments assigned by mpirun to\n'
                + '# the executable with some MPI-1 implementations:\n'
                + '# we use $@ to forward arguments passed to this script'
                + ' to the executable files.\n\n')

        e.write('MPI_RANK=`'
                + self.package.get_runcase_script('runcase_mpi_rank')
                + ' $@`\n')

        app_id = 0
        nr = 0
        test_pf = 'if [ $MPI_RANK -lt '
        test_sf = ' ] ; then\n'

        for d in self.syr_domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_command()
            e.write('  cd ' + s_args[0] + '\n')
            e.write('  ' + s_args[1] + s_args[2] + ' $@\n')
            if app_id == 0:
                test_pf = 'el' + test_pf
            app_id += 1

        for d in self.domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_command(app_id=app_id)
            e.write('  cd ' + s_args[0] + '\n')
            e.write('  ' + s_args[1] + s_args[2] + ' $@\n')
            if app_id == 0:
                test_pf = 'el' + test_pf
            app_id += 1

        for d in self.py_domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_command()
            e.write('  cd ' + s_args[0] + '\n')
            e.write('  ' + s_args[1] + s_args[2] + ' $@\n')
            if app_id == 0:
                test_pf = 'el' + test_pf
            app_id += 1

        e.write('fi\n'
                + 'CS_RET=$?\n'
                + 'exit $CS_RET\n')

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
        if len(self.domains) + len(self.syr_domains) + len(self.py_domains) > 1:
            if mpi_env.mpmd & cs_exec_environment.MPI_MPMD_mpiexec:
                mpiexec_mpmd = True
            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_configfile:
                mpiexec_mpmd = True

        # Start assembling command

        mpi_cmd = ''
        mpi_cmd_exe = ''
        mpi_cmd_args = ''
        if n_procs > 1 and mpi_env.mpiexec != None:
            mpi_cmd = mpi_env.mpiexec
            if mpi_env.mpiexec_opts != None:
                mpi_cmd += ' ' + mpi_env.mpiexec_opts
            if mpiexec_mpmd == False:
                if mpi_env.mpiexec_n != None:
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

        # Case with only one cs_solver instance possibly under MPI

        if len(self.domains) == 1 and len(self.syr_domains) == 0 \
           and len(self.py_domains) == 0:

            s_args = self.domains[0].solver_command()

            s.write('cd ' + s_args[0] + '\n\n')
            cs_exec_environment.write_script_comment(s, 'Run solver.\n')
            s.write(mpi_cmd + s_args[1] + mpi_cmd_args + s_args[2])
            s.write(' ' + cs_exec_environment.get_script_positional_args() +
                    '\n')

        # General case

        else:

            if mpi_env.mpmd & cs_exec_environment.MPI_MPMD_mpiexec:

                if mpi_env.mpiexec_separator != None:
                    mpi_cmd += mpi_env.mpiexec_separator + ' '

                e_path = self.generate_solver_mpmd_mpiexec(n_procs,
                                                           mpi_env)

            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_configfile:

                if mpi_env.type == 'BGQ_MPI':
                    e_path = self.generate_solver_mpmd_configfile_bgq(n_procs,
                                                                      mpi_env)
                    if mpi_env.mpiexec == 'srun':
                        mpi_cmd += '--launcher-opts=\'--mapping ' + e_path + '\' '
                    else:
                        mpi_cmd += '--mapping ' + e_path + ' '
                    if mpi_env.mpiexec_separator != None:
                        mpi_cmd += mpi_env.mpiexec_separator + ' '
                    mpi_cmd += self.package_compute.get_solver()

                elif mpi_env.mpiexec == 'srun':
                    e_path = self.generate_solver_mpmd_configfile_srun(n_procs,
                                                                       mpi_env)
                    mpi_cmd += '--multi-prog ' + e_path

                else:
                    e_path = self.generate_solver_mpmd_configfile(n_procs,
                                                                  mpi_env)
                    mpi_cmd += '-configfile ' + e_path

                e_path = ''

            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_script:

                if mpi_env.mpiexec_separator != None:
                    mpi_cmd += mpi_env.mpiexec_separator + ' '

                e_path = self.generate_solver_mpmd_script(n_procs, mpi_env)

            else:
                raise RunCaseError(' No allowed MPI MPMD mode defined.\n')

            s.write(mpi_cmd + e_path + '\n')

    #---------------------------------------------------------------------------

    def fsi_script_body(self, s):
        """
        Generate script body for Code_Aster FSI coupling.
        """

        template = """\
%(salomeenv)s;
export PATH=%(exec_dir)s/bin:$PATH
appli=%(exec_dir)s/appli
unset PYTHONHOME
port_log=`pwd`/.salome_port.log
$appli/salome start -t --ns-port-log=$port_log
$appli/salome shell -- driver -e fsi_yacs_scheme.xml
$appli/salome kill `cat $port_log`

"""
        s.write(template % {'salomeenv': self.package.config.salome_env,
                            'exec_dir': self.exec_dir})

    #---------------------------------------------------------------------------

    def generate_fsi_scripts(self, n_procs, mpi_env):
        """
        # Generate scripts for Code_Aster domains and YACS wrappers
        # Script per Code_Aster domain generated
        """

        # Generate simple solver command script

        s_path = os.path.join(self.exec_dir, 'cfd_by_yacs')
        s = open(s_path, 'w')

        cs_exec_environment.write_shell_shebang(s)

        self.solver_script_body(n_procs-1, mpi_env, s)

        cs_exec_environment.write_export_env(s, 'CS_RET',
                                             cs_exec_environment.get_script_return_code())

        if sys.platform.startswith('win'):
            s.write('\nexit %CS_RET%\n')
        else:
            s.write('\nexit $CS_RET\n')
        s.close()

        # Generate additional scripts

        for d in self.ast_domains:
            d.generate_script()

        self.generate_yacs_wrappers()

        oldmode = (os.stat(s_path)).st_mode
        newmode = oldmode | (stat.S_IXUSR)
        os.chmod(s_path, newmode)

        return s_path

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

        if n_procs == None:
            n_procs = 0
            for d in (self.syr_domains + self.py_domains):
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

        # Add path for SALOME in case of coupling with YACS

        if self.ast_domains:
            salome_libdir = os.path.join(self.package_compute.get_dir("libdir"),
                                         'salome')
            cs_exec_environment.write_prepend_path(s,
                                                   'LD_LIBRARY_PATH',
                                                   salome_libdir)

        # Handle rcfile and environment modules if used

        rcfile = cs_exec_environment.get_rcfile(self.package_compute)

        if rcfile or self.package_compute.config.env_modules != "no":
            cs_exec_environment.write_script_comment(s, \
               'Load environment if this script is run directly.\n')
            s.write('if test "$CS_ENVIRONMENT_SET" != "true" ; then\n')
            if self.package_compute.config.env_modules != "no":
                s.write('  module purge\n')
                for m in self.package_compute.config.env_modules.strip().split():
                    s.write('  module load ' + m + '\n')
            if rcfile:
                s.write('  source ' + rcfile + '\n')
            s.write('fi\n\n')

        # Add paths for plugins or dynamic library dependencies

        plugin_lib_dirs, plugin_pythonpath_dirs, plugin_env_vars \
            = self.package_compute.config.get_run_environment_dependencies()
        for d in plugin_lib_dirs:
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

        # Boot MPI daemons if necessary

        if mpi_env.gen_hostsfile != None:
            cs_exec_environment.write_script_comment(s, 'Generate hostsfile.\n')
            s.write(mpi_env.gen_hostsfile + ' || exit $?\n\n')

        if n_procs > 1 and mpi_env.mpiboot != None:
            cs_exec_environment.write_script_comment(s, 'Boot MPI daemons.\n')
            s.write(mpi_env.mpiboot + ' || exit $?\n\n')

        # Generate script body

        if self.ast_domains:
            self.generate_fsi_scheme()
            self.fsi_script_body(s)
            self.generate_fsi_scripts(n_procs, mpi_env)

        else:
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
                if src_tmp_name == None:
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
        p = os.path.join(self.exec_dir, 'status.exceeded_time_limit')

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

        if self.run_id == None:
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

        os.chdir(self.exec_dir)

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
                        " ****************************************\n" \
                        "  Compiling user subroutines and linking\n" \
                        " ****************************************\n\n"
                    sys.stdout.write(msg)
                    sys.stdout.flush()
                d.compile_and_link()
                if len(d.error) > 0:
                    self.error = d.error

        if len(self.error) > 0:
            self.update_scripts_tmp('preparing', 'failed', self.error)
            err_str = ' Error in ' + self.error + ' stage.\n\n'
            sys.stderr.write(err_str)
            return 1

        # Setup data
        #===========

        sys.stdout.write('\n'
                         ' ****************************\n'
                         '  Preparing calculation data\n'
                         ' ****************************\n\n')
        sys.stdout.flush()

        for d in ( self.domains + self.syr_domains + self.ast_domains \
                 + self.py_domains):
            d.prepare_data()
            if len(d.error) > 0:
                self.error = d.error

        # Output coupling parameters for staging

        if self.coupling_parameters:
            s = open(os.path.join(self.result_dir,
                                  'coupling_parameters.py'), 'w')
            s.write(self.coupling_parameters)
            s.close()

        # Rename temporary file to indicate new status

        if len(self.error) == 0:
            status = 'prepared'
        else:
            status = 'failed'

        self.update_scripts_tmp('preparing', status, self.error)

        # Standard or error exit

        if len(self.error) > 0:
            err_str = ' Error in ' + self.error + ' stage.\n\n'
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
            if os.path.isfile(os.path.join(d.exec_dir,
                                           d.package_compute.solver)):
               d.solver_path = os.path.join('.',
                                            d.package_compute.solver)

        for d in self.syr_domains:
            d.solver_path = os.path.join('.', 'syrthes')

        for d in self.py_domains:
            d.solver_path = os.path.join('.', 'python')

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
        if len(self.domains) == 1 and len(self.syr_domains) == 0 \
           and len(self.py_domains) == 0:
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

        # Transfer parameters MPI parameters from user scripts here.

        if len(self.domains) == 1 and len(self.syr_domains) == 0 \
           and len(self.py_domains) == 0:
            d = self.domains[0]
            if d.user_locals:
                m = 'define_mpi_environment'
                if m in d.user_locals.keys():
                    eval(m + '(exec_env.mpi_env)', locals(), d.user_locals)
                    del d.user_locals[m]

        # Compute number of processors.

        n_procs_tot = self.distribute_procs(exec_env.resources.n_procs)

        if n_procs_tot > 1:
            exec_env.resources.n_procs = n_procs_tot

        self.print_procs_distribution()

        # Preprocessing
        #==============

        sys.stdout.write('\n'
                         ' ***************************\n'
                         '  Preprocessing calculation\n'
                         ' ***************************\n\n')
        sys.stdout.flush()

        self.summary_init(exec_env)

        for d in (self.domains + self.syr_domains + self.py_domains):
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

    def run_solver(self, run_id = None):

        """
        Run solver proper.
        """

        if not self.exec_solver:
            return 0

        os.chdir(self.exec_dir)

        # Indicate status using temporary file for SALOME.

        self.update_scripts_tmp('ready', 'running')

        sys.stdout.write('\n'
                         ' **********************\n'
                         '  Starting calculation\n'
                         ' **********************\n\n')
        sys.stdout.flush()

        # Maximum remaining time for PBS or similar batch system.

        b = cs_exec_environment.batch_info()
        max_time = b.get_remaining_time()
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

        name = self.package.code_name

        if retcode != 0:
            self.error = 'solver'
            err_str = \
                ' solver script exited with status ' \
                + str(retcode) + '.\n\n'
            sys.stderr.write(err_str)

        if self.error == 'solver':

            if len(self.syr_domains) > 0:
                err_str = \
                    'Error running the coupled calculation.\n\n' \
                    'Either ' + name + ' or SYRTHES may have failed.\n\n' \
                    'Check ' + name + ' log (listing) and SYRTHES log (syrthes.log)\n' \
                    'for details, as well as error* files.\n\n'
            elif len(self.py_domains) > 0:
                err_str = \
                    'Error running the coupled calculation.\n\n' \
                    'Either ' + name + ' or Python script may have failed.\n\n' \
                    'Check ' + name + ' log (listing) and Python log (python.log)\n' \
                    'for details, as well as error* files.\n\n'
            else:
                err_str = \
                    'Error running the calculation.\n\n' \
                    'Check ' + name + ' log (listing) and error* files for details.\n\n'
            sys.stderr.write(err_str)

            # Update error status for domains

            for d in (self.domains + self.syr_domains + self.py_domains):
                d.error = self.error

        # Indicate status using temporary file for SALOME.

        if retcode == 0:
            status = 'finished'
        else:
            status = 'failed'

        self.update_scripts_tmp(('running',), status, self.error)

        return retcode

    #---------------------------------------------------------------------------

    def save_results(self,
                     run_id = None):

        """
        Save calcultation results from execution directory to result
        directory.
        """

        os.chdir(self.exec_dir)

        self.update_scripts_tmp(('ready', 'finished'), 'saving')

        # Now save results

        sys.stdout.write('\n'
                         ' *****************************\n'
                         '  Post-calculation operations\n'
                         ' *****************************\n\n')
        sys.stdout.flush()

        self.summary_finalize()
        self.copy_log('summary')

        n_domains = len(self.domains) + len(self.syr_domains) + len(self.py_domains)
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
        for d in (self.domains + self.syr_domains + self.py_domains):
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
            mpiexec_options=None,
            scratchdir = None,
            run_id = None,
            force_id = False,
            stages = None):

        """
        Main script.
        """

        # Define run stages if not provided

        if not stages:
            stages = {'prepare_data':True,
                      'initialize':True,
                      'run_solver':True,
                      'save_results':True}

        # Define scratch directory
        # priority: argument, environment variable, preference setting.

        if scratchdir == None:
            scratchdir = os.getenv('CS_SCRATCHDIR')

        if scratchdir == None:

            # Read the possible config files

            if sys.platform.startswith('win'):
                username = os.getenv('USERNAME')
            else:
                username = os.getenv('USER')

            config = configparser.ConfigParser({'user':username})
            config.read(self.package.get_configfiles())

            # Determine default execution directory if not forced;
            # If the case is already in a sub-directory of the execution
            # directory determined in the configuration file, run in place.

            if config.has_option('run', 'scratchdir'):
                scratchdir = os.path.expanduser(config.get('run', 'scratchdir'))
                scratchdir = os.path.realpath(os.path.expandvars(scratchdir))
                if os.path.realpath(self.case_dir).find(scratchdir) == 0:
                    scratchdir = None

        if scratchdir != None:
            self.exec_prefix = os.path.join(scratchdir, self.package.scratchdir)

        # Define MPI execution options
        # priority: argument, environment variable, preference setting, defaults.

        if mpiexec_options == None:
            mpiexec_options = os.getenv('CS_MPIEXEC_OPTIONS')

        # If preparation stage is not requested, it must have been done
        # previously, and the id must be forced.

        if not stages['prepare_data']:
            force_id = True

        # Run id and associated directories

        self.set_run_id(run_id)

        # If preparation stage is missing, force it
        if stages['initialize'] and not stages['prepare_data']:
            self.define_exec_dir()
            if not os.path.isdir(self.exec_dir):
                stages['prepare_data'] = True

        # Set result copy mode

        self.set_result_dir(force_id)

        # Set working directory
        # (nonlocal filesystem, reachable by all the processes)

        self.set_exec_dir(force_id)

        # Greeting message.

        msg = \
            '\n' \
            + '                      ' + self.package.code_name + '\n' \
            + '                      ************\n' \
            + '\n' \
            + ' Version:   ' + self.package.version + '\n' \
            + ' Path:      ' + self.package.get_dir('exec_prefix') + '\n'
        if self.package_compute != self.package:
            msg += '   compute: ' + self.package_compute.get_dir('exec_prefix') + '\n\n'
        else:
            msg += '\n'
        msg += ' Result directory:\n' \
               + '   ' +  str(self.result_dir) + '\n\n'

        if self.exec_dir != self.result_dir:
            msg += ' Working directory (to be periodically cleaned):\n' \
                + '   ' +  str(self.exec_dir) + '\n'
            self.add_exec_dir_stamp()

        msg += '\n'

        sys.stdout.write(msg)
        sys.stdout.flush()

        # Remove possible status from previous run

        self.update_scripts_tmp(('failed,',
                                 'exceeded_time_limit'), None)
        # Now run

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

        if len(self.error) > 0:
            check_stage = {'preprocess':'preprocessing',
                           'solver':'calculation'}
            if self.error in check_stage:
                error_stage = check_stage[self.error]
            else:
                error_stage = self.error
            err_str = ' Error in ' + error_stage + ' stage.\n\n'
            sys.stderr.write(err_str)
            return 1
        else:
            return 0

    #---------------------------------------------------------------------------

    def suggest_id(self,
                   run_id_prefix = None,
                   run_id_suffix = None):

        """
        Suggest run id.
        """

        now = datetime.datetime.now()
        run_id_base = now.strftime('%Y%m%d-%H%M')

        r = os.path.join(self.case_dir, 'RESU')

        if len(self.domains) + len(self.syr_domains) \
         + len(self.ast_domains) + len(self.py_domains) > 1:
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
                run_id = run_id_suffix + run_id

            result_dir = os.path.join(r, run_id)

            if os.path.isdir(result_dir):
                j += 1
            else:
                break

        return run_id

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
