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
import os
import os.path
import sys
import stat

import cs_config
import cs_config_build
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
        cmd += ' -help'
        run_command(cmd)
        sys.exit(0)
    else:
        homard_options=' -v'
        cmd += ' -Saturne_Script ' + saturne_script + ' ' + case_dir
        cmd += ' -Pilotage_Adaptation ' + adaptation + homard_options

    if run_command(cmd) != 0:
        sys.exit(0)

#===============================================================================
# Main class
#===============================================================================

class case:
    """Base class from which classes handling running case should inherit."""

    #---------------------------------------------------------------------------

    def __init__(self,
                 case_dir,
                 domains,
                 syr_domains = None,
                 results_by_suffix = True,
                 exec_preprocess = True,
                 exec_partition = True,
                 exec_solver = True):

        # Names, directories, and files in case structure

        (self.study_dir, self.name) = os.path.split(case_dir)

        self.case_dir = case_dir

        self.data_dir = os.path.join(self.case_dir, 'DATA')
        self.result_dir = os.path.join(self.case_dir, 'RESU')
        self.src_dir = os.path.join(self.case_dir, 'SRC')
        self.script_dir = os.path.join(self.case_dir, 'SCRIPTS')

        self.results_by_suffix = results_by_suffix
        self.results_suffix = None

        self.mesh_dir = os.path.join(self.study_dir, 'MESH')

        # Associate case domains and set case directory

        if type(domains) == tuple or  type(domains) == list:
            self.domains = domains
        else:
            self.domains = (domains,)

        if len(self.domains) == 1:
            self.domains[0].set_case_dir(self.case_dir)
        else:
            d_tag = 0
            for d in self.domains:
                d_tag += 1
                d.set_tag(d_tag)
                d.set_case_dir(self.case_dir)

        # Syrthes coupling

        if syr_domains == None:
            self.syr_domains = ()
        elif type(syr_domains) == tuple or  type(syr_domains) == list:
                self.syr_domains = syr_domains
        else:
            self.syr_domains = (syr_domains,)

        if len(self.syr_domains) == 1:
            self.syr_domains[0].set_case_dir(self.case_dir)
        else:
            d_tag = 0
            for sd in self.syr_domains:
                d_tag += 1
                sd.set_tag(d_tag)
                sd.set_case_dir(self.case_dir)

        # Working directory

        self.exec_dir = None

        # Execution and debugging options

        self.exec_preprocess = exec_preprocess
        self.exec_partition = exec_partition
        self.exec_solver = exec_solver

        self.exec_prefix = None

        # Date or other name

        self.date = None
        self.suffix = None

        # Error reporting
        self.error = ''

    #---------------------------------------------------------------------------

    def print_procs_distribution(self):

        """
        Print info on the processor distribution.
        """

        # Print process info

        if len(self.domains) == 1:
            if self.domains[0].n_procs > 1:
                msg = ' Parallel Code_Saturne on ' \
                    + str(self.domains[0].n_procs) + ' processes.\n'
            else:
                msg = ' Single processor Code_Saturne simulation.\n'
            sys.stdout.write(msg)
        else:
            for d in self.domains:
                msg = ' Code_Saturne domain ' + str(d.tag) + ' on ' \
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
                msg = ' SYRTHES domain ' + str(d.tag) + ' on ' \
                    + str(d.n_procs) + ' processes.\n'
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

        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':
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
            ' Error:\n' \
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
                / (n_procs_tot - n_fixed_procs)

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
                              / (len(np_list) - n_fixed_apps)))

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

        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':
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

    def define_exec_dir(self, exec_basename):
        """
        Create execution directory.
        """

        if self.exec_prefix != None:
            names = os.path.split(self.case_dir)
            exec_dir_name = os.path.split(names[0])[1] + '.' + names[1]
            exec_dir_name += '.' + exec_basename
            self.exec_dir = os.path.join(self.exec_prefix, exec_dir_name)
        else:
            self.exec_dir = os.path.join(self.case_dir, 'RESU', exec_basename)

    #---------------------------------------------------------------------------

    def set_exec_dir(self, path):
        """
        Set execution directory.
        """

        # Set execution directory

        self.exec_dir = path

        for d in self.domains:
            d.set_exec_dir(self.exec_dir)
        for d in self.syr_domains:
            d.set_exec_dir(self.exec_dir)

    #---------------------------------------------------------------------------

    def mk_exec_dir(self, exec_basename):
        """
        Create execution directory.
        """

        # Define execution directory name

        self.define_exec_dir(exec_basename)

        # Check that directory does not exist (unless it is also
        # the results directory, which may already have been created,
        # but should be empty at this stage.

        if self.result_dir != self.exec_dir:

            if os.path.isdir(self.exec_dir):
                err_str = \
                    '\nWorking directory: ' + self.exec_dir \
                    + ' already exists.\n' \
                    + 'Calculation will not be run.\n'
                raise RunCaseError(err_str)

            else:
                os.makedirs(self.exec_dir)

        elif (  len(os.listdir(self.exec_dir))
              > len(self.domains) + len(self.syr_domains)):
            err_str = \
                '\nWorking/results directory: ' + self.exec_dir \
                + ' not empty.\n' \
                + 'Calculation will not be run.\n'
            raise RunCaseError(err_str)

        # Set execution directory

        self.set_exec_dir(self.exec_dir)

    #---------------------------------------------------------------------------

    def set_result_dir(self, name):
        """
        If suffix = true, add suffix to all names in result dir.
        Otherwise, create subdirectory
        """

        if self.results_by_suffix == True:
            self.results_suffix = name

        else:
            self.result_dir = os.path.join(self.case_dir, 'RESU', name)

        if not os.path.isdir(self.result_dir):
            os.mkdir(self.result_dir)

        for d in self.domains:
            d.set_result_dir(name, self.results_by_suffix)

        for d in self.syr_domains:
            d.set_result_dir(name, self.results_by_suffix)

    #---------------------------------------------------------------------------

    def summary_init(self, exec_env):

        """
        Build summary start.
        """

        s_path = os.path.join(self.exec_dir, 'summary')
        s = open(s_path, 'w')

        preprocessor = os.path.join(cs_config.dirs.bindir, 'cs_preprocess')
        partitioner = os.path.join(cs_config.dirs.bindir, 'cs_partition')
        if not os.path.isfile(partitioner):
            partitioner = ''
        solver = os.path.join(self.exec_dir, 'cs_solver')
        if not os.path.isfile(solver):
            solver = os.path.join(cs_config.dirs.bindir, 'cs_solver')

        r = exec_env.resources

        date = (datetime.datetime.now()).strftime("%A %B %d %H:%M:%S CEST %Y")
        t_uname = os.uname()
        s_uname = ''
        for t in t_uname:
            s_uname = s_uname + t + ' '
        n_procs = r.n_procs
        if n_procs == None:
            n_procs = 1

        dhline = '========================================================\n'
        hline =  '  ----------------------------------------------------\n'

        s.write(dhline)
        s.write('  Start time       : ' + date + '\n')
        s.write(hline)
        s.write('    Solver          : ' + cs_config.dirs.exec_prefix + '\n')
        s.write('    Preprocessor    : ' + preprocessor + '\n')
        s.write('    Partitioner     : ' + partitioner + '\n')
        s.write(hline)
        s.write('    SYRTHES         : ' + cs_config.dirs.syrthes_prefix + '\n')
        if homard_prefix != None:
            s.write('    HOMARD          : ' + homard_prefix + '\n')
        s.write(hline)
        s.write('    MPI path        : ' + cs_config.mpi_lib.bindir + '\n')
        if len(cs_config.mpi_lib.type) > 0:
            s.write('    MPI type        : ' + cs_config.mpi_lib.type + '\n')
        s.write(hline)
        s.write('    User            : ' + exec_env.user  + '\n')
        s.write(dhline)
        s.write('    Machine         : ' + s_uname  + '\n')
        s.write('    N Procs         : ' + str(n_procs)  + '\n')
        if r.manager == None and r.hosts_list != None:
            s.write('    Processors      :')
            for p in r.hosts_list:
                s.write(' ' + p)
            s.write('\n')
        s.write(dhline)
        s.write('    Case dir.       : ' + self.name + '\n')
        s.write('    Exec. dir.      : ' + self.exec_dir + '\n')
        if self.exec_solver:
            s.write('    Executable      : ' + solver + '\n')
            s.write(hline)

        s.close()

    #---------------------------------------------------------------------------

    def summary_finalize(self):

        """
        Output summary body.
        """

        date = (datetime.datetime.now()).strftime("%A %B %d %H:%M:%S CEST %Y")

        if self.exec_preprocess:
            preprocess = "yes"
            if self.error == 'preprocess':
                preprocess = "failed"
        else:
            preprocess = "no"

        if self.exec_partition:
            partition = "yes"
            if self.error == 'partition':
                partition = "failed"
        else:
            partition = "no"

        if self.exec_solver:
            solver = "yes"
            if self.error == 'solver':
                solver = "failed"
        else:
            solver = "no"

        dhline = '========================================================\n'
        hline =  '  ----------------------------------------------------\n'

        s_path = os.path.join(self.exec_dir, 'summary')
        s = open(s_path, 'a')

        s.write('    Preprocessing   : ' + preprocess + '\n')
        s.write('    Partitioning    : ' + partition + '\n')
        s.write('    Calculation     : ' + solver + '\n')

        s.write(hline)
        s.write('  Finish time      : ' + date + '\n')
        s.write(dhline)

        s.close()

    #---------------------------------------------------------------------------

    def copy_result(self, name):
        """
        Retrieve result from the execution directory
        """

        if os.path.isabs(name):
            src = name
        else:
            src = os.path.join(self.exec_dir, name)

        dest = os.path.join(self.result_dir, os.path.basename(name))

        if (self.results_suffix != None):
            dest += '.' + self.results_suffix

        # Copy single file

        if os.path.isfile(src) and src != dest:
            shutil.copy2(src, dest)

    #---------------------------------------------------------------------------

    def generate_partition_script(self, d, exec_env):
        """
        Generate partitioning script.
        """

        # If n_procs not already given by environment, determine it

        n_procs = exec_env.resources.n_procs
        mpi_env = exec_env.mpi_env

        if n_procs != None:
            if n_procs < d.partition_n_procs:
                ' Error:\n' \
                    '   The current partitioning scheme requires at least ' \
                    + str(n_procs_min) + 'processes,\n' \
                    + '   while the execution environment provides only ' \
                    + str(n_procs) + '.\n' \
                    + '   You may either allocate more processes or try to\n' \
                    + '   oversubscribe by forcing the number of processes\n' \
                    + '   in the toplevel script.'
                raise RunCaseError(err_str)

        # Initialize simple command script

        s_path = os.path.join(d.exec_dir, 'partition.sh')

        s = open(s_path, 'w')

        s.write('#!/bin/sh\n\n')

        # Add MPI directories to PATH if in nonstandard path

        s.write('# Export paths here if necessary or recommended.\n')
        if len(cs_config.mpi_lib.bindir) > 0:
            s.write('export PATH='+ cs_config.mpi_lib.bindir + ':$PATH\n')
        if len(cs_config.mpi_lib.libdir) > 0:
            s.write('export LD_LIBRARY_PATH='+ cs_config.mpi_lib.libdir \
                        + ':$LD_LIBRARY_PATH\n')
        s.write('\n')

        # Boot MPI daemons if necessary

        if mpi_env.gen_hostsfile != None:
            s.write('# Generate hostsfile.\n')
            s.write(mpi_env.gen_hostsfile + ' || exit $?\n\n')

        if n_procs > 1 and mpi_env.mpiboot != None:
            s.write('# Boot MPI daemons.\n')
            s.write(mpi_env.mpiboot + ' || exit $?\n\n')

        # Start assembling command

        mpi_cmd = ''
        mpi_cmd_exe = ''
        mpi_cmd_args = ''
        if d.partition_n_procs > 1 and mpi_env.mpiexec != None:
            mpi_cmd = mpi_env.mpiexec
            if mpi_env.mpiexec_n != None:
                mpi_cmd += mpi_env.mpiexec_n + str(n_procs)
            mpi_cmd += ' '
            if mpi_env.mpiexec_exe != None:
                mpi_cmd += mpi_env.mpiexec_exe + ' '
            if mpi_env.mpiexec_args != None:
                mpi_cmd_args = mpi_env.mpiexec_args + ' '

        p_args = d.partitioner_args()

        s.write('cd ' + p_args[0] + '\n\n')
        s.write('# Run partitioner.\n')
        s.write(mpi_cmd + p_args[1] + mpi_cmd_args + p_args[2] + ' $@\n')

        # Obtain return value (or sum thereof)

        s.write('\nCS_RET=$?\n')

        # Halt MPI daemons if necessary

        if n_procs > 1 and mpi_env.mpihalt != None:
            s.write('\n# Halt MPI daemons.\n')
            s.write(mpi_env.mpihalt + '\n\n')

        if mpi_env.del_hostsfile != None:
            s.write('# Remove hostsfile.\n')
            s.write(mpi_env.del_hostsfile + '\n\n')

        s.write('\nexit $CS_RET\n\n')
        s.close()

        oldmode = (os.stat(s_path)).st_mode
        newmode = oldmode | (stat.S_IXUSR)
        os.chmod(s_path, newmode)

        return s_path

    #---------------------------------------------------------------------------

    def solver_script_path(self):
        """
        Determine name of solver script file.
        """
        return os.path.join(self.exec_dir, 'run_solver.sh')

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_mpiexec(self, n_procs, mpi_env):
        """
        Generate MPMD mpiexec command.
        """

        cmd = ''

        app_id = 0

        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':
                s_args = d.solver_args(app_id=app_id)
                if len(cmd) > 0:
                    cmd += ' : '
                cmd += '-n ' + str(d.n_procs) + ' -wdir ' + s_args[0] \
                    + ' ' + s_args[1] + s_args[2]
                app_id += 1

        for d in self.domains:
            s_args = d.solver_args(app_id=app_id)
            if len(cmd) > 0:
                cmd += ' : '
            cmd += '-n ' + str(d.n_procs) + ' -wdir ' + s_args[0] \
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

        e.write('# MPMD configuration file for mpiexec\n')

        app_id = 0

        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':
                s_args = d.solver_args(app_id=app_id)
                cmd = '-n ' + str(d.n_procs) + ' -wdir ' + s_args[0] \
                    + ' ' + s_args[1] + s_args[2] + '\n'
                e.write(cmd)
                app_id += 1

        for d in self.domains:
            s_args = d.solver_args(app_id=app_id)
            cmd = '-n ' + str(d.n_procs) + ' -wdir ' + s_args[0] \
                + ' ' + s_args[1] + s_args[2] + '\n'
            e.write(cmd)
            app_id += 1

        e.close()

        return e_path

    #---------------------------------------------------------------------------

    def generate_solver_mpmd_script(self, n_procs, mpi_env):
        """
        Generate MPMD dispatch file.
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_exec.sh')
        e = open(e_path, 'w')

        e.write('#!/bin/sh\n\n')
        e.write('# Make sure to transmit possible additional '
                + 'arguments assigned by mpirun to\n'
                + '# the executable with some MPI-1 implementations:\n'
                + '# we use $@ to forward arguments passed to this script'
                + ' to the executable files.\n\n')

        e.write('MPI_RANK=`'
                + cs_config.dirs.pkgdatadir
                + '/runcase_mpi_rank $@`\n')

        app_id = 0
        nr = 0
        test_pf = 'if [ $MPI_RANK -lt '
        test_sf = ' ] ; then\n'

        for d in self.syr_domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_args(app_id=app_id)
            e.write('  cd ' + s_args[0] + '\n')
            e.write('  ' + s_args[1] + s_args[2] + ' $@ > listsyr 2>&1\n')
            if app_id == 0:
                test_pf = 'el' + test_pf
            app_id += 1

        for d in self.domains:
            nr += d.n_procs
            e.write(test_pf + str(nr) + test_sf)
            s_args = d.solver_args(app_id=app_id)
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

    def generate_solver_execve_launcher(self, n_procs, mpi_env):
        """
        Generate execve-based launcher (for machines such as Blue Gene/L).
        """

        e_path = os.path.join(self.exec_dir, 'mpmd_exec.c')
        o_path = os.path.join(self.exec_dir, 'mpmd_exec')
        e = open(e_path, 'w')

        e.write('#include <stdio.h>\n'
                '#include <stdlib.h>\n'
                '#include <unistd.h>\n\n'
                'int\n'
                'main(const int argc, char *const argv[], char *const envp[])\n'
                '{\n')

        e.write('  /* Get MPI rank before MPI_Init() is called (as MPI_Init\n'
                '     may be called only once, by the called program). */\n\n'
                '  int rank = -1;\n\n')

        e.write('#if defined(__blrts__)  /* IBM Blue Gene/L, pid() = rank */\n'
                '  rank = (int)getpid();\n'
                '#endif\n\n')

        e.write('  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL) /* OpenMPI */\n'
                '    rank = atoi(getenv("OMPI_MCA_ns_nds_vpid"));\n'
                '  else if (getenv("OMPI_COMM_WORLD_RANK") != NULL)\n'
                '    rank = atoi(getenv("OMPI_COMM_WORLD_RANK"));\n'
                '  else if (getenv("SLURM_PROCID") != NULL) /* under SLURM */\n'
                '    rank = atoi(getenv("SLURM_PROCID"));\n'
                '  else if (getenv("PMI_RANK") != NULL) /* MPICH 2 */\n'
                '    rank = atoi(getenv("PMI_RANK"));\n'
                '  else if (getenv("LAMRANK") != NULL) /* LAM-MPI */\n'
                '    rank = atoi(getenv("LAMRANK"));\n'
                '  else if (getenv("GMPI_ID") != NULL) /* MPICH-GM */\n'
                '    rank = atoi(getenv("GMPI_ID"));\n\n')

        e.write('  if (rank == -1)\n'
                '    return EXIT_FAILURE;\n\n')

        e.write('  /* Launch appropriate executable */\n\n')

        app_id = 0
        nr = 0

        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':

                nr += d.n_procs
                s_args = self.domains[0].solver_args()
                arg0 = os.path.split(s_args[1])[1]

                e.write('  else if (rank < ' + str(nr) + ') {\n')
                e.write('    const char *filename = "' + s_args[0] + '";\n')

                a_s = '    char *const argv[] = {"' + arg0 + '"'
                for arg in s_args[2].split(' ')[1:]:
                    a_s += ', "' + arg + '"'
                a_s += ', (char *)NULL};\n'
                e.write(a_s)
                e.write('    FILE *fp;\n')
                e.write('    chdir("' + s_args[0] + '");\n')
                e.write('    freopen("listsyr", "w", stdout);\n'
                        '    dup2(fileno(fp), fileno(stderr));\n')
                e.write('    execve(filename, argv, envp);\n'
                        '  }\n')

                app_id += 1

        for d in self.domains:

            nr += d.n_procs
            s_args = self.domains[0].solver_args()
            arg0 = os.path.split(s_args[1])[1]

            e.write('  else if (rank < ' + str(nr) + ') {\n')
            e.write('    const char *filename = "' + s_args[0] + '";\n')

            a_s = '    char *const argv[] = {"' + arg0 + '"'
            for arg in s_args[2].split(' ')[1:]:
                a_s += ', "' + arg + '"'
            a_s += ', (char *)NULL};\n'
            e.write(a_s)
            e.write('    chdir("' + s_args[0] + '");\n')
            e.write('    execve(filename, argv, envp);\n'
                    '  }\n')

            app_id += 1

        e.write('\n'
                '  /* We should not arrive here. */\n'
                '  return EXIT_FAILURE;\n'
                '}\n')

        e.close()

        cmd = cs_config_build.build.cc + ' -o ' + o_path + ' -g ' + e_path

        sys.stdout.write('\n'
                         'Generating MPMD launcher:\n\n')

        retcode = cs_exec_environment.run_command(cmd, echo = True)

        sys.stdout.write('\n')

        if retcode != 0:
            raise RunCaseError(' Error generating MPMD launcher.\n')

        return o_path

    #---------------------------------------------------------------------------

    def generate_solver_script(self, exec_env):
        """
        Generate localexec file.
        """

        # If n_procs not already given by environment, determine it

        n_procs = exec_env.resources.n_procs
        mpi_env = exec_env.mpi_env

        if n_procs == None:
            n_procs = 0
            for d in self.syr_domains:
                if d.coupling_mode == 'MPI':
                    n_procs += d.n_procs

        n_mpi_syr = 0
        for d in self.syr_domains:
            if d.coupling_mode == 'MPI':
                n_mpi_syr += 1

        # Determine if an MPMD syntax (mpiexec variant) will be used

        mpiexec_mpmd = False
        if len(self.domains) > 1:
            if mpi_env.mpmd & cs_exec_environment.MPI_MPMD_mpiexec:
                mpiexec_mpmd = True
            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_configfile:
                mpiexec_mpmd = True

        # Avoid mpiexec variant with SYRTHES as stdout must be redirected;

        if n_mpi_syr > 0:
            mpiexec_mpmd = False
            mpi_env.unset_mpmd_mode(cs_exec_environment.MPI_MPMD_mpiexec)
            mpi_env.unset_mpmd_mode(cs_exec_environment.MPI_MPMD_configfile)

        # Initialize simple solver command script

        s_path = self.solver_script_path()
        s = open(s_path, 'w')

        s.write('#!/bin/sh\n\n')

        # Add detection and handling of SALOME YACS module if run from
        # this environment.

        yacs_test = \
"""
# Detect and handle running under SALOME YACS module.
YACS_ARG=
if test "$SALOME_CONTAINERNAME" != "" -a "$CFDRUN_ROOT_DIR" != "" ; then
  YACS_ARG="--yacs-module=${CFDRUN_ROOT_DIR}"/lib/salome/libCFD_RunExelib.so
fi
"""
        s.write(yacs_test + '\n')

        # Add MPI directories to PATH if in nonstandard path

        s.write('# Export paths here if necessary or recommended.\n')
        if len(cs_config.mpi_lib.bindir) > 0:
            s.write('export PATH='+ cs_config.mpi_lib.bindir + ':$PATH\n')
        if len(cs_config.mpi_lib.libdir) > 0:
            s.write('export LD_LIBRARY_PATH='+ cs_config.mpi_lib.libdir \
                        + ':$LD_LIBRARY_PATH\n')
        s.write('\n')

        # Boot MPI daemons if necessary

        if mpi_env.gen_hostsfile != None:
            s.write('# Generate hostsfile.\n')
            s.write(mpi_env.gen_hostsfile + ' || exit $?\n\n')

        if n_procs > 1 and mpi_env.mpiboot != None:
            s.write('# Boot MPI daemons.\n')
            s.write(mpi_env.mpiboot + ' || exit $?\n\n')

        # Start assembling command

        mpi_cmd = ''
        mpi_cmd_exe = ''
        mpi_cmd_args = ''
        if n_procs > 1 and mpi_env.mpiexec != None:
            mpi_cmd = mpi_env.mpiexec
            if mpiexec_mpmd == False:
                if mpi_env.mpiexec_n != None:
                    mpi_cmd += mpi_env.mpiexec_n + str(n_procs)
                mpi_cmd += ' '
                if mpi_env.mpiexec_exe != None:
                    mpi_cmd += mpi_env.mpiexec_exe + ' '
                if mpi_env.mpiexec_args != None:
                    mpi_cmd_args = mpi_env.mpiexec_args + ' '
            else:
                mpi_cmd += ' '

        # Case with only one cs_solver instance possibly under MPI

        if len(self.domains) == 1 and n_mpi_syr == 0:

            if len(self.syr_domains) == 0:

                s_args = self.domains[0].solver_args()

                s.write('cd ' + s_args[0] + '\n\n')
                s.write('# Run solver.\n')
                s.write(mpi_cmd + s_args[1] + mpi_cmd_args + s_args[2]
                        + ' $YACS_ARGS' + ' $@\n')

            else: # coupling through sockets

                s.write('CS_PORT=35623\n')

                syr_id = 0
                for d in self.syr_domains:
                    s_args = d.solver_args(host_port='localhost:$CS_PORT')
                    s.write('cd ' + s_args[0] + '\n')
                    s.write(s_args[1] + s_args[2] + ' $@ > '
                            + s_args[0] + '/listsyr 2>&1 &\n')
                    s.write('SYR_PID' + str(syr_id) + '=$!\n')
                    syr_id += 1

                s_args = self.domains[0].solver_args(syr_port='$CS_PORT')
                s.write('cd ' + s_args[0] + '\n')
                s.write(mpi_cmd + s_args[1] + mpi_cmd_args + s_args[2] + ' $@\n')
                s.write('CS_PID=$!\n')

                s.write('wait $CS_PID ; CS_RET=$?\n')

                syr_id = 0
                for d in self.syr_domains:
                    s.write('wait SYR_PID' + str(syr_id)
                            + 'SYR_RET=$? ; ((CS_RET=CS_RET+SYR_RET))\n')
                    syr_id += 1

        # General case

        else:

            if mpi_env.mpmd & cs_exec_environment.MPI_MPMD_mpiexec:

                e_path = self.generate_solver_mpmd_mpiexec(n_procs,
                                                           mpi_env)

            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_configfile:

                e_path = self.generate_solver_mpmd_configfile(n_procs,
                                                              mpi_env)
                e_path = '-configfile ' + e_path

            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_script:
                e_path = self.generate_solver_mpmd_script(n_procs, mpi_env)

            elif mpi_env.mpmd & cs_exec_environment.MPI_MPMD_execve:
                e_path = self.generate_solver_execve_launcher(n_procs, mpi_env)

            else:
                raise RunCaseError(' No allowed MPI MPMD mode defined.\n')

            s.write(mpi_cmd + e_path + '\n')

        # Obtain return value (or sum thereof)

        s.write('\nCS_RET=$?\n')

        # Halt MPI daemons if necessary

        if n_procs > 1 and mpi_env.mpihalt != None:
            s.write('\n# Halt MPI daemons.\n')
            s.write(mpi_env.mpihalt + '\n\n')

        if mpi_env.del_hostsfile != None:
            s.write('# Remove hostsfile.\n')
            s.write(mpi_env.del_hostsfile + '\n\n')

        s.write('\nexit $CS_RET\n\n')
        s.close()

        oldmode = (os.stat(s_path)).st_mode
        newmode = oldmode | (stat.S_IXUSR)
        os.chmod(s_path, newmode)

        return s_path

    #---------------------------------------------------------------------------

    def update_scripts_tmp(self, src, dest):

        """
        Create a stamp file in the scripts directory, rename it, or destroy it.
        """

        # Create a temporary file for SALOME (equivalent to "ficstp")

        src_tmp_name = None
        dest_tmp_name = None

        if src != None:
            if type(src) == tuple:
                for s in src:
                    p = os.path.join(self.script_dir, s + '.' + self.suffix)
                    if os.path.isfile(p):
                        src_tmp_name = p

            else:
                p = os.path.join(self.script_dir, src + '.' + self.suffix)
                if os.path.isfile(p):
                    src_tmp_name = p

        try:
            if dest != None:
                dest_tmp_name = os.path.join(self.script_dir,
                                             dest + '.' + self.suffix)
                if src_tmp_name == None:
                    scripts_tmp = open(dest_tmp_name, 'w')
                    scripts_tmp.write(self.exec_dir + '\n')
                    scripts_tmp.close()
                else:
                    os.rename(src_tmp_name, dest_tmp_name)

            else:
                os.remove(src_tmp_name)

        except Exception:
            pass

    #---------------------------------------------------------------------------

    def prepare_data(self,
                     n_procs = None,
                     hosts_list = None,
                     mpi_environment = None,
                     suffix = None):

        """
        Prepare data for calculation.
        """

        # General values

        now = datetime.datetime.now()

        self.date = now.strftime('%m%d%H%M')

        if suffix == None:
            self.suffix = self.date
        else:
            self.suffix = suffix

        for d in self.domains:

            d.exec_preprocess = self.exec_preprocess
            d.exec_partition = self.exec_partition
            d.exec_solver = self.exec_solver

            if d.adaptation:
                adaptation(d.adaptation, saturne_script, self.case_dir)

        if len(self.syr_domains) > 0:
            coupling_mode = self.syr_domains[0].coupling_mode
            for d in self.syr_domains:
                if d.coupling_mode != coupling_mode:
                    err_str = 'This script can only handle SYRTHES couplings ' \
                        + 'using the same coupling mode.\n'
                    raise RunCaseError(err_str)
            for d in self.domains:
                if d.solcom != False:
                    err_str = 'SYRTHES coupling is not compatible ' \
                        + 'with SOLCOM-type meshes.'
                    raise RunCaseError(err_str)

        # Now that all domains are defined, set result copy mode

        self.set_result_dir(self.suffix)

        # Before creating or generating file, create stage 'marker' file.

        self.update_scripts_tmp(None, 'preparing')

        # Create working directory (reachable by all the processors)

        self.mk_exec_dir(self.suffix)

        # Copy script before changing to the working directory
        # (as after that, the relative path will not be up to date).

        self.copy_result(sys.argv[0])

        os.chdir(self.exec_dir)

        # Determine execution environment.

        exec_env = cs_exec_environment.exec_environment(self.exec_dir,
                                                        n_procs,
                                                        hosts_list)

        # Set user MPI environment if required.

        if mpi_environment != None:
            exec_env.mpi_env = mpi_environment

        # Compute number of processors

        n_procs_tot = self.distribute_procs(exec_env.resources.n_procs)

        if n_procs_tot > 1:
            exec_env.resources.n_procs = n_procs_tot

        # Greeting message.

        msg = \
            '\n' \
            + '                      Code_Saturne is running\n' \
            + '                      ***********************\n' \
            + '\n' \
            + ' Working directory (to be periodically cleaned):\n' \
            + '   ' +  str(self.exec_dir) + '\n\n' \
            + ' Kernel version:  ' + cs_config.package.version + '\n' \
            + ' Preprocessor:    ' + os.path.join(cs_config.dirs.bindir,
                                                  'cs_preprocess') + '\n\n'
        sys.stdout.write(msg)

        self.print_procs_distribution()

        # Compile user subroutines if necessary.

        if self.exec_solver == True:

            need_compile = False

            for d in self.domains:
                if d.check_model_consistency() != 0:
                    raise RunCaseError('Incompatible model options.')
                if d.needs_compile() == True:
                    need_compile = True

            if need_compile == True or len(self.syr_domains) > 0:

                msg = \
                    " ****************************************\n" \
                    "  Compiling user subroutines and linking\n" \
                    " ****************************************\n\n"
                sys.stdout.write(msg)

                for d in self.domains:
                    d.compile_and_link()

                for d in self.syr_domains:
                    d.compile_and_link()

        # Setup data
        #===========

        sys.stdout.write('\n'
                         ' ****************************\n'
                         '  Preparing calculation data\n'
                         ' ****************************\n\n')

        if self.exec_preprocess:
            for d in self.domains:
                d.copy_preprocessor_data()
        else:
            for d in self.domains:
                d.copy_preprocessor_output_data()

        if self.exec_solver:

            for d in self.domains:
                d.copy_solver_data()

            for d in self.syr_domains:
                d.copy_data()

        sys.stdout.write('\n'
                         ' ***************************\n'
                         '  Preprocessing calculation\n'
                         ' ***************************\n\n')

        self.summary_init(exec_env)

        if self.exec_preprocess:
            for d in self.domains:
                d.run_preprocessor()
                if len(d.error) > 0:
                    self.error = d.error

        if self.exec_partition:
            for d in self.domains:
                d.check_partitioner()
                if d.partition_n_procs > 1:
                    p_path = self.generate_partition_script(d, exec_env)
                else:
                    d.run_partitioner()
                if len(d.error) > 0:
                    self.error = d.error

        s_path = self.generate_solver_script(exec_env)

        # Rename temporary file to indicate new status

        if len(self.error) == 0:
            status = 'ready'
        else:
            status = 'failed'

        self.update_scripts_tmp('preparing', status)

        # Standard or error exit

        if len(self.error) > 0:
            stage = {'preprocess':'preprocessing',
                     'partition':'partitioning'}
            if self.error in stage:
                error_stage = stage[self.error]
            else:
                error_stage = self.error
            err_str = ' Error in ' + error_stage + ' stage.\n\n'
            sys.stderr.write(err_str)
            return 1
        else:
            return 0

    #---------------------------------------------------------------------------

    def run_partitioner(self, suffix = None):

        """
        Run partitioner.
        """

        if not self.exec_partition:
            return 0

        # Check if parallel partitioning is required

        n_partitionings = 0

        for d in self.domains:
            d.check_partitioner()
            if d.partition_n_procs > 1:
                n_partitionings += 1

        if n_partitionings == 0:
            return 0

        # Now prepare for parallel partitionings

        retcode = 0

        if suffix != None:
            self.suffix = suffix

        if self.exec_dir == None:
            self.define_exec_dir(self.suffix)
            self.set_exec_dir(self.exec_dir)

        os.chdir(self.exec_dir)

        # Indicate status using temporary file for SALOME.

        self.update_scripts_tmp('ready', 'partitioning')

        sys.stdout.write('\n'
                         ' **********************\n'
                         '  Parallel partitioning\n'
                         ' **********************\n\n')

        # Now run the partitionings

        for d in self.domains:

            d.check_partitioner()
            if d.partition_n_procs > 1:

                s_path = os.path.join(d.exec_dir, 'partition.sh')
                retcode = cs_exec_environment.run_command(s_path)

                if retcode != 0:

                    err_str = \
                        'Partitioner script:\n' + s_path \
                        + '\nexited with status ' + str(retcode) + '\n' \
                        + 'Check the partition.log file for details.\n\n'
                    sys.stderr.write(err_str)

                    # Update error codes

                    d.error = 'partition'
                    self.error = d.error

                    break

        # Indicate status using temporary file for SALOME.

        if len(self.error) == 0:
            status = 'ready'
        else:
            status = 'failed'

        self.update_scripts_tmp('partitioning', 'status')

        return retcode

    #---------------------------------------------------------------------------

    def run_solver(self, suffix = None):

        """
        Run solver proper.
        """

        if not self.exec_solver:
            return 0

        if suffix != None:
            self.suffix = suffix

        if self.exec_dir == None:
            self.define_exec_dir(self.suffix)
            self.set_exec_dir(self.exec_dir)

        os.chdir(self.exec_dir)

        # Indicate status using temporary file for SALOME.

        if self.domains[0].logging_args == '--log 0':
            self.update_scripts_tmp('ready', 'runningstd')
        else:
            self.update_scripts_tmp('ready', 'runningext')

        sys.stdout.write('\n'
                         ' **********************\n'
                         '  Starting calculation\n'
                         ' **********************\n\n')

        # Maximum remaining time for PBS or similar batch system.

        b = cs_exec_environment.batch_info()
        max_time = b.get_remaining_time()
        if max_time != None:
            os.putenv('CS_MAXTIME', max_time)

        # Now run the calculation

        s_path = self.solver_script_path()

        retcode = cs_exec_environment.run_command(s_path)

        # Update error codes

        if retcode != 0:
            self.error = 'solver'
            err_str = \
                'Code_Saturne solver script exited with status ' \
                + str(retcode) + '.\n\n'
            sys.stderr.write(err_str)

        if self.error == 'solver':
            if len(self.syr_domains) > 0:
                err_str = \
                    'Error running the coupled calculation.\n\n' \
                    'Either Code_Saturne or SYRTHES may have failed.\n\n' \
                    'Check Code_Saturne log (listing) and SYRTHES log (listsyr)\n' \
                    'for details, as well as error* files.\n\n'
            else:
                err_str = \
                    'Error running the calculation.\n\n' \
                    'Check Code_Saturne log (listing) and error* files for details.\n\n'
            sys.stderr.write(err_str)

        # Indicate status using temporary file for SALOME.

        if retcode == 0:
            status = 'finished'
        else:
            status = 'failed'

        self.update_scripts_tmp(('runningstd', 'runningext'), status)

        return retcode

    #---------------------------------------------------------------------------

    def save_results(self,
                     suffix = None):

        """
        Save calcultation results from execution directory to result
        directory.
        """

        if suffix != None:
            self.suffix = suffix

        if self.exec_dir == None:
            self.define_exec_dir(self.suffix)
            self.set_exec_dir(self.exec_dir)

        self.set_result_dir(self.suffix)

        os.chdir(self.exec_dir)

        self.update_scripts_tmp(('ready', 'failed', 'finished'), 'saving')

        # Now save results

        sys.stdout.write('\n'
                         ' ****************************\n'
                         '  Saving calculation results\n'
                         ' ****************************\n\n')

        if self.exec_preprocess:
            for d in self.domains:
                d.copy_preprocessor_results()

        if self.exec_partition:
            for d in self.domains:
                d.copy_partition_results()

        if self.exec_solver:
            for d in self.domains:
                d.copy_solver_results(self.date)

            for d in self.syr_domains:
                d.copy_results()

        self.summary_finalize()

        self.copy_result('summary')

        # Remove the Salome temporary file

        self.update_scripts_tmp('saving', None)

    #---------------------------------------------------------------------------

    def run(self,
            n_procs = None,
            hosts_list = None,
            mpi_environment = None,
            exec_prefix = None,
            suffix = None,
            prepare_data = True,
            run_solver = True,
            save_results = True):

        """
        Main script.
        """

        tmpdir = None

        if exec_prefix != None:
            tmpdir = exec_prefix
        else:
            # Read the possible config files

            if sys.platform == 'win32' or sys.platform == 'win64':
                username = os.getenv('USERNAME')
            else:
                username = os.getenv('USER')

            config = ConfigParser.ConfigParser({'user':username})
            config.read([os.path.expanduser('~/.code_saturne.cfg'),
                         os.path.join(cs_config.dirs.sysconfdir,
                                      'code_saturne.cfg')])

            if config.has_option('run', 'tmpdir'):
                tmpdir = os.path.expanduser(config.get('run', 'tmpdir'))
                tmpdir = os.path.expandvars(tmpdir)

        if tmpdir != None:
            self.exec_prefix = os.path.join(tmpdir, 'tmp_Saturne')

        if suffix != None:
            self.suffix = suffix

        try:
            retcode = 0
            if prepare_data == True:
                retcode = self.prepare_data(n_procs,
                                            hosts_list,
                                            mpi_environment)
                if retcode == 0:
                    self.run_partitioner()
            if run_solver == True:
                self.run_solver()

            if save_results == True:
                self.save_results()

        finally:
            if self.suffix != None:
                self.update_scripts_tmp(('preparing',
                                         'ready',
                                         'runningstd',
                                         'runningext',
                                         'finished',
                                         'failed'), None)

        # Standard or error exit

        if len(self.error) > 0:
            stage = {'preprocess':'preprocessing',
                     'partition':'partitioning',
                     'solver':'calculation'}
            if self.error in stage:
                error_stage = stage[self.error]
            else:
                error_stage = self.error
            err_str = ' Error in ' + error_stage + ' stage.\n\n'
            sys.stderr.write(err_str)
            return 1
        else:
            return 0

    #---------------------------------------------------------------------------

