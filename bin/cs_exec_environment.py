#!/usr/bin/env python
#
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009  EDF
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

import datetime
import fnmatch
import os.path
import sys
import tempfile

python_version = sys.version[:3]

try:
    import subprocess        # Python 2.4+
    have_subprocess = True
except ImportError:
    import popen2            # Python 2.3-
    have_subprocess = False

from optparse import OptionParser

from cs_config import mpi_lib

#===============================================================================
# Utility functions
#===============================================================================

def abs_exec_path(path):
    """
    Find an executable in the system path.
    """

    abs_path = None

    if os.path.isabs(path):
        return path

    else:
        try:
            for d in os.getenv('PATH').split():
                f = os.path.join(d, path)
                if os.path.isfile(f):
                    return f
        except Exception:
            pass

    return None

#-------------------------------------------------------------------------------

def run_command(cmd, echo = False, stdout = sys.stdout, stderr = sys.stderr):
    """
    Run a command.
    """
    if echo == True:
        stdout.write(cmd + '\n')

    if have_subprocess == True:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        returncode = p.returncode

    else:
        p = popen2.Popen3(cmd, capturestderr=True)
        returncode = p.wait()
        output = (p.fromchild.read(), p.childerr.read())

    if len(output[0]) > 0:
        stdout.write(output[0])
    if len(output[1]) > 0:
        stderr.write(output[1])

    return returncode

#-------------------------------------------------------------------------------

def get_command_output(cmd):
    """
    Run a command and return it's standard output.
    """
    if have_subprocess == True:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        if p.returncode != 0:
            sys.stderr.write(output[1] + '\n')
            return ''
        else:
            return output[0]

    else:
        p = popen2.Popen3(cmd, capturestderr=True)
        returncode = p.wait()
        if returncode != 0:
            sys.stderr.write(p.childerr.read())
            return ''
        else:
            return p.fromchild.read()

#-------------------------------------------------------------------------------

class batch_info:

    #---------------------------------------------------------------------------

    def __init__(self):

        """
        Get batch system information.
        """

        self.batch_type = None
        self.job_file = None
        self.job_name = None
        self.job_id = None
        self.queue = None

        # Check for specific batch environments

        s = os.getenv('LSB_JOBID') # LSF
        if s != None:
            self.batch_type = 'LSF'
            self.job_file = os.getenv('LSB_JOBFILENAME')
            self.job_name = os.getenv('LSB_JOBNAME')
            self.job_id = os.getenv('LSB_BATCH_JID')
            self.queue = os.getenv('LSB_QUEUE')

        if self.batch_type == None:
            s = os.getenv('PBS_JOBID') # PBS
            if s != None:
                self.batch_type = 'PBS'
                self.job_name = os.getenv('PBS_JOBNAME')
                self.job_id = os.getenv('PBS_JOBID')
                self.queue = os.getenv('PBS_QUEUE')

        if self.batch_type == None:
            s = os.getenv('LOADL_JOB_NAME') # LoadLeveler
            if s != None:
                self.batch_type = 'LOADL'

        if self.batch_type == None:
            s = os.getenv('SGE_TASKID') # Sun Grid Engine
            if s != None:
                self.batch_type = 'SGE'
                self.job_name = os.getenv('JOB_NAME')
                self.job_id = os.getenv('JOB_ID')
                self.queue = os.getenv('QUEUE')

    #---------------------------------------------------------------------------

    def get_remaining_time(self):

        """
        Get remaining time if available from batch system.
        """

        rtime = None

        if self.batch_type == 'PBS':
            cmd = "qstat -r $PBS_JOBID | grep $PBS_JOBID" \
                + " | sed -e's/ \{1,\}/ /g' | cut -d ' ' -f 9"
            rtime = get_command_output(cmd)

        return rtime

#-------------------------------------------------------------------------------

class resource_info(batch_info):

    #---------------------------------------------------------------------------

    def __init__(self, n_procs = None, hosts_list = None):

        """
        Get execution resources information.
        """

        batch_info.__init__(self)

        self.manager = None
        self.n_procs = None
        self.n_nodes = None

        # If obtained from an environment variable, express
        # the hosts file using a shell variable rather than
        # an absolute name (for use in generated scripts).

        self.hosts_file = None
        self.hosts_list = None

        # Check for resource manager and eventual hostsfile

        s = os.getenv('SLURM_NPROCS')
        if s != None:
            self.manager = 'SLURM'
            self.n_procs = int(s)
            s = os.getenv('SLURM_NNODES')
            if s != None:
                self.n_nodes = int(s)

        if self.manager == None and self.batch_type == 'LSF':
            s = os.getenv('LSB_HOSTS')
            self.manager = 'LSF'
            if s != None:
                self.hosts_list = s.split(' ')

        if self.manager == None and self.batch_type == 'LOADL':
            s = os.getenv('LOADL_PROCESSOR_LIST')
            self.manager = 'LOADL'
            if s != None:
                self.hosts_list = s.split(' ')

        if self.manager == None and self.batch_type == 'PBS':
            s = os.getenv('PBS_NODEFILE')
            if s != None:
                self.manager = 'PBS'
                self.n_procs_from_hosts_file(s)
                self.hosts_file = '$PBS_NODEFILE'

        if self.manager == None and self.batch_type == 'SGE':
            s = os.getenv('NSLOTS')
            if s != None:
                self.n_procs = int(s)
            s = os.getenv('NHOSTS')
            if s != None:
                self.n_nodes = int(s)
            s = os.getenv('TMPDIR')
            if s != None:
                s += '/machines'
                if os.path.isfile(s):
                    s.manager = 'SGE'
                    self.n_procs_from_hosts_file(s)
                    self.hosts_file = '$TMPDIR/machines'

        # Set an optional list of hosts if we are not running under
        # a resource manager.

        if hosts_list != None:
            if self.manager != None:
                sys.stderr.write('Warning:\n'
                                 + '   Host list will be ignored because a'
                                 + ' resource manager (' + self.manager
                                 + ') is in use.\n\n')
            else:
                self.resources.hosts_list = hosts_list

        # Determine number of processors from hosts file or list

        if self.n_procs == None:
            if self.hosts_file != None:
                self.n_procs = 0
                f = open(self.hosts_file, 'r')
                for line in f:
                    self.n_procs += 1
                f.close()
            elif self.hosts_list != None:
                self.n_procs = len(self.hosts_list)

        # Check and possibly set number of processes

        if n_procs != None:
            if self.n_procs != None:
                if self.n_procs != n_procs:
                    sys.stderr.write('Warning:\n'
                                     +'   Will use ' + str(self.n_procs)
                                     + ' processes while resource manager ('
                                     + self.resources + ')\n   allows for '
                                     + str(self.resources.n_procs) + '.\n\n')
            self.n_procs = n_procs

    #---------------------------------------------------------------------------

    def n_procs_from_hosts_file(self, hosts_file):

        """
        Compute number of hosts from a hostsfile.
        """

        self.n_procs = 0
        f = open(hosts_file, 'r')
        for line in f:
            self.n_procs += 1
        f.close()

    #---------------------------------------------------------------------------

    def get_hosts_file(self, wdir = None):
        """
        Returns the name of the hostsfile associated with the
        resource manager. A hostsfile is built from a hosts
        list if necessary.
        """

        hosts_file = self.hosts_file

        if self.hosts_file == None:

            if self.hosts_list != None:
                if wdir != None:
                    hosts_file = os.path.join(wdir, 'hostsfile')
                else:
                    hosts_file = 'hostsfile'
                f = open(hosts_file, 'w')
                # If number of procs not specified, determine it by list
                if self.n_procs == None or self.n_procs < 1:
                    n_procs = 0
                    for host in self.hosts_list:
                        f.write(host + '\n')
                        n_procs += 1
                # If the number of procs is known, use only beginning of list
                # if it contains more entries than the required number, or
                # loop through list as many time as necessary to reach
                # prescribed proc count
                else:
                    proc_count = 0
                    while proc_count < self.n_procs:
                        for host in self.hosts_list:
                            if proc_count < self.n_procs:
                                f.write(host + '\n')
                                proc_count += 1
                f.close()

        return hosts_file

#-------------------------------------------------------------------------------
# MPI environments and associated commands
#-------------------------------------------------------------------------------

MPI_MPMD_mpiexec = (1<<0)
MPI_MPMD_script =  (1<<1)
MPI_MPMD_execve =  (1<<2)

class mpi_environment:

    def __init__(self, resource_info=None):
        """
        Returns MPI environment info.
        """

        # Note that self.mpiexec will usually be ' -n ' if present;
        # blanks are used to separate from the surrounding arguments,
        # but in the case of srun, which uses a -n<n_procs> instead
        # of -n <n_procs> syntax, setting it to ' -n' will be enough.

        self.type = mpi_lib.type
        self.bindir = mpi_lib.bindir

        self.gen_hostsfile = None
        self.del_hostsfile = None
        self.mpiboot = None
        self.mpihalt = None
        self.mpiexec = None
        self.mpiexec_n = None
        self.mpiexec_exe = None
        self.mpiexec_args = None
        self.mpmd = MPI_MPMD_mpiexec

        self.info_cmds = None

        # Initialize based on known MPI types, or default.

        init_method = self.__init_other__

        if len(mpi_lib.type) > 0:
            mpi_env_by_type = {'MPICH2':self.__init_mpich2__,
                               'MPICH1':self.__init_mpich1__,
                               'OpenMPI':self.__init_openmpi__,
                               'LAM_MPI':self.__init_lam__,
                               'BGL_MPI':self.__init_bgl__,
                               'BGP_MPI':self.__init_bgp__,
                               'HP_MPI':self.__init_hp_mpi__,
                               'MPIBULL2':self.__init_mpibull2__}
            if mpi_lib.type in mpi_env_by_type:
                init_method = mpi_env_by_type[mpi_lib.type]

        p = os.getenv('PATH').split(':')
        if len(mpi_lib.bindir) > 0:
            p.insert(0, mpi_lib.bindir)

        init_method(p, resource_info)

    #---------------------------------------------------------------------------

    def __init_mpich2__(self, p, resource_info=None):

        """
        Initialize for MPICH2 environment.
        """

        # Determine base executable paths

        launcher_names = ['mpiexec', 'mpirun']

        for name in launcher_names:
            for dir in p:
                absname = os.path.join(dir, name)
                if os.path.isfile(absname):
                    if dir == mpi_lib.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        for dir in p:
            if os.path.isfile(os.path.join(dir, 'mpdboot')):
                if dir == mpi_lib.bindir:
                    self.mpiboot = os.path.join(dir, 'mpdboot')
                    self.mpihalt = os.path.join(dir, 'mpdallexit')
                else:
                    self.mpiboot = 'mpdboot'
                    self.mpihalt = 'mpdallexit'
                break

        # Note that the self.mpiboot and self.mpihalt commands may be
        # too simplistic in certain cases. For example, it would be
        # interesting to run mpdringtest first, and depending on its
        # exit code (0 if ring is up, 255 if not up), run mpdboot.
        # The same, mpdlistjobs could be run before mpdallexit to ensure
        # that we do not kill an mpdring if another job is running.

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        if launcher_base[:7] == 'mpiexec':
            self.mpmd = MPI_MPMD_mpiexec
            self.mpiexec_n = ' -n '
        elif launcher_base[:7] == 'mpirun':
            self.mpiexec_n = ' -np '
            self.mpmd = MPI_MPMD_script

        # Other options to add

        # Resource manager info

        if resource_info != None:
            if resource_info.manager == 'SLURM':
                # This requires linking with SLURM's implementation
                # of the PMI library.
                self.mpiexec = 'srun'
                self.mpiexec_n = ' -n'
                self.mpmd = MPI_MPMD_script
                self.mpiboot = None
                self.mpihalt = None
            elif resource_info.manager == 'PBS':
                # Convert PBS to MPD format (based on MPICH2 documentation)
                # before MPI boot.
                self.gen_hostsfile = 'sort $PBS_NODEFILE | uniq -C ' \
                    + '| awk \'{ printf("%s:%s", $2, $1); }\' > ./mpd.nodes'
                self.del_hostsfile = 'rm -f ./mpd.nodes'
                if self.mpiboot != None:
                    self.mpiboot = pbs_to_mpd \
                        + self.mpiboot + ' --file ./mpd.nodes\n'
            else:
                hostsfile = resource_info.get_hosts_file(self)
                if hostsfile != None and self.mpiboot != None:
                    self.mpiboot += ' --file ' + hostsfile

        # Info commands

        self.info_cmds = ['mpich2version']

    #---------------------------------------------------------------------------

    def __init_mpich1__(self, p, resource_info=None):

        """
        Initialize for MPICH1 environment.
        """

        # Determine base executable paths

        launcher_names = ['mpirun.mpich',
                          'mpirun.mpich-mpd',
                          'mpirun.mpich-shmem',
                          'mpirun']

        for name in launcher_names:
            for dir in p:
                absname = os.path.join(dir, name)
                if os.path.isfile(absname):
                    if dir == mpi_lib.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        # Determine processor count and MPMD handling

        self.mpiexec_n = ' -np '
        self.mpmd = MPI_MPMD_script

        # MPD daemons must be launched using the ch_p4mpd device,
        # but may not be automated as easily as with MPICH2's mpdboot,
        # so it is not handled here, and must be managed explicitely
        # by the user or environment.

        # Other options to add

        # Resource manager info

        if resource_info != None:
            hostsfile = resource_info.get_hosts_file(self)
            if hostsfile != None:
                self.mpiexec += ' -machinefile ' + hostsfile

        # Info commands

        self.info_cmds = ['mpichversion']

    #---------------------------------------------------------------------------

    def __init_openmpi__(self, p, resource_info=None):
        """
        Initialize for OpenMPI environment.
        """

        # Determine base executable paths

        launcher_names = ['mpiexec.openmpi', 'mpirun.openmpi',
                          'mpiexec', 'mpirun']

        for name in launcher_names:
            for dir in p:
                absname = os.path.join(dir, name)
                if os.path.isfile(absname):
                    if dir == mpi_lib.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        self.mpiexec_n = ' -n '
        if launcher_base[:7] == 'mpiexec':
            self.mpmd = MPI_MPMD_mpiexec
        elif launcher_base[:7] == 'mpirun':
            self.mpmd = MPI_MPMD_script

        # Other options to add

        # Resource manager info (SLURM and PBS Pro/Torque are
        # handled automatically by OpenMPI, so there is
        # nothing to do in this case).

        if resource_info != None:
            if not resource_info.manager in ['SLURM', 'PBS']:
                hostsfile = resource_info.get_hosts_file()
                if hostsfile != None:
                    self.mpiexec += ' --machinefile ' + hostsfile

        # Info commands

        self.info_cmds = ['ompi_info -a']

    #---------------------------------------------------------------------------

    def __init_lam__(self, p, resource_info=None):
        """
        Initialize for LAM/MPI environment.
        """

        # Determine base executable paths

        launcher_names = ['mpiexec.lam', 'mpirun.lam', 'mpiexec', 'mpirun']

        for name in launcher_names:
            for dir in p:
                absname = os.path.join(dir, name)
                if os.path.isfile(absname):
                    if dir == mpi_lib.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        for dir in p:
            if os.path.isfile(os.path.join(dir, 'lamboot')):
                if dir == mpi_lib.bindir:
                    self.mpiboot = os.path.join(dir, 'lamboot')
                    self.mpihalt = os.path.join(dir, 'lamhalt')
                else:
                    self.mpiboot = 'lamboot'
                    self.mpihalt = 'lamhalt'
                break

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        if launcher_base[:7] == 'mpiexec':
            self.mpiexec_n = ' -n '
            self.mpmd = MPI_MPMD_mpiexec
        elif launcher_base[:7] == 'mpirun':
            self.mpiexec_n = ' -np '
            self.mpmd = MPI_MPMD_script

        # Determine options to add

        if self.mpiboot != None:
            self.mpiboot += ' -v'
            self.mpihalt += ' -v'

        # Resource manager info

        if resource_info != None:
            hostsfile = resource_info.get_hosts_file(self)
            if hostsfile != None and self.mpiboot != None:
                self.mpiboot += ' ' + hostsfile
                self.mpihalt += ' ' + hostsfile

        # Info commands

        self.info_cmds = ['laminfo -all']

    #---------------------------------------------------------------------------

    def __init_bgl__(self, p, resource_info=None):

        """
        Initialize for Blue Gene/L environment.
        """

        # Set base executable path

        self.mpiexec = '/bgl/BlueLight/ppcfloor/bglsys/mpi/mpirun'

        # Determine processor count and MPMD handling

        self.mpiexec_n = None
        self.mpmd = MPI_MPMD_execve

        # Other options to add

        self.mpiexec_exe = '-exe'
        self.mpiexec_args = '-args'

        # Info commands

    #---------------------------------------------------------------------------

    def __init_bgp__(self, p, resource_info=None):

        """
        Initialize for Blue Gene/P environment.
        """

        # Set base executable path

        self.mpiexec = 'mpiexec'

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        self.mpiexec_n = None
        self.mpmd = MPI_MPMD_mpiexec

        # Other options to add

        # Info commands

    #---------------------------------------------------------------------------

    def __init_hp_mpi__(self, p, resource_info=None):
        """
        Initialize for HP MPI environment.
        """

        # Determine base executable paths

        self.mpiexec = 'mpirun'

        if not os.path.isabs(self.mpiexec):
            for dir in p:
                absname = os.path.join(dir, self.mpiexec)
                if os.path.isfile(absname):
                    if dir == mpi_lib.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = launcher_name
                    break

        # Determine processor count and MPMD handling

        self.mpmd = MPI_MPMD_script
        self.mpiexec_n = '-np'

        # Determine options to add

        # If resource manager is used, add options

        if resource_info != None:
            if resource_info.manager == 'SLURM':
                self.mpiexec += ' -srun'
                self.mpiexec_n = None
            elif resource_info.manager == 'LSF':
                self.mpiexec += ' -lsb_hosts'

        # Info commands

    #---------------------------------------------------------------------------

    def __init_mpibull2__(self, p, resource_info=None):
        """
        Initialize for MPIBULL2 environment.
        """

        self.__init_mpich2__(p, resource_info)

        # On a Bull Novascale machine using MPIBULL2 (based on
        # MPICH2), mpdboot and mpdallexit commands exist, as well as mpiexec
        # (which connects to the mpd ring), but the SLURM srun launcher or
        # mpibull2-launch meta-launcher (which can integrate directly to
        # several resource managers using prun, srun, orterun, or mprun)
        # should be simpler to use.

        # The SLURM configuration is slightly different from that of MPICH2.

        if resource_info != None:
            if resource_info.manager == 'SLURM':
                self.mpiexec = 'srun'
                self.mpmd = MPI_MPMD_script
                self.mpiexec_n = None
                self.mpiboot = None
                self.mpihalt = None
            elif resource_info.manager != None:
                err_str = 'Resource manager type ' + resource_info.manager \
                    + ' options not handled yet for MPIBULL2.'
                raise ValueError, err_str

        # Info commands

        self.info_cmds = ['mpibull2-version']

    #---------------------------------------------------------------------------

    def __init_other__(self, resource_info=None):
        """
        Initialize MPI environment info for environments not handled
        in one the the previous cases.
        """

        # If possible, select launcher based on resource manager or
        # specific systems (would be better done in derived classes,
        # but these are cases on systems we have used in the past
        # but do not currently have access to).

        if os.uname()[0] == 'OSF1':
            if abs_exec_path('prun') != None:
                self.mpiexec = 'prun'
                self.mpiexec_n = ' -n '

        elif os.uname()[0] == 'AIX':
            if abs_exec_path('poe') != None:
                self.mpiexec = 'poe'
                self.mpiexec_n = None

    #---------------------------------------------------------------------------

    def info(self):
        """
        Outputs MPI environment information in known cases.
        """
        output = ''

        if self.info_cmds != None:
            for cmd in self.info_cmds:
                output += get_command_output(cmd) + '\n'

        return output

#-------------------------------------------------------------------------------
# Execution environment (including MPI, OpenMP, ...)
#-------------------------------------------------------------------------------

class exec_environment:

    #---------------------------------------------------------------------------

    def __init__(self, wdir=None, n_procs=None, hosts_list=None):
        """
        Returns Execution environment.
        """

        if sys.platform == 'win32' or sys.platform == 'win64':
            self.user = os.getenv('USERNAME')
        else:
            self.user = os.getenv('USER')

        self.wdir = wdir
        if self.wdir == None:
            self.wdir = os.getcwd()

        self.resources = resource_info(n_procs, hosts_list)

        self.mpi_env = None

        # Associate default launcher and associated options based on
        # known MPI types, use default otherwise

        self.mpi_env = mpi_environment(self.resources)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    e = exec_environment()

    mpi_env = mpi_environment()

    print 'mpi_env.bindir =        ', mpi_env.bindir
    print 'mpi_env.mpiexec =       ', mpi_env.mpiexec
    print 'mpi_env.mpiexec_args =  ', mpi_env.mpiexec_args
    print 'mpi_env.mpiexec_exe =   ', mpi_env.mpiexec_exe
    print 'mpi_env.mpiexec_n =     ', mpi_env.mpiexec_n
    print 'mpi_env.gen_hostsfile = ', mpi_env.gen_hostsfile
    print 'mpi_env.del_hostsfile = ', mpi_env.del_hostsfile
    print 'mpi_env.mpiboot =       ', mpi_env.mpiboot
    print 'mpi_env.mpihalt =       ', mpi_env.mpihalt
    print 'mpi_env.info_cmds =     ', mpi_env.info_cmds
    print 'mpi_env.mpmd =          ', mpi_env.mpmd
    print 'mpi_env.type =          ', mpi_env.type

    print mpi_env.info()

