#!/usr/bin/env python
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009-2011  EDF
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

from cs_config import mpi_lib, env_modules

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

        # As a workaround for a bug in which the standard output an error
        # are "lost" (observed in an apparently random manner, with Python 2.4),
        # we only add the stdout and stderr keywords if they are non-default.

        kwargs = {}
        if (stdout != sys.stdout):
            kwargs['stdout'] = stdout
        if (stderr != sys.stderr):
            kwargs['stderr'] = stderr

        p = subprocess.Popen(cmd, shell=True, **kwargs)
        p.communicate()
        returncode = p.returncode

    else:

        p = popen2.Popen3(cmd, capturestderr=True)
        returncode = p.wait()
        output = (p.fromchild.read(), p.childerr.read())
        if len(output[0]) > 0:
            stdout.write(output[0])
        if len(output[1]) > 0:
            stderr.write(output[1])
        p.fromchild.close()
        p.childerr.close()

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

def set_modules():
    """
    Set environment modules if present.
    """
    if env_modules.modules == "no" or not os.environ.has_key('MODULESHOME'):
        return

    cmd_prefix = os.path.join(os.environ['MODULESHOME'], 'bin', 'modulecmd')

    if have_subprocess == True:
        cmds = ['purge']
        for m in env_modules.modules.strip().split():
            cmds.append('load ' + m)
        for cmd in cmds:
            (output, error) \
                = subprocess.Popen([cmd_prefix, 'python'] + cmd.split(),
                                   stdout=subprocess.PIPE).communicate()
            exec output

    else:
        cmds = ['purge']
        for m in env_modules.modules.strip().split():
            cmds.append('load ' + m)
        for cmd in cmds:
            p = os.popen(cmd_prefix + ' python ' + cmd).read()
            exec p

#-------------------------------------------------------------------------------

class batch_info:

    #---------------------------------------------------------------------------

    def __init__(self):

        """
        Get batch system information.
        """

        self.batch_type = None
        self.submit_dir = None
        self.job_file = None
        self.job_name = None
        self.job_id = None
        self.queue = None

        # Check for specific batch environments

        s = os.getenv('LSB_JOBID') # LSF
        if s != None:
            self.batch_type = 'LSF'
            self.submit_dir = os.getenv('LS_SUBCWDIR')
            self.job_file = os.getenv('LSB_JOBFILENAME')
            self.job_name = os.getenv('LSB_JOBNAME')
            self.job_id = os.getenv('LSB_BATCH_JID')
            self.queue = os.getenv('LSB_QUEUE')

        if self.batch_type == None:
            s = os.getenv('PBS_JOBID') # PBS
            if s != None:
                self.batch_type = 'PBS'
                self.submit_dir = os.getenv('PBS_O_WORKDIR')
                self.job_name = os.getenv('PBS_JOBNAME')
                self.job_id = os.getenv('PBS_JOBID')
                self.queue = os.getenv('PBS_QUEUE')

        if self.batch_type == None:
            s = os.getenv('LOADL_JOB_NAME') # LoadLeveler
            if s != None:
                self.batch_type = 'LOADL'
                self.submit_dir = os.getenv('LOADL_STEP_INITDIR')
                self.job_name = os.getenv('LOADL_JOB_NAME')
                self.job_id = os.getenv('LOADL_STEP_ID')
                self.queue = os.getenv('LOADL_STEP_CLASS')

        if self.batch_type == None:
            s = os.getenv('SGE_TASK_ID') # Sun Grid Engine
            if s != None:
                self.batch_type = 'SGE'
                self.submit_dir = os.getenv('SGE_O_WORKDIR')
                self.job_name = os.getenv('JOB_NAME')
                self.job_id = os.getenv('JOB_ID')
                self.queue = os.getenv('QUEUE')

        if self.batch_type == None:
            s = os.getenv('SLURM_JOBID') # SLURM
            if s != None:
                self.batch_type = 'SLURM'
                self.submit_dir = os.getenv('SLURM_SUBMIT_DIR')
                self.job_name = os.getenv('SLURM_JOB_NAME')
                self.job_id = os.getenv('SLURM_JOBID')
                self.queue = os.getenv('SLURM_PARTITION')

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

        # Check for resource manager

        # Test for SLURM (Simple Linux Utility for Resource Management).

        s = os.getenv('SLURM_NPROCS')
        if s != None:
            self.manager = 'SLURM'
            self.n_procs = int(s)
            s = os.getenv('SLURM_NNODES')
            if s != None:
                self.n_nodes = int(s)
        else:
            s = os.getenv('SLURM_NNODES')
            if s != None:
                self.manager = 'SLURM'
                self.n_nodes = int(s)
                s = os.getenv('SLURM_TASKS_PER_NODE')
                if s != None:
                    # Syntax may be similar to SLURM_TASKS_PER_NODE=2(x3),1"
                    # indicating three nodes will each execute 2 tasks and
                    # the  fourth node will execute 1 task.
                    self.n_procs = 0
                    for s0 in s.split(','):
                        i = s0.find('(')
                        if i > -1:
                            self.n_procs += int(s0[0:i])*int(s0[i+2:-1])
                        else:
                            self.n_procs += int(s0)
                else:
                    self.n_procs = self.n_nodes

        # Test for Platform LSF.

        if self.manager == None and self.batch_type == 'LSF':
            self.manager = 'LSF'
            self.n_procs = 0
            self.n_nodes = 0
            s = os.getenv('LSB_MCPU_HOSTS')
            if s != None:
                mcpu_list = s.split(' ')
                self.n_nodes = len(mcpu_list)/2
                for i in range(self.n_nodes):
                    self.n_procs += int(mcpu_list[i*2 + 1])
            else:
                s = os.getenv('LSB_HOSTS')
                if s != None:
                    hl = s.split(' ')
                    self.n_procs_from_hosts_list(hl, True)

        # Test for IBM LoadLeveler.

        if self.manager == None and self.batch_type == 'LOADL':
            s = os.getenv('LOADL_TOTAL_TASKS')
            if s == None:
                s = os.getenv('LOADL_BG_SIZE')
            if s != None:
                self.manager = 'LOADL'
                self.n_procs = int(s)
            else:
                s = os.getenv('LOADL_PROCESSOR_LIST')
                if s != None:
                    self.manager = 'LOADL'
                    hl = s.strip().split(' ')
                    self.n_procs_from_hosts_list(hl, True)
            s = os.getenv('LOADL_HOSTFILE')
            if s != None:
                self.manager = 'LOADL'
                self.hosts_file = '$LOADL_HOSTFILE'

        # Test for TORQUE or PBS Pro.

        if self.manager == None and self.batch_type == 'PBS':
            s = os.getenv('PBS_NODEFILE')
            if s != None:
                self.manager = 'PBS'
                self.hosts_file = '$PBS_NODEFILE'

        # Test for Oracle Grid Engine.

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
                    self.hosts_file = '$TMPDIR/machines'
            else:
                s = os.getenv('PE_HOSTFILE')
                if s != None:
                    self.hosts_file = '$PE_HOSTFILE'

        # Set an optional list of hosts if we are not running under
        # a resource manager.

        if hosts_list != None:
            if self.manager != None:
                sys.stderr.write('Warning:\n'
                                 + '   Host list will be ignored because a'
                                 + ' resource manager (' + self.manager
                                 + ') is in use.\n\n')
            else:
                self.hosts_list = hosts_list

        # Determine number of processors from hosts file or list

        if self.n_procs == None:
            if self.hosts_file != None:
                self.n_procs_from_hosts_file(self.hosts_file)
            elif self.hosts_list != None:
                self.n_procs_from_hosts_list(self.hosts_list)

        # Check and possibly set number of processes

        if n_procs != None:
            if self.n_procs != None:
                if self.n_procs != n_procs:
                    sys.stderr.write('Warning:\n'
                                     +'   Will use ' + str(self.n_procs)
                                     + ' processes while resource manager ('
                                     + self.manager + ')\n   allows for '
                                     + str(n_procs) + '.\n\n')
            self.n_procs = n_procs

    #---------------------------------------------------------------------------

    def n_procs_from_hosts_file(self, hosts_file):

        """
        Compute number of hosts from a hostsfile.
        """

        self.n_procs = 0
        if hosts_file[0] == '$':
           path = os.getenv(hosts_file[1:])
        else:
           path = hosts_file
        f = open(path, 'r')
        for line in f:
            self.n_procs += 1
        f.close()

    #---------------------------------------------------------------------------

    def n_procs_from_hosts_list(self, hosts_list, is_copy=False):

        """
        Determine number of processors and nodes from hosts list.
        """

        self.n_procs = len(hosts_list)
        self.n_nodes = 1

        # If the hosts list is not already a copy, build one so
        # that sorting will not alter the original list.

        if is_copy == False:
            hl = []
            for s in hosts_list:
                hl.append(s)
        else:
            hl = hosts_list

        hl.sort()

        for i in range(self.n_procs - 1):
            if hl[i] != hl[i+1]:
                self.n_nodes += 1

    #---------------------------------------------------------------------------

    def get_hosts_list(self):

        """
        Get execution resources information.
        """

        hosts_list = None

        # Hosts list may already have been defined by constructor

        if self.hosts_list != None:
            hosts_list = self.hosts_list

        # Check for resource manager and eventual hostsfile

        elif self.manager == 'SLURM':
            s = os.getenv('SLURM_TASKS_PER_NODE')
            hosts_count = None
            if s != None:
                # Syntax may be similar to SLURM_TASKS_PER_NODE=2(x3),1"
                # indicating three nodes will each execute 2 tasks and
                # the  fourth node will execute 1 task.
                hosts_count = []
                for s0 in s.split(','):
                    i = s0.find('(')
                    if i > -1:
                        for j in range(int(s0[i+2:-1])):
                            hosts_count.append(int(s0[0:i]))
                    else:
                        hosts_count.append(int(s0))
            s = os.getenv('SLURM_NODELIST')
            if s != None:
                hosts_list = []
                l = 0
                # List uses a compact representation
                for s0 in s.split(','):
                    i = s0.find('[')
                    if i > -1:
                        basename = s0[0:i]
                        for s1 in s0[i+1:-1].split(','):
                            s2 = s1.split('-')
                            if len(s2) > 1:
                                fmt_s = basename + '%0' + str(len(s2[0])) + 'd'
                                for j in range(int(s2[0]), int(s2[1])+1):
                                    host_name = fmt_s % j
                                    if hosts_count:
                                        for k in range(hosts_count[l]):
                                            hosts_list.append(host_name)
                                        l += 1
                                    else:
                                        hosts_list.append(host_name)
                            else:
                                if hosts_count:
                                    for k in range(hosts_count[l]):
                                        hosts_list.append(s2[0])
                                    l += 1
                                else:
                                    hosts_list.append(s2[0])
                    else:
                        if hosts_count:
                            for k in range(hosts_count[l]):
                                hosts_list.append(s0)
                            l += 1
                        else:
                            hosts_list.append(s0)
            else:
                hosts_list = get_command_output('srun hostname -s').split()
                hosts_list.sort()

        elif self.manager == 'LSF':
            s = os.getenv('LSB_MCPU_HOSTS')
            if s != None:
                mcpu_list = s.split(' ')
                hosts_list = []
                for i in range(len(mcpu_list)/2):
                    host = mcpu_list[i*2]
                    count = int(mcpu_list[i*2 + 1])
                    for j in range(count):
                        hosts_list.append(host)
            else:
                s = os.getenv('LSB_HOSTS')
                if s != None:
                    hosts_list = s.split(' ')

        elif self.manager == 'LOADL':
            hosts_list = []
            s = os.getenv('LOADL_PROCESSOR_LIST')
            if s != None:
                hosts_list = s.split(' ')

        return hosts_list

    #---------------------------------------------------------------------------

    def get_hosts_file(self, wdir = None):
        """
        Returns the name of the hostsfile associated with the
        resource manager. A hostsfile is built from a hosts
        list if necessary.
        """

        hosts_file = self.hosts_file

        if self.hosts_file == None:

            hosts_list = self.get_hosts_list()

            if hosts_list != None:
                if wdir != None:
                    hosts_file = os.path.join(wdir, 'hostsfile')
                else:
                    hosts_file = 'hostsfile'
                f = open(hosts_file, 'w')
                # If number of procs not specified, determine it by list
                if self.n_procs == None or self.n_procs < 1:
                    n_procs = 0
                    for host in hosts_list:
                        f.write(host + '\n')
                        n_procs += 1
                # If the number of procs is known, use only beginning of list
                # if it contains more entries than the required number, or
                # loop through list as many time as necessary to reach
                # prescribed proc count
                else:
                    proc_count = 0
                    while proc_count < self.n_procs:
                        for host in hosts_list:
                            if proc_count < self.n_procs:
                                f.write(host + '\n')
                                proc_count += 1
                f.close()

            self.hosts_file = hosts_file

        return hosts_file

#-------------------------------------------------------------------------------
# MPI environments and associated commands
#-------------------------------------------------------------------------------

MPI_MPMD_none       = 0
MPI_MPMD_mpiexec    = (1<<0) # mpiexec colon-separated syntax
MPI_MPMD_configfile = (1<<1) # mpiexec -configfile syntax
MPI_MPMD_script     = (1<<2)
MPI_MPMD_execve     = (1<<3)

class mpi_environment:

    def __init__(self, resource_info=None, wdir = None):
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
        self.mpmd = MPI_MPMD_none

        self.info_cmds = None

        # Initialize based on known MPI types, or default.

        init_method = self.__init_other__

        if len(self.type) > 0:
            mpi_env_by_type = {'MPICH2':self.__init_mpich2__,
                               'MPICH1':self.__init_mpich1__,
                               'OpenMPI':self.__init_openmpi__,
                               'LAM_MPI':self.__init_lam__,
                               'BGL_MPI':self.__init_bgl__,
                               'BGP_MPI':self.__init_bgp__,
                               'HP_MPI':self.__init_hp_mpi__,
                               'MPIBULL2':self.__init_mpibull2__}
            if self.type in mpi_env_by_type:
                init_method = mpi_env_by_type[self.type]

        p = os.getenv('PATH').split(':')
        if len(self.bindir) > 0:
            p.insert(0, self.bindir)

        init_method(p, resource_info, wdir)

    #---------------------------------------------------------------------------

    def __init_mpich2__(self, p, resource_info=None, wdir = None):

        """
        Initialize for MPICH2 environment.

        MPICH2 allows for 4 different process managers, all or some
        of which may be built depending on installation options:

        - HYDRA is the default starting with MPICH2-1.3, and is available
          in MPICH2-1.2. It natively uses existing daemons on the system
          such as ssh, SLURM, PBS, etc.

        - MPD is the traditional process manager, which consists of a
          ring of daemons. A hostsfile may be defined to start this
          ring (a machinefile may still be used for mpiexec for
          finer control, though we do not use it here).
          We try to test if such a ring is already running, and
          start and stop it if this is not the case.

        - spmd may be used both on Windows and Linux. It consists
          of independent daemons, so if a hostsfile is used, it
          must be passed to mpiexec (we do not attempt to start
          or stop the daemons here if this manager is used).

        - gforker is a simple manager that creates all processes on
          a single machine (the equivalent seem possible with HYDRA
          using "mpiexec -bootstrap fork")
        """

        # Determine base executable paths

        # Also try to determine if mpiexec is that of MPD or of another
        # process manager, using the knowledge that MPD's mpiexec
        # is a Python script, while other MPICH2 mpiexec's are binary.
        # We could have a false positive if another wrapper is used,
        # but we still need to find mpdboot.
        # In a similar fashion, we may determine if we are using the
        # Hydra process manager, as the 'mpiexec --help' will contain
        # a 'Hydra' string.

        # Also, mpirun is a wrapper to mpdboot + mpiexec + mpiallexec,
        # so it does not require running mpdboot and mpdallexit separately.

        launcher_names = ['mpiexec.mpich2', 'mpiexec',
                          'mpiexec.hydra', 'mpiexec.mpd', 'mpiexec.gforker',
                          'mpirun.mpich2', 'mpirun']
        pm = ''

        for d in p:
            for name in launcher_names:
                absname = os.path.join(d, name)
                if os.path.isfile(absname):
                    # Try to determine launcher type
                    basename = os.path.basename(name)
                    if basename in ['mpiexec.mpich2', 'mpiexec',
                                    'mpirun.mpich2', 'mpirun']:
                        info = get_command_outputs(absname)
                        if info.find('Hydra') > -1:
                            pm = 'hydra'
                        elif info.find(' mpd ') > -1:
                            pm = 'mpd'
                        elif info.find('-usize') > -1:
                            pm = 'gforker'
                    elif basename == 'mpiexec.hydra':
                        pm = 'hydra'
                    elif basename == 'mpiexec.mpd':
                        pm = 'mpd'
                    elif basename == 'mpiexec.gforker':
                        pm = 'gforker'
                    # Set launcher name
                    if d == self.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        # Determine if MPD should be handled
        # (if we are using a root MPD, no need for setup)

        if pm == 'mpd' and basename[:6] != 'mpirun':

            mpd_setup = True
            s = os.getenv('MPD_USE_ROOT_MPD')
            if s != None and int(s) != 0:
                mpd_setup = False

            # If a setup seems necessary, check paths
            if mpd_setup:
                if os.path.isfile(os.path.join(d, 'mpdboot')):
                    if d == self.bindir:
                        self.mpiboot = os.path.join(d, 'mpdboot')
                        self.mpihalt = os.path.join(d, 'mpdallexit')
                        mpdtrace = os.path.join(d, 'mpdtrace')
                        mpdlistjobs = os.path.join(d, 'mpdlistjobs')
                    else:
                        self.mpiboot = 'mpdboot'
                        self.mpihalt = 'mpdallexit'
                        mpdtrace = 'mpdtrace'
                        mpdlistjobs = 'mpdlistjobs'

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        if launcher_base[:7] == 'mpiexec':
            self.mpmd = MPI_MPMD_mpiexec | MPI_MPMD_configfile | MPI_MPMD_script
            self.mpiexec_n = ' -n '
        elif launcher_base[:6] == 'mpirun':
            self.mpiexec_n = ' -np '
            self.mpmd = MPI_MPMD_script

        # Other options to add

        # Resource manager info

        rm = None
        if resource_info != None:
            rm = resource_info.manager

        if pm == 'mpd':
            if rm == 'SLURM':
                # This requires linking with SLURM's implementation
                # of the PMI library.
                self.mpiexec = 'srun'
                self.mpiexec_n = ' -n'
                self.mpmd = MPI_MPMD_script
                self.mpiboot = None
                self.mpihalt = None
            elif rm == 'PBS':
                # Convert PBS to MPD format (based on MPICH2 documentation)
                # before MPI boot.
                if self.mpiboot != None:
                    self.gen_hostsfile = 'sort $PBS_NODEFILE | uniq -C ' \
                        + '| awk \'{ printf("%s:%s", $2, $1); }\' > ./mpd.nodes'
                    self.del_hostsfile = 'rm -f ./mpd.nodes'
                    self.mpiboot += ' --file=./mpd.nodes'
            else:
                hostsfile = resource_info.get_hosts_file(wdir)
                if hostsfile != None:
                    self.mpiboot += ' --file=' + hostsfile

        elif pm == 'hydra':
            # Nothing to do for resource managers directly handled by Hydra
            if rm not in ['PBS', 'LOADL', 'LSF', 'SGE', 'SLURM']:
                hostsfile = resource_info.get_hosts_file(wdir)
                if hostsfile != None:
                    self.mpiexec += ' -f ' + hostsfile

        elif pm == 'gforker':
            hosts = False
            hostslist = resource_info.get_hosts_list()
            if hostslist != None:
                hosts = True
            else:
                hostsfile = resource_info.get_hosts_file(wdir)
                if hostsfile != None:
                    hosts = True
            if hosts == True:
                sys.stderr.write('Warning:\n'
                                 + '   Hosts list will be ignored by'
                                 + ' MPICH2 gforker program manager.\n\n')

        # Finalize mpiboot and mpihalt commands.
        # We use 'mpdtrace' to determine if a ring is already running,
        # and mpdlistjobs to determine if other jobs are still running.
        # This means that a hostsfile will be ignored if an MPD ring
        # is already running, but will avoid killing other running jobs.

        if self.mpiboot != None:
            self.mpiboot = \
                mpdtrace + ' > /dev/null 2>&1\n' \
                + 'if test $? != 0 ; then ' + self.mpiboot + ' ; fi'
            self.mpihalt = \
                'listjobs=`' + mpdlistjobs + ' | wc -l`\n' \
                + 'if test $listjobs = 0 ; then ' + self.mpihalt + ' ; fi'

        # Info commands

        self.info_cmds = ['mpich2version']

    #---------------------------------------------------------------------------

    def __init_mpich1__(self, p, resource_info=None, wdir = None):

        """
        Initialize for MPICH1 environment.
        """

        # Determine base executable paths

        launcher_names = ['mpirun.mpich',
                          'mpirun.mpich-mpd',
                          'mpirun.mpich-shmem',
                          'mpirun']

        for d in p:
            for name in launcher_names:
                absname = os.path.join(d, name)
                if os.path.isfile(absname):
                    if d == self.bindir:
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
            hostsfile = resource_info.get_hosts_file(wdir)
            if hostsfile != None:
                self.mpiexec += ' -machinefile ' + hostsfile

        # Info commands

        self.info_cmds = ['mpichversion']

    #---------------------------------------------------------------------------

    def __init_openmpi__(self, p, resource_info=None, wdir = None):
        """
        Initialize for OpenMPI environment.
        """

        # Determine base executable paths

        launcher_names = ['mpiexec.openmpi', 'mpirun.openmpi',
                          'mpiexec', 'mpirun']
        info_name = ''

        for d in p:
            for name in launcher_names:
                absname = os.path.join(d, name)
                if os.path.isfile(absname):
                    if d == self.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    info_name = os.path.join(d, 'ompi_info')
                    break
            if self.mpiexec != None:
                break

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        self.mpiexec_n = ' -n '
        if launcher_base[:7] == 'mpiexec':
            self.mpmd = MPI_MPMD_mpiexec | MPI_MPMD_script
        elif launcher_base[:7] == 'mpirun':
            self.mpmd = MPI_MPMD_script

        # Other options to add

        # Detect if resource manager is known by this Open MPI build

        if resource_info != None:
            known_manager = False
            if resource_info.manager == 'PBS':
                known_manager = True
            elif os.path.isfile(info_name):
                rc_mca_by_type = {'SLURM':' slurm ',
                                  'LSF':' lsf ',
                                  'LOADL':' loadleveler ',
                                  'SGE':' gridengine '}
                if resource_info.manager in rc_mca_by_type:
                    info = get_command_output(info_name)
                    if info.find(rc_mca_by_type[resource_info.manager]) > -1:
                        known_manager = True
            if known_manager == False:
                hostsfile = resource_info.get_hosts_file(wdir)
                if hostsfile != None:
                    self.mpiexec += ' --machinefile ' + hostsfile

        # Info commands

        self.info_cmds = ['ompi_info -a']

    #---------------------------------------------------------------------------

    def __init_lam__(self, p, resource_info=None, wdir = None):
        """
        Initialize for LAM/MPI environment.
        """

        # Determine base executable paths

        launcher_names = ['mpiexec.lam', 'mpirun.lam', 'mpiexec', 'mpirun']

        for d in p:
            for name in launcher_names:
                absname = os.path.join(d, name)
                if os.path.isfile(absname):
                    if d == self.bindir:
                        self.mpiexec = absname
                    else:
                        self.mpiexec = name
                    break
            if self.mpiexec != None:
                break

        if os.path.isfile(os.path.join(d, 'lamboot')):
            if d == self.bindir:
                self.mpiboot = os.path.join(d, 'lamboot')
                self.mpihalt = os.path.join(d, 'lamhalt')
            else:
                self.mpiboot = 'lamboot'
                self.mpihalt = 'lamhalt'

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        if launcher_base[:7] == 'mpiexec':
            self.mpiexec_n = ' -n '
            self.mpmd = MPI_MPMD_mpiexec | MPI_MPMD_script
        elif launcher_base[:7] == 'mpirun':
            self.mpiexec_n = ' -np '
            self.mpmd = MPI_MPMD_script

        # Determine options to add

        if self.mpiboot != None:
            self.mpiboot += ' -v'
            self.mpihalt += ' -v'

        # Resource manager info

        if resource_info != None:
            hostsfile = resource_info.get_hosts_file(wdir)
            if hostsfile != None and self.mpiboot != None:
                self.mpiboot += ' ' + hostsfile
                self.mpihalt += ' ' + hostsfile

        # Finalize mpiboot and mpihalt commands

        if self.mpiboot != None:
            self.mpiboot += '|| exit $?'

        # Info commands

        self.info_cmds = ['laminfo -all']

    #---------------------------------------------------------------------------

    def __init_bgl__(self, p, resource_info=None, wdir = None):

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

    def __init_bgp__(self, p, resource_info=None, wdir = None):

        """
        Initialize for Blue Gene/P environment.
        """

        # Set base executable path

        self.mpiexec = 'mpiexec'

        # Determine processor count and MPMD handling

        launcher_base = os.path.basename(self.mpiexec)

        self.mpiexec_n = None
        self.mpmd = MPI_MPMD_configfile

        # Other options to add

        # Info commands

    #---------------------------------------------------------------------------

    def __init_hp_mpi__(self, p, resource_info=None, wdir = None):
        """
        Initialize for HP MPI environment.

        The last version of HP MPI is version 2.3, released early 2009.
        HP MPI was then acquired by Platform MPI (formerly Scali MPI),
        which merged Scali MPI and HP MPI in Platform MPI 8.0.
        Platform MPI still seems to use the mpirun launcher syntax.
        """

        # Determine base executable paths

        self.mpiexec = 'mpirun'

        for d in p:
            if not os.path.isabs(self.mpiexec):
                absname = os.path.join(d, self.mpiexec)
                if os.path.isfile(absname):
                    if d == self.bindir:
                        self.mpiexec = absname
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

    def __init_mpibull2__(self, p, resource_info=None, wdir = None):
        """
        Initialize for MPIBULL2 environment.
        """

        self.__init_mpich2__(p, resource_info)

        # On a Bull Novascale machine using MPIBULL2 (based on MPICH2),
        # mpdboot and mpdallexit commands and related mpiexec may be found,
        # but the SLURM srun launcher or mpibull2-launch meta-launcher
        # (which can integrate directly to several resource managers
        # using prun, srun, orterun, or mprun) should be simpler to use.

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
                raise ValueError(err_str)

        # Info commands

        self.info_cmds = ['mpibull2-version']

    #---------------------------------------------------------------------------

    def __init_other__(self, p, resource_info=None, wdir = None):
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

        self.mpmd = MPI_MPMD_script

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

    #---------------------------------------------------------------------------

    def unset_mpmd_mode(self, mpi_mpmd_mode):
        """
        Unset mask allowing a given mpmd mode.
        """

        if self.mpmd & mpi_mpmd_mode:
            self.mpmd = self.mpmd ^ mpi_mpmd_mode

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

        self.mpi_env = mpi_environment(self.resources, wdir)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    e = exec_environment()

    mpi_env = mpi_environment()

    print('mpi_env.bindir =        ', mpi_env.bindir)
    print('mpi_env.mpiexec =       ', mpi_env.mpiexec)
    print('mpi_env.mpiexec_args =  ', mpi_env.mpiexec_args)
    print('mpi_env.mpiexec_exe =   ', mpi_env.mpiexec_exe)
    print('mpi_env.mpiexec_n =     ', mpi_env.mpiexec_n)
    print('mpi_env.gen_hostsfile = ', mpi_env.gen_hostsfile)
    print('mpi_env.del_hostsfile = ', mpi_env.del_hostsfile)
    print('mpi_env.mpiboot =       ', mpi_env.mpiboot)
    print('mpi_env.mpihalt =       ', mpi_env.mpihalt)
    print('mpi_env.info_cmds =     ', mpi_env.info_cmds)
    print('mpi_env.mpmd =          ', mpi_env.mpmd)
    print('mpi_env.type =          ', mpi_env.type)

    print(mpi_env.info())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
