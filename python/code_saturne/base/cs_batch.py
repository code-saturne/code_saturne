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

import sys
import os, os.path, shutil, sys, string, types, re, subprocess

#===============================================================================
# Utility functions
#===============================================================================

def parse_wall_time_slurm(s):
    """
    Parse SLURM wall time
    """
    t = None

    wt0 = s.split('-')
    if len(wt0) == 2:
        th = int(wt0[0])*3600*24
        wt = wt0[1].split(':')
    else:
        th = 0
        wt = wt0[0].split(':')
    if len(wt) == 3:
        t = th + int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
    elif len(wt) == 2:
        if len(wt0) == 2:
            t = th + int(wt[0])*3600 + int(wt[1])*60
        else:
            t = th + int(wt[0])*60 + int(wt[1])
    elif len(wt) == 1:
        if len(wt0) == 2:
            t = th + int(wt[0])*3600
        else:
            try:
                t = th + int(wt[0])*60
            except Exception:
                pass

    return t

#-------------------------------------------------------------------------------

def generate_header(batch_template=None, job_name=None, package=None):
    """
    Generate batch header lines based on batch template configuration
    """

    lines = []

    if not package:
        from code_saturne.base import cs_package
        package = cs_package.package()

    if not batch_template:
        config = configparser.ConfigParser()
        config.read(package.get_configfiles())
        if config.has_option('install', 'batch'):
            batch_template = config.get('install', 'batch')

    if not batch_template:
        return lines

    if not job_name:
        job_name = package.name + ':' + os.path.basename(os.getcwd())

    # For some systems, name lengths may be limited.
    # With SLURM 18.08, the max. length seems to be 37 characters.
    job_name = job_name[:38]

    if not os.path.isabs(batch_template):
        batch_template = os.path.join(package.get_batchdir(),
                                      'batch.' + batch_template)

    fdt = open(batch_template, 'r')

    import re, string
    kwd1 = re.compile('nameandcase')

    # Determine or build default names if required

    for line in fdt:
        line = line.rstrip()
        line = re.sub(kwd1, job_name, line)
        lines.append(line)

    fdt.close()

    return lines

#===============================================================================
# Class used to manage batch directives
#===============================================================================

#-------------------------------------------------------------------------------

class batch:
    """
    Parse and modify special batch options in text file lines
    """

    #---------------------------------------------------------------------------

    def __init__(self, package, install_config=None):
        """
        Constructor.
        """

        self.rm_type = None

        self.submit_cmd = ''

        self.params = {}

        self.params['job_name'] = None
        self.params['job_nodes'] = None
        self.params['job_ppn'] = None
        self.params['job_procs'] = None
        self.params['job_threads'] = None
        self.params['job_walltime'] = None
        self.params['job_class'] = None
        self.params['job_account'] = None
        self.params['job_wckey'] = None

        # Are we using a resource manager ?

        self.rm_type = None
        self.rm_template = None
        self.submit_command = ''

        submit_command = None

        if install_config:
            self.rm_template = install_config['batch']
            self.rm_help = install_config['batch_help']

        elif package:
            config = configparser.ConfigParser()
            config.read(package.get_configfiles())

            if config.has_option('install', 'batch'):
                self.rm_template = config.get('install', 'batch')
            if config.has_option('install', 'submit_command'):
                submit_command = config.get('install', 'submit_command')

        if self.rm_template:
            if os.path.isabs(self.rm_template):
                i = self.rm_template.rfind(".")
                if i > -1:
                    self.rm_template = self.rm_template[i+1:]

            if self.rm_template[0:5] == 'SLURM':
                self.rm_type = 'SLURM'
                self.submit_command = 'sbatch'
                for k in ('SBATCH_WCKEY', 'SLURM_WCKEY'):
                    if k in os.environ:
                        v = os.environ[k]
                        if v:
                            self.params['job_wckey'] = v
            elif self.rm_template[0:3] == 'CCC':
                self.rm_type = 'CCC'
                self.submit_command = 'ccc_msub'
            elif self.rm_template[0:3] == 'LSF':
                self.rm_type = 'LSF'
                self.submit_command = 'bsub'
            elif self.rm_template[0:3] == 'OAR':
                self.rm_type = 'OAR'
                self.submit_command = 'oarsub'
            elif self.rm_template[0:3] == 'PBS':
                self.rm_type = 'PBS'
                self.submit_command = 'qsub'
            elif self.rm_template[0:3] == 'SGE':
                self.rm_type = 'SGE'
                self.submit_command = 'qsub'
            else:
                i = self.rm_template.rfind(".")
                if i > -1:
                    self.rm_type = self.rm_template[i+1:]
                else:
                    self.rm_type = os.path.basename(rm_template)

        if submit_command:
            self.submit_command = submit_command

    #---------------------------------------------------------------------------

    def __get_command_output__(self, cmd):
        """
        Run a command and return it's standard output.
        """
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        lines = []
        while True:
            l = p.stdout.readline()
            lines.append(l.strip())
            if len(l) == 0 and p.poll() != None:
                break
        output = p.communicate()

        if p.returncode == 0:
            return lines
        else:
            print(output[1])
            return []

    #---------------------------------------------------------------------------

    def __pre_parse__(self, s):
        """
        Pre-parse batch line
        """
        r = ' '
        i = s.find('#')
        if i > -1:
            s = s[:i]
        s = r.join(s.split())

        return s

    #---------------------------------------------------------------------------

    def __parse_lines_env_vars__(self, lines):
        """
        Parse environment variables in batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            j = batch_lines[i].find('#')
            if j > -1:
                toks = batch_lines[i][:j].split()
            else:
                toks = batch_lines[i].split()
            if len(toks) > 1:
                var = None
                val = None
                if toks[0] in ('set', 'export'):
                    k = toks[1].find('=')
                    if k > 1:
                        var = toks[1][0:k]
                        val = toks[1][k+1:]
                elif toks[0] in ('setenv'):
                    if len(toks) > 2:
                        var = toks[1]
                        val = toks[2]
                if var == 'OMP_NUM_THREADS':
                    try:
                        self.params['job_threads'] = int(val)
                    except Exception:
                        pass

    #---------------------------------------------------------------------------

    def __update_lines_env_vars__(self, lines):
        """
        Update environment variables in batch file lines
        """
        if lines:
            batch_lines = lines

            for i in range(len(batch_lines)):
                j = batch_lines[i].find('#')
                if j > -1:
                    toks = batch_lines[i][:j].split()
                else:
                    toks = batch_lines[i].split()
                if len(toks) > 1:
                    var = None
                    val = None
                    if toks[0] in ('set', 'export'):
                        k = toks[1].find('=')
                        if k > 1:
                            var = toks[1][0:k]
                            val = toks[1][k+1:]
                    elif toks[0] in ('setenv'):
                        if len(toks) > 2:
                            var = toks[1]
                            val = toks[2]
                    if var == 'OMP_NUM_THREADS' and self.params['job_threads']:
                        s_threads = str(self.params['job_threads'])
                        if toks[0] in ('set', 'export'):
                            s = toks[0] + ' ' + var + '=' + s_threads
                        elif toks[0] in ('setenv'):
                            s = toks[0] + ' ' + var + ' ' + s_threads
                        if j > 1:
                            s += ' ' + batch_lines[i][j:]
                        batch_lines[i] = s

    #---------------------------------------------------------------------------

    def __parse_lines_slurm__(self, lines):
        """
        Parse SLURM batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:7] == '#SBATCH':
                batch_args = self.__pre_parse__(batch_lines[i][7:])
                if batch_args[0:2] == '--':
                    tok = batch_args.split('=')
                    if len(tok) < 2:
                        continue
                    kw = tok[0] + '='
                    val = tok[1].split(',')[0].strip()
                elif batch_args[0] == '-':
                    kw = batch_args[0:2]
                    val = batch_args[2:].split(',')[0].strip()
                else:
                    continue
                if kw == '--job-name=' or kw == '-J':
                    self.params['job_name'] = val
                elif kw == '--ntasks=' or kw == '-n':
                    self.params['job_procs'] = val
                elif kw == '--nodes=' or kw == '-N':
                    self.params['job_nodes'] = val
                elif kw == '--ntasks-per-node=':
                    self.params['job_ppn'] = val
                elif kw == '--cpus-per-task=':
                    self.params['job_threads'] = val
                elif kw == '--time=' or kw == '-t':
                    self.params['job_walltime'] = parse_wall_time_slurm(val)
                elif kw == '--partition=' or kw == '-p':
                    self.params['job_class'] = val
                elif kw == '--account=' or kw == '-A':
                    self.params['job_account'] = val
                elif kw == '--wckey=':
                    self.params['job_wckey'] = val

    #---------------------------------------------------------------------------

    def __update_lines_slurm__(self, lines):
        """
        Update the SLURM batch file from dictionary self.params.
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:7] == '#SBATCH':
                batch_args = self.__pre_parse__(batch_lines[i][7:])
                if batch_args[0:2] == '--':
                    tok = batch_args.split('=')
                    if len(tok) < 2:
                        continue
                    kw = tok[0] + '='
                    val = tok[1].split(',')[0].strip()
                elif batch_args[0] == '-':
                    kw = batch_args[0:2]
                    val = batch_args[2:].split(',')[0].strip()
                if kw == '--job-name=' or kw == '-J':
                    val = str(self.params['job_name'])
                elif kw == '--ntasks=' or kw == '-n':
                    val = str(self.params['job_procs'])
                elif kw == '--nodes=' or kw == '-N':
                    val = str(self.params['job_nodes'])
                elif kw == '--ntasks-per-node=':
                    val = self.params['job_ppn']
                elif kw == '--cpus-per-task=':
                    val = self.params['job_threads']
                elif kw == '--time=' or kw == '-t':
                    wt = self.params['job_walltime']
                    if wt > 86400: # 3600*24
                        val = '%d-%d:%02d:%02d' % (wt//86400,
                                                   (wt%86400)//3600,
                                                   (wt%3600)//60,
                                                   wt%60)
                    else:
                        val = '%d:%02d:%02d' % (wt//3600,
                                                (wt%3600)//60,
                                                wt%60)
                elif kw == '--partition=' or kw == '-p':
                    val = self.params['job_class']
                elif kw == '--account=' or kw == '-A':
                    val = self.params['job_account']
                elif kw == '--wckey=':
                    val = self.params['job_wckey']
                else:
                    continue
                batch_lines[i] = '#SBATCH ' + kw + str(val)

    #---------------------------------------------------------------------------

    def __parse_lines_ccc__(self, lines):
        """
        Parse CCC (CCRT) batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.__pre_parse__(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-r':
                    self.params['job_name'] = val
                elif kw == '-n':
                    self.params['job_procs'] = int(val)
                elif kw == '-N':
                    self.params['job_nodes'] = int(val)
                elif kw == '-T':
                    self.params['job_walltime'] = int(val)
                elif kw == '-q':
                        self.params['job_class'] = val

    #---------------------------------------------------------------------------

    def __update_lines_ccc__(self, lines):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#MSUB':
                batch_args = self.__pre_parse__(batch_lines[i][5:])
                tok = batch_args.split()
                if len(tok) < 2:
                    continue
                kw = tok[0]
                if kw == '-r':
                    val = str(self.params['job_name'])
                elif kw == '-n':
                    val = str(self.params['job_procs'])
                elif kw == '-N':
                    val = str(self.params['job_nodes'])
                elif kw == '-T':
                    val = str(self.params['job_walltime'])
                elif kw == '-q':
                    val = self.params['job_class']
                else:
                    continue
                batch_lines[i] = '#MSUB ' + kw + ' ' + str(val)

    #---------------------------------------------------------------------------

    def __parse_lines_lsf__(self, lines):
        """
        Parse Platform LSF batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.__pre_parse__(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                val = tok[1].split(',')[0].strip()
                if kw == '-J':
                    self.params['job_name'] = val
                elif kw == '-n':
                    self.params['job_procs'] = int(val)
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = val.split(':')
                    if len(wt) == 1:
                        self.params['job_walltime'] = int(wt[0])*60
                    elif len(wt) == 2:
                        self.params['job_walltime'] \
                            = int(wt[0])*3600 + int(wt[1])*60
                elif kw == '-q':
                        self.params['job_class'] = val

    #---------------------------------------------------------------------------

    def __update_lines_lsf__(self, lines):
        """
        Update the Platform LSF batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:5] == '#BSUB':
                batch_args = self.__pre_parse__(batch_lines[i][5:])
                tok = batch_args.split()
                kw = tok[0]
                if kw == '-J':
                    val = str(self.params['job_name'])
                elif kw == '-n':
                    val = str(self.params['job_procs'])
                elif kw == '-W' or kw == '-wt' or kw == '-We':
                    wt = self.params['job_walltime']
                    val = '%d:%02d' % (wt//3600, (wt%3600)//60)
                elif kw == '-q':
                    val = self.params['job_class']
                else:
                    continue
                batch_lines[i] = '#BSUB ' + kw + ' ' + str(val)

    #---------------------------------------------------------------------------

    def __parse_lines_pbs__(self, lines):
        """
        Parse PBS batch file lines
        """
        batch_lines = lines

        # TODO: specialize for PBS Professional and TORQUE (OpenPBS has not been
        # maintained since 2004, so we do not support it).
        # The "-l nodes=N:ppn=P" syntax is common to all PBS variants,
        # but PBS Pro considers the syntax depecated, and prefers its
        # own "-l select=N:ncpus=P:mpiprocs=P" syntax.
        # We do not have access to a PBS Professional system, but according to
        # its documentation, it has commands such as "pbs-report" or "pbs_probe"
        # which are not part of TORQUE, while the latter has "pbsnodelist" or
        # #pbs-config". The presence of either could help determine which
        # system is available.

        for i in range(len(batch_lines)):
            if batch_lines[i][0:4] == '#PBS':
                batch_args = ' ' + self.__pre_parse__(batch_lines[i][4:])
                index = batch_args.rfind(' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        self.params['job_name'] = arg.split()[1]
                    elif arg[0:9] == '-l nodes=':
                        arg_tmp = arg[9:].split(':')
                        self.params['job_nodes'] = arg_tmp[0]
                        for s in arg_tmp[1:]:
                            j = s.find('ppn=')
                            if j > -1:
                                self.params['job_ppn'] \
                                    = s[j:].split('=')[1]
                    elif arg[0:10] == '-l select=':
                        arg_tmp = arg[10:].split(':')
                        self.params['job_nodes'] = arg_tmp[0]
                        for s in arg_tmp[1:]:
                            j = s.find('ncpus=')
                            if j > -1:
                                self.params['job_ppn'] \
                                    = s[j:].split('=')[1]
                    elif arg[0:12] == '-l walltime=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.params['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.params['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.params['job_walltime'] \
                                = int(wt[0])
                    elif arg[0:2] == '-q':
                            self.params['job_class'] = arg.split()[1]
                    index = batch_args.rfind(' -')

    #---------------------------------------------------------------------------

    def __update_lines_pbs__(self, lines):
        """
        Update the PBS batch file from dictionary self.params.
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:4] == '#PBS':
                ch = ''
                batch_args = ' ' + self.__pre_parse__(batch_lines[i][4:])
                index = batch_args.rfind(' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.params['job_name'] + ch
                    elif arg[0:9] == '-l nodes=':
                        arg_tmp = arg[9:].split(':')
                        ch1 = ' -l nodes=' + self.params['job_nodes']
                        for s in arg_tmp[1:]:
                            j = s.find('ppn=')
                            if j > -1:
                                ch1 += ':' + s[0:j] \
                                       + 'ppn=' + self.params['job_ppn']
                            else:
                                ch1 += ':' + s
                        ch = ch1 + ch
                    elif arg[0:10] == '-l select=':
                        arg_tmp = arg[10:].split(':')
                        ch1 = ' -l select=' + self.params['job_nodes']
                        for s in arg_tmp[1:]:
                            j = s.find('ncpus=')
                            if j > -1:
                                ch1 += ':' + s[0:j] \
                                       + 'ncpus=' + self.params['job_ppn']
                            else:
                                ch1 += ':' + s
                        ch = ch1 + ch
                    elif arg[0:12] == '-l walltime=':
                        wt = self.params['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt//3600,
                                                 (wt%3600)//60,
                                                 wt%60)
                        ch = ' -l walltime=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.params['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = batch_args.rfind(' -')
                ch = '#PBS' + ch
                batch_lines[i] = ch

    #---------------------------------------------------------------------------

    def __parse_lines_sge__(self, lines):
        """
        Parse Sun Grid Engine batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                batch_args = ' ' + self.__pre_parse__(batch_lines[i][2:])
                index = batch_args.rfind(' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        self.params['job_name'] = arg.split()[1]
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            self.params['job_procs'] = arg_tmp[2]
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = (arg.split('=')[1]).split(':')
                        if len(wt) == 3:
                            self.params['job_walltime'] \
                                = int(wt[0])*3600 + int(wt[1])*60 + int(wt[2])
                        elif len(wt) == 2:
                            self.params['job_walltime'] \
                                = int(wt[0])*60 + int(wt[1])
                        elif len(wt) == 1:
                            self.params['job_walltime'] = int(wt[0])
                    elif arg[0:2] == '-q':
                        self.params['job_class'] = arg.split()[1]
                    index = batch_args.rfind(' -')

    #---------------------------------------------------------------------------

    def __update_lines_sge__(self, lines):
        """
        Update the Sun Grid Engine batch file lines
        """
        batch_lines = lines

        for i in range(len(batch_lines)):
            if batch_lines[i][0:2] == '#$':
                ch = ''
                batch_args = ' ' + self.__pre_parse__(batch_lines[i][2:])
                index = batch_args.rfind(' -')
                while index > -1:
                    arg = batch_args[index+1:]
                    batch_args = batch_args[0:index]
                    if arg[0:2] == '-N':
                        ch = ' -N ' + self.params['job_name'] + ch
                    elif arg[0:3] == '-pe':
                        try:
                            arg_tmp = arg[3:].split(' ')
                            ch = ' -pe ' + arg_tmp[1] + ' ' \
                                + str(self.params['job_procs']) + ch
                        except Exception:
                            pass
                    elif arg[0:8] == '-l h_rt=':
                        wt = self.params['job_walltime']
                        s_wt = '%d:%02d:%02d' % (wt//3600,
                                                 (wt%3600)//60,
                                                 wt%60)
                        ch = ' -l h_rt=' + s_wt + ch
                    elif arg[0:2] == '-q':
                        ch = ' -q ' + self.params['job_class'] + ch
                    else:
                        ch = ' ' + arg + ch
                    index = batch_args.rfind(' -')
                    ch = '#$' + ch
                    batch_lines[i] = ch

    #---------------------------------------------------------------------------

    def parse_lines(self, lines):
        """
        Fill self.params reading the job file.
        """

        if self.rm_type:
            if self.rm_type == 'SLURM':
                self.__parse_lines_slurm__(lines)
            elif self.rm_type == 'CCC':
                self.__parse_lines_ccc__(lines)
            elif self.rm_type == 'LSF':
                self.__parse_lines_lsf__(lines)
            elif self.rm_type == 'PBS':
                self.__parse_lines_pbs__(lines)
            elif self.rm_type == 'SGE':
                self.__parse_lines_sge__(lines)

        self.__parse_lines_env_vars__(lines)

    #---------------------------------------------------------------------------

    def update_lines(self, lines, keyword=None):
        """
        Update the batch file from reading dictionary self.params.
        If a keyword is given, its presence is checked.
        """

        if keyword != None:
            if keyword not in self.params:
                msg = str(keyword) + " is not in list " \
                    + str(list(self.params.keys())) + os.linesep
                raise ValueError(msg)

        if self.rm_type:
            if self.rm_type == 'SLURM':
                self.__update_lines_slurm__(lines)
            elif self.rm_type == 'CCC':
                self.__update_lines_ccc__(lines)
            elif self.rm_type == 'LSF':
                self.__update_lines_lsf__(lines)
            elif self.rm_type == 'PBS':
                self.__update_lines_pbs__(lines)
            elif self.rm_type == 'SGE':
                self.__update_lines_sge__(lines)

        self.__update_lines_env_vars__(lines)

    #---------------------------------------------------------------------------

    def submit_command_prefix(self):
        """
        Return the command prefix used to submit a job.
        """

        return self.submit_command + ' '

    #---------------------------------------------------------------------------

    def get_class_list(self):
        """
        Return the list of available classes
        """

        class_list = []

        try:

            rm_type = self.rm_type

            if rm_type == 'SLURM':
                output = self.__get_command_output__('sinfo -s')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        name = l.split(' ')[0]
                        if name[-1:] == '*':
                            name = name[:-1]
                        class_list.append(name)

            elif rm_type == 'CCC':
                output = self.__get_command_output__('class')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        class_list.append(l.split(' ')[0])

            elif rm_type == 'LSF':
                output = self.__get_command_output__('bqueues')
                ignore = True
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        class_list.append(l.split(' ')[0])

            elif rm_type == 'PBS':
                output = self.__get_command_output__('qstat -q')
                ignore = True
                for l in output:
                    if l[0:3] == '---':
                        ignore = not ignore
                    elif ignore == False:
                        class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'SGE':
                output = self.__get_command_output__('qconf -sc')
                for l in output:
                    if l[0:1] != '#':
                        class_list.append(l.split(' ')[0])

        except Exception:
            import traceback
            exc_info = sys.exc_info()
            bt = traceback.format_exception(*exc_info)
            print(bt)
            pass

        return class_list

    #---------------------------------------------------------------------------

    def get_help_text(self, package=None):
        """
        Read batch help text based on configuration
        """

        text = None
        help_path = None

        if not package:
            from code_saturne.base import cs_package
            package = cs_package.package()

        config = configparser.ConfigParser()
        config.read(package.get_configfiles())
        if config.has_option('install', 'batch_help'):
            help_path = config.get('install', 'batch_help')
            if not os.path.isabs(help_path):
                i = help_path.rfind(".")
                if i > -1:
                    help_path = 'batch_help.' + help_path
                help_path = os.path.join(package.get_batchdir(),
                                         help_path)

        elif self.rm_template != None:
            help_path = os.path.join(package.get_batchdir(),
                                     'batch_help.' + self.rm_template)
            if not os.path.isfile(help_path):
                help_path = os.path.join(package.get_batchdir(),
                                         'batch_help.' + self.rm_type)

        if help_path != None:
            try:
                fdt = open(help_path, 'r')
                text = fdt.read()
                fdt.close()
            except Exception:
                print('help file: ' + help_path + ' not present or readable.')
                pass

        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

