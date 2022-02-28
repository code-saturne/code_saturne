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

import sys, os
import configparser

from code_saturne.base.cs_exec_environment import separate_args, assemble_args, \
    enquote_arg, get_command_single_value, update_command_single_value, \
    update_command_no_value

#===============================================================================
# Class used to manage runcase files
#===============================================================================

class runcase(object):

    def __init__(self,
                 path,
                 package=None,
                 submit=False,
                 job_header=None,
                 prologue=None,
                 epilogue=None):
        """
        Initialize runcase info object.
        """

        self.path = path
        self.package = package

        if submit:
            self.build_template(job_header=job_header,
                                prologue=prologue, epilogue=epilogue)
        else:
            f = open(self.path, mode = 'r')
            self.lines = f.readlines()
            f.close()

            for i in range(len(self.lines)):
                self.lines[i] = self.lines[i].rstrip()

            self.get_run_command()

    #---------------------------------------------------------------------------

    def save(self):
        """
        Save runcase.
        """

        f = open(self.path, mode = 'w')
        for line in self.lines:
            if line.startswith('export PATH='):
                if self.package != None:
                    line = "export PATH=" + self.package.get_dir("bindir")  + ":$PATH"
            f.write(line + '\n')
        f.close()

    #---------------------------------------------------------------------------

    def get_run_command(self):
        """
        Determine the name of the main command of the runcase, and the associated
        line in the file; this allows mixing code_saturne and neptune_cfd
        cases in the same study.
        """
        # Read the runcase script from the Repository.

        self.cmd_name = None
        self.run_cmd_line_id = -1

        for i in range(len(self.lines) - 1, -1, -1):

            line = self.lines[i]

            # Skip comment and empty lines

            if len(line) == 0:
                continue
            if line[0] == '#' or line[0:4] == 'rem ':
                continue
            j = line.find('#')
            if j > -1:
                line = line[0:j]

            args = separate_args(line.strip())
            if args.count('run') == 1:
                if args.index('run') == 1: # "<package_name> run"
                    for name in ('code_saturne', 'neptune_cfd'):
                        if sys.platform.startswith('win'):
                            test_name = name
                            if args[0].find(test_name) == 0:
                                self.cmd_name = name
                                self.run_cmd_line_id = i
                                return

                        test_name = '\\' + name
                        if args[0].find(test_name) == 0:
                            self.cmd_name = name
                            self.run_cmd_line_id = i
                            return
                        elif args[0] == name:
                            self.cmd_name = name
                            self.run_cmd_line_id = i
                            return

        # We should have exited before reaching this.

        err_str = "Error: unable to determine the name of the script for " + self.path + os.linesep
        sys.stderr.write(err_str)

    #---------------------------------------------------------------------------

    def get_run_args(self):
        """
        Get the run command and arguments, as a list
        """

        return separate_args(self.lines[self.run_cmd_line_id])

    #---------------------------------------------------------------------------

    def build_template(self, job_header=None,
                       prologue=None, epilogue=None):
        """
        Build batch file template
        """

        import os, stat
        from code_saturne.base.cs_exec_environment import append_shell_shebang, \
            append_script_comment

        if not self.package:
            from code_saturne.base import cs_package
            self.package = cs_package.package()

        # Use alternate wrapper if configured

        wrapper_postfix = None

        config = configparser.ConfigParser()
        config.read(self.package.get_configfiles())
        if config.has_option('install', 'wrapper_postfix'):
            wrapper_postfix = config.get('install', 'wrapper_postfix')

        self.lines = []

        append_shell_shebang(self.lines)

        # Add batch system info and user prologue if necessary

        if job_header:
            for line in job_header.split(os.linesep):
                self.lines.append(line)
            self.lines.append('')

        # Ensure switch to script directory

        append_script_comment(self.lines,
                              'Set working directory:' + os.linesep)
        self.lines.append('cd ' + enquote_arg(os.path.dirname(self.path)))
        self.lines.append('')

        if prologue:
            append_script_comment(self.lines,
                                  'User prologue:' + os.linesep)
            for line in prologue.split(os.linesep):
                self.lines.append(line)
            self.lines.append('')

        # Add command to execute.

        append_script_comment(self.lines, 'Run command:' + os.linesep)

        if wrapper_postfix is None:
            exec_path = os.path.join(self.package.get_dir("bindir"), self.package.name)
        else:
            exec_path = '\\' + self.package.name + wrapper_postfix

        run_cmd = enquote_arg(exec_path) + ' run'
        self.cmd_name = self.package.name
        self.run_cmd_line_id = len(self.lines)
        self.lines.append(run_cmd)
        self.lines.append('')

        if epilogue:
            append_script_comment(self.lines,
                                  'User epilogue:' + os.linesep)
            for line in epilogue.split(os.linesep):
                self.lines.append(line)
            self.lines.append('')

        self.save()

        st   = os.stat(self.path)
        mode = st[stat.ST_MODE]
        os.chmod(self.path, mode | stat.S_IEXEC)

    #---------------------------------------------------------------------------

    def get_compute_build(self):
        """
        Get the compute-build option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--compute-build',
                                         '--compute-build='))

    #---------------------------------------------------------------------------

    def get_parameters(self):
        """
        Get the parameters option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--param', '--param=', '-p'))

    #---------------------------------------------------------------------------

    def get_nprocs(self):
        """
        Get the nprocs option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--nprocs', '--nprocs=', '-n'))

    #---------------------------------------------------------------------------

    def get_nthreads(self):
        """
        Get the nthreads option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--threads-per-task',
                                         '--threads-per-task=',
                                         '-nt'))

    #---------------------------------------------------------------------------

    def get_run_id(self):
        """
        Get the run id, id_prefix, and id_suffix options in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        run_id = get_command_single_value(args, ('--id', '--id='))
        run_id_prefix = get_command_single_value(args, ('--id-prefix',
                                                        '--id-prefix='))
        run_id_suffix = get_command_single_value(args, ('--id-suffix',
                                                        '--id-suffix='))

        return run_id, run_id_prefix, run_id_suffix

    #---------------------------------------------------------------------------

    def get_run_stage(self, stage):
        """
        Return True if a given stage is specified in the run command,
        False otherwise
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        s_arg = '--' + stage
        if s_arg in args:
            return True

        return False

    #---------------------------------------------------------------------------

    def run_conf_sections(self,
                          resource_name='job_defaults',
                          batch_template=None):
        """
        Build "run_conf" sections from an existing runcase.
        """

        sections = {}

        param = self.get_parameters()

        setup_dict = {}
        if param != 'setup.xml':
            setup_dict['param'] = param
        if len(setup_dict) > 0:
            sections['setup'] = setup_dict

        run_dict = {}
        compute_build = self.get_compute_build()
        if compute_build:
            run_dict['compute_build'] = compute_build
        run_id, run_id_prefix, run_id_suffix = self.get_run_id()
        if run_id:
            run_dict['id'] = run_id
        if run_id_prefix:
            run_dict['id_prefix'] = run_id_prefix
        if run_id_suffix:
            run_dict['id_suffix'] = run_id_suffix
        stage_map = {'stage': 'stage',
                     'initialize': 'preprocess',
                     'execute': 'compute',
                     'finalize': 'finalize'}
        for stage in stage_map.keys():
            if self.get_run_stage(stage):
                run_dict[stage_map[stage]] = True
        if len(run_dict) > 0:
            sections['run'] = run_dict

        resource_dict = {}
        nprocs = self.get_nprocs()
        if nprocs:
            resource_dict['n_procs'] = nprocs
        nthreads = self.get_nthreads()
        if nthreads:
            resource_dict['n_threads'] = nthreads

        # Handle batch info
        if batch_template:
            from code_saturne.base import cs_batch
            batch_src = cs_batch.batch(self.package)
            batch_src.parse_lines(self.lines)

            topdir, scriptdir = os.path.split(os.path.split(self.path)[0])
            if scriptdir == 'SCRIPTS':
                studydir, casedir = os.path.split(topdir)
                studydir = os.path.split(studydir)[1]
            else:
                casedir = ''
                studydir = scriptdir
            job_name = studydir.lower() + casedir.lower()

            dst_header = cs_batch.generate_header(batch_template=batch_template,
                                                  job_name=job_name,
                                                  package=self.package)
            batch_dst = cs_batch.batch(self.package)
            batch_dst.parse_lines(dst_header)
            for k in batch_src.params.keys():
                if batch_src.params[k] != None:
                    batch_dst.params[k] = batch_src.params[k]
            batch_dst.update_lines(dst_header)

            br = os.linesep

            job_header = ''
            for i, l in enumerate(dst_header):
                if i == 0:
                    job_header = l
                else:
                    job_header += br + l

            resource_dict['job_header'] = job_header

        if len(resource_dict) > 0:
            sections[resource_name] = resource_dict

        return sections

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
