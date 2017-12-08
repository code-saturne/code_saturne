#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

try:
    from code_saturne.cs_exec_environment import separate_args, assemble_args, \
        enquote_arg, get_command_single_value, update_command_single_value
except exception:
    from cs_exec_environment import separate_args, assemble_args, \
        enquote_arg, get_command_single_value, update_command_single_value

#===============================================================================
# Class used to manage runcase files
#===============================================================================

class runcase(object):

    def __init__(self,
                 path,
                 package=None,
                 create_if_missing=True,
                 rebuild=False,
                 study_name=None,
                 case_name=None,
                 ignore_batch=False):
        """
        Initialize runcase info object.
        """

        self.path = path
        self.package = package

        try:
            f = open(self.path, mode = 'r')
            self.lines = f.readlines()
            f.close()
            if rebuild:
                self.get_run_command()
                command_line = self.lines[self.run_cmd_line_id]
                self.build_template(package,
                                    study_name,
                                    case_name,
                                    ignore_batch)
                self.lines[self.run_cmd_line_id] = command_line

        except IOError:
            if create_if_missing or rebuild:
                self.build_template(package, study_name, case_name, ignore_batch)
            else:
                print("Error: can not open or read %s\n" % self.path)
                sys.exit(1)

        for i in range(len(self.lines)):
            self.lines[i] = self.lines[i].rstrip()

        self.get_run_command()

    #---------------------------------------------------------------------------

    def save(self, path=None):
        """
        Save runcase, optionally to a different path.
        """

        if (path):
            self.path = path

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
        line in the file; this allows mixing Code_Saturne and NEPTUNE_CFD
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

            args = separate_args(line.rstrip())
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

        # We should have exited before reaching this.

        err_str = "Error: unable to determine the name of the script for " + self.path
        raise ValueError(err_str)

    #---------------------------------------------------------------------------

    def get_run_args(self):
        """
        Get the run command and arguments, as a list
        """

        return separate_args(self.lines[self.run_cmd_line_id])

    #---------------------------------------------------------------------------

    def set_run_args(self, args):
        """
        Set the run command and arguments from a list
        """

        self.lines[self.run_cmd_line_id] = assemble_args(args)

    #---------------------------------------------------------------------------

    def build_template(self, package=None, study_name=None, case_name=None,
                       ignore_batch=False):
        """
        Build batch file template
        """

        import os, stat
        from cs_exec_environment import append_shell_shebang, \
            append_script_comment, prepend_path_command

        if not package:
            import cs_package
            package = cs_package.package()

        self.lines = []

        append_shell_shebang(self.lines)

        # Add batch system info if necessary

        batch_template = None
        config = configparser.ConfigParser()
        config.read(package.get_configfiles())

        if not ignore_batch and config.has_option('install', 'batch'):

            batch_template = config.get('install', 'batch')

            if not os.path.isabs(batch_template):
                batch_template = os.path.join(package.get_batchdir(),
                                             'batch.' + batch_template)

            fdt = open(batch_template, 'r')

            import re, string
            kwd1 = re.compile('nameandcase')

            # Determine or build default names if required

            if not case_name or not study_name:

                topdir, scriptdir = os.path.split(os.path.split(self.path)[0])
                if scriptdir == 'SCRIPTS':
                    studydir, casedir = os.path.split(topdir)
                    studydir = os.path.split(studydir)[1]
                else:
                    casedir = ''
                    studydir = scriptdir

                if not case_name:
                    if casedir:
                        case_name = casedir
                    else:
                        case_name = ''
                if not study_name:
                    study_name = studydir

            studycasename = study_name.lower() + case_name.lower()

            # For some systems, names are limited to 15 caracters
            studycasename = studycasename[:15]

            for line in fdt:
                line = line.rstrip()
                line = re.sub(kwd1, studycasename, line)
                self.lines.append(line)

            fdt.close()

        # Add command to execute.

        append_script_comment(self.lines, 'Ensure the correct command is found:')

        self.lines.append(prepend_path_command('PATH',
                                               package.get_dir("bindir")))
        self.lines.append('')
        append_script_comment(self.lines, 'Run command:\n')
        # On Linux systems, add a backslash to prevent aliases
        if sys.platform.startswith('win'):
            run_cmd = ''
        else:
            run_cmd = '\\'
        run_cmd += package.name + ' run'
        self.run_cmd_line_id = len(self.lines)
        self.lines.append(run_cmd)

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

    def set_compute_build(self, parameters):
        """
        Set the compute-build option in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = update_command_single_value(separate_args(line),
                                           ('--compute-build',
                                            '--compute-build='),
                                           enquote_arg(parameters))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

    #---------------------------------------------------------------------------

    def get_coupling(self):
        """
        Get the coupling option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--coupling',
                                         '--coupling='))

    #---------------------------------------------------------------------------

    def set_coupling(self, coupling):
        """
        Set the coupling option in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = update_command_single_value(separate_args(line),
                                           ('--coupling',),
                                           enquote_arg(coupling))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

    #---------------------------------------------------------------------------

    def get_parameters(self):
        """
        Get the parameters option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--param', '--param=', '-p'))

    #---------------------------------------------------------------------------

    def set_parameters(self, parameters):
        """
        Set the parameters option in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = update_command_single_value(separate_args(line),
                                           ('--param', '--param=', '-p'),
                                           enquote_arg(parameters))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

    #---------------------------------------------------------------------------

    def get_nprocs(self):
        """
        Get the nprocs option in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        return get_command_single_value(args,
                                        ('--nprocs', '--nprocs=', '-n'))

    #---------------------------------------------------------------------------

    def set_nprocs(self, parameters):
        """
        Set the nprocs option in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = update_command_single_value(separate_args(line),
                                           ('--nprocs', '--nprocs=', '-n'),
                                           enquote_arg(parameters))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

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

    def set_nthreads(self, parameters):
        """
        Set the nthreads option in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = update_command_single_value(separate_args(line),
                                           ('--threads-per-task',
                                            '--threads-per-task=',
                                            '-nt'),
                                           enquote_arg(parameters))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

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

    def set_run_id(self, run_id=None, run_id_prefix=None, run_id_suffix=None):
        """
        Set the run id, id_prefix, and id_suffix options in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = separate_args(line)
        if run_id != None:
            args = update_command_single_value(args,
                                               ('--id', '--id='),
                                               enquote_arg(run_id))
        if run_id_prefix != None:
            args = update_command_single_value(args,
                                               ('--id-prefix', '--id-prefix='),
                                               enquote_arg(run_id_prefix))
        if run_id_suffix != None:
            args = update_command_single_value(args,
                                               ('--id-suffix', '--id-suffix='),
                                               enquote_arg(run_id_suffix))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
