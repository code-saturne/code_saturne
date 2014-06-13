#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

import sys

from cs_exec_environment import separate_args, assemble_args, enquote_arg, \
    get_command_single_value, update_command_single_value

#===============================================================================
# Class used to manage runcase files
#===============================================================================

class runcase(object):

    def __init__(self, path):
        """
        Initialize runcase info object.
        """
        self.path = path

        try:
            f = open(self.path, mode = 'r')
        except IOError:
            print("Error: can not open %s\n" % self.path)
            sys.exit(1)

        self.lines = f.readlines()
        f.close()

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
            if line[0] == '#':
                continue
            j = line.find('#')
            if j > -1:
                line = line[0:j]

            args = separate_args(line)
            if args.count('run') == 1:
                if args.index('run') == 1: # "<package_name> run"
                    for name in ('code_saturne', 'neptune_cfd'):
                        if not sys.platform.startswith('win'):
                            test_name = '\\' + name
                        else:
                            test_name = name
                        if args[0].find(test_name) == 0:
                            self.cmd_name = name
                            self.run_cmd_line_id = i
                            return

        # We should have exited before reaching this.

        err_str = "Error: unable to determine the name of the script for " + self.path
        raise ValueError(err_str)

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

    def get_run_id(self):
        """
        Get the run id, id_prefix, and id_suffix options in the run command
        """

        args = separate_args(self.lines[self.run_cmd_line_id])

        run_id = get_command_single_value(args, ('--id',))
        run_id_prefix = get_command_single_value(args, ('--id-prefix',))
        run_id_suffix = get_command_single_value(args, ('--id-suffix',))

        return run_id, run_id_prefix, run_id_suffix

    #---------------------------------------------------------------------------

    def set_run_id(self, run_id=None, run_id_prefix=None, run_id_suffix=None):
        """
        Set the run id, id_prefix, and id_suffix options in the run command
        """

        line = self.lines[self.run_cmd_line_id]

        args = separate_args(line)
        if run_id:
            args = update_command_single_value(args,
                                               ('--id',),
                                               enquote_arg(run_id))
        if run_id_prefix:
            args = update_command_single_value(args,
                                               ('--id-prefix',),
                                               enquote_arg(run_id_prefix))
        if run_id_suffix:
            args = update_command_single_value(args,
                                               ('--id-suffix',),
                                               enquote_arg(run_id_suffix))

        self.lines[self.run_cmd_line_id] = assemble_args(args)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
