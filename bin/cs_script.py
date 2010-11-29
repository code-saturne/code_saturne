#!/usr/bin/env python
#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009-2010 EDF
#
#   The Code_Saturne Preprocessor is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Preprocessor is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import sys

#-------------------------------------------------------------------------------

class master_script:

    def __init__(self, argv = None, package = None):

        self.command = argv
        self.package = package

        self.commands = {'check_consistency':self.check_consistency,
                         'check_mesh':self.check_mesh,
                         'compile':self.compile,
                         'config':self.config,
                         'create':self.create,
                         'gui':self.gui,
                         'info':self.info,
                         'plot_probes':self.plot_probes}

    def execute(self):

        help_commands = ("help", "--help", "-h")
        vers_commands = ("--version", "-v")

        if (len(self.command) < 1):
            self.usage()
            sys.exit(0)

        command = self.command[0]

        if command in help_commands:
            self.usage()
            sys.exit(0)

        if command in self.commands:
            options = self.command[1:]
            self.commands[command](options)
        else:
            self.usage()
            sys.exit(1)

    def usage(self):

        usage = \
            """Usage: %(prog)s <topic>

Topics:
  help
  check_consistency
  check_mesh
  compile
  config
  create
  gui
  info
  plot_probes

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})

    def check_consistency(self, options = None):
        import cs_check_consistency
        cs_check_consistency.main(options, self.package)

    def check_mesh(self, options = None):
        import cs_check_mesh
        cs_check_mesh.main(options, self.package)

    def compile(self, options = None):
        import cs_compile
        cs_compile.main(options, self.package)

    def config(self, options = None):
        import cs_config
        cs_config.main(options, self.package)

    def create(self, options = None):
        import cs_create
        cs_create.main(options, self.package)

    def gui(self, options = None):
        import cs_gui
        cs_gui.main(options, self.package)

    def info(self, options = None):
        import cs_info
        cs_info.main(options, self.package)

    def plot_probes(self, options = None):
        import cs_plot_probes
        cs_plot_probes.main(options, self.package)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
