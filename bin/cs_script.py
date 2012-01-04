#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

#-------------------------------------------------------------------------------

class master_script:

    def __init__(self, argv = None, package = None):

        self.command = argv
        self.package = package

        self.commands = {'autovnv':self.autovnv,
                         'compile':self.compile,
                         'config':self.config,
                         'create':self.create,
                         'gui':self.gui,
                         'info':self.info,
                         'run':self.run,
                         'salome':self.salome}

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
            return self.commands[command](options)
        else:
            self.usage()
            sys.exit(1)

    def usage(self):

        usage = \
            """Usage: %(prog)s <topic>

Topics:
  help
  autovnv
  compile
  config
  create
  gui
  info
  run
  salome

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})

    def autovnv(self, options = None):
        import cs_autovnv
        return cs_autovnv.main(options, self.package)

    def compile(self, options = None):
        import cs_compile
        return cs_compile.main(options, self.package)

    def config(self, options = None):
        import cs_config
        return cs_config.main(options, self.package)

    def create(self, options = None):
        import cs_create
        return cs_create.main(options, self.package)

    def gui(self, options = None):
        import cs_gui
        return cs_gui.main(options, self.package)

    def info(self, options = None):
        import cs_info
        return cs_info.main(options, self.package)

    def run(self, options = None):
        import cs_run
        return cs_run.main(options, self.package)

    def salome(self, options = None):
        import cs_salome
        return cs_salome.main(options, self.package)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
