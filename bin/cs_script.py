#!/usr/bin/env python

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

import os
import sys

#-------------------------------------------------------------------------------

class master_script:

    def __init__(self, argv = None, package = None):

        self.command = argv
        self.package = package

        self.commands = {'studymanager':self.studymanager,
                         'smgr':self.studymanager,
                         'bdiff':self.bdiff,
                         'bdump':self.bdump,
                         'compile':self.compile,
                         'config':self.config,
                         'create':self.create,
                         'gui':self.gui,
                         'studymanagergui':self.studymanager_gui,
                         'smgrgui':self.studymanager_gui,
                         'trackcvg':self.trackcvg,
                         'info':self.info,
                         'run':self.run,
                         'salome':self.salome,
                         'submit':self.submit}

        if package != None:
            sys.path.insert(1, package.get_dir('pythondir'))

    def execute(self):

        help_commands = ("help", "--help", "-h")
        vers_commands = ("--version", "-v")

        if len(self.command) < 1:
            # On Windows platform, we have two executables frozen by cx_freeze:
            # a .com designed for a console mode (it can also be accessed
            # without the extension) and a .exe designed for a GUI mode.
            # At the moment, running the main script does not launch the GUI
            # mode. Therefore, we force the 'gui' option in this case.
            # This can be changed later on.
            if sys.platform.startswith('win') and \
                    os.path.basename(sys.argv[0]).endswith('exe'):
                command = 'gui'
            else:
                self.usage()
                sys.exit(0)
        else:
            command = self.command[0]

        if command in help_commands:
            subcommand = False
            if len(self.command) > 1:
                subcommand = True
                command = self.command[1]
                if self.command[1] in self.commands:
                    options = self.command[1:] + ['--help']
                    return self.commands[command](options)
            if not subcommand:
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
  studymanager
  smgr
  bdiff
  bdump
  compile
  config
  create
  gui
  studymanagergui
  smgrgui
  trackcvg
  info
  run
  salome
  submit

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})

    def studymanager(self, options = None):
        import cs_studymanager
        return cs_studymanager.main(options, self.package)

    def bdiff(self, options = None):
        import cs_bdiff
        return cs_bdiff.main(options, self.package)

    def bdump(self, options = None):
        import cs_bdump
        return cs_bdump.main(options, self.package)

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

    def studymanager_gui(self, options = None):
        import cs_studymanager_gui
        return cs_studymanager_gui.main(options, self.package)

    def trackcvg(self, options = None):
        import cs_trackcvg
        return cs_trackcvg.main(options, self.package)

    def info(self, options = None):
        import cs_info
        return cs_info.main(options, self.package)

    def run(self, options = None):
        import cs_run
        return cs_run.main(options, self.package)

    def salome(self, options = None):
        import cs_salome
        return cs_salome.main(options, self.package)

    def submit(self, options = None):
        import cs_submit
        return cs_submit.main(options, self.package)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
