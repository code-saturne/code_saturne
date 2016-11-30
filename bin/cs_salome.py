#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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

"""
This module describes the script used to launch SALOME with the CFDSTUDY module.

This module defines the following functions:
- process_cmd_line
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys

from optparse import OptionParser

import cs_config
import cs_exec_environment
from cs_config import config

#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------


def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    (options, args) = parser.parse_args(argv)

    return


#-------------------------------------------------------------------------------
# Launch SALOME platform with CFDSTUDY module
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    template = """\
%(salomeenv)s;
CFDSTUDY_ROOT_DIR=%(prefix)s;
PYTHONPATH=%(pythondir)s/salome:%(pkgpythondir)s${PYTHONPATH:+:$PYTHONPATH};
export CFDSTUDY_ROOT_DIR PYTHONPATH;
%(runsalome)s --modules=%(modules)s
"""

    cfg = config()

    if cfg.have_salome == "no":
        sys.stderr.write("SALOME is not available in this installation.\n")
        sys.exit(1)

    # Skipped modules (version 6.3.0): YACS,JOBMANAGER,HOMARD,OPENTURNS
    default_modules = "GEOM,SMESH,MED,CFDSTUDY,PARAVIS"

    run_cmd = cfg.salome_run
    if not run_cmd:
        run_cmd = "${KERNEL_ROOT_DIR}/bin/salome/envSalome.py python ${KERNEL_ROOT_DIR}/bin/salome/runSalome.py"

    path = pkg.get_dir('pkgpythondir')
    if pkg.name == 'neptune_cfd':
        cspath = os.path.join(pkg.get_dir('csdir'), "lib",
                                                   "python" + sys.version[:3],
                                                   "site-packages", "code_saturne")
        path = path+":"+cspath

        # Test if EOS modules could be imported
        cfg = cs_config.config()
        if cfg.libs['eos'].have == "yes":
            eosprefix = cfg.libs['eos'].prefix
            try:
                from distutils import sysconfig
                eospath = os.path.join(sysconfig.get_python_lib(0, 0, prefix=eosprefix), 'eos')
            except Exception:
                eospath = ''
            if eospath:
                if os.path.isdir(eospath) and not eospath in sys.path:
                    path = path+":"+eospath

    cmd = template % {'salomeenv': cfg.salome_env,
                      'prefix': pkg.get_dir('prefix'),
                      'pythondir': pkg.get_dir('pythondir'),
                      'pkgpythondir': path,
                      'runsalome': run_cmd,
                      'modules': default_modules}

    process_cmd_line(argv, pkg)

    retcode = cs_exec_environment.run_command(cmd,
                                              stdout=None,
                                              stderr=None)


if __name__ == "__main__":
    main(sys.argv[1:], None)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
