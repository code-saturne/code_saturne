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

"""
This module describes the script used to update a code_saturne case
links and paths, as well as user functions examples.

This module defines the following functions:
- process_cmd_line
- set_executable
- unset_executable
- update_case
- main
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, shutil, stat
import types, string, re, fnmatch
from optparse import OptionParser

from code_saturne.base import cs_exec_environment
from code_saturne.base import cs_runcase
from code_saturne.base import cs_run_conf
from code_saturne.base.cs_case import get_case_dir
from code_saturne.base.cs_create import set_executable, unset_executable
from code_saturne.base.cs_create import create_local_launcher

#-------------------------------------------------------------------------------
# Process the passed command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog create [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-c", "--case", dest="case_names", type="string",
                      metavar="<case>", action="append",
                      help="case to update")

    parser.add_option("-q", "--quiet",
                      action="store_const", const=0, dest="verbose",
                      help="do not output any information")

    parser.add_option("-v", "--verbose",
                      action="store_const", const=2, dest="verbose",
                      help="dump study creation parameters")

    parser.set_defaults(case_names=[])
    parser.set_defaults(verbose=1)

    (options, args) = parser.parse_args(argv)

    return (options, args)

#-------------------------------------------------------------------------------
# Copy a directory without raising an error if the destination allready exists
#
# If ref is set to true, ignore the directory if the destination does not
# exist, so that if the case was created with a "--noref" option, we  can
# avoid adding references or examples in an update.
# In this case, previous files are also removed, as they may be obsolete.
#-------------------------------------------------------------------------------

def copy_directory(src, dest, ref=False):

    if not os.path.isdir(dest):
        if ref == False:
            shutil.copytree(src, dest)

    else:
        if ref:
            for itm in os.listdir(dest):
                os.remove(os.path.join(dest, itm))

        for itm in os.listdir(src):
            shutil.copy2(os.path.join(src, itm), os.path.join(dest, itm))

#-------------------------------------------------------------------------------
# Update the case paths and examples/reference files
#-------------------------------------------------------------------------------

def update_case(options, pkg):

    topdir = os.getcwd()
    study_name = os.path.basename(os.getcwd())

    i_c = cs_run_conf.get_install_config_info(pkg)
    resource_name = cs_run_conf.get_resource_name(i_c)

    for case in options.case_names:
        os.chdir(topdir)

        if case == ".":
            case, staging_dir = get_case_dir()
            if not case:
                sys.stderr.write("  o Skipping '%s', which  does not seem "
                                 "to be a case directory\n" % topdir)
                continue
            casename = os.path.basename(case)
        else:
            casename = case

        if options.verbose > 0:
            sys.stdout.write("  o Updating case '%s' paths...\n" % casename)

        datadir = os.path.join(pkg.get_dir("pkgdatadir"))

        os.chdir(case)

        # Write a local wrapper to main command
        data = 'DATA'
        if not os.path.isdir(data):
            os.mkdir(data)

        dataref_distpath = os.path.join(datadir, 'user')
        user = os.path.join(data, 'REFERENCE')

        # Only update user_scripts reference, not data
        # (we should try to deprecate the copying of reference data
        # or use the GUI to align it with the active options)
        if os.path.exists(user):
            abs_f = os.path.join(datadir, 'data', 'user', 'cs_user_scripts.py')
            shutil.copy(abs_f, user)
            unset_executable(user)

        for s in ("SaturneGUI", "NeptuneGUI"):
            old_gui_script = os.path.join(data, s)
            if os.path.isfile(old_gui_script):
                os.remove(old_gui_script)

        # Rebuild launch script
        create_local_launcher(pkg, data)

        # User source files directory
        src = 'SRC'
        if not os.path.isdir(src):
            os.mkdir(src)

        user_ref_distpath = os.path.join(datadir, 'user_sources')
        for srcdir in ('REFERENCE', 'EXAMPLES', 'EXAMPLES_neptune_cfd'):
            if os.path.isdir(os.path.join(user_ref_distpath, srcdir)):
                copy_directory(os.path.join(user_ref_distpath, srcdir),
                               os.path.join(src, srcdir),
                               True)

            unset_executable(os.path.join(src, srcdir))

        # Results directory (only one for all instances)

        resu = 'RESU'
        if not os.path.isdir(resu):
            os.mkdir(resu)

        # Script directory (only one for all instances)

        run_conf_file = os.path.join(topdir, case, 'DATA', 'run.cfg')
        batch_file = os.path.join(topdir, case, 'SCRIPTS', 'runcase')
        if sys.platform.startswith('win'):
            batch_file = batch_file + '.bat'

        run_conf = cs_run_conf.run_conf(run_conf_file,
                                        package=pkg,
                                        rebuild=True)

        if os.path.isfile(batch_file):
            runcase = cs_runcase.runcase(batch_file,
                                         package=pkg)
            sections = runcase.run_conf_sections(resource_name=resource_name,
                                                 batch_template=i_c['batch'])
            for sn in sections:
                if not sn in run_conf.sections:
                    run_conf.sections[sn] = {}
                for kw in sections[sn]:
                    run_conf.sections[sn][kw] = sections[sn][kw]
            os.remove(batch_file)
            scripts_dir = os.path.join(topdir, case, 'SCRIPTS')
            try:
                os.rmdir(scripts_dir)
            except Exception:
                pass

        run_conf.save()

#-------------------------------------------------------------------------------
# Main function
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    welcome = """\
%(name)s %(vers)s case update
"""
    opts, args = process_cmd_line(argv, pkg)

    if opts.case_names == []:
        if len(args) > 0:
            opts.case_names = args
        else:
            # Value is set to "." instead of "" to avoid chdir errors due
            # to non existant paths
            opts.case_names = ["."]

    if opts.verbose > 0:
        sys.stdout.write(welcome % {'name':pkg.name, 'vers':pkg.version})

    update_case(opts, pkg)

if __name__=="__main__":
    main(sys.argv[1:], None)
