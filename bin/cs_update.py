#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module describes the script used to update a Code_Saturne case
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

from __future__ import print_function

import os, sys, shutil, stat
import types, string, re, fnmatch
from optparse import OptionParser
try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

from code_saturne import cs_exec_environment
from code_saturne import cs_runcase

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
# Assign executable mode (chmod +x) to a file
#-------------------------------------------------------------------------------

def set_executable(filename):
    """
    Give executable permission to a given file.
    Equivalent to `chmod +x` shell function.
    """

    st   = os.stat(filename)
    mode = st[stat.ST_MODE] | stat.S_IXUSR
    if mode & stat.S_IRGRP:
        mode = mode | stat.S_IXGRP
    if mode & stat.S_IROTH:
        mode = mode | stat.S_IXOTH
    os.chmod(filename, mode)

    return

#-------------------------------------------------------------------------------
# Remove executable mode (chmod -x) from a file or files inside a directory
#-------------------------------------------------------------------------------

def unset_executable(path):
    """
    Remove executable permission from a given file or files inside a directory.
    Equivalent to `chmod -x` shell function.
    """

    if os.path.isfile(path):
        try:
            st   = os.stat(path)
            mode = st[stat.ST_MODE] | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            mode = mode ^ (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            os.chmod(path, mode)
        except Exception:
            pass

    elif os.path.isdir(path):
        l = os.listdir(path)
        for p in l:
            unset_executable(p)

    return

#-------------------------------------------------------------------------------
# Copy a directory without raising an error if the destination allready exists
#-------------------------------------------------------------------------------

def copy_directory(src, dest):

    if not os.path.isdir(dest):
        shutil.copytree(src, dest)

    else:
        for itm in os.listdir(src):
            shutil.copy2(os.path.join(src, itm), os.path.join(dest, itm))

    return

#-------------------------------------------------------------------------------
# Update the case paths and examples/reference files
#-------------------------------------------------------------------------------

def update_case(options, pkg):

    repbase = os.getcwd()
    study_name=os.path.basename(os.getcwd())
    for case in options.case_names:
        os.chdir(repbase)

        if case == ".":
            casename = os.path.split(repbase)[-1]
        else:
            casename = case

        if options.verbose > 0:
            sys.stdout.write("  o Updating case '%s' paths...\n" % casename)

        datadir = os.path.join(pkg.get_dir("pkgdatadir"))

        os.chdir(case)

        # Write a wrapper for GUI launching
        data = 'DATA'
        if not os.path.isdir(data):
            os.mkdir(data)

        guiscript = os.path.join(data, pkg.guiname)

        fd = open(guiscript, 'w')
        cs_exec_environment.write_shell_shebang(fd)

        cs_exec_environment.write_script_comment(fd,
            'Ensure the correct command is found:\n')
        cs_exec_environment.write_prepend_path(fd, 'PATH',
                                               pkg.get_dir("bindir"))
        fd.write('\n')
        cs_exec_environment.write_script_comment(fd, 'Run command:\n')
        # On Linux systems, add a backslash to prevent aliases
        if sys.platform.startswith('linux'): fd.write('\\')
        fd.write(pkg.name + ' gui ' +
                 cs_exec_environment.get_script_positional_args() + '\n')

        fd.close()

        set_executable(guiscript)

        # User source files directory
        src = 'SRC'
        if not os.path.isdir(src):
            os.mkdir(src)

        user_distpath = os.path.join(datadir, 'user')
        user_examples_distpath = os.path.join(datadir, 'user_examples')

        user = os.path.join(src, 'REFERENCE')
        user_examples = os.path.join(src, 'EXAMPLES')

        copy_directory(user_distpath, user)
        copy_directory(user_examples_distpath, user_examples)

        add_datadirs = []
        if pkg.name == 'neptune_cfd' :
            add_datadirs.append(os.path.join(pkg.get_dir("datadir"),
                                             pkg.name))

        for d in add_datadirs:
            user_distpath = os.path.join(d, 'user')
            if os.path.isdir(user_distpath):
                copy_directory(user_distpath, user)

            user_examples_distpath = os.path.join(d, 'user_examples')
            if os.path.isdir(user_examples_distpath):
                copy_directory(user_examples_distpath, user_examples)


        unset_executable(user)
        unset_executable(user_examples)

        # Results directory (only one for all instances)

        resu = 'RESU'
        if not os.path.isdir(resu):
            os.mkdir(resu)

        # Script directory (only one for all instances)

        scripts = 'SCRIPTS'
        if not os.path.isdir(scripts):
            os.mkdir(scripts)

        batch_file = os.path.join(repbase, case, scripts, 'runcase')
        if sys.platform.startswith('win'):
            batch_file = batch_file + '.bat'

        # Add info from parent in case of copy

        runcase = cs_runcase.runcase(batch_file,
                                     package=pkg,
                                     rebuild=True,
                                     study_name=study_name,
                                     case_name=case)
        runcase.save()



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
