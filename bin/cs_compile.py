#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import fnmatch
import os
import sys
import tempfile

from optparse import OptionParser

from cs_config import build_syrthes

from cs_exec_environment import run_command, set_modules

#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string,
             it is split to create a list of arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-t", "--test", dest="test_mode",
                      action="store_true",
                      help="test only, discard compilation result")

    parser.add_option("-f", "--force", dest="force_link",
                      action="store_true",
                      help="force link, even with no source files")

    parser.add_option("-k", "--keep-going", dest="keep_going",
                      action="store_true",
                      help="continue even if errors are encountered")

    parser.add_option("-s", "--source", dest="src_dir", type="string",
                      metavar="<src_dir>",
                      help="choose source file directory")

    parser.add_option("-d", "--dest", dest="dest_dir", type="string",
                      metavar="<dest_dir>",
                      help="choose executable file directory")

    parser.add_option("--version", dest="version", type="string",
                      metavar="<version>",
                      help="select installed version")

    parser.add_option("--opt-libs", dest="opt_libs", type="string",
                      metavar="<libs>",
                      help="optional libraries")

    parser.add_option("--syrthes", dest="link_syrthes",
                      action="store_true",
                      help="SYRTHES link")

    parser.set_defaults(test_mode=False)
    parser.set_defaults(force_link=False)
    parser.set_defaults(keep_going=False)
    parser.set_defaults(src_dir=os.getcwd())
    parser.set_defaults(dest_dir=os.getcwd())
    parser.set_defaults(version="")
    parser.set_defaults(opt_libs="")
    parser.set_defaults(log_name=None)
    parser.set_defaults(link_syrthes=False)

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)

    src_dir = os.path.expanduser(options.src_dir)
    src_dir = os.path.expandvars(src_dir)
    src_dir = os.path.abspath(src_dir)
    if not os.path.isdir(src_dir):
        sys.stderr.write('Error: ' + src_dir + ' is not a directory')
        sys.exit(1)

    dest_dir = os.path.expanduser(options.dest_dir)
    dest_dir = os.path.expandvars(dest_dir)
    dest_dir = os.path.abspath(dest_dir)
    if not os.path.isdir(dest_dir):
        sys.stderr.write('Error: ' + dest_dir + ' is not a directory')
        sys.exit(1)

    return options.test_mode, options.force_link, options.keep_going, \
           src_dir, dest_dir, options.version, options.opt_libs, \
           options.link_syrthes

#-------------------------------------------------------------------------------

def so_dirs_path(flags, rpath):
    """
    Assemble path for shared libraries in nonstandard directories.
    """
    retval = ""
    first = True

    args = flags.split(" ")

    for arg in args:
        if arg[0:2] == '-L' and arg[0:10] != '-L/usr/lib' and arg[0:6] != '-L/lib':
            if first == True:
                retval = " " + rpath + ":" + arg[2:]
                first = False
            else:
                retval = retval + ":" + arg[2:]

    return retval

#-------------------------------------------------------------------------------

def compile_and_link(pkg, srcdir, destdir, optlibs,
                     force_link = False, keep_going = False,
                     stdout = sys.stdout, stderr = sys.stderr):
    """
    Compilation and link function.
    """
    retval = 0

    # Determine executable name

    exec_name = pkg.solver
    if destdir != None:
        exec_name = os.path.join(destdir, exec_name)

    # Change to temporary directory

    call_dir = os.getcwd()
    temp_dir = tempfile.mkdtemp(suffix=".cs_compile")
    os.chdir(temp_dir)

    # Find files to compile in source path

    dir_files = os.listdir(srcdir)

    c_files = fnmatch.filter(dir_files, '*.c')
    h_files = fnmatch.filter(dir_files, '*.h')
    cxx_files = fnmatch.filter(dir_files, '*.cxx') + fnmatch.filter(dir_files, '*.cpp')
    hxx_files = fnmatch.filter(dir_files, '*.hxx') + fnmatch.filter(dir_files, '*.hpp')
    f_files = fnmatch.filter(dir_files, '*.[fF]90')

    for f in c_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.cc
        if len(h_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " -I" + pkg.pkgincludedir
        cmd = cmd + " " + pkg.cppflags
        cmd = cmd + " " + pkg.cflags
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    for f in cxx_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.cxx
        if len(hxx_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " -I" + pkg.pkgincludedir
        cmd = cmd + " " + pkg.cppflags
        cmd = cmd + " " + pkg.cxxflags
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    user_mod_name = 'user_modules.f90'
    if user_mod_name in f_files:
        f_files.remove(user_mod_name)
        f_files.insert(0, user_mod_name)

    for f in f_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.fc
        cmd = cmd + " -I" + srcdir
        if pkg.fcmodinclude != "-I":
            cmd += " " + pkg.fcmodinclude + srcdir
        cmd = cmd + " -I" + pkg.fcmoddir
        if pkg.fcmodinclude != "-I":
            cmd += " " + pkg.fcmodinclude + pkg.fcmoddir
        cmd = cmd + " " + pkg.fcflags
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    if retval == 0 and (force_link or (len(c_files) + len(cxx_files) + len(f_files)) > 0):
        cmd = pkg.ld
        cmd = cmd + " -o " + exec_name
        if (len(c_files) + len(cxx_files) + len(f_files)) > 0:
          cmd = cmd + " *.o"
        cmd = cmd + " -L" + pkg.libdir
        if optlibs != None:
            if len(optlibs) > 0:
                cmd = cmd + " " + optlibs
        cmd = cmd + " " + pkg.ldflags + " " + pkg.libs
        cmd = cmd + " " + pkg.deplibs
        if pkg.rpath != "":
            cmd = cmd + " " + so_dirs_path(cmd, pkg.rpath)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    # Cleanup

    for f in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, f))

    # Return to original directory

    os.chdir(call_dir)
    os.rmdir(temp_dir)

    return retval

#-------------------------------------------------------------------------------

def compile_and_link_syrthes(pkg, srcdir, destdir,
                             stdout = sys.stdout, stderr = sys.stderr):
    """
    Compilation and link function.
    """
    retval = 0

    exec_name = "syrthes"
    if destdir != None:
        exec_name = os.path.join(destdir, exec_name)

    # Change to temporary directory

    call_dir = os.getcwd()
    temp_dir = tempfile.mkdtemp(suffix=".cs_compile_syrthes")
    os.chdir(temp_dir)

    # Find files to compile in source path

    if srcdir != None:
        dir_files = os.listdir(srcdir)
    else:
        dir_files = []

    c_files = fnmatch.filter(dir_files, '*.c')
    h_files = fnmatch.filter(dir_files, '*.h')
    f_files = fnmatch.filter(dir_files, '*.[fF]')

    for f in c_files:
        cmd = build_syrthes.cc
        if len(h_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " " + build_syrthes.cppflags
        cmd = cmd + " " + build_syrthes.cflags
        cmd = cmd + " " + pkg.ple_cppflags
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    for f in f_files:
        cmd = build_syrthes.fc
        if len(h_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " " + build_syrthes.cppflags
        cmd = cmd + " " + build_syrthes.fcflags
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    if retval == 0:
        # Link with Code_Saturne C compiler
        cmd = pkg.cc
        cmd = cmd + " -o " + exec_name
        if (len(f_files)) > 0:
          cmd = cmd + " *.o"
        cmd = cmd + " -L" + pkg.libdir + " " + pkg.ldflags
        cmd = cmd + " -lsyrcs"
        cmd = cmd + " " + build_syrthes.ldflags + " " + build_syrthes.libs
        if pkg.rpath != "":
            cmd = cmd + " " + so_dirs_path(cmd, pkg.rpath)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    # Cleanup

    for f in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, f))

    # Return to original directory

    os.chdir(call_dir)
    os.rmdir(temp_dir)

    return retval

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    test_mode, force_link, keep_going, src_dir, dest_dir, \
        version, opt_libs, link_syrthes = process_cmd_line(argv, pkg)

    if (version):
        pkg = pkg.get_alternate_version(version)

    set_modules(pkg) # Set environment modules if present

    if test_mode == True:
        dest_dir = None

    if link_syrthes == True:
        retcode = compile_and_link_syrthes(pkg, src_dir, dest_dir)
    else:
        retcode = compile_and_link(pkg, src_dir, dest_dir, opt_libs,
                                   force_link, keep_going)

    sys.exit(retcode)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
