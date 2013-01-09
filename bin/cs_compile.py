#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import fnmatch
import os
import sys
import tempfile

from optparse import OptionParser

from cs_exec_environment import run_command

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

    parser.add_option("--cflags", dest="cflags", type="string",
                      metavar="<cflags>",
                      help="additional C compiler and preprocessor flags")

    parser.add_option("--cxxflags", dest="cxxflags", type="string",
                      metavar="<cxxflags>",
                      help="additional C++ compiler and preprocessor flags")

    parser.add_option("--fcflags", dest="fcflags", type="string",
                      metavar="<fcflags>",
                      help="additional Fortran compiler flags")

    parser.add_option("--libs", dest="libs", type="string",
                      metavar="<libs>",
                      help="additional libraries")

    parser.set_defaults(test_mode=False)
    parser.set_defaults(force_link=False)
    parser.set_defaults(keep_going=False)
    parser.set_defaults(src_dir=os.getcwd())
    parser.set_defaults(dest_dir=os.getcwd())
    parser.set_defaults(version="")
    parser.set_defaults(cflags=None)
    parser.set_defaults(cxxflags=None)
    parser.set_defaults(fccflags=None)
    parser.set_defaults(libs=None)
    parser.set_defaults(log_name=None)

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
           src_dir, dest_dir, options.version, options.cflags, \
           options.cxxflags, options.fcflags, options.libs

#-------------------------------------------------------------------------------

def so_dirs_path(flags, pkg):
    """
    Assemble path for shared libraries in nonstandard directories.
    """
    retval = " " + pkg.rpath
    count = 0

    pkg_lib = os.path.join(pkg.get_dir('libdir'), pkg.name)
    if os.path.isdir(pkg_lib):
        retval = retval + ":" + pkg_lib
        count += 1

    args = flags.split(" ")

    for arg in args:
        if arg[0:2] == '-L' and arg[0:10] != '-L/usr/lib' and arg[0:6] != '-L/lib':
            retval = retval + ":" + arg[2:]
            count += 1

    if count == 0:
        retval = ""

    return retval

#-------------------------------------------------------------------------------

def compile_and_link(pkg, srcdir, destdir,
                     opt_cflags=None, opt_cxxflags=None, opt_fcflags=None,
                     opt_libs=None, force_link=False, keep_going=False,
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

    # Special handling for some linkers (such as Mac OS X), for which
    # no multiple definitions are allowable in static mode;
    # in this case, extract archive, then overwrite with user files.

    p_libs = pkg.get_flags('libs')
    if pkg.special_user_link == 'ar_x':
        if force_link or (len(c_files) + len(cxx_files) + len(f_files)) > 0:
            i = p_libs.find(' ')
            if (i > 0):
                lib0 = os.path.join(pkg.get_dir('libdir'), 'lib' + p_libs[2:i] + '.a')
                p_libs = p_libs[i+1:]
            else:
                lib0 = os.path.join(pkg.get_dir('libdir'), 'lib' + p_libs[2:] + '.a')
                p_libs = ''
            cmd = 'ar x ' + lib0
            if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
                retval = 1

    # Compile files

    for f in c_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.get_compiler('cc')
        if opt_cflags != None:
            cmd = cmd + " " + opt_cflags
        if len(h_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " -I" + pkg.get_dir('pkgincludedir')
        if f == 'cs_base.c':
            cmd = cmd + ' -DLOCALEDIR=\\"' + pkg.get_dir('localedir') \
                      + '\\" -DPKGDATADIR=\\"' + pkg.get_dir('pkgdatadir') + '\\"'
        cmd = cmd + " " + pkg.get_flags('cppflags')
        cmd = cmd + " " + pkg.get_flags('cflags')
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    for f in cxx_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.get_compiler('cxx')
        if opt_cxxflags != None:
            cmd = cmd + " " + opt_cxxflags
        if len(hxx_files) > 0:
            cmd = cmd + " -I" + srcdir
        cmd = cmd + " -I" + pkg.get_dir('pkgincludedir')
        cmd = cmd + " " + pkg.get_flags('cppflags')
        cmd = cmd + " " + pkg.get_flags('cxxflags')
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    user_mod_name = 'cs_user_modules.f90'
    if user_mod_name in f_files:
        f_files.remove(user_mod_name)
        f_files.insert(0, user_mod_name)

    for f in f_files:
        if (retval != 0 and not keep_going):
            break
        cmd = pkg.get_compiler('fc')
        if opt_fcflags != None:
            cmd = cmd + " " + opt_fcflags
        cmd = cmd + " -I" + srcdir
        if pkg.fcmodinclude != "-I":
            cmd += " " + pkg.fcmodinclude + srcdir
        cmd = cmd + " -I" + pkg.get_flags('fcmoddir')
        if pkg.fcmodinclude != "-I":
            cmd += " " + pkg.fcmodinclude + pkg.get_flags('fcmoddir')
        cmd = cmd + " " + pkg.get_flags('fcflags')
        cmd = cmd + " -c " + os.path.join(srcdir, f)
        if run_command(cmd, echo=True, stdout=stdout, stderr=stderr) != 0:
            retval = 1

    if retval == 0 and (force_link or (len(c_files) + len(cxx_files) + len(f_files)) > 0):
        cmd = pkg.get_compiler('ld')
        cmd = cmd + " -o " + exec_name
        if (len(c_files) + len(cxx_files) + len(f_files)) > 0:
          cmd = cmd + " *.o"
        cmd = cmd + " -L" + pkg.get_dir('libdir')
        if opt_libs != None:
            if len(opt_libs) > 0:
                cmd = cmd + " " + opt_libs
        cmd = cmd + " " + pkg.get_flags('ldflags') + " " + p_libs
        cmd = cmd + " " + pkg.get_flags('deplibs')
        if pkg.rpath != "":
            cmd = cmd + " " + so_dirs_path(cmd, pkg)
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

    from cs_exec_environment import set_modules, source_rcfile

    test_mode, force_link, keep_going, src_dir, dest_dir, \
        version, cflags, cxxflags, fcflags, libs = process_cmd_line(argv, pkg)

    if (version):
        pkg = pkg.get_alternate_version(version)

    set_modules(pkg)    # Set environment modules if present
    source_rcfile(pkg)  # Source rcfile if defined

    if test_mode == True:
        dest_dir = None

    retcode = compile_and_link(pkg, src_dir, dest_dir, cflags, cxxflags,
                               fcflags, libs, force_link, keep_going)

    sys.exit(retcode)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
