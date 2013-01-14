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

def get_compiler(pkg, compiler):

    # First, handle the standard "non relocatable" case
    if pkg.config.features['relocatable'] == "no":
        return pkg.config.compilers[compiler]

    # On Windows, compilers are installed in "bindir" for DLL search
    if sys.platform.startswith("win"):
        bindir = os.path.dirname(pkg.scriptdir)
        return os.path.join(bindir, pkg.config.compilers[compiler])

    # On Linux systems, one assumes compilers are installed in the system
    else:
        return pkg.config.compilers[compiler]

#-------------------------------------------------------------------------------

def get_flags(pkg, flag):

    cmd_line = []

    # Build the command line, and split possible multiple arguments
    # in lists, assuming at this step no white-space in paths or names
    for lib in pkg.config.deplibs:
        if (pkg.config.libs[lib].have == "yes"
            and (not pkg.config.libs[lib].dynamic_load)):
            cmd_line += pkg.config.libs[lib].flags[flag].split()

    # Add CPPFLAGS and LDFLAGS information for the current package
    if flag == 'cppflags':
        cmd_line.insert(0, "-I" + pkg.get_dir("pkgincludedir"))
    elif flag == 'ldflags':
        cmd_line.insert(0, "-L" + pkg.get_dir("libdir"))

    # Add CPPFLAGS information when PLE is "internal"
    if pkg.config.libs['ple'].variant == "internal" and flag == 'cppflags':
        cmd_line.insert(0, "-I" + pkg.get_dir("includedir"))

    # On Windows, flags must be adapted so as to handle the relocation
    # of system headers (together with the compiler)
    if sys.platform.startswith("win"):
        for i in range(len(cmd_line)):
            s = cmd_line[i]
            for p in ["-I/mingw", "-Ic:/mingw"]:
                s_tmp = s.replace(p, "-I" + pkg.get_dir("prefix"))
            for p in ["-L/mingw", "-Lc:/mingw"]:
                s_tmp = s.replace(p, "-L" + pkg.get_dir("prefix"))
            cmd_line[i] = s_tmp

    return cmd_line

#-------------------------------------------------------------------------------

def so_dirs_path(flags, pkg):
    """
    Assemble path for shared libraries in nonstandard directories.
    """
    retval = pkg.config.rpath.split(" ")
    count = 0

    pkg_lib = os.path.join(pkg.get_dir('libdir'), pkg.name)
    if os.path.isdir(pkg_lib):
        retval[1] +=  ":" + pkg_lib
        count += 1

    if type(flags) == str:
        args = flags.split(" ")
    else:
        args = flags

    for arg in args:
        if arg[0:2] == '-L' and arg[0:10] != '-L/usr/lib' and arg[0:6] != '-L/lib':
            retval[1] += ":" + arg[2:]
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

    p_libs = get_flags(pkg, 'libs')
    if pkg.config.special_user_link == 'ar_x':
        if force_link or (len(c_files) + len(cxx_files) + len(f_files)) > 0:
            lib0 = os.path.join(pkg.get_dir('libdir'),
                                'lib' + p_libs[0][2:] + '.a')
            p_libs = [lib0] + p_libs[1:]
            cmd = ['ar', 't', lib0]
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1

    # Compile files

    for f in c_files:
        if (retval != 0 and not keep_going):
            break
        cmd = [get_compiler(pkg, 'cc')]
        if opt_cflags != None:
            cmd = cmd + opt_cflags.split()
        if len(h_files) > 0:
            cmd = cmd + ["-I", srcdir]
        cmd.append('-DHAVE_CONFIG_H')
        if f == 'cs_base.c':
            cmd = cmd + ['-DLOCALEDIR=\\"' + pkg.get_dir('localedir') + '\\"', \
                         '-DPKGDATADIR=\\"' + pkg.get_dir('pkgdatadir') + '\\"']
        cmd = cmd + get_flags(pkg, 'cppflags')
        cmd = cmd + pkg.config.flags['cflags'].split()
        cmd = cmd + ["-c", os.path.join(srcdir, f)]
        if run_command(cmd, pkg=pkg, echo=True,
                       stdout=stdout, stderr=stderr) != 0:
            retval = 1

    for f in cxx_files:
        if (retval != 0 and not keep_going):
            break
        cmd = [get_compiler(pkg, 'cxx')]
        if opt_cxxflags != None:
            cmd = cmd + opt_cxxflags.split()
        if len(hxx_files) > 0:
            cmd = cmd + ["-I", srcdir]
        cmd = cmd + get_flags(pkg, 'cppflags')
        cmd = cmd + pkg.config.flags['cxxflags'].split()
        cmd = cmd + ["-c", os.path.join(srcdir, f)]
        if run_command(cmd, pkg=pkg, echo=True,
                       stdout=stdout, stderr=stderr) != 0:
            retval = 1

    user_mod_name = 'cs_user_modules.f90'
    if user_mod_name in f_files:
        f_files.remove(user_mod_name)
        f_files.insert(0, user_mod_name)

    for f in f_files:
        if (retval != 0 and not keep_going):
            break
        cmd = [get_compiler(pkg, 'fc')]
        if opt_fcflags != None:
            cmd = cmd + opt_fcflags.split()
        cmd = cmd + ["-I", srcdir]
        if pkg.config.fcmodinclude != "-I":
            cmd = cmd + [pkg.config.fcmodinclude, srcdir]
        cmd = cmd + ["-I", pkg.get_dir('pkgincludedir')]
        if pkg.config.fcmodinclude != "-I":
            cmd = cmd + [pkg.config.fcmodinclude, pkg.get_dir('pkgincludedir')]
        cmd = cmd + pkg.config.flags['fcflags'].split()
        cmd = cmd + ["-c", os.path.join(srcdir, f)]
        if run_command(cmd, pkg=pkg, echo=True,
                       stdout=stdout, stderr=stderr) != 0:
            retval = 1

    if retval == 0 and (force_link or (len(c_files) + len(cxx_files) + len(f_files)) > 0):
        cmd = [get_compiler(pkg, 'ld')]
        cmd = cmd + ["-o", exec_name]
        if (len(c_files) + len(cxx_files) + len(f_files)) > 0:
            dir_files = os.listdir(temp_dir)
            o_files = fnmatch.filter(dir_files, '*.o')
            cmd = cmd + o_files
        if opt_libs != None:
            if len(opt_libs) > 0:
                cmd.append(opt_libs)
        cmd = cmd + get_flags(pkg, 'ldflags')
        cmd = cmd + p_libs
        if pkg.config.rpath != "":
            cmd += so_dirs_path(cmd, pkg)
        if run_command(cmd, pkg=pkg, echo=True,
                       stdout=stdout, stderr=stderr) != 0:
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
