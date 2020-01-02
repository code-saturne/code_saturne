#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

try:
    from code_saturne.cs_exec_environment import run_command, separate_args
except Exception:
    from cs_exec_environment import run_command, separate_args

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

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog compile [options]"

    parser = OptionParser(usage=usage)

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

def files_to_compile(src_dir):
    """
    Return files to compile in source path (empty list if none are found)
    """

    dir_files = os.listdir(src_dir)

    src_files = (fnmatch.filter(dir_files, '*.c')
                 + fnmatch.filter(dir_files, '*.cxx')
                 + fnmatch.filter(dir_files, '*.cpp')
                 + fnmatch.filter(dir_files, '*.[fF]90'))

    return src_files

#===============================================================================
# Class used to manage compilation
#===============================================================================

class cs_compile(object):

    def __init__(self,
                 package=None):
        """
        Initialize compiler object.
        """
        self.pkg = package

    #---------------------------------------------------------------------------

    def get_compiler(self, compiler):
        """
        Determine the compiler path for a given compiler type.
        """

        # First, handle the standard "non relocatable" case
        if self.pkg.config.features['relocatable'] == "no":
            return self.pkg.config.compilers[compiler]

        # On Windows, compilers are installed in "bindir" for DLL search
        if sys.platform.startswith("win"):
            return os.path.join(self.pkg.get_dir('bindir'),
                                self.pkg.config.compilers[compiler])

        # On Linux systems, one assumes compilers are installed in the system
        else:
            return self.pkg.config.compilers[compiler]

    #---------------------------------------------------------------------------

    def flags_relocation(self, flag, cmd_line):
        """
        Adjust flags for relocation on some systems.
        """

        # On MinGW hosts, flags must be adapted so as to handle the relocation
        # of system headers (together with the compiler)
        # ...unless the compilation is done during the build stage
        if sys.platform.startswith("win"):
            for i in range(len(cmd_line)):
                s = cmd_line[i]
                for p in ["-I/mingw", "-Ic:/mingw"]:
                    s_tmp = s.replace(p, "-I" + self.pkg.get_dir("prefix"))
                for p in ["-L/mingw", "-Lc:/mingw"]:
                    s_tmp = s.replace(p, "-L" + self.pkg.get_dir("prefix"))
                cmd_line[i] = s_tmp

    #---------------------------------------------------------------------------

    def get_pkg_path_flags(self, flag):
        """
        Determine compilation flags for a given flag type.
        """

        flags = []

        # Add CPPFLAGS and LDFLAGS information for the current package
        if flag == 'cppflags':
            pkgincludedir = self.pkg.get_dir("pkgincludedir")
            flags.insert(0, "-I" + pkgincludedir)
            if self.pkg.config.libs['ple'].variant == "internal":
                flags.insert(0, "-I" + self.pkg.get_dir("includedir"))
            if self.pkg.name != os.path.basename(pkgincludedir):
                d = os.path.join(self.pkg.get_dir("includedir"),
                                 self.pkg.name)
                if os.path.isdir(d):
                    flags.insert(0, "-I" + d)

        elif flag == 'ldflags':
            flags.insert(0, "-L" + self.pkg.get_dir("libdir"))
            # Add library paths which may be indirectly required
            re_lib_dirs = self.pkg.config.get_compile_dependency_paths()
            for d in re_lib_dirs:
                flags.insert(0, "-L" + d)

        return flags

    #---------------------------------------------------------------------------

    def get_flags(self, flag):
        """
        Determine compilation flags for a given flag type.
        """

        cmd_line = self.get_pkg_path_flags(flag)

        # Build the command line, and split possible multiple arguments in lists.
        for lib in self.pkg.config.deplibs:
            if (self.pkg.config.libs[lib].have == "yes"
                and (not self.pkg.config.libs[lib].dynamic_load)):
                cmd_line += separate_args(self.pkg.config.libs[lib].flags[flag])

        # Specific handling of low-level libraries, which should come last,
        # such as -lm and -lpthread:
        # If -lm appears multiple times, only add it at the end of the
        # libraries, so that fast versions of the library may appear first

        if flag == 'libs':
            for lib in ['-lpthread', '-lm']:
                n = cmd_line.count(lib)
                if n > 0:
                    for i in range(n):
                        cmd_line.remove(lib)
                    cmd_line.append(lib)

        # Adapt for relocation on build stage on some systems

        self.flags_relocation(flag, cmd_line)

        return cmd_line

    #---------------------------------------------------------------------------

    def get_ar_lib_dir(self):
        """
        Determine directory containing library in archive mode.
        """

        return self.pkg.get_dir('libdir')

    #---------------------------------------------------------------------------

    def so_dirs_path(self, flags):
        """
        Assemble path for shared libraries in nonstandard directories.
        """
        retval = separate_args(self.pkg.config.rpath)
        i = len(retval) - 1
        if i < 0:
            return
        count = 0

        pkg_lib = os.path.join(self.pkg.get_dir('libdir'), self.pkg.name)
        if os.path.isdir(pkg_lib):
            retval[i] +=  ":" + pkg_lib
            count += 1

        if type(flags) == str:
            args = separate_args(flags)
        else:
            args = flags

        for arg in args:
            if (arg[0:2] == '-L' and arg[0:10] != '-L/usr/lib'
                and arg[0:6] != '-L/lib'):
                retval[i] += ":" + arg[2:]
                count += 1

        if count == 0:
            retval = ""

        return retval

    #---------------------------------------------------------------------------

    def obj_name(self, f):
        """
        Determine object file name.
        """

        o_name = os.path.basename(f)
        i = o_name.rfind(".")
        if i > -1:
            o_name = o_name[:i] + '.o'

        return o_name

    #-------------------------------------------------------------------------------

    def compile_src(self, src_list=None,
                    opt_cflags=None, opt_cxxflags=None, opt_fcflags=None,
                    keep_going=False,
                    stdout=sys.stdout, stderr=sys.stderr):
        """
        Compilation function.
        """
        retval = 0

        # Short names

        pkg = self.pkg

        # Find files to compile in source path

        c_files = fnmatch.filter(src_list, '*.c')
        h_files = fnmatch.filter(src_list, '*.h')
        cxx_files = fnmatch.filter(src_list, '*.cxx') + fnmatch.filter(src_list, '*.cpp')
        hxx_files = fnmatch.filter(src_list, '*.hxx') + fnmatch.filter(src_list, '*.hpp')
        f_files = fnmatch.filter(src_list, '*.[fF]90')
        o_files = fnmatch.filter(src_list, '*.o')

        # Determine additional include directories

        c_include_dirs = []
        for f in h_files:
            c_include_dirs.append(os.path.dirname(f))
        c_include_dirs = sorted(set(c_include_dirs))

        cxx_include_dirs = []
        for f in hxx_files:
            cxx_include_dirs.append(os.path.dirname(f))
        cxx_include_dirs = sorted(set(cxx_include_dirs))

        f_include_dirs = []
        for f in f_files:
            f_include_dirs.append(os.path.dirname(f))
        f_include_dirs = sorted(set(cxx_include_dirs))

        # Compile files

        for f in c_files:
            if (retval != 0 and not keep_going):
                break
            cmd = [self.get_compiler('cc')]
            if opt_cflags != None:
                cmd += separate_args(opt_cflags)
            for d in c_include_dirs:
                cmd += ["-I", d]
            cmd.append('-DHAVE_CONFIG_H')
            if os.path.basename(f) == 'cs_base.c':
                cmd += ['-DLOCALEDIR=\\"' + pkg.get_dir('localedir') + '\\"', \
                        '-DPKGDATADIR=\\"' + pkg.get_dir('pkgdatadir') + '\\"']
            cmd += self.get_flags('cppflags')
            cmd += separate_args(pkg.config.flags['cflags'])
            cmd += ["-c", f]
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1
            o_files.append(self.obj_name(f))

        for f in cxx_files:
            if (retval != 0 and not keep_going):
                break
            cmd = [self.get_compiler('cxx')]
            if opt_cxxflags != None:
                cmd += separate_args(opt_cxxflags)
            for d in cxx_include_dirs:
                cmd += ["-I", d]
            for d in c_include_dirs:
                cmd += ["-I", d]
            cmd.append('-DHAVE_CONFIG_H')
            cmd += self.get_flags('cppflags')
            cmd += separate_args(pkg.config.flags['cxxflags'])
            cmd += ["-c", f]
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1
            o_files.append(self.obj_name(f))

        for f in f_files:
            if (retval != 0 and not keep_going):
                break
            cmd = [self.get_compiler('fc')]
            f_base = os.path.basename(f)
            o_name = self.obj_name(f)
            if f_base == 'cs_user_boundary_conditions.f90':
                o_name = "cs_f_user_boundary_conditions.o"
                cmd += ["-o", o_name]
            if f_base == 'cs_user_parameters.f90':
                o_name = "cs_f_user_parameters.o"
                cmd += ["-o", o_name]
            if f_base == 'cs_user_extra_operations.f90':
                o_name = "cs_f_user_extra_operations.o"
                cmd += ["-o", o_name]
            if f_base == 'cs_user_initialization.f90':
                o_name = "cs_f_user_initialization.o"
                cmd += ["-o", o_name]
            if f_base == 'cs_user_physical_properties.f90':
                o_name = "cs_f_user_physical_properties.o"
                cmd += ["-o", o_name]
            if f_base == 'cs_user_porosity.f90':
                o_name = "cs_f_user_porosity.o"
                cmd += ["-o", o_name]
            if opt_fcflags != None:
                cmd += separate_args(opt_fcflags)
            for d in f_include_dirs:
                cmd += ["-I", d]
            if pkg.config.fcmodinclude != "-I":
                cmd += [pkg.config.fcmodinclude, srcdir]
            cmd += ["-I", pkg.get_dir('pkgincludedir')]
            if pkg.config.fcmodinclude != "-I":
                cmd += [pkg.config.fcmodinclude, pkg.get_dir('pkgincludedir')]
            cmd += separate_args(pkg.config.flags['fcflags'])
            cmd += ["-c", f]
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1
            o_files.append(o_name)

        return retval, o_files

    #---------------------------------------------------------------------------

    def link_obj(self, exec_name, obj_files=None, opt_libs=None,
                 stdout=sys.stdout, stderr=sys.stderr):
        """
        Link function.
        """
        retval = 0

        # Short names

        pkg = self.pkg

        # Directories

        call_dir = None
        temp_dir = None

        o_files = obj_files

        p_libs = self.get_flags('libs')

        # Special handling for some linkers (such as Mac OS X), for which
        # no multiple definitions are allowable in static mode;
        # in this case, extract archive, then overwrite with user files.

        if pkg.config.special_user_link == 'ar_x':

            call_dir = os.getcwd()
            temp_dir = tempfile.mkdtemp(suffix=".cs_link")
            os.chdir(temp_dir)

            lib0 = os.path.join(self.get_ar_lib_dir(),
                                'lib' + p_libs[0][2:] + '.a')
            p_libs = p_libs[1:]
            cmd = ['ar', 'x', lib0]
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1

            if obj_files:
                import shutil
                for f in obj_files:
                    if os.path.isabs(f):
                        f_src = f
                    else:
                        f_src = os.path.join(call_dir, f)
                    shutil.copy2(f_src, temp_dir)

            dir_files = os.listdir(os.getcwd())
            o_files = fnmatch.filter(dir_files, '*.o')

        # Prepare link command

        cmd = [self.get_compiler('ld')]
        cmd += ["-o", exec_name]
        if o_files:
            cmd += o_files

        # If present, address sanitizer needs to come first

        if '-lasan' in p_libs:
            p_libs.remove('-lasan')
            cmd += ['-lasan']

        if os.path.basename(exec_name) in self.pkg.config.exec_libs:
            cmd += [self.pkg.config.exec_libs[os.path.basename(exec_name)]]

        if opt_libs != None:
            if len(opt_libs) > 0:
                cmd += separate_args(opt_libs)
        cmd += self.get_flags('ldflags')
        cmd += p_libs
        if pkg.config.rpath != "":
            cmd += self.so_dirs_path(cmd)

        # Clean paths in link flags
        cmd_new = []
        for c in cmd:
            if c[:2] == '-L':
                c = '-L' + os.path.normpath(c[2:])
                if cmd_new.count(c) < 1:
                    cmd_new.append(c)
            else:
                cmd_new.append(c)
        cmd = cmd_new

        # Call linker

        if retval == 0:
            if run_command(cmd, pkg=pkg, echo=True,
                           stdout=stdout, stderr=stderr) != 0:
                retval = 1

        # Cleanup for special cases

        if temp_dir:
            if not os.path.isabs(exec_name):
                import shutil
                shutil.copy2(exec_name, os.path.join(call_dir, exec_name))
            for f in os.listdir(temp_dir):
                os.remove(os.path.join(temp_dir, f))
            os.chdir(call_dir)
            os.rmdir(temp_dir)

        return retval

    #---------------------------------------------------------------------------

    def compile_and_link(self, base_name, srcdir, destdir=None, src_list=None,
                         opt_cflags=None, opt_cxxflags=None, opt_fcflags=None,
                         opt_libs=None, force_link=False, keep_going=False,
                         stdout=sys.stdout, stderr=sys.stderr):
        """
        Compilation and link function.
        """
        retval = 0

        # Determine executable name

        exec_name = base_name
        if destdir != None:
            exec_name = os.path.join(destdir, exec_name)

        # Change to temporary directory

        call_dir = os.getcwd()
        temp_dir = tempfile.mkdtemp(suffix=".cs_compile")
        os.chdir(temp_dir)

        # Find files to compile in source path (special case
        # for user modules which must be compiled first)

        dir_files = os.listdir(srcdir)
        user_mod_name = 'cs_user_modules.f90'
        if user_mod_name in dir_files:
            dir_files.remove(user_mod_name)
            dir_files.insert(0, user_mod_name)

        src_list = []
        for f in dir_files:
            src_list.append(os.path.join(srcdir, f))

        retval, obj_list = self.compile_src(src_list,
                                            opt_cflags, opt_cxxflags, opt_fcflags,
                                            keep_going, stdout, stderr)

        if retval == 0 and (force_link or len(obj_list)) > 0:
            retval = self.link_obj(exec_name, obj_files=obj_list,
                                   opt_libs=opt_libs,
                                   stdout=stdout, stderr=stderr)

        # Cleanup

        for f in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, f))

        # Return to original directory

        os.chdir(call_dir)
        os.rmdir(temp_dir)

        return retval

#---------------------------------------------------------------------------

def compile_and_link(pkg, base_name, srcdir, destdir=None,
                     opt_cflags=None, opt_cxxflags=None, opt_fcflags=None,
                     opt_libs=None, force_link=False, keep_going=False,
                     stdout=sys.stdout, stderr=sys.stderr):
    """
    Compilation and link function.
    """
    c = cs_compile(pkg)

    retcode = c.compile_and_link(base_name,
                                 srcdir,
                                 destdir=destdir,
                                 opt_cflags=opt_cflags,
                                 opt_cxxflags=opt_cxxflags,
                                 opt_fcflags=opt_fcflags,
                                 opt_libs=opt_libs,
                                 force_link=force_link,
                                 keep_going=keep_going,
                                 stdout=stdout,
                                 stderr=stderr)

    return retcode

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    try:
        from code_saturne.cs_exec_environment import set_modules, source_rcfile
    except Exception:
        from cs_exec_environment import set_modules, source_rcfile

    test_mode, force_link, keep_going, src_dir, dest_dir, \
        version, cflags, cxxflags, fcflags, libs = process_cmd_line(argv, pkg)

    if (version):
        pkg = pkg.get_alternate_version(version)

    set_modules(pkg)    # Set environment modules if present
    source_rcfile(pkg)  # Source rcfile if defined

    if test_mode == True:
        dest_dir = None

    retcode = compile_and_link(pkg,
                               pkg.solver,
                               src_dir,
                               destdir=dest_dir,
                               opt_cflags=cflags,
                               opt_cxxflags=cxxflags,
                               opt_fcflags=fcflags,
                               opt_libs=libs,
                               force_link=force_link,
                               keep_going=keep_going)

    sys.exit(retcode)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Retrieve package information (name, version, installation dirs, ...)

    from code_saturne.cs_package import package
    pkg = package()

    # Create an instance of the main script
    main(sys.argv[1:], pkg)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
