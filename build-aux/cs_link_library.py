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

import sys, os.path
import argparse
import subprocess, fnmatch

#-------------------------------------------------------------------------------

# Utility functions for build system to assemble shared or archive libraries.
# Note that some settings/options adapted to various systems are also
# defined in `config/cs_auto_flags.sh`.

# These functions have beed tested on Linux, and may need to be adapted to
# various systems using other linker configurations and options.
# In this case, adding functions dedicated to various cases is recommended.

# These functions play a role similar to the Libtool scripts previously built.
# All the built-in knowledge of Libtool is dropped here, and may need to be
# re-learned here, but most of that knowledge regarded obsolete systems (or
# systems not encountered in a computational environment), and some aspects
# of Libtool's automation (especially the fact that incorrect paths in .la
# files could not be ignored) caused too frequent issues.

#===============================================================================
# Utility functions
#===============================================================================

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def parse_cmd_line(argv):
    """
    Process the passed command line arguments.
    """

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    description="Link object or library files into a shared library."
    parser = argparse.ArgumentParser(description=description)

    ar_group = parser.add_mutually_exclusive_group()

    parser.add_argument("--linker", dest="linker", type=str, default='',
                        metavar="<linker>",
                        help="define linker")

    parser.add_argument("--version", dest="version", type=str, default='',
                        metavar="<version>",
                        help="define library version number")

    ar_group.add_argument("--archive", dest="archive", default=False,
                          action="store_true",
                          help="build archive library")

    h = "include every object in subsequent archives (for shared library build)"

    ar_group.add_argument("--whole-archive-start", dest="whole_archives",
                          nargs="*", help=h)

    parser.add_argument("--whole-archive-end", default=False,
                        action="store_true",
                        help="ends list beginning with --whole-archive-start")

    h =   "if 'yes', use standard library search paths; " \
        + "if 'no', only used specified paths, " \
        + "if 'try', check without paths first."

    parser.add_argument("--std-search-paths", dest="stdlib", type=str,
                        default='yes',
                        metavar="<yes/no/try>",
                        help=h)

    parser.add_argument("--echo", dest="echo", default=False,
                        action="store_true",
                        help="echo commands run")

    parser.add_argument("-o", dest="output_name", type=str,
                        metavar="<output_name>",
                        help="output library name")

    options, ext = parser.parse_known_args(argv)

    return options, ext

#-------------------------------------------------------------------------------
# Call external command, with checking and simple logging for errors
#-------------------------------------------------------------------------------

def run_command(cmd, echo=False):
    """
    Call external command, with checking and simple logging for errors
    """

    if echo:
        print(' '.join(cmd))
    p = subprocess.Popen(cmd)
    p.communicate()
    if p.returncode != 0 and not echo:
        print(sys.argv[0], ":", cmd, "returned", p.returncode)

    return p.returncode

#-------------------------------------------------------------------------------
# Build (usually static) library archive
#-------------------------------------------------------------------------------

def build_archive(output,
                  archives,
                  objects,
                  ar_options=[], ranlib_options=[],
                  echo=False):
    """
    Build an archive given the archives and objects list.
    """

    # Add objects from libraries

    if archives:
        objdir = ('.tmp_objs')
        if not os.path.isdir(objdir):
            os.mkdir(objdir)
        wd = os.getcwd()
        os.chdir(objdir)
        if os.path.isabs(output):
            p_output = output
        else:
            p_output = os.path.join('..', output)
        for a in archives:
            if os.path.isabs(a):
                f = a
            else:
                f = os.path.join('..', a)
            cmd = ['ar', 'x', f]

            retcode = run_command(cmd, echo)

            dir_files = os.listdir('.')
            o_files = fnmatch.filter(dir_files, '*.o')

            cmd = ['ar', 'cr'] + ar_options + [p_output] + o_files

            retcode = run_command(cmd, echo)

            for f in dir_files:
                os.remove(f)

        os.chdir(wd)
        os.rmdir(objdir)

    # Add objects provided directly

    if objects:
        cmd = ['ar', 'cr'] + ar_options + [output] + objects
        retcode = run_command(cmd, echo)

    # Now generate index

    cmd = ['ranlib'] + ranlib_options + [output]
    retcode = run_command(cmd, echo)

    return retcode

#-------------------------------------------------------------------------------
# Build shared library archive.
#-------------------------------------------------------------------------------

def build_shared_library(linker,
                         output,
                         version,
                         archives,
                         objects,
                         other=[],
                         echo=False,
                         stdlib='yes'):
    """
    Build an archive given the archives and objects list.
    """

    if not linker:
        print("No linker specified for output library ", output)
        return 1

    o_name = os.path.basename(output)
    if version:
        f, e = os.path.splitext(o_name)
        o_name_v = f + '-' + version + e
        f, e = os.path.splitext(output)
        output_v = f + '-' + version + e
    else:
        o_name_v = o_name
        output_v = output

    # TODO:
    # Some flags should already be provided by the caller through
    # the command-line arguments. More things could be moved to
    # this file (or the Python scripts in general), and later
    # removed from `config/cs_auto_flags.sh`.

    cmd = [linker, "-o", output_v]

    if stdlib != 'yes':
        cmd.append("-nostdlib")

    cmd += ["-Wl,-soname", "-Wl," + o_name_v]

    # Add objects from libraries

    if archives:
        cmd.append("-Wl,-whole-archive")
        for a in archives:
            cmd.append(a)
        cmd.append("-Wl,-no-whole-archive")

    # Add external objects and archives provided directly

    if objects:
        for o in objects:
            cmd.append(o)

    # Pass all other options to linker;
    # Convert libtool-like "-R' syntax to rpath.

    for o in other:
        if o in ("-lptscotch", "-lscotch"):
            continue
        if o[:2] == '-R':
            cmd += ["-Wl,-rpath", "-Wl,"+o[2:]]
        else:
            cmd.append(o)

    # Now call linker

    retcode = run_command(cmd, echo)

    # Try again allowing standard libraries if failed.
    # This occurs for example for the ParaView Catalyst adaptor plugin
    # using GCC, where adding GCC's own `crti.o` and `crtbeginS.o`
    # files is an alternative solution (used by Libtool).

    if retcode != 0 and stdlib == 'try':
        cmd.remove("-nostdlib")
        print()
        print("Retry, allowing standard librairies search path for link:")
        retcode = run_command(cmd, echo)

    # Add symbolic link in case of version number.

    if output_v != output and retcode == 0:
        try:
            if os.path.islink(output):
                os.remove(output)
            os.symlink(o_name_v, output)
        except Exception:
            fmt = "Failed to create sybolic link from {0} to {1}."
            print(fmt.format(output, o_name_v))
            return 1

    return retcode

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys, os.path

    retcode = 0

    options, ext = parse_cmd_line(sys.argv[1:])

    objects = []
    archives = []
    other = []

    if options.echo:  # Add blank line for readability in echo mode.
        print()

    for s in ext:
        f, e = os.path.splitext(s)
        if e in ('.o', '.obj'):
            objects.append(s)
        elif e == '.a':
            archives.append(s)
        else:
            if s == '-Wl,--end-group':  # Eliminate empty start/end group sequences
                p = other.pop()
                if p != '-Wl,--start-group':
                    other.apend(p)
                    other.apend(s)
            else:
                other.append(s)

    # Clean duplicates:
    # - keep last instance for objects and libraries
    # - keep first instance for paths

    objects.reverse()
    l_new = []
    for o in objects:
        if l_new.count(o) > 0:
            continue
        l_new.append(o)
    objects = l_new
    objects.reverse()

    archives.reverse()
    l_new = []
    for o in archives:
        if l_new.count(o) > 0:
            continue
        l_new.append(o)
    archives = l_new
    archives.reverse()

    other.reverse()
    l_new = []
    for o in other:
        l_new.append(o)
    other = l_new
    other.reverse()

    l_new = []
    for o in other:
        if o[:2] == '-L':
            o = '-L' + os.path.normpath(o[2:])
            if l_new.count(o) > 0:
                continue
        if o[:10] == '-Wl,-rpath':
            if l_new.count(o[10:]) > 0:
                continue
        l_new.append(o)
    other = l_new

    # Build archive or dynamic library
    # Remark: For dynamic libraries, archives not in whole_archives are
    #         assumed to be external dependencies, so grouped with "other".
    #         For archive libraries, dependencies should not be provided,
    #         so specifying whole archives is not necessary, and archives

    if options.archive:
        retcode = build_archive(options.output_name,
                                archives, objects,
                                ar_options=[], ranlib_options=[],
                                echo=options.echo)

    else:
        retcode = build_shared_library(options.linker,
                                       options.output_name,
                                       options.version,
                                       options.whole_archives,
                                       objects,
                                       archives + other,
                                       echo=options.echo,
                                       stdlib=options.stdlib)

    if options.echo:  # Add blank line for readability in echo mode.
        print()

    sys.exit(retcode)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
