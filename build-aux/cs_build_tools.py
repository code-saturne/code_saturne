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

import fnmatch
import os
import sys
import tempfile

#-------------------------------------------------------------------------------

#===============================================================================
# Functions
#===============================================================================

def get_ld_default_search_path():
    """
    On Linux systems, try to determine the search path used by
    the dynamic linker (so as to allow filtering commands to avoid
    repeating paths already in that path).
    """

    import subprocess

    path = ""
    try:
        cmd = ['/sbin/ldconfig', '-v']
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        output = p.communicate()

        for l in output[0].splitlines():
            if l[:1] == '/':
                idx = l.find(':')
                if idx > -1:
                    path += l[:idx+1]
                else:
                    path += l

    except Exception:
        pass

    return path[:-1]

def clean_lib_search_flags(flags):
    """
    Clean library search flags, to normalize paths and remove
    options which are not search flags.
    """

    n_flags = []
    for f in flags:
        prf = ""
        if f[0:2] == '-L':
            prf = f[0:2]
            p = f[2:]
        else:
            p = f
            if f[0:2] == '-l':
                # Also remove names which do not seem to match library names
                # (should not be needed but helps remove options such as
                # "-loopopt=0" which appears in the autoconf'ed FCLIBS using
                # the Intel ifx Fortran.
                if f.find("=") > 0:   # Non-library option
                    continue
        if os.path.isfile(p) or os.path.isdir(p):
            p = os.path.normpath(p)
        n_flags.append(prf + p)

    c_flags = []
    while n_flags:
        f = n_flags.pop(0)
        if n_flags.count(f) < 1:
            if f[0:2] == '-L':
                if f[2:] in ['/lib', '/usr/lib']:
                    continue
            c_flags.append(f)

    return c_flags

    sys.exit(0)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    tool = sys.argv[1]

    if tool == 'clean-lib-flags':
        c_flags = clean_lib_search_flags(sys.argv[2:])
        m = ''
        for f in c_flags:
            m += ' ' + f
        print(m[1:])

    elif tool == 'get-ld-default-search-path':
        print(get_ld_default_search_path())

    else:
        print(sys.argv[0], ": unknown command: ", tool)
        sys.exit(1)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
