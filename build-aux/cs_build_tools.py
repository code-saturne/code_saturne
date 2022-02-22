#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2021 EDF S.A.
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

def clean_lib_search_flags(flags):

    n_flags = []
    for f in flags:
        prf = ""
        if f[0:2] == '-L':
            prf = f[0:2]
            p = f[2:]
        else:
            p = f
        if os.path.isfile(p) or os.path.isdir(p):
            p = os.path.normpath(p)
        n_flags.append(prf + p)

    c_flags = []
    while n_flags:
        f = n_flags.pop(0)
        if f not in ('-L/lib', '-L/usr/lib') and n_flags.count(f) < 1:
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

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
