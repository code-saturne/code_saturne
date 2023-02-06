#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

import sys
import os.path

#===============================================================================
# Utility functions
#===============================================================================

#-------------------------------------------------------------------------------

def pythonpath_filter(l):
    """
    Build a minimal string for the PYTHONPATH shell variable containing the
    files passed as arguments.
    """

    p = []
    for d in sys.path[1:]:
        required = False
        for f in l:
            if os.path.isfile(os.path.join(d, f)):
                required = True
            if required:
                p.append(d)

    s = ''
    for d in p:
        if d.find(' ') > -1:
            s += '\"' + d + '\":'
        else:
            s += d + ':'
    if len(s) > 0:
        s = s[:-1]
    print(s)

#-------------------------------------------------------------------------------

def ld_library_path_filter(l):
    """
    Build a minimal string for the LD_LIBRARY_PATH shell variable containing
    the file patterns passed as arguments.
    """

    e = []
    k = 'LD_LIBRARY_PATH'
    if not k in os.environ:
        return

    import re

    patterns = []
    for f in l:
        patterns.append(f)

    e = os.environ[k].split(':')

    p = []
    for d in e:
        if not os.path.isdir(d):
            continue
        required = False
        for f in os.listdir(d):
            found = []
            for pattern in patterns:
                if re.search(pattern, f):
                    found.append(pattern)
            if found:
                required = True
                for pattern in found:
                    patterns.remove(pattern)
        if required:
            p.append(d)

    s = ''
    for d in p:
        if d.find(' ') > -1:
            s += '\"' + d + '\":'
        else:
            s += d + ':'
    if len(s) > 0:
        s = s[:-1]
    print(s)

#-------------------------------------------------------------------------------

def salome_paraview_ld_add_path_filter(l):
    """
    Build a minimal string for the LD_LIBRARY_PATH shell variable additions
    for ParaView Catalyst from SALOME distribution.
    """

    e = []
    k = 'LD_LIBRARY_PATH'
    if not k in os.environ:
        return

    # 'cgns' and 'TTK' do not seem to be needed at least in a simple case,
    # but might be with more complex scripts.
    keep_list = ('catalyst', 'cgns', 'embree', 'hdf5', 'openturns',
                 'openVKL', 'ospray', 'rkCommon', 'TTK')

    e = os.environ[k].split(':')

    p = []
    for d in e:
        if not os.path.isdir(d):
            continue
        for c in keep_list:
            if d.find(c) > -1:
                p.append(d)
                break

    s = ''
    for d in p:
        if d.find(' ') > -1:
            s += '\"' + d + '\":'
        else:
            s += d + ':'
    if len(s) > 0:
        s = s[:-1]
    print(s)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    # Simple parse for command line

    cmd = ''
    if len(sys.argv) == 1:
        cmd = '--help'
    elif sys.argv[1] in ['-h', '--help']:
        cmd = '--help'
    else:
        cmd = sys.argv[1]

    # Call tools

    if cmd == 'pythonpath_filter':
        pythonpath_filter(sys.argv[2:])

    elif cmd == 'ld_library_path_filter':
        ld_library_path_filter(sys.argv[2:])

    elif cmd == 'salome_paraview_ld_add_path_filter':
        salome_paraview_ld_add_path_filter(sys.argv[2:])

    else:
        usage = \
            """Usage: %(prog)s <tool> [arguments]

Tools:
  pythonpath_filter
  ld_library_path_filter
  salome_paraview_ld_add_path_filter

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})
        sys.exit(1)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
