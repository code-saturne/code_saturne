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

import sys
import os.path

#===============================================================================
# Install filter for user-defined functions.
#===============================================================================

# Reads from standard input, output to standard output.
# Takes version as argument.

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    # Simple parse for command line

    cmd = ''
    version = 'VERS'

    if len(sys.argv) > 1:
        version = sys.argv[1]

    lines = sys.stdin.readlines()

    p_skip = False
    p_space = False

    for l in lines:

        l = l.rstrip()

        if l == '':
            if not (p_skip and p_space):
                print()
            p_skip = False
            p_space = True
            continue

        l_s = l.lstrip()
        if l_s[:12] == "#pragma weak" or l_s[:5] == "/*! [":
            # Remove "pragma weak" attributes for functions
            p_skip = True
            continue

        elif l_s[:5] == "/*! [":
            # remove snippet
            p_skip = True
            continue

        elif l_s[:38] == "! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE":
            # Fortran only
            p_skip = True
            continue

        else:
            # Version header
            if l_s == "/* VERS */":
                l = "/* code_saturne version " + version + " */"

            # All lines
            p_skip = False
            p_space = False

        print(l)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
