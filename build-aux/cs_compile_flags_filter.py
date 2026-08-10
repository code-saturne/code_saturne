#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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
import subprocess, fnmatch, platform

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
# Filter functions
#-------------------------------------------------------------------------------

def filter_ginkgo(args):
    """
    Filter compile flags for compilation including Ginkgo library,
    which makes heavy use of shadowed variables.
    """

    for i, a in enumerate(args):
        if a == "-Werror=shadow":
            args[i] = "-Wshadow"

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    flags = sys.argv[1:]
    if len(sys.argv) > 1:
        filter_name = sys.argv[1]
        flags = sys.argv[2:]
        if filter_name == "ginkgo":
            filter_ginkgo(flags)

    sep = " "
    sys.stdout.write(sep.join(flags))

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
