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

import os

#===============================================================================
# Local functions
#===============================================================================

#===============================================================================
# Defining parameters for a calculation domain
#===============================================================================

def domain_prepare_data_add(domain):
    """
    Additional steps to prepare data
    (called in data preparation stage, between copy of files
    in DATA and copy of link of restart files as defined by domain).
    """

    return

#-------------------------------------------------------------------------------

def domain_copy_results_add(domain):
    """
    Additional steps to copy results or cleanup execution directory
    (called at beginning of results copy stage).
    """

    return

#-------------------------------------------------------------------------------

def define_domain_parameters(domain):
    """
    Define domain execution parameters.
    (called in data preparation stage, just after copy of base setup files
    in execution directory, and before compiling user-defined functions).
    """

    # Note: when data is already staged
    #----------------------------------

    # When a case is prepared (staged) then run or submitted in a separate
    # step (such as when using the `code_saturne submit` command
    # or running under the Study Manager), this function is called again,
    # As most settings are done in the first pass, the `domain.data_is_staged`
    # attribute can be checked to avoid duplicating some operations.

    # Reusing output from previous runs
    #----------------------------------

    # By default, upon restart, if the previous mesh is available in the
    # checkpoint, preprocessing is skipped. To allow additional preprocessing,
    # set:
    #   domain.preprocess_on_restart = True

    # To use the output of a previous Preprocessor (mesh import) run, set:
    #   domain.mesh_input = 'RESU/<run_id>/mesh_input'

    # To use the mesh output of a previous solver run, set:
    #   domain.mesh_input = 'RESU/<run_id>/mesh_output'

    # To reuse a previous partitioning, set:
    #   domain.partition_input = 'RESU/<run_id>/partition'

    # To continue a calculation from a previous checkpoint, set:
    #  domain.restart_input = 'RESU/<run_id>/checkpoint'

    # To continue a calculation from a previous checkpoint with
    # a different mesh, it may be necessary to provide a mesh_input
    # matching that checkoint if is not contained, so set:
    #  domain.restart_mesh_input = 'RESU/<run_id>/restart_mesh_input'

    domain.preprocess_on_restart = False
    domain.mesh_input = None
    domain.partition_input = None
    domain.restart_input = None
    domain.restart_mesh_input = None

    # Defining meshes to import (only if domain.mesh_input = None)
    #-------------------------------------------------------------

    # A case-specific mesh directory can be defined through domain.mesh_dir.

    # domain.mesh_dir = None

    # Mesh names can be given as absolute path names; otherwise, file are
    # searched in the following directories, in the order:
    #   - the   case mesh database directory, domain.mesh_dir (if defined)
    #   - the  study mesh database directory, ../MESH directory (if present)
    #   - the   user mesh database directory, see ~/.code_saturne.cfg
    #   - the global mesh database directory, see $prefix/etc/code_saturne.cfg

    # Meshes should be defined as a list of strings defining mesh
    # file names. If preprocessor options must be used for some meshes,
    # a tuple containing a file name followed by additional options
    # may be used instead of a simple file name, for example:
    #   domain.meshes = ['part1.des',
    #                    ('part2_med', '--format', 'med',
    #                     '--num', '2', '--reorient'),
    #                    '~/meshdatabase/part3.unv']

    # domain.meshes = None

    # Logging arguments
    #------------------

    # Command-line arguments useful for logging, or determining the calculation
    # type may be defined here, for example:
    #   domain.logging_args = '--logp'

    # domain.logging_args = None

    # Solver options
    #---------------

    # Running the solver may be deactivated by setting
    #   domain.exec_solver = False

    # The type of run may be determined using domain.solver_args.
    # For example:
    #   domain.solver_args = '--quality'
    # allows activation of elementary mesh quality criteria output, while:
    #   domain.solver_args = '--benchmark'
    # allows running  basic linear algebra operation benchmarks.
    # To run the solver's preprocessing stage only (mesh joining, smoothing,
    # and other modifications), use:
    #   domain.solver_args = '--preprocess'

    # domain.exec_solver = True
    # domain.solver_args = None

    # Compile and build options
    #--------------------------

    # Additionnal compiler flags may be passed to the C, C++, or Fortran
    # compilers, and libraries may be added, in case linking of user
    # functions against external libraries is needed.

    # Note that compiler flags will be added before the default flags;
    # this helps ensure added search paths have priority, but also implies
    # that user optimization options may be superceded by the default ones.

    domain.compile_cflags = None
    domain.compile_cxxflags = None
    domain.compile_fcflags = None
    domain.compile_libs = None

    # import pprint
    # pprint.pprint(domain.__dict__)

    return

#-------------------------------------------------------------------------------
