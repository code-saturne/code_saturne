<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
-->

\page cs_ug_udf User-defined functions

[TOC]

User-defined functions {#cs_ug_udf_intro}
======================

In addition to the settings defined with a `setup.xml` file and usually
defined using the GUI, both basic and more advanced settings and models
can be defined using *user-defined functions*, called at specific points
in the code's execution.

As mentioned in the [case structure description](@ref sec_prg_run_udf),
and [step by step preparation](@ref sec_prg_stepbystepcalculation) documentation,
all C, C++, and Fortran files located directly in a case's `SRC`
subdirectory will be automatically compiled and linked with the solver during
a run's [initialization](@ref sec_prg_exec_stages) stage.

This section provides basic information regarding the use of such functions.

User-defined functions usually allow for more complete or detailed
settings than those defined through the GUI, and can be combined
with them. All settings defined through the GUI (except for the activation
of non-default physical or turbulence models) can be overridden
in the matching user-defined functions, so the recommended
practice is to define all that is possible using the GUI, and
user-defined functions for the rest.

It is possible, though not recommended, to define a computation solely
through user-defined functions, without using an XML setup or needing
the GUI. This is not recommended for general use, as consistency checking
and upgrading to a later version require much more effort, but _might_
be useful for application specific or embedded models.

Function and file names
-----------------------

In most cases, a given user function reference may be found in the C or
Fortran file of the same name.

When both C and Fortran versions of a given function of function series exist,
the C version may be prefixed by `cs`, while the fortran version may be
prefixed by `cs_f`. The Fortran version is usually called first, so the
C function has the "last word". For example, `cs_user_extra_operations`
may be found in the `cs_user_extra_operations.c` file, while
`cs_f_user_extra_operations` is in `cs_user_extra_operations.f90`.

Callback functions  {#cs_ug_udf_intro_callback}
------------------

For advanced settings, many settings involve passing pointers to
[callback](https://en.wikipedia.org/wiki/Callback_%28computer_programming%29)
functions which will be be called at specific call points.

Callback functions need to match a specific profile (i.e. type
and number of arguments) based on their role, but can otherwise
be named in any manner, so choosing a self_explanatory name
is recommended.

In most of this user documentation both *callback function*
and *function pointer* terms may be used interchangeably.
To be precise, callbacks are assigned using function pointers.

Examples
--------

Many illustrations are provided in the
[examples](@ref cs_user_examples) section of the documentation.
The matching files are also found under the
`src/cs_user_examples` directory, in the source tree, or
under `${install_prefix}/code_saturne/user_sources/EXAMPLES` in
most installations).

When multiple examples are provided, example file names are
defined by appending the name of the matching example
to the matching reference file name. For example, general
examples for user source terms are provided in
`cs_user_source_terms-base.c`, while examples specific to advanced
turbulence model modification are found in
`cs_user_source_terms-turbulence.c`

The user is encouraged to check which examples are available, and to study
those that are relevant to a given setup.

Comparing template and example files with a graphical file comparison tool
should help the user highlight the matching sections from the examples,
so it is recommended as good practice for those not already very familiar
with those user functions. Code blocks from examples can of
course be cherry-picked and adapted to the user's needs. indiscriminately
copying complete files to use only a fraction of the examples therein
is considered bad practice, as it leads to confusing code, and is a form of
[cargo-cult programming](https://en.wikipedia.org/wiki/Cargo_cult_programming).

Base user-defined functions {#cs_ug_udf_common}
===========================

The following sections lists some of the most frequently
used user-defined functions, and their associated roles.
It is not by any means exhaustive (for a full list, inspect
all files in the source tree's `src/user` directory, or
`${install_prefix}/code_saturne/user_sources/REFERENCE` in
most installations).

Functions called during the computation setup
---------------------------------------------

The following functions are called before the resolution stage,
and before even the mesh is read. As a consequence, access to mesh
elements and counts, and any field or other array value defined
on the mesh, is not available in these functions
(though some settings may involve setting pointers to
[callback](@ref cs_ug_udf_intro_callback)
functions which will be able to access those elements.

### For all physical models

- \ref cs_user_model (in \ref cs_user_parameters.c)

  Allows defining user scalars (species), variances, or activating
  a specific physical model (by setting \ref cs_physical_model_type_t
  entry values in the \ref cs_glob_physical_model_flag array).

  It is called before all other physical or numerical
  oriented user functions (only system settings and mesh definitions are called
  earlier), so that the the variable and property fields implied by those
  models are instanciated and available in the following user function call points.

  The equivalent Fortran function is named \ref usppmo.

- \ref cs_user_zones

  Allows defining which zones (based on mesh groups or geometric criteria)
  will be used for the computation. This allows advanced definitions,
  such as time-evolving zones based un user callback functions, in
  addition to the basic definitions provided by the GUI..

  It is called before all physical or numerical oriented user functions.

- \ref cs_user_parameters (in \ref cs_user_parameters.c)

  Allows defining most general settings, such as reference physical properties
  model and numerical settings  for main variable and property fields, etc.
  By default, settings should go here.

  It is called before all other physical or numerical
  oriented user functions (only system settings and mesh definitions are called
  earlier), so that the the variable and property fields implied by those
  models are instanciated and available in the following user function call points.

  The equivalent Fortran function is named \ref usppmo.

- \ref cs_user_postprocess_writers, \ref cs_user_postprocess_meshes,
  and \ref cs_user_postprocess_probes (in \ref cs_user_postprocess.c)

  May be used to define or modify postprocessing extracts using the
  supported output formats, using the
  [mesh and writer](@ref cs_ug_postprocess_intro) concepts.

- \ref cs_user_boundary_conditions_setup (in \ref cs_user_boundary_conditions.c)

  May be used to define advanced boundary conditions using zone-based
  definitions.

- \ref cs_user_finalize_setup (in \ref cs_user_parameters.c)

  May be used for additionl definitions, or as an alternative or extension
  to \ref cs_user_boundary_conditions_setup.

### For specific models

- \ref cs_user_combustion (in \ref cs_user_parameters.f90)

  Allows specifying calculation options specific to the selected
  combustion model (*gas*, *pulverized coal*, or *heavy fuel*).

- \ref cs_user_atmospheric_model.f90

  Contains several user subroutines used to define atmospheric
  model settings such as ground properties and 1-d atmospheric profiles.

- \ref cs_user_lagr_model

   Allows defining physical, numerical and post-processing options
   for the Lagrangian model (dispersed phase).

Functions called before time stepping
-------------------------------------

- \ref cs_user_mesh_modify

  Allows performing various mesh modifications during the preprocessing stage.

- \ref cs_user_initialization

  Allows setting the initial values of variables and properties.
  In case of computation restart, the values from the restart files are loaded before
  this function is called, so can be either further modified or left alone.

- \ref cs_user_extra_operations_initialize (in \ref cs_user_extra_operations.c)

  Called just after \ref cs_user_initialization (with only a few updates
  for some specific models in between), this function is used for more
  general operation than simply initializing variables.

  Since it is placed in the same file as \ref cs_user_extra_operations, it may
  be used to initialize user variables or structures local to that file
  before time stepping.

Functions called during time stepping
-------------------------------------

### For all physical models

- \ref cs_user_physical_properties

  Allows defining variable physical property (such as fluid density, viscosity ...)
  values. It is called at each time step to allow for updating the relevant fields.

- \ref cs_user_boundary_conditions

  Allows defining complex boundary conditions.

- \ref cs_user_source_terms

  Allows defining complex source terms.

- \ref cs_user_postprocess_values (in \ref cs_user_postprocess.c)

  May be used to output locally-computed volume or surface values,
  such as formulas involving fields, or for fine grained association
  of given fields with different writers.

- \ref cs_user_extra_operations

  Allows defining variable physical property (such as fluid density, viscosity ...)
  values. It is called at each time step to allow for updating the relevant fields.

### For specific models

- \ref cs_user_lagr_boundary_conditions and \ref cs_user_lagr_volume_conditions

  Allow defining and modifying boundary conditions and volume injections for
  the Lagrangian particles (dispersed phase).

Functions called during after time stepping
-------------------------------------------

- \ref cs_user_extra_operations_finalize (in \ref cs_user_extra_operations.c)

  Called just after the time stepping/resolution stage, this function
  allows handling operations required only at the end of the computation (such
  as some specific post-processing extracts), and possibly cleaning up and freeing
  structures used in the main \ref cs_user_extra_operations function.

User-defined Fortran modules
----------------------------

When compiling user sources in a case's `SRC` directory, the order
of compilation is not based on any dependency check. This is not
an issue for additional user C or C++ code, but can be an issue
for Fortran code with user-defined modules.

If a file named `cs_user_modules.f90` is present, it
will be compiled before any other Fortran file. So if needed,
user-defined modules should be defined in that file, to ensure they
are available in other user subroutines.

Main variables and structures
=============================

The main [variables and structures reference](@ref cs_var_dico)
provides descriptions and recommendations for variables
needed in user-defined functions.
