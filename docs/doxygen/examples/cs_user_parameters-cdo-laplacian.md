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

\page cs_cdo_laplacian Minimalist example to solve a Laplacian using CDO schemes

[TOC]

<!--
    References used in this page
-->

[Bonel14]: https://hal.archives-ouvertes.fr/tel-01116527

Introduction {#sec_parameters_cdo_laplacian_intro}
============

Here is a minimalist example to solve a scalar-valued Laplacian
problem (an isotropic diffusion equation). Since one considers the CDO
framework to solve this problem, this enables a quick start with this
framework.  This simple example is performed in three steps:
* Define a mesh and its boundary zones (done in the GUI)
* Activate the CDO and create the equation to solve (done in the
  function \ref cs_user_model of the user source file \ref
  cs_user_parameters.cpp)
* Define the equation to solve and set the boundary conditions (done
  in the function \ref cs_user_finalize_setup of the user source file \ref
  cs_user_parameters.cpp)

For a very beginner, one strongly advises to read the page [Case directory structure](@ref cs_ug_case_structure).

In a terminal, write the following syntax for the creation of a new
study with a new case:

```
code_saturne create --study STUDY CASE_NAME1
```

which creates a study directory `STUDY` with case sub-directory
`CASE_NAME1`. If no case name is given, a default case directory
called `CASE1` is created.


First step: Preprocessing {#sec_cdo_laplacian_prepro}
=======================

One assumes that the mesh is a Cartesian mesh generated using the GUI
but other ways are also possible. A more detail description of the
different preprocessing possibilities is available [here](@ref cs_ug_mesh_prepare)

\anchor gui_cdo_laplacian_cartesian_mesh
\image html gui_cartesian_mesh.png "GUI: Cartesian mesh definition"

The generated Cartesian is built with 6 mesh face groups collecting
boundary faces and named `X0`, `X1`, `Y0`, `Y1`, `Z0` and `Z1`. These
mesh face groups are used to define the two boundary zones needed for
this example.  One proceeds as detailed in the GUI screenshot.

![GUI: Definition of boundary zones](gui_boundary_zone_cartesian.png)


Second step: Activate CDO and add a user-defined equation {#sec_cdo_laplacian_init_model}
=======================

The second step corresponds to the edition of the user source file
named **cs_user_parameters.c** and especially the function \ref
cs_user_model

Click on the icon displayed just below in the GUI toolbar which
corresponds to "Open the source file editor" <img src="src_editor-icon.png" width="60px">

Then right-click on the file named `REFERENCE/cs_user_parameters.c`
and select **Copy to SRC**. Then, you can choose your favorite file
editor to add in the function \ref cs_user_model the following lines.

\snippet cs_user_parameters-cdo-laplacian.c param_cdo_laplacian_init

This first activates the CDO module (please refer to \ref
cs_user_parameters_h_cdo_activation for more details) and then add a
new scalar equation called `Laplacian` with an unknown named
`potential`. This will be the name of the associated variable field. A
default boundary condition is also defined and it corresponds to an
homogeneous Neumann boundary condition. If no other boundary condition
is defined, then all the boundary faces will be associated to this
default definition.


Third step: Define the equation to solve {#sec_cdo_laplacian_finalize}
=======================

The last step corresponds to the modification of the function \ref
cs_user_finalize_setup in the file cs_user_parameters.c

\snippet cs_user_parameters-cdo-laplacian.c param_cdo_laplacian_finalize

After having retrieved the structure \ref cs_equation_param_t
associated to the equation to solve, one first add a diffusion term
which is associated to the default property named `unity`. Then, one
defines two Dirichlet boundary conditions: a first one on the left
side (`X0`) with a value equal to `0` and a second one on the right
side (`X1`) with a value equal to `1`.

Last step: Run the computation
=======================

For a very beginner, one strongly advises to read the page [Running a calculation](@ref cs_ug_run_computation).

In a terminal, write the following syntax when your current directory
corresponds to one of the sub-directories of a case.

```
code_saturne run
```


To go beyond
=======================

In order to change the numerical settings related to an equation call
the function \ref cs_equation_param_set inside the user function named
\ref cs_user_parameters in the file \ref cs_user_parameters.c

Here are some examples of numerical settings:

\snippet cs_user_parameters-cdo-condif.c param_cdo_numerics

This relies on  a `key` `value` principle. The available keys are listed [here](\ref cs_equation_key_t)



For the reader willing to get a better understanding of the
mathematical concepts underpinning the CDO schemes, one refers to the
[PhD thesis entitled *Compatible Discrete Operator schemes on
polyhedral meshes for elliptic and Stokes equations*][Bonel14] \cite Bonel14
