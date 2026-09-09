<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

\page cs_ug_time_settings Time settings

Time Settings
=============

The **Time Settings** page allows users to define the time step strategy and the
overall physical duration of the computation.

\anchor gui_time_step
\image html gui_time_step.png "Time step settings"

## Time step options

- **Constant**: uses a reference time step that remains constant throughout the
  computation and is identical for all cells.
- **Time varying** (adaptive): the time step is initialized with the reference
  value and automatically evolves during the computation while satisfying the
  specified maximum CFL and Fourier numbers. Additional parameters allow users
  to control the maximum time step variation and scaling factors.
- **Steady** (local time step): intended **only for steady-state computations**.
  In this mode, the time step may vary both in time and from one cell to
  another. It is initialized with the reference value and then evolves while
  satisfying the specified maximum CFL and Fourier numbers.

## Velocity-pressure algorithm

- **SIMPLEC** (default)
- **Inner iterations**

## Stopping criterion

- **Number of time steps**
- **Physical time (s)**
- **Additional time steps**
- **Additional physical time (s)**

The distinction between **Number of time steps** and **Additional time steps**
(and similarly between **Physical time** and **Additional physical time**) only
applies when restarting a computation.

For example, if a computation has already completed 100 iterations and is then
restarted, users can either set **Additional time steps** to 100 or increase
**Number of time steps** to 200 in order to perform 100 more iterations.
Otherwise, no additional iterations will be executed.

Subpage
=======

- \subpage start_restart

<!-- ======================================================================= -->

\page start_restart Start/Restart

Checkpoint/Restart
==================

This page indicates if a computation resumes from a previous one.

\image html gui_restart.png

## Restart from checkpoint

- **Off**: start a new computation without restarting from a checkpoint.
- **On**: restart the computation from a user-selected checkpoint directory
  located in a `<case>/RESU/<run_id>` folder.
- **Automatic**: automatically select a checkpoint directory from the available
  `<run_id>` directories located in the `<case>/RESU/` folder.

## Determine restart behavior relative to the mesh: 

- **Unmodified** (disable preprocessing): read mesh directly from the restart
  checkpoint directory, with no additional preprocessing.
- **Different mesh** (interpolate): execute standard import or generation and
preprocessing steps, interpolating field data from the restart checkpoint mesh.
- **Rebuild same mesh**: Do not use mesh in restart checkpoint directory even if
  present, importing and/or reprocessing mesh as for initial run. This can be
  useful when saving the mesh modified by preprocessing was disabled in the
  previous run (presumably to save disk space or I/O time for large runs), and
  applying the same preprocessing steps will rebuild the mesh matching the other
  restart files.
- **Automatic** (unmodified if present): same as **Unmodified** if
  `restart/mesh_input.csm` is present, **Rebuild same mesh** otherwise

In complex cases where the `restart/mesh_input.csm` file does not match the
other files in the restart directory, but a matching `mesh_input.csm` can be
located or generated, placing it in the matching checkpoint directory post-hoc
is suggested. An alternative is also to use the `cs_restart_map_set_mesh_input`
user_defined function to specify the path of the file which should be used,
combined with the **Different Mesh** option above (as the path defined through
the user function will override the default one).

## Calculation on frozen dynamic

When this option is enabled, the thermal evolution and species transport are
computed while keeping the velocity, pressure, and turbulence fields constant.

## Advanced options

- **Read auxiliary restart file**: when disabled, additional restart data are
  not read from the auxiliary restart file, including the time step, reference
  point, mass fluxes, boundary condition coefficients, extrapolated source
  terms, time moments, and fields associated with specific physical models.

- **Frequency of restart checkpoints**: controls how often restart files are
  generated. Available options are:
  - **Never**
  - **At the end of the computation**
  - **Four restart checkpoints** (default)
  - **User-defined frequency**
