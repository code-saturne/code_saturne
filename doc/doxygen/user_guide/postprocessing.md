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

\page cs_ug_postprocessing Postprocessing

# Output Control tab

Use this tab to define how often information is written to the listing (log) file.
The log file contains information related to convergence, variable values, and
the time evolution of the computation.

## Log Frequency

- **No output**: no information is written to the log file.
- **Output every *n* time steps**: write information to the log file at a
  user-defined frequency.
- **Automatic (decreasing frequency)**: the output frequency is automatically
  reduced as the number of iterations increases (default setting).


# Writer Tab

Writers are objects used to configure post-processing outputs for subsequent
visualization. **A writer must be selected to access and modify its settings.**

The **Writer** tab allows users to create, edit, and delete writers. For each
writer, it is possible to define its name, output format (EnSight, MED, CGNS,
Catalyst, or Histogram), and the directory in which result files will be stored.

## Frequency

Define the output frequency for the selected writer:

- **No periodic output**: output is generated only when explicitly requested
  (for example, at the beginning or end of the calculation).
- **Output every *n* time steps**: write results at a user-defined time-step
  interval.
- **Output every *x* seconds**: write results at a user-defined physical-time
  interval.
- **Output using formula**: control the output frequency using a user-defined
  expression.

Two additional check boxes allow output to be generated at the start and/or at
the end of the calculation.

## Time Dependency

Define the time-dependent behavior of meshes associated with the selected
writer:

- **Fixed mesh**: all associated meshes remain unchanged throughout the
  calculation.
- **Transient coordinates**: the mesh connectivity remains unchanged, but node
  coordinates may vary over time (deforming mesh).
- **Transient connectivity**: the mesh topology and the number of elements may
  vary over time.

## Options

- **Separate sub-writer for each mesh**: generate separate result files for
  each mesh associated with the selected writer.
- **File type**: binary native, binary big-endian, text.
- **Polygons**: display, discard, subdivide.
- **Polyhedra**: display, discard, subdivide.

The **Polygons** and **Polyhedra** options are useful when exporting results to
formats or visualization tools that do not support polygonal or polyhedral
elements. In such cases, these elements can be discarded or automatically
subdivided into supported element types.

# Mesh Tab

The **Mesh** tab allows users to create, edit, and delete mesh parts. For each
mesh part, it is possible to define a name, select its type (cells, interior
faces, boundary faces, volume zones, or boundary zones), and specify a
selection criterion. **A mesh part must be selected to access and modify its
settings.**

## Variables

Select the variables to be exported on the selected mesh.

- **Auto** (default): automatically export the standard set of variables
  available for the selected mesh.

## Associated Writers

Associate the selected mesh part with one or more previously defined writers.
The mesh data will then be exported according to the settings of the associated
writers.

# Monitoring Points Tab

The **Monitoring Points** tab allows users to create, import, edit, and delete
monitoring points (probes). Monitoring points can be used to record the
evolution of variables at specific locations during the calculation.

It is also possible to define the output frequency and select the format of the
result files.

Two additional options are available:

- **Limit to one probe per cell**: if several probes are located within the
  same cell, only one probe is retained.
- **Interpolation**: interpolate variable values at the exact probe locations
  instead of using the values of the containing cells.

# Supages

- \subpage calculator
- \subpage additional_user_arrays
- \subpage time_averages
- \subpage volume_solution_control
- \subpage surface_solution_control
- \subpage profiles
- \subpage balance_by_zone

<!-- ======================================================================= -->

\page calculator Calculator

# User-defined calculator functions

The calculator subpage allows users to define custom variables, fields, or
post-processing quantities by combining existing variables through mathematical
expressions, without modifying the source code.

<!-- ======================================================================= -->

\page additional_user_arrays Additional user arrays

# User-defined fields

The user-defined fields subpage allows users to create additional scalar,
vector, or tensor fields that can be initialized, transported, post-processed,
or used in models and boundary conditions without modifying the source code.

<!-- ======================================================================= -->

\page time_averages Time averages

# Time Averages Tab

The **Time Averages** tab allows users to define time-averaged quantities that
are accumulated during the calculation. These averages can be based on one or
more variables and are computed starting from a user-defined time step or
physical time.

For each time average, users can specify:
- A name for the averaged quantity.
- The starting time step or physical time from which averaging begins.
- Whether the average should be restarted from a previous calculation.
- The variables to be included in the averaging process.

<!-- ======================================================================= -->

\page volume_solution_control Volume solution control

# Volume Solution Control Tab

The **Volume Solution Control** tab allows users to select which volume fields
(solution variables, physical properties, and auxiliary quantities) are written
to the listing, post-processing, and monitoring-point output files during the
simulation.

## Iterative process error estimators

TODO

\page surface_solution_control Surface solution control

# Surface Solution Control Tab

The **Surface Solution Control** tab allows users to select which surface fields
(wall stresses, wall distance and thermal quantities) are written to the output
files associated with the corresponding writers.

<!-- ======================================================================= -->

\page surface_solution_control Surface solution control

<!-- ======================================================================= -->

\page profiles Profiles

# Definition of 1D profiles

The **Definition of 1D Profiles** tab allows users to define profiles of selected
variables along user-defined line segments. For each profile, it is possible to
specify:
- A file name for the output.
- The output frequency.
- The line definition with the mathematical expression editor.
- The number of sampling points along the line.
- The volume fields (solution variables, physical properties, and auxiliary
  quantities) to be written to the profile output file.
**A profile must be selected to access and modify its settings.**

As in the **probes** tab, two additional options are available:
- **Limit to one probe per cell**: if several probes are located within the
  same cell, only one probe is retained.
- **Interpolation**: interpolate variable values at the exact probe locations
  instead of using the values of the containing cells.

<!-- ======================================================================= -->

\page balance_by_zone Balance by zone

## Pressure Drop Definition

The **Pressure Drop Definition** tab allows users to define regions based on
selection criteria where pressure-drop balances are evaluated during the
simulation. These balances can be used to analyze pressure losses across
selected volume zones or user-defined regions.
 
## Scalar Balance
 
The **Scalar Balance** tab allows users to define regions based on selection
criteria where balances are computed for selected scalar quantities.
