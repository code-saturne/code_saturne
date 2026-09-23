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

[TOC]

General concepts {#cs_ug_postprocess_intro}
================

Logging and post-processing in code_saturne is handled so as to allow
combining the ease of automatic outputs for the main variables and more advanced
user-defined values or extracts.

Both the GUI and several functions present in \ref cs_user_postprocess.cpp allow
for managing post-processing output.

Writers and meshes
------------------

The main concepts are those of \em writers and \em meshes, which must be
associated to produce outputs.

- A \em writer combines the definition of an output type, frequency, path, and name.
  One or more \em writers can be defined using the GUI and the
  \ref cs_user_postprocess_writers user function.

- A \em mesh is based on a subset of the the computational mesh, or point
  sets such as particles or probe sets. One or more \em meshes can be defined
  using the GUI and the \ref cs_user_postprocess_meshes user function.

The combination of writers and meshes allows generating chronological outputs
in *EnSight*, *MED*, or *CGNS* format, as well as in-situ visualization
using [ParaView Catalyst](https://www.paraview.org/in-situ) or
ensemble data output to [Melissa](https://melissa-sa.github.io).

Probe sets and profiles
-----------------------

Probe sets allow defining a number of probes at user-defined locations.
The values of selected variables at these points may be output to specific
writers, usually *time plots* (CSV or regular text files with one column
per probe position and one line line output time step).
These are typically used:

- To monitor the computation's convergence at chosen points.
- To easily extract values at points for which measured data may be available.

Profiles are defined as a special type of probe set, for which a curvilinear
coordinate (the local profile coordinate) is defined, and are typically used not
with *time plots*, but with simple *plots* (also in CSV or text format)
where one file is associated to a single time step, and columns correspond
to different variables, and lines to different positions along the profile.

Probe sets and profiles are handled as a special category of post-processing
meshes, so they can be also be used in combination with most other types
of writers.

Output variable selection
-------------------------

For the main computed values defined as *fields*, the GUI allows  defining
whether a given field with be associated using the default output meshes,
the default probe set, and the main log.

This can also be handled using the \ref post_vis and \ref log field key values.

Also, the "auto variables" flag may be assigned to additional post-processing meshes
so that output of fields behaves in a manner similar to that of the
default post-processing mesh in the same category (for example the main volume
mesh for a volume extract).

To output additional values, whether simply restricting output to a small subset of field
or adding values based on specific computations or formulas, see the
\ref cs_user_postprocess_values user function and associated examples
(in [Definition of the variables to post-process](@ref cs_user_postprocess_h_var_p)).

Output Control
==============

In the GUI, the matching tab allows setting global output-related options.

## Log Frequency

The logging output frequency may be managed using the GUI, as shown below.
The first 10 time steps are always logged.
When running many time steps, it recommended to use a logging period
greater than 1 (20 or 100 are often sufficient), as logging each time step
with even a moderate verbosity can lead to large `run_solver.log` files.

\anchor fig_gui_output_log
![Logging parameters](gui_output_log.png)

Note that logging verbosity is also influenced by the verbosity used for solved
variables.

- **No output**: no information is written to the log file.
- **Output every *n* time steps**: write information to the log file at a
  user-defined frequency.
- **Automatic (decreasing frequency)**: the output frequency is automatically
  reduced as the number of iterations increases (default setting).

If desired, finer-grained control is possible through user-defined functions,
using the `cs_log_iteration_set_interval` or `cs_log_iteration_set_automatic`
functions, and the associted control structure may be accessed using the
`cs_log_iteration_get_time_control` function.

Writer
======

Writers are objects used to configure post-processing outputs for subsequent
visualization. in the GUI, **A writer must be selected to access and modify its settings**,
such as in the following example:

\anchor fig_gui_output_writers
![Management of postprocessing writers](gui_output_writers.png)

The **Writer** tab allows users to create, edit, and delete writers. For each
writer, it is possible to define its name, output format (EnSight, MED, CGNS,
Catalyst, or Histogram for example), and the directory in which result
files will be stored.

\remark

When a format is not available on the code_saturne build associated
with the GUI, its name appears in **red** in the list of available formats:
this allows warning the user that it may not be available, but allows setting up
a computation which will run on a different build including support for that
format.

\remark

The default writer cannot be "undefined", though it may be renamed and its
settings modified. A writer with no associated mesh is inactive, so removing
the default writer from the list of writers associated
to all meshes (see following section) deactivates the default output.

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

## In_situ output with ParaView Catalyst

For the proper use of a Catalyst writer, see the detailed
[step-by step instructions and recommendations](@ref cs_ug_catalyst).

Note that it is currently possible to define only a single Catalyst writer,
at least when using legacy Catalyst (i.e. Catalyst1).

Mesh
====

Postprocessing meshes can be handled in the **Mesh** tab, which allows
defining, editing, and deletin postprocessing meshes. For each such mesh
(usually an extract of the computational mesh), it is possible to define a name,
select its associated type (cells, interior faces, boundary faces), and specify a
selection criterion or associated zone .

In the GUI, **a mesh part must be selected to access and modify its
settings.**

\anchor fig_gui_output_meshes
![Management of postprocessing meshes](gui_output_meshes.png)

Postprocessing meshes are initially associated with the default writer
when defined through the GUI.

\remark

Note that the default meshes cannot be "undefined", though they may be renamed
and their settings modified. To avoid building these meshes, simply remove
all associated writers.

## Variables

Select the variables to be exported on the selected mesh.

- **Auto** (default): automatically export the standard set of variables
  available for the selected mesh.

When automatic export is deactivated, finer control on associated output
may be obtained using the `cs_user_postprocess_values` function.

## Associated Writers

Associate the selected mesh part with one or more previously defined writers.
The mesh data will then be exported according to the settings of the associated
writers.

## User-defined functions

In order to allow the user to add an output format to the main output,
or to add a mesh to the default output, the lists of standard and user
meshes and writers are not separated. Negative numbers are reserved for
the non-user items. For instance, the mesh numbers -1
(\ref CS_POST_MESH_VOLUME) and -2 (\ref CS_POST_MESH_BOUNDARY) correspond
respectively to the global mesh and to boundary faces, generated by default,
and the writer -1 (\ref CS_POST_WRITER_DEFAULT) corresponds to the
default post-processing writer.

The user chooses the numbers corresponding to the post-processing
meshes and writers he wants to create. These numbers must be positive
integers. It is for example possible to associate a user mesh with the default
post-processing writer (-1), or to add outputs regarding the boundary
faces (-2) associated with a user writer.

For safety, the output frequency and the possibility to modify the
post-processing meshes are associated with the writers rather than
with the meshes. This logic avoids unwanted generation of
inconsistent post-processing outputs. For instance, ParaView or EnSight would not
be able to correctly read a case in which one field is output to a given part
every 10 time steps, while another field is output to the same part
every 200 time steps. If some fields should be output using different
frequencies, using separate writers allow maintaining consistent output sets.

Monitoring Points
=================

The **Monitoring Points** tab allows users to create, import, edit, and delete
monitoring points (probes). Monitoring points can be used to record the
evolution of variables at specific locations during the calculation.

The GUI currently only allows handling a single (default) probe set
for time plots. Probe coordinates can be entered directly or read
from a CSV file.
It is also possible to define the output frequency and select the format of the
result files.

Two additional options are available:

- **Limit to one probe per cell**: if several probes are located within the
  same cell, only one probe is retained.
- **Interpolation**: interpolate variable values at the exact probe locations
  instead of using the values of the containing cells.

The \ref cs_user_postprocess_probes function may be used to define additional
probe sets (including probes lying on boundary faces),

Subpages
========

- \subpage calculator "User-defined calculator functions"
- \subpage additional_user_arrays "User-defined fields"
- \subpage time_averages "Time averages"
- \subpage volume_solution_control "Volume solution control"
- \subpage surface_solution_control "Surface solution control"
- \subpage profiles "Profile"
- \subpage balance_by_zone "Balance by zone"

<!-- ======================================================================= -->

\page calculator Calculator

# User-defined calculator functions {#cs_ug_postprocess_calculator}

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

It is possible to define *time moments* for existing variables (fields) or
user-defined formulas, including both means and variances. This is limited to
"visualizable" output, so for example although the variance of a vector
(which is a tensor due to the covariance terms) may be computed, the variance
of a tensor may not (though that of specific components may always be computed).

This computation uses recurrence formulas, so at any given time steps,
the mean or variance field values are always updated and accessible as a regular
field. The computation of time moments may start at a user-defined time step of
physical time so as not to include transient initialization data. In this case, the
associated field values before actually starting the moment updates is 0.

Time moments may be defined in the GUI, as shown below, or using
the \ref cs_user_time_moments function in \ref cs_user_parameters.cpp
(see [Time moment related options examples](@ref cs_user_parameters_h_cs_user_moments)).
Variances are not yet accessible through the GUI.

\anchor fig_gui_time_averages
![Management of time averaged variables](gui_time_averages.png)

The **Time Averages** GUI tab allows users to define time-averaged quantities that
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

# Volume Solution Control

The **Volume Solution Control** tab allows users to select which volume fields
(solution variables, physical properties, and auxiliary quantities) are written
to the `run_solver.log`, post-processing, and monitoring-point output files during
the simulation.

The activation of error estimators for the standard Navier-Stokes
computation may be also be managed in this tab.

\page surface_solution_control Surface solution control

# Surface Solution Control

The **Surface Solution Control** tab allows users to select which surface fields
(wall stresses, wall distance and thermal quantities) are written to the output
files associated with the corresponding writers.

<!-- ======================================================================= -->

\page surface_solution_control Surface solution control

<!-- ======================================================================= -->

\page profiles Profiles

Profiles are based on probe sets whith an associated a curvilinear coordinate
definition.

In the GUI, they are handled in a specific manner, as they are
automatically associated with default simple plot writers (CSV or text), and the
associated variables may be selected specifically.

\anchor fig_gui_output_profiles
![Management of 1D solution profiles](gui_output_profiles.png)

This presentation difference is mainly due to historic reasons (it predates
the unification of probes, profiles, and other writers), not to a difference
in the concepts used (so the logic may seem more consistent with user-defined
functions).

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
