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

\page cs_ug_catalyst In situ postprocessing with ParaView Catalyst

[TOC]

Basic principles
================

It is possible to can generate postprocessing output directly from a running
code_saturne simulation using
[ParaView Catalyst](https://kitware.github.io/paraview-catalyst),
provided support for this tool is specified at build time (see
[installation documentation](@ref autotoc_md0)).

To define in-situ postprocessing output, the following steps are required:

- Produce an initial output using a post-hoc output format (e.g. EnSight, MED, or CGNS).
- Define your visualization in ParaView.
- Export a catalyst Python script from ParaView.
  * Modify that script if necessary or desired.
- Copy the script to a case's DATA directory.
- Define an output writer using the "catalyst" format (or change the format
  of an existing one).
- Run the computation.

Step by step example {#cs_ug_catalyst_step_by_step}
====================

Let us consider that after an initial code_saturne run (at least over a few
time steps), a visualization pipeline is set up using ParaView, such as
that below (using the code_saturne Tee_Junction tutorial, and based on
the default EnSight format output):

Preparing the Catalyst script
-----------------------------

\anchor fig_ug_catalyst_01_initial_pipeline
<img src="ug_catalyst_01_initial_pipeline.png" width="90%" alt="">
<div class="caption">Initial Visualization</div></div>

\anchor ug_catalyst_02_rename_input
<img src="ug_catalyst_02_rename_input.png" width="90%" alt="">
<div class="caption">Renaming the pipeline input</div></div>

\anchor fig_ug_catalyst_03_local_path
<img src="ug_catalyst_03_local_path.png" width="90%" alt="">
<div class="caption">Prune absolute path</div></div>

\anchor fig_ug_catalyst_04_add_ghost_cells_1
<img src="ug_catalyst_04_add_ghost_cells_1.png" width="90%" alt="">
<div class="caption">Insert ghost cells filter</div></div>

\anchor fig_ug_catalyst_05_add_ghost_cells_2
<img src="ug_catalyst_05_add_ghost_cells_2.png" width="30%" alt="">
<div class="caption">Insert ghost cells dialog</div></div>

\anchor fig_ug_catalyst_06__extractor_1
<img src="ug_catalyst_06_extractor_1.png" width="90%" alt="">
<div class="caption">Define extractor</div></div>

\anchor fig_ug_catalyst_07__extractor_2
<img src="ug_catalyst_07_extractor_2.png" width="90%" alt="">
<div class="caption">Extractor dialog</div></div>

\anchor fig_ug_catalyst_08_save_state_1
<img src="ug_catalyst_08_save_state_1.png" width="90%" alt="">
<div class="caption">Save Catalyst state</div></div>

\anchor fig_ug_catalyst_09_save_state_2
<img src="ug_catalyst_09_save_state_2.png" width="80%" alt="">
<div class="caption">Catalyst state file dialog</div></div>

\anchor fig_ug_catalyst_10_save_state_3
<img src="" width="30%" alt="">
<div class="caption"></div></div>
![Catalyst state settings dialog](ug_catalyst_10_save_state_3.png)

Associating the Catalyst script with code_saturne
-------------------------------------------------

\anchor fig_ug_catalyst_20_alternative_cleanup
<img src="ug_catalyst_20_alternative_cleanup.png" width="80%" alt="">
<div class="caption">Direct edit of Catalyst script</div></div>

\anchor fig_ug_catalyst_30_add_cs_writer
<img src="ug_catalyst_30_add_cs_writer.png" width="85%" alt="">
<div class="caption">Add code_saturne Catalyst writer</div></div>

\anchor fig_ug_catalyst_31_associate_cs_mesh
<img src="ug_catalyst_31_associate_cs_mesh.png" width="85%" alt="">
<div class="caption">Associate postprocessing mesh to Catalyst</div></div>

Generated output
----------------

\anchor fig_ug_catalyst_40_stream_lines_50
<img src="ug_catalyst_stream_lines1_000050.png" width="90%" alt="">
<div class="caption">Catalyst output at 50 time steps</div></div>

\anchor fig_ug_catalyst_41_stream_lines_200
<img src="ug_catalyst_stream_lines1_000200.png" width="90%" alt="">
<div class="caption">Catalyst output at 200 time steps</div></div>

Additional resources
====================

The top-level [**ParaView Catalyst web page**](https://kitware.github.io/paraview-catalyst/)
contains links to most the the resources described below, as well as descriptions
of many use cases and example.

A detailed description on Catalyst workflows and instrumentation is provided in
a [Kitware blog post](https://www.kitware.com/real-time-insight-with-paraview-catalyst-a-hands-on-guide-part-1-the-basics/). This example also shows how it is possible to
define a Catalyst script from scratch, without a prior computation, and generate
simple statistics. In fact, using this approach should allow using many VTK and
other postprocessing libraries accessible from Python.

In this case, remember that in MPI computations, code_saturne data is distributed
across ranks so Python MPI operations may be necessary to ensure that the generated
output accounts for this correctly.

All Kitware [blog posts related to Catalyst](https://www.kitware.com/tag/catalyst/) and
the [ParaView in-situ support Channel](https://discourse.paraview.org/c/in-situ-support/8)
can also be of great interest.

The [ParaView Catalyst guide](https://kitware.github.io/paraview-catalyst/guide/concepts.html) contains some sections which are of interest mostly
for developers and maintainers of sofware (such as code_saturne) implementing
Catalyst output (such as the *Getting Started* section), with other sections
(especially *Using ParaView to Create the ParaView Catalyst Script* section)
providing additional examples on Catalyst script generation and use.

The [Catalyst section of the ParaView User's guide](https://docs.paraview.org/en/latest/Catalyst/index.html) also provides detailed information and
[debugging tips](https://docs.paraview.org/en/latest/Catalyst/debugging.html#).
