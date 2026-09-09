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
