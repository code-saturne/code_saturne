<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

\page cs_ug_troubleshooting What if a calculation fails ?

[TOC]

What if a calculation fails ? {#sec_ug_troubleshhoting}
=============================

The first things to do in case of a failure of the calculation are
the following:

- First of all, check for messages in the terminal from which the code was
  launched, or if run under a batch system, in the job error and/or standard
  output files.

  * These should either contain messages explaining the error, or (most often)
    indicate which other files should be checked for detailed messages.

  * For example, when errors such as out of bounds arrays accesses are detected
    by the Fortran runtime, they will appear here instead before the code_saturne
    error detection and logging mechanism can catch them.

- If errors occur when importing a mesh, the file to check will be the
  `preprocessor.log`.

- If they occur during the computation, browse the indicated log files
  (especially `run_solver.log`) and error files if present (`error*`).

  * In general, the content of the main `error` file is also present at the
    end of the `run_solver.log` file.

  * For parallel runs, if `error_r*` files are present, they may contain more
    useful information than the main `error` file or `run_solver.log`.

    - When a fatal error occurs on a given MPI rank _p_, the matching process
      waits for a few seconds before generating an `error_r`_p_ file and
      interrupting the whole computation (by calling `MPI_Abort`), except for
       rank 0; which generates an `error_r` file immediately.

      So if no error had occurred on rank 0, it will simply report that is has
      been interrupted by the environment (initiated by the process generating
      the error), and the useful information will be in that rank's error file.
      If the initial error also occurs on rank 0 (such as an incorrect setting
      being detected), the only error file generated will be `error`, for better
      readability.

      This approach was chosen so as to try to generate as few as possible
      redundant `error_r*` files, while still assuming that in the worst case,
      it is preferable to have extra error files than not have the required ones
      (for example if the delay is not adapted to a given machine).

- In most cases, the `error_r*` files should contain a stack trace, indicating at
  least approximately where the error occurred, and what type of error occurred.

- In some cases (such as incorrect boundary condition definitions), a
  postprocessing output (`postprocessing/error*`) output may be generated, so as
  to allow vissualization of the error location.

- Check that the version of the code is the expected one and that it is
  effectively accessible on the machine.

- Simplify the case: remove as many user-defined functions as possible to better
  identify the source of the failure (for example, the simulation should be
  checked with constant physical properties before applying variable values).

- If the code has not been used on this computer before, and errors seem to
  be related to the system, the case should be transferred to a machine that the
  user is more familiar with so as to avoid any cause associated with the
  operating system or with the installation of the code.

When the computation starts to run but diverges {#sec_ug_troubleshhoting_diverg}
-----------------------------------------------

If the computation starts to run but diverges, the following points
should be checked:

- If a _runaway_ (.i.e. diverging) computation error is reported, the code crashes
  at the first iteration, or field values in the log seem unexpected, increase
  the level of verbosity for the resolution of the various equations (in the GUI,
  in the _Numerical parameters/Equation parameters_ section, and re-run the case.
  This should provide more info on convergence an possibly help identify the
  source of the problem.

- Check the Courant number in the logs. If it is significantly higher than 1,
  and increases from one time step to the next, reducing the time step value
  (or reference time step) may help.

- Check that the initial state is “reasonable”: it is not essential to prescribe
  an initial solution close to the solution that is searched for, however an
  non-physical initial estimate of the solution may lead to a failure of the
  computation.

  * For example: initial velocity oriented in the opposite direction,
    initial turbulent variables incoherent between each other, too large initial
    density gradients, ...).

  * Some difficulties can be associated with the choice
    of the initial state (they are normally characterized by a transient that can
    be evacuated through the outlet of the domain by setting a time step value
    corresponding to a large Courant number during the first few time steps
    (for example: the standard value of the time step may be multiplied by a
    factor of 10).

- If one suspects a difficulty associated with the mesh quality, or with meshes
  with strong non-orthogonality (such as tetrahedral meshes), one may use the
  gradient computations based on least-squared with extended cell neighbors
  (at least _opposite adjacent cell centers_).

  * The most robust (but most costly) option is to use the full extended
    cell neighbors.
  * Using a Green-Gauss based gradient (iterative or with least-squares based face
    values) is usually preferred, but pure least-squared gradients usually
    exhibit smoother behavior (better respecting th maximum principle) and may
    be more robust.

- With tetrahedral meshes, extruding the mesh to obtain at least one layer of
  prisms at the outlet, or using an imposed pressure outlet rather than a
  free outlet condition can improve behavior.

The following steps may also be taken if the above are not sufficient to
solve the issue:

- If one observes oscillations on the advected variables, it is advised to
  modify the convective scheme (in the worst, most unstable case, a first-order
  upwind scheme should be used for all the advected variables).

- If one observes that the pressure iterative solver takes a long time to
  converge or that overshoots of the velocity occur locally where the quality
  of the mesh may be questionable, it is advised to under-relax the pressure
  increments (using 0.7 for the relaxation of pressure increase for example).

- If one observes that the convergence is difficult on a low quality mesh, the
  stability of a RANS computation with a turbulent viscosity model
  (k-epsilon et k-omega) may be improved by switching off the flux reconstruction
  for turbulence variables (k, epsilon, or omega).

- If spurious velocities appear in a region where the pressure gradient is
  discontinuous because of the presence of a stratification or of a head loss
  zone, make sure the "improved pressure interpolation" option
  (`iphydr = 1`) is used (it is activated by default).

Safe mode parameters {#sec_ug_troubleshhoting_safe_mode}
====================

An example of a “safe mode” set of parameters is provided hereafter:
with this choice, the priority is given to robustness and low computational cost
over accuracy:

- Turbulence model: k-epsilon with linear production.
- Convective scheme: SOLU for the velocity, UPWIND for the turbulence and
  the scalars.
- No flux reconstruction for the turbulence variables (k and epsilon).
- Threshold for the convergence of the iterative solvers: 10<sup>-5</sup>.
- Least-squares gradient reconstruction with extended neighborhood.
- Pressure increments under-relaxation (0.7 instead of 1).

If the computation runs correctly with these options but not with default
options, improving the mesh quality at least locally may be necessary
to enable running with default, "higher-precision" options.

Some users also report success initializing a computation with “safe”
options (also deactivating flux reconstruction an using an upwind convective
scheme for all variables), and restoring the higher precision options
progressively after 50 or more iterations.
