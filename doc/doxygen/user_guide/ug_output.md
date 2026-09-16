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

\page cs_ug_output Managing and analyzing the output

[TOC]

Analyzing the output
====================

Checking the convergence
------------------------

Checking the convergence is difficult to automate, but code_saturne
provides several tools to help manage this.

### Main log file

In the *run_solver.log* file for each run's output, a _convergence_ section similar
to the example below is available:

```
   ** INFORMATION ON CONVERGENCE
      --------------------------
-----------------------------------------------------------------------------
   Variable    Rhs norm      N_iter  Norm. residual      Drift  Time residual
-----------------------------------------------------------------------------
c  Velocity     0.43576E+00      16   0.43716E-02   0.28615E-07 0.74561E-02
c  Velocity[X]                                      0.17295E-07
c  Velocity[Y]                                      0.44914E-08
c  Velocity[Z]                                      0.68285E-08
c  Pressure     0.55189E-02      47   0.31730E-05   0.10737E-03 0.61715E-05
c  TurbEner     0.15629E-03      13   0.18499E-05   0.17727E-15 0.18405E-02
c  Dissip       0.33146E-04      13   0.83187E-06   0.16668E-16 0.31349E-02
c  TempC        0.10728E+07      13   0.12667E-02   0.13931E-02 0.27664E-02
-----------------------------------------------------------------------------
```

For each solved variable, it provides:
- The number of iterations required for the linear solvers (which depends on the
  chosen solver type, and may be converted to an "equivalent" cost estimation
  in some cases).
- The normalized residual.
- The *time drift*;\n
  For a given variable \f$ \varia \f$ this is usually the following term:

  \f[
    \int_\Omega \left| \der{\varia}{t} \right|^2 \Delta t \dd \Omega / \int_\Omega \dd \Omega
  \f]

  For the pressure variable it is computed in a different manner.

- The *time residual*;\n
  For a given variable \f$ \varia \f$ this is the normalized L<sup>2</sup> unsteady term:

  \f[
    \sqrt{\int_\Omega \left| \der{\varia}{t} \right|^2 \dd \Omega / \int_\Omega \left| \varia \right|^2 \dd \Omega}
  \f]

### Residuals file

A file named *residuals.csv* is usually produced, containing the residuals for solved
variables (as defined above), in an easy to plot CSV format.

### Time plots

It is recommended to place some *probes* at selected points of interest, so as to activate
time plots of the main variable values at those points.

This allows both checking how values evolve over time in selected points and if that
evolution is regular or "noisy".
