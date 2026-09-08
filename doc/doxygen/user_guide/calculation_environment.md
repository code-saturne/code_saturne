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

\page cs_ug_calculation_environment Calculation environment

Main page
=========

\image html gui_calc_env_directories.png

This page displays the study name and path, as well as the case name and its
associated subdirectories (DATA, SRC, and RESU). A different case can be
selected using the folder selection dialog box.

Subpages
========

# Notebook

Notebook parameters are user-defined variables used to parameterize a
code_saturne case and share values between the GUI and user code.

\image html gui_calc_env_notebook.png

The Notebook subpage allows users to create, modify, and delete notebook
variables. It is also possible to import a list of variables from a CSV file.

In addition to the variable name and value, several options are available:
- OpenTurns Variable: indicates that the variable can be used as an input or
  output parameter in OpenTurns workflows for uncertainty quantification and
  sensitivity analyses.
- Editable: allows the notebook value to be modified during the execution of
  the case.
- Read at Restart: initializes the notebook value from a previous computation
  when restarting a simulation.
- Print to Default Log File: prints the selected variables to
  `run_solver.log` at each time step.

# Time tables

Time tables are used to define time-dependent input data, allowing boundary
conditions, source terms, physical properties, or user parameters to vary
during a simulation according to tabulated values.

\image html gui_calc_env_time_tables.png

The Time Tables subpage allows users to import and manage time tables from
CSV, TXT, or DAT files.
