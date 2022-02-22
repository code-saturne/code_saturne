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

\page cs_ug_quick_start Quick start

We assume in this section that the user has the calculation data files
(calculation setup) at his disposal, or has already prepared it following
for instance the step-by-step guidance provided in the code_saturne tutorials.
The steps described below are intended to provide a way to quickly run through the
Graphical User Interface (GUI).

Setting up your environment {#sec_prg_environement_cs}
===========================

It is recommended before running code_saturne to define an alias to the
`code_saturne` script, for example:

```
alias cs='${install_prefix}/bin/code_saturne'
```

where '${install_prefix} is the base directory where code_saturne and its components
have been installed.

More detailed instructions are provided in the dedicated
[setting up your environment](@ref cs_ug_user_settings) section.

Preparing a case
----------------

The second thing is to prepare the computation directories. For instance, the
`T_JUNCTION` study directory, containing a single calculation directory
`CASE1`, will be created by typing the command (see [sec_prg_cscreate]):
```
code_saturne create -s T_JUNCTION
```
The mesh files should be copied in the `MESH` directory (though they may also be
selected from another directory  (see [sec:prg_stepbystepcalculation])
and the C and Fortran user files necessary for the calculation in the directory
`CASE1/SRC`.  Finally, the calculation data file `setup.xml` managed  by the GUI
should be copied to the `CASE1/DATA` directory.
Once these steps are completed, the user should go in the directory `CASE1/DATA` and run
```./code_saturne gui  setup.xml``` to load the calculation file into the interface.
A window similar to the one below will appear.

\anchor fig_3_e1
![GUI main page](gui_case_dir.png "GUI Main Page")

Click on the "Run computation" button in the toolbar, as shown
[below](@ref fig_43_e1).
After having chosen the number of processors,
press "start calculation" to run the calculation.

\anchor fig_43_e1
![Prepare execution](gui_prepare_execution.png)

If no problem arises, the simulation results can be found in the `CASE1/RESU`
directory and be read directly by ParaView or EnSight in
`CASE1/RESU/<YYYYMMDD-hhmm>/postprocessing`.
The main calculation log can be found in the file `<YYYYMMDD-hhmm>/run_solver.log`.

Troubleshooting
---------------

If the calculation does not run properly, the user is advised to check the
following points in `CASE1/RESU/<YYYYMMDD-hhmm>`

- if the calculation stops in the preprocessor, the user should check for
  error messages in the `preprocessor*.log` file;
- if the problem is related to boundary conditions, the user should visualize
  the `error.ensight` file with ParaView or EnSight;
- if the calculation stops in the code_saturne solver, the user should look for
  messages at the end of the files `run_solver.log` and `error*`.
  In addition, the user can track the following keywords in the log
  (these are specific error signals):
  * `SIGFPE`: a floating point exception occurred. It happens when there is a
              division by 0, when the calculation did not converge,
              or when a real number reached a value over *10<sup>300</sup>*.
              Depending on the host architecture and build options,
              this type of exception may be caught or ignored.
  * `SIGSEGV`: a memory error such as a segmentation violation occurred.
               An array may have exceeded its allocated memory size and a
               memory location in use was overwritten.

In order to easily find and fix the problem, it is also recommended to use a debug
build of code_saturne (see the installation documentation), possibly in combination
with the use of the Valgrind tool or a debugger. The use of such a tool can be
specified in the GUI in the advanced options of the item "Prepare batch calculation",
or without the GUI, in the `cs_user_scripts.py` file (which can be found in
`DATA/REFERENCE` and should be copied in `DATA`, see the
[step by step calculation instructions](@ref sec_prg_stepbystepcalculation)).
