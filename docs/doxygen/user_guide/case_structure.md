<!--
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

\page cs_case_structure Case directory structure

[TOC]

Introduction {#cs_case_structure_intro}
============

This page describes the directory structure used by code_saturne to
handle case data, model and numerical parameters, and run options.

To create a case, either the GUI or the```code_saturne create``` command
can be used. As usual for all code_saturne commands,
the ```code_saturne create --help```
will list the available options. More details are provided in the
dedicated [case generator](@ref #sec_prg_cscreate) section.

Organization and set-up of code_saturne computations is based on several concepts:

* A **study** is a group of related *cases*, usually sharing a common mesh or family of meshes;
* A **case** represents a given computation, and represents the base level of a simulation;
* A **run** represents a given computational step. Successive runs in a given *case* may
  represent simple restarts of the same setup to compute additional time steps, may
  involve setup modifications, or both.

Standard directory hierarchy {#case_structure_standard_hierarchy}
----------------------------

Studies, cases, and runs are usually organized in a standardized directory structure, allowing
a good traceability of computations while trying to minimize duplication of large files. While
the execution of the solver itself could uses a 'flat' structure, with prescribed file names,
the GUI and high level `code_saturne` commands assume the standard structure is used.

The standard architecture for the simulation studies is:

An optional study directory containing:

* A `MESH` directory containing the required mesh(es)
* A `POST` directory for post-processing scripts (not used directly by the code)
* One or several calculation directories

Every calculation directory contains:

* A `DATA` directory for the setup data
  (xml file from the GUI, input profiles, thermo-chemical data, ...).
* A `SRC` directory for user-defined functions (C or Fortran)
* A `RESU` directory for the computation results

To improve the calculation traceability, the files and directories
sent to `RESU` after a calculation are  placed in a sub-directory
named after that run's `id`, which is by default based on the run date
and time, using the format: *YYYYMMDD-hhmm*.
It is also possible to force a specific run id, using the `--id`
option of `code_saturne run`.

Below are typical contents of a case directory *Case1* in a study *Study*:

![Example study and case directory structure](cs_directory_structure.svg)

### Coupled computation hierarchy {#case_structure_coupled_hierarchy}

For coupled calculations, whether with code_saturne itself or Syrthes, each coupled
calculation domain is defined by its own directory (bearing the same name as the
domain), but results are placed in a `RESU COUPLING` directory, with a sub-directory
for each run, itself containing one sub-directory per coupled domain.

Coupled cases are run through the standard the code_saturne run command, but
require a coupling parameters file (`coupling parameters.py`) specified using
the `--coupling` option. The run command must be called from the top-level (`Study`)
directory, so an additional `Study/run.cfg` file is also used in this case.
Note that case-local scripts (such as `Study/Domain_1/DATA/run.cfg`)
are still used by the master script to determine which parameter file to use.
So in the coupled case, calculation results would not be placed in
`Study/Domain_1/RESU/YYYYMMDD-hhmm`, but in `Study/RESU_COUPLING/YYYYMMDD-hhmm/Domain_1`,
with the summary file being directly placed in `Study/RESU_COUPLING/YYYYMMDD-hhmm`
(as it references all coupled domains).

The following example illustrates a coupled case with one code_saturne domain (named *Fluid*)
and one Syrthes domain (named *Solid*):

![Example study directory structure with Syrthes coupling](cs_directory_structure_with_syrthes.svg)

Files Copied and referenced by a run {#case_structure_run_copy}
------------------------------------

When running a case, a new entry in `RESU` (or `RESU_COUPLING`) is generated,
based on a *run_id*, which can either be given, or automatically generated (in which
case it is based on the *YYYYMMDD-hhmm* (year/month/day-hour/minute) format, possibly
extended by `_1`, `_2,` ... if multiple runs are launched during the same minute.

Once a computation is started or submitted, it should be possible to modify
files in the case setup without interfering with the running or pending computation,
and for "quality control" considerations, it is often useful to keep a trace of
files used, so the following rules apply when a computation run is prepared:

* Files directly in `DATA` are copied to `RESU/<run_id>`, except for `run.cfg`
  (which is not copied directly but may be used to generate `runcase`)
  if the currently active XML file is not named `setup.xml`, a symbolic link
  named `setup.xml` is added to that file, as the solver expects that name.
* Files directly in `SRC` are copied to `RESU/src`
* Sub-directories of `DATA` and `SRC` are ignored
* for large files or directories referenced in `setup.xml` or `user_scripts.py`, such
  as `mesh_input*`, `checkpoint`, or `partition_input`, symbolic links are used in
  the run directory rather than a full copy. The link may have a different name: for
  example, a `checkpoint` from a previous run is linked as `restart` for the new run.

In most cases, the solver is run using `RESU/<run_id>` as its work directory.
For coupled cases, `RESU/<run_id>/<domain>` is used for each domain, so as to avoid
write conflicts.

### Reserved file and directory names {#case_structure_reserved_names}

As can be seen in the examples above, some file and directory names have special roles,
so should not be used as specific user-defined inputs or outputs:

Name             |  Type                 | Role
-----------------|-----------------------|----------------------------------------------------
DATA             | base input directory  | data and setup definition files for the computation
SRC              | base input directory  | user-defined sources for the computation
run.cfg          | input file            | definition of run resources and job allocation
setup.xml        | input file            | main computational setup
user_scripts.py  | input file            | additional computational setup
RESU             | output directory      | directory in which `<run_id>` run outputs are generated
RESU_COUPLING    | output directory      | same as `RESU`, for coupled cases
runcase          | generated script      | generated script for job submission
run_solver       | generated script      | generated low-level script for computation execution
mesh_input       | input directory       | directory of imported meshes
mesh_input.csm   | input file            | imported mesh
monitoring       | output directory      | directory for probe history
postprocessing   | output directory      | directory for post-processing and visualization output (EnSight Gold, MED, CGNS)
partition_input  | input directory       | link to `partition_output` from previous run
partition_output | output directory      | directory for partition maps
restart          | input directory       | link to `checkpoint` from previous run
checkpoint       | output directory      | directory for checkpoint files
src              | input directory       | directory for compiled source files (in run directory)
compile.log      | output file           | compilation log for files in `src`
performance.log  | output file           | summary of performance data
preprocessor.log | ouptut file           | mesh import log
run_solver.log   | output file           | main computation log
setup.log        | output file           | setup options log
listing          | output file           | symbolic link to run_solver.log (legacy name)
residuals.csv    | output file           | residuals per iteration
timer_stats.csv  | output file           | main timings per iteration
summary          | output file           | global execution and environment summary

### Temporary execution directory {#case_structure_scratchdir}

In special cases, it is possible to define a separate *scratch* execution directory,
by setting a 'CS_SCRATCHDIR` environment variable or defining such a directory in the
general `code_saturne.cfg` settings.
In this case, results are copied from the scratch directory to
the run output directory during the `finalize` stage of a computation.

This is recommended only if the compute environment includes different
file-systems, some better suited to data storage, others to intensive I/O.
If this is not the case, there is no point in running in a scratch directory
rather than the results directory, as this incurs additional file copies.

If the `CS_SCRATCHDIR` environment variable  is defined, its value has priority
over that defined in the preference file, so if necessary, it is possible to define
a setting specific to a given run using this mechanism.

\warning in case of an error, the temporary directories are not deleted
after a calculation, so that they may be used for debugging. They may then
accumulate and lead to loss of usable disk space.
It is therefore essential to remove them regularly.

Case generator {#sec_prg_cscreate}
--------------

The `code_saturne create` case generation command  automatically creates
a study or case directory according to the typical architecture and copies
the required files.

The syntax is briefly described here:

```
code_saturne create --study STUDY CASE_NAME1
```
creates a study directory `STUDY` with case sub-directory
`CASE_NAME1`. If no case name is given, a default case directory called
`CASE1` is created. While:

```
code_saturne create --case Flow3 --case Flow4
```
executed in the `STUDY` directory adds the case directories `Flow3` and `Flow4`.
Whenever multiple cases are created simultaneously, it is assumed they may be
coupled, so top-level `run.cfg` and `coupling_parameters.py` files and a
`RESU_COUPLING` directory are also created.

In each case's `DATA` directory, reference (minimal) `setup.xml` and
`run.cfg` files are generated.

Unless the `--noref` option is used, under `DATA`, a `REFERENCE` sub-directory
containing a `cs_user_scripts.py` advanced settings template and
examples of thermochemical data files used for pulverized coal combustion,
gas combustion, electric arcs, or a meteorological profile.
The files to be actually used for the calculation must be copied directly in
the `DATA` directory and its name may either be unchanged, or be referenced using
the GUI or using the [usppmo](@ref usppmo) user subroutine.
In same manner, under the `SRC` directory, a sub-directory named `REFERENCE`
containing all the available user-defined function templates and a
the sub-directory named `EXAMPLES`  containing multiple examples are copied.

As a rule of thumb, all files in `DATA` or `SRC` except for the
`code_saturne` script are copied for use during code execution,
but subdirectories are not.
