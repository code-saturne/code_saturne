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

\page cs_ug_case_structure Case directory structure

[TOC]

Introduction {#cs_case_structure_intro}
============

This page describes the directory structure used by code_saturne to
handle case data, model and numerical parameters, and run options.

To create a case, either the GUI or the ```code_saturne create``` command
can be used. As usual for all code_saturne commands,
the ```code_saturne create --help```
will list the available options. More details are provided in the
dedicated [case generator](@ref sec_prg_cscreate) section.

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

Coupled cases are run through the standard the code_saturne run command, which
must be called from the top-level (`Study`) directory, so an additional
`Study/run.cfg` file containing coupling metadata is also used in this case.
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
general configuration ([code_saturne.cfg](@ref cs_user_configuration_file})) settings.

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
coupled, so a top-level `run.cfg and a `RESU_COUPLING` directory are
also created.

In each case's `DATA` directory, reference (minimal) `setup.xml` and
`run.cfg` files are generated.

If the `--copy-ref` option is used, under `DATA`, a `REFERENCE` sub-directory
is created, containing a `cs_user_scripts.py` advanced settings template and
examples of thermochemical data files used for pulverized coal combustion,
gas combustion, electric arcs, or a meteorological profile.
These files are also always available in the installation directory, usually
in `${install_prefix}/share/code_saturne/data/user`.
The files to be actually used for the calculation must be copied directly in
the `DATA` directory and its name may either be unchanged, or be referenced using
the GUI or using the [cs_user_model](@ref cs_user_model) user function.
In same manner, under the `SRC` directory, a sub-directory named `REFERENCE`
containing all the available user-defined function templates and a
the sub-directory named `EXAMPLES`  containing multiple examples are copied.

As a rule of thumb, all files in `DATA` or `SRC` except for the
`code_saturne` script are copied for use during code execution,
but subdirectories are not.

Using the GUI and user-defined functions {#sec_prg_run_gui_udf}
----------------------------------------

A Graphical User Interface (GUI) is available with code_saturne.  This tool
creates or reads an XML file according to a specific code_saturne schema which
is then interpreted by the main script and by the Solver.

The GUI manages calculation parameters, standard initialization values and
boundary conditions, most available specific physical models (coal and gas combustion,
atmospheric flows, Lagrangian module, electrical model, compressible model and radiative
transfers).

Using the GUI is optional, but highly recommended. Each setting or definition
that can be specified through the GUI can also be specified in the user-defined sources.

The GUI and user-defined functions are designed to be used in combination:
it is generally preferable to use the GUI for as many settings
as possible, and resort to user-defined functions only for more complex
settings which cannot be done through the GUI. This may also include
settings with many elements that can be better defined using programmatic
loops. As a general rule, the most concise and easily verifiable approach
should be used.

In general, user functions and subroutines are called after the GUI-defined
settings for the relevant settings are loaded, so that when a given parameter
is specified both in the interface and in a user-defined function or subroutine,
the value in the user function has priority, or rather has the last word.

\warning
There are a few limitations to the changes that can be made between the GUI and
the user routines, related to which variables are solved. In particular, it is
not possible to activate a specific physical or turbulence model in the GUI and
activate a conflicting one in use functions (for example, specifying the use
of a <em>k-ε</em> model in the GUI and change it to
<em>R<sub>ij</sub>-ε</em> \ref cs_user_model.

For example, in order to set the boundary conditions of a calculation
corresponding to a channel flow with a given inlet velocity profile, the recommended
practice is to:

- Using the GUI:
  - Set the boundary conditions corresponding to the wall and the output.
  - Set a dummy boundary condition for the inlet (uniform velocity for instance)
    so as to define the appropriate zone.
- With user-defined functions:
  - set the proper velocity profile at inlet in \ref cs_user_boundary_conditions.f90.
    The dummy velocity entered in the GUI will not be taken into account as it is
    superceded by this definition (but should appear as the initial value
    in the corresponding arrays).

The GUI is launched with the `./code_saturne` command in a case's
`DATA` directory. The first step is then to load an existing parameter file (in
order to modify it) or to create a new one. By default, the assumed file name
is `setup.xml`, and changing it is not recommended (though many setting files
from older versions using various names may be encountered).

The settings available for a typical calculation are the following:

- Calculation environment: case path info,
  definition of notebook (parametric) variables.

- Mesh: definition of the mesh file(s),
  mesh preprocessing options, and mesh checking mode.

- Calculation features: choice of physical model, ALE mobile mesh features,
      turbulence model, thermal model, coupling with SYRTHES...

- Fluid properties: reference pressure, fluid characteristics, gravity.
  It is also possible to write user laws for the density, the viscosity,
  the specific heat and the thermal conductivity in the interface through
  the use of a formulae generator.

- Volume zones: variables initialization, and definition of
  the zones where to apply head losses or source terms.

- Boundary zones: definition of the boundary conditions for
  each variable. The colors of the boundary faces may be read
  directly from a `preprocessor.log*` files created by the Preprocessor
  or a `run_solver.log` file from a previous solver run.

- Time settings: time stepping scheme, number of time steps,
   management of calculation restart from a previous run.

- Numerical parameters: advanced parameters
  for the numerical solution of the equations.

- Postprocessing: visualizable output settings, time averages,
  probe sets and 1-d profile definitions.

- Performance settings: advanced parallel computing settings
  (such as partitioning and, IO options).

### User-defined function templates and examples {#sec_prg_run_udf}

Reference user-defined functions and subroutines may be found in
the `SRC/REFERENCE` subdirectory of a given case, if it was created
with the `--copy-ref` option. Otherwise, they may always be found
in the code's installation directory, usually under
`${install_prefix}/share/code_saturne/user_sources/REFERENCE`.

In a similar manner, examples may be found in
the `SRC/EXAMPLES` subdirectory of a given case if it was created
with the `--copy-ref` option, and may always be found
in the code's installation directory, usually under
`${install_prefix}/share/code_saturne/user_sources/EXAMPLES`.

Note that all C, C++, and Fortran files present directly under a case's
`SRC` directory will be used when running, while those in subdirectories
will be ignored. To use a given user-defined function, it should be
copied from the reference to `SRC` and adapted, possibly using code snippets
from the examples. To temporarily deactivate a given source file,
a recommended practice is to create a `SRC/STASH` subdirectory
nd move them to that subdirectory (rather than renaming them).

The GUI also includes a tool which can help manage and edit user-defined
functions, and check their validity (i.e. correct compilation).

### Upgrading to a newer code_saturne version

Note that when upgrading to a new code_saturne version, the GUI can
automatically update the XML file (and in the rare case where somelements cannot
be updated, a warning will be issued). Whereas although an effort is made
not to break user-defined functions too often, those functions
are guaranteed to be "stable" only within a same release series.

So for example functions written for v6.0.0 need not be changed in
bug-fix release 6.0.4, but should be at least verified and possibly updated
when moving to a release from the 6.1.* series for example.

The easier upgrade mechanism using the GUI is one of the main reasons for which
defining as many settings as possible using the GUI and keeping user-defined functions
to the minimum required is so strongly encouraged.

Run configuration file (run.cfg) {#sec_prg_run_cfg}
----------------------

For a given case, various execution-related settings can be defined
using a file named `run.cfg` in a case's `DATA` sub-directory (or in
a coupling's main directory). This file uses a format similar to classical
[`.ini`](https://docs.python.org/3/library/configparser.html#supported-ini-file-structure)
style file format, with some special section types being handled differently

A section named *section-name* is denoted by a line starting with `[section-name]`.
Key-value pairs in a section are defined using simple `key = value` or `key: value`
statements. Value definitions continuing over multiple lines must be indented by
at least one character.

Note that for key values that can take *true* or *false* values,
either `true`, `yes`, and `1` can be used for *true*, and
either `false`, `no`, and `0` can be used for *false*. Any case (capitalization)
combination can be used for the key value.

Also, as this file may be modified automatically by the code_saturne
scripts and GUI, is is recommended to place matching commments before
section and key definitions, so thet may be written in the correct
place when the file is regenerated.

As an extension of the common [`.ini`] file format, an alternative way of
defining key-value pairs in a given section is to define a section named
*section-name:key*. In this case, all lines inside that section are
associated to the value (except for initial and final empty lines).
This avoids indentation requirements using multiline entries and generally
keeps things more readable. Such sections may optionally be closed by
an empty `[]` section declaration. This is only useful if comments
must be added before a following section, as they are implicitly closed
when a following section declaration or the end of file is reached.

The relevant sections and associated keywords are described below:

### [setup]
<!-- -->

Optional section relative to associated setup. Allowed keywords are:

* `parameters`

  Name of the parameters file (default: setup.xml).

* `coupled_domains`

  List of domains that should be coupled, separated by colons (_:_).
  When present (for a coupled run's top-level configuration file),
  a section named after each domain's (transformed to lowercase)
  may also be present to define additional options for that case.

As the recommended `setup.xml` is used by default if not specified here but
present in the directory structure, this section is optional except for
top-level coupled case `run.cfg` files, and useful only for compatibility
with older cases containing multiple setup files (which is not recommended).
An absolute path may also be used, but is usually not recommended as the case
structure is then not self-contained.

### [job_defaults]
<!-- -->

This section defines defaults when no associated
[${resource_name}](@ref case_structure_run_conf_resource_options)
or [${batch_type}] section is present. The same key-value
pairs may be used.

### [run]
<!-- -->

Optional (recommended) section relative to run stages and other aspects.

In case of multiple available builds (such as when standard and debug builds
are available), the compute build may be defined here:

* `compute_build`

   name or path of alternate compute build; if not set, the
   install settings or defaults apply; <br/>

The run *ids* (defining the results directory names in the `RESU` or
`RESU_COUPLING` directory) can be specified using the following keys!

* `id`

  id assigned to run in results directory;

* `id_prefix`

  prefix to assign to automatically-generated run id's
  (i.e. in absence of `id` value);

* `id_suffix`

  suffix to assign to automatically-generated run id's
  (i.e. in absence of `id` value);

* `force`

  if *true*, the computation is allowed when a result directory with
  the same id is already present and the stage step is required;
  by default, (*false*), an error is returned and the computation is not run,
  to avoid overwriting existing data; <br/>

The following keywords allow determining which stages which should be run.
By default, all steps are executed, unless some steps are specified, in which
case the specified steps and those in between are run; if a single step
is specified, all steps up to that one are assumed.

* `stage`

  *true* or *false* to indicate whether the staging step should be run; must
  be *true* unless a result directory with the same id has already  been staged.

* `initialize`

  *true* or *false* to indicate whether the preprocessing step should be run;

* `compute`

  *true* or *false* to indicate whether the compute step should be run;

* `finalize`

  *true* or *false* to indicate whether the finalization step should be run.

### [${resource_name}] {#case_structure_run_conf_resource_options}
<!-- -->

A section defining the requested options specific to compute environment
(and associated resource) can be defined by using the resource's name
which can be configured using the `resource_name` keyword in the
`[install]` section of the global install or user
[configuration](@ref cs_user_configuration_file) (with the system rather
than user setting being recommended).

In this documentation section, ${resource_name} should be replaced by
the actual active resource name. If no resource name is specified
in the main code_saturne (or user) install configuration, the name
of the configured batch system type (in lowercase) is used instead.
If this is not available either, `${resource_name}` finally defaults
to `job_defaults`.

Optional (recommended) section relative to job defaults. Allowed keywords are:

* `n_procs`

   Number of MPI processes for computation (default: 1).

   This option overrides the values determined automatically through the resource
   manager (batch system) when both are present. This may be useful in case of
   an incorrect automatic determination of the number of MPI ranks on some
   systems (or for advanced uses such as debugging a case on a number of MPI ranks
   which is not a multiple of the number of processes per node available using the
   batch system), but should otherwise only be defined in the absence of a
   resource manager.

* `n_threads`

  Number of OpenMP threads for computation (default: 1).

* `time_limit`

  Time limit for computation, in seconds; unlimited if < 0 (default).
  When running under a resource manager (batch system), the actual limit
  will usually be lower.

When a batch system is configured, associated batch settings may be given
using one of several keywords:

* `job_parameters`

  List of parameters which should be passed to the resource manager-specific
  command (for example `sbatch`, `llsubmit`, 'qsub`, ... depending on system);

* `job_header`

   Job jeader that should be inserted at the beginning of the generated and
   submitted `runcase` or `run_solver` scripts.

   Using the special [[${resource_name}:job_header]] section type instead
   is recommended, as it avoids indentation requirements.

* `job_header_file`

   Defines the path to a file that contains the job header to insert
   (otherwise as above). Either an  absolute or relative (to `run.cfg`)
   path may be used.

* `jobmanager`

    This option is reserved for future use with the SALOME platform's
    JOBMANAGER tool, but is not yet available.

If more than one of these options are defined, the priority, from highest to
lowest, is as follows: `job_parameters`, `job_header`, `job_header_file`.

### [${resource_name}:${key}] {#case_structure_run_conf_section_key}
<!-- -->

Sections of this type are used to define key values associated to the
*${resource_name}* section that may spread over multiple lines.
All lines (except empty initial and final lines) are used as the key value.
The main usage is to store batch job headers, with the following key:

* `job_header`

  The associated lines are inserted in the generated and submitted `runcase`
  file.

Additional resource:key combinations allow inserting additional snippets
in the generated scripts, and may be useful mostly to define or modify
additional environment variables. The associated key names are:

* `run_prologue`

  The associated entry is inserted before the active run steps are executed.

* `run_epilogue`

  The associated entry is inserted after the active run steps are executed.

* `compute_prologue`

  The associated entry is inserted in the generated `run_solver` script,
  before the main solver execution, and is restricted to the computation
  environment; it is thus usually preferred to `run_prologue` when both
  could be used.

* `compute_epilogue`

  The associated entry is inserted in the generated `run_solver` script,
  after the main solver execution, and is restricted to the computation
  environment; it is thus usually preferred to `run_epilogue` when both
  could be used.

* `debug_args`

   Allows running the solver through a debugger using a wrapper to
   allow simple handling of options and arguments (see
   [section in developer guide](@ref code_saturne_debug_launch_set)
   for more details).

* `tool_args`

   Allows running the solver through a tool whose command-line
   arguments can be passed before those of the solver executable.
   This is similar to `debug_args`, but with no transformation of the
   given arguments.

* `mpi_tool_args`

   Allows running the solver through a tool whose command-line
   arguments can be passed before those of the MPI launch command.
   This is similar to `tool_args`, but prefixes the MPI command
   when it is present, so should be used for tools designed to
   integrate with MPI, such as parallel debuggers.

### [${coupled_case_name}] {#case_structure_coupling_options}
<!-- -->

In case of code coupling, for each domain, a section whose name is
based on the domain name may be present.
The section name should always be in lowercase (per file formmat
specifications) even if the domain name is not.

* `solver`

  Defines the solver type; currently allowed names (case-independent) are:
  `code_saturne`, `neptune_cfd`, `SYRTHES`, `CATHARE`, `python_code`.
  Additional allowed or required keywords may depend on the solver
  type.

* `domain`

  Directory name (with exact capitalization) associated to the given
  domain. By default, this should be the same as the domain name.

* `n_procs_weight`

  How many MPI ranks will be assigned to this domain will be based
  on the ratio of this weight relative to the total `n_procs_weight`
  of all coupled domains and the total number of ranks assigned
  to the coupled computation.

  The weight to assign to each domain may be estimated based on the
  relative domain sizes and associated computational cost, so
  as to balance the load as well as possible. Checking performance log
  coupling timings may help improving the load balance based on previous
  runs (when the coupling communication time represents the largest part
  of the coupling exchange cost, this can usually be interpreted as including
  time waiting for other domains, so more resources should be allocated
  to domains with lower communication time, and less to those with higher
  communication time.

* `n_procs_min`

  Minimum number of MPI ranks assigned to this domain. By default,
  this value is 1. This setting may be useful if the weight-based computation
  could lead to an insufficient number of assigned ranks for some resource
  configurations, for example due to rounding.

* `n_procs_max`

  Maximmum number of MPI ranks assigned to this domain. This may
  be useful if the computational tool associated to a given domain
  is not parallel or is expected not to scale well beyond a given number
  of MPI ranks.

* `opt`

  For Syrthes domains, additional options (for example, postprocessing with
  `-v ens` or `-v med`).

* `param`

  For Syrthes domains, name of associated parameters file.

* `cathare_case_file`

  For CATHARE domains, name of the associated dataset file.

* `neptune_cfd_domain`

  For CATHARE domains, name of the computational domain assigned
  to the false neptune_cfd instance which actually wraps CATHARE.

* `script`

  For Python-based solver domains, name of the main matching
  Python script.

* `command_line`

  For Python-based solver domains, name of the associated
  command-line arguments.

### [paths] {#case_structure_run_conf_paths}
<!-- -->

Sections of this type are used to define paths in specific cases,
and should only appear in `run.cfg` files generated by a first
run stage.

* `case`

  Specifies the base case directory. When using a standard directory structure,
  this is not needed, but is used to determine the parent case when a
  separate results top directory has been specified (i.e. using the `--dest`
  run and submit option).

* `top_results_directory`

  Specifies the base case directory. This is useful to determine whether a
  separate results top directory has been specified (i.e. using the `--dest`
  run and submit option), or if a temporary execution directory is being used.

