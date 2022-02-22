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

\page cs_ug_run_computation Running a calculation

[TOC]

Step by step preparation {#sec_prg_stepbystepcalculation}
========================

This paragraph summarizes the different steps which are necessary to
prepare and run a standard case:

* Check the version of code_saturne set for use in the environment variables
  (`code_saturne info --version`). If it does not correspond to
  the desired version, update the user environment or aliases to use the
  required version, logging out of the session and in again if necessary
  (cf. [setting up you environment](@ref sec_prg_environement_cs)).
  <br/><br/>

* Prepare the different directories using the GUI or
  [`code_saturne create`](@ref sec_prg_cscreate) command .
  <br/><br/>

* It is recommended to place the mesh(es) in the `MESH` directory,
  but they may be selected from other directories, either with the Graphical User
  Interface (GUI) or the `cs_user_scripts.py` file (see below).
  Make sure they are in a [supported format](@ref sec_prg_meshes). There can be
  several meshes in case of mesh joining, code coupling, or simply parametric
  studies.
  <br/><br/>

* Switch to the `DATA` directory and start the GUI using the following command:
  `./code_saturne gui`
  <br/><br/>

* If not using the GUI or when requiring additional advanced settings, copy the
  `DATA/REFERENCE/cs_user_scripts.py` file to `DATA` and edit it, so that the
  correct run options and paths may be set.
  Just as with user-defined functions described later, settings defined in this
  file have priority over those defined through the GUI;
  <br/><br/>

* Before running a full simulation, it is recommended to check the mesh using
  one of the mesh import or preprocessing [execution modes](@ref sec_prg_executionmodes);
  <br/><br/>

* Define the computation settings, preferably using the GUI.
  <br/><br/>

* If needed, place in the `DATA` directory the different external data
  files which may be necessary;
  <br/><br/>

* Place the necessary user source files in the `SRC` directory
  - those files can be copied from the references in `SRC/REFERENCE`, and are
    documented in the *files* section of this documentation;
  - it is recommended that only settings which cannot be defined through the GUI
    be defined through these files; Note that GUI settings and user-defined
    functions can be combined, with user-defined settings having priority.
  - code snippets can often be cherry-picked and adapted from the many
    available [user examples](@ref cs_user_examples);
    avoid copying unneeded sections from examples, as they make the code harder
    to read and maintain;
  - all files with a .c, .f90, or. C++ extension in `SRC` will be compiled
    and linked with the code at execution;
  - avoid placing non-user files (modified or not) from the main code_saturne source
    tree in `SRC`, unless specifically provided by the core development team as a
    temporary patch for an issue, as this can cause more or less subtle issues,
    and removes any semblance of quality control.
  <br/><br/>

* Execution parameters for the current system may be defined either
  using the GUI or editing the relevant [`run.cfg`](@ref sec_prg_run_cfg)
  file sections;
  <br/><br/>

* Run the calculation and analyze the results
  <br/><br/>

* Remember to occasionally, purge temporary or unneeded files (in `RESU/<run\_id>`)
  (or in the `<scratch>/<run_id>`
  [temporary execution directory](@ref case_structure_scratchdir) if set).

Execution modes {#sec_prg_executionmodes}
===============

As explained before, code_saturne is composed of several main modules, including
the *Preprocessor* and the *Solver*. The Preprocessor reads the meshes.
The resulting data is transferred to the Solver through specific
files, usually named `mesh_input.csm` or having the `.csm`
extension, or placed in a directory named `mesh_input` when multiple meshes
are imported. Files simply named `mesh_input` or `mesh_output` from
older versions may also be used.

Since the Preprocessor does not run in parallel and may require a
large amount of memory, the various execution modes allow minimizing
its use or running it separately. This also allows preprocessing meshes
using code_saturne builds with support for specific file formats, and
using them on another machine with a minimal build not supporting those
same formats.

Currently, the code_saturne GUI allows for 4 different execution modes
under the ``Mesh`` page:

* __mesh import__: the Preprocessor is run to transform one or more meshes into
  an internal `mesh_input.csm` file (or `mesh_input` directory in case of
  multiple meshes). The Solver is not run;

* __mesh preprocessing__: the Solver is run in preprocessing mode, so as to
  handle all mesh modification operations, such as joining, periodicity,
  smoothing, *etc.* If a `mesh_input.csm` file or `mesh_input` directory is
  provided, it is used directly; otherwise, the Preprocessor is run first
  to import meshes;

* __mesh quality criteria__: similar to preprocessing, with the addition of
  mesh quality criteria computation, and post-processing output of those
  criteria. Some additional mesh consistency checks are also run;

* __standard__: this includes preprocessing, followed by a standard computation.

At a lower level, the launch scripts allow specifically
choosing which modules to run, either through the GUI or through the
`cs_user_scripts.py` file:

* If an input file or directory is defined (which may be either a `mesh_input`
  or `mesh_input.csm` from a previous Preprocessor run or a `mesh_output.csm`
  from a previous Solver run, the script will copy or link it to
  the execution directory, and the Preprocessor will not be run again;

* If `domain.exec_kernel = False` in `cs_user_scripts.py`, the Solver will not
  be run. This is useful when only the mesh import stage is required,
  so as to emulate the GUI's *mesh import* only mode;

* When running a follow-up computation (restarted from a previous computation),
  preprocessing in general is avoided, and the `mesh_input.csm` file
  from the previous checkpoint directory is used.
  This behavior can be modified in the GUI (under "Mesh/Execution mode",
  un-checking "Use unmodified checkpoint mesh in case of restart",
  or setting `domain.preprocess_on_restart = True` in `cs_user_scripts.py`).
  This may be useful in 2 cases:
  - when the pre-processed mesh has not been saved from a previous run
    (usually so as to save disk space)
  - when the computation restarted from used a different mesh than the
    current one (to which data will be interpolated).

In a similar manner, the Solver accepts several command-line options relative to
execution mode, notably `domain.solver_args = '--preprocess` or `--quality`,
restricting the run to the preprocessing stages, or preprocessing stages augmented
by mesh quality criteria computation. Whenever the preprocessing stages defined
lead to an effective mesh modification, a `mesh_output.csm` file is produced, which
can be used directly as an input for a successive calculation, unless
deactivated using the GUI (unckecking "Save mesh if modified by preprocessing")
or using the [cs_user_mesh_save](@ref cs_user_mesh_save) user-defined function.

To allow preprocessing in multiple passes, all defined preprocessing
operations may be allowed on previously preprocessed meshes, as described
above.

It is encouraged to separate the preprocessing and calculation runs, as
this not only speeds up calculations, but also ensures that the mesh is
identical, regardless of the architecture or number of processors it is run
on. Indeed, when running the same preprocessing stages such as mesh joining
on a different machine or a different number of processors, very minor
floating-point truncation errors may lead to very slightly different
preprocessed meshes. The GUI option to "Use unmodified checkpoint mesh
in case of restart" encourages this usage.

Note also that mesh partitioning is done directly by the Solver.
Depending on the partitioning algorithm used, a partition map
(`partition_output/domain_number_*`) may be output,
allowing the use of the same partitioning in future calculations.
By default, this file is output when using graph-based partitioners, which may
use randomization and do not guarantee a reproducible output, and is not output
when using a deterministic space-filling curve based partitioning.

If the code was built only with a serial partitioning library,
graph-based partitioning may best be run in a serial preprocessing stage.
In some cases, serial partitioning might also provide better partitioning
quality than parallel partitioning, so if both are available, comparing
the performance of the code may be worthwhile, at least for calculations
expected to run for many iterations.

Execution steps {#sec_prg_exec_stages}
---------------

Each calculation (whether run from the GUI or using `code_saturne run`)
is split into the following execution steps:

* __initialize__

  * __stage__
    - Create `RESU/<run-id>` directory
    - Copy data from `DATA` to `RESU/<run-id>`
    - Copy user sources from `SRC` to `RESU/<run-id>/src`, and compile them.

  * __initialize__
    - Import meshes using Preprocessor (serial) tool to `mesh_input*`
    - For Syrthes couplings, partition Syrthes mesh.
    - Generate `run_solver` mini-script for execution stage.

* __compute__
  - Execute `run_solver` mini-script:
    sets module or environment variables, runs MPI or serial command;

* __finalize__
  - Cleanup only when previous stages successful: remove executable and
    possibly other user-specified files.

By default, all steps are executed. If some stages are specified, by adding
`--stage`, `--initialize`, `--compute`, and/or `--finalize` to the
`code_saturne run` options, only steps between the first and last one
specified are executed.

### Job submission on cluster {#sec_prg_exec_stages_hpc}

The *initialize* step itself can be split into two sub-steps, so that when
submitting a run to a batch system (using the GUI or the `code_saturne submit`
command), *stage* is run interactively, and the following steps are executed:

* run the __stage__ step immediately;
* switch to `RESU/<run-id>` directory;
* generate `runcase` script based on `run.cfg` parameters, containing
  batch directives and a `code_sature run --id <run_id>` command removing the
  (already done) `--stage` step;
* submit the `runcase` script

This allows errors in data paths (such as incorrect restart paths) or user-defined
sources to be detected immediately, before the job is submitted, and
also modifying the base setup without impacting an already submitted job.
On some HPC systems, the compilers may also be available only on the front-end nodes,
so this also avoids possible issues related to trying to compile user-defined
functions and re-link the code from a compute node.

Computation restart {#sec_prg_restart}
===================

As some computations may require long run times it is possible to
run a computation in multiple steps using a checkpoint/restart mechanism.

By default, a running computation produces a `checkpoint` directory
containing several files:

- `mesh.csm` contains the preprocessed mesh
- `main.csc` contains the values of the main fields.
- `auxiliary.csc` contains additional field values, such as mass fluxes,
   which are needed for a perfectly continuous restart, but not
   absolutely necessary.
- possible additional files based on the active physical models
  (`radiative_transfer`, `lagrangian`, `lagrangian stats`, ...).

This directory can be used as `restart` directory for a subsequent
computation (this can be easily defined in the GUI or in the
additional `cs_user_scripts.py` user scripts).

If the code cannot find one or several elements of data required for the
calculation restart, default values are then used.
If the number of faces has been modified (for instance in case of
modification of the mesh merging or of periodicity), reading
the auxiliary restart file should be deactivated (in the GUI
or setting \ref cs_glob_restart_auxiliary->read_auxiliary to `false`).

Checkpoint files are binary, but their contents can be queried
and dumpled using the `code_saturne bdump` command,
and compared using `code_saturne bdiff`.

When running in parallel, data is read and written in in partitioning-independent
manner, based on element global numbers, so restarting can be transparently done
on a different number of processes, or using a different partitioning.
When the mesh is assembled from multiple files, the global element
numbers are assigned in sequence, so the order of assembly should not be
modified between standard restarts.

It is possible to restart a computation that was run using a different
mesh. In that case, the original `mesh.csc` file must also be provided
along with the restart directory (in most cases, it should be already be
present in the matching `checkpoint` directory).

Checkpointing of the mesh and various files may usually be deactivated
globally, or using parameters specific to each file.

Note also that when unchanged between succeeding computations, the `mesh.csm`
file is linked (in the Unix/Linux/POSIX) sense rather than copied: when
as long as checkpoint directories are present on the sale file system,
the files are shared rather than actually copied, to save space.

Interactive modification of selected parameters {#sec_prg_control_file}
===============================================

During a calculation, it is possible to interactively modify the time step or time value
limit defined in the setup. To do so, a file named `control_file` must be placed in the
execution directory.
The existence of this file is checked at the beginning of each time step.

To change the maximum number of time steps, this file must contain a line
indicating the value of the new limit number of time steps. If this new limit
has already been reached, code_saturne will stop properly at the end of the current
time step (the results and checkpoint/restart files will be written correctly).
This procedure allows the user to stop a calculation in a clean and interactive
way whenever needed.

The `control_file` may also contain a few other commands, allowing the user to
force checkpointing or post-processing at a given time step or physical time, or
to force an update of log files.
The following commands are available (using the common notations "`< >`" to
indicate a required argument, "`[ ]`" to indicate an optional argument).

<table>
<caption id="control_file_commands">control_file syntax</caption>
<tr><th> command                          <th> arguments
<tr><td> max_time_step                    <td> <time_step_number>
<tr><td> max_time_value                   <td> <time_value>
<tr><td> max_wall_time                    <td> <wall_time>
<tr><td>                                  <td>
<tr><td> checkpoint_time_step             <td> <time_step_number>
<tr><td> checkpoint_time_value            <td> <time_value>
<tr><td> checkpoint_wall_time             <td> <wall_clock_time>
<tr><td>                                  <td>
<tr><td> checkpoint_time_step_interval    <td> <time_step_interval>
<tr><td> checkpoint_time_value_interval   <td> <time_interval>
<tr><td> checkpoint_wall_time_interval    <td> <wall_time_interval>
<tr><td>                                  <td>
<tr><td> control_file_wtime_interval      <td> <wall_time_interval>
<tr><td>                                  <td>
<tr><td> flush                            <td> [time_step_number]
<tr><td>                                  <td>
<tr><td> notebook_set                     <td> <parameter_name> <value>
<tr><td>                                  <td>
<tr><td> postprocess_time_step            <td> <time_step_number> [writer_id]
<tr><td> postprocess_time_value           <td> <time_step_value> [writer_id]
<tr><td>                                  <td>
<tr><td> time_step_limit                  <td> <time_step_count>
</table>

The `time_step_limit` differs from the `max_time_step` command,
in the sense that it allows reducing the maximum number of time steps,
but not increasing it. Also, in the case of a restart, it refers to the
number of additional time steps, not to the number of absolute time steps.
It is used mainly by the Studymanager component.

Note that for the `postprocess_time_*` options, the last argument
(`writer_id`) is optional. If not defined, or 0, postprocessing
is activated for all writers; if specified, only the writer with the specified
id is affected. Also, postprocessing output by one ore more writers at a
future time step may be cancelled using the negative value of that time step.

For the `flush` option, the time step is also optional. If not
specified, logs and time plots are updated at the beginning of the
next time step. Also, if the `control_file` is empty (such as
when created by the `touch control_file` command on Unix/Linux
systems, a `flush` request for the next time step.

Multiple entries may be defined in this file, with one line per entry.

Environment variables {#sec_env_var}
=====================

Setting some environment variables allows modifying code_saturne's
default behavior. The environment variables relevant to and specific to code_saturne
are described here:

General environment variables relevant to code_saturne {#sec_env_var_std}
------------------------------------------------------

Variable             | Role
---------------------|------------------------------------------------------------
`OMP_NUM_THREADS`    | Set the number of OpenMP threads if supported (`OMP_NUM_THREADS=1` for single-threaded behavior). This variable is normally set by the run scripts and GUI when OpenMP support is available.
`LD_LIBRARY_PATH`    | Search path for locate shared libraries. Needed in many cases to locate libraries outside the default system paths
`LD_PRELOAD`         | List of dynamic libraries to load before all others. This allows experimenting with other versions of these libraries without recompiling the code, or to load debugging libraries. This may be risky in case of incompatible library versions.
`PATH`               | Search path for executable

Note that the path management variables may be modified by the code_saturne scripts, but also
by other scripts, wrappers, and environment modules when present (which is the case
on many HPC systems).

To determine which shared libraries are used by an executable file (and whether they are found using the currently loaded environment), use
the following command: `ldd <executable_path>`.

Environment variables specific to the code_saturne environment {#sec_env_var_cs}
--------------------------------------------------------------

Variable               | Role
-----------------------|------------------------------------------------------------
`CS_SCRATCHDIR`        | Allows defining the execution directory (see [temporary directory](@ref case_structure_scratchdir)), overriding the default path or settings from the global or user `code_saturne.cfg`.
`CS_MEM_LOG`           | Allows defining a file name in which memory management based on the [BFT_MALLOC](@ref BFT_MALLOC), [BFT_REALLOC](@ref BFT_REALLOC), and [BFT_FREE](@ref BFT_FREE) is logged (useful to check for some memory leaks).
`CS_MPIEXEC_OPTIONS`   | This variable allows defining extra arguments to be passed to the MPI execution command by the run scripts.  If this option is defined, it will have priority over the value defined in the preferences file (or by computed defaults), so if necessary, it is possible to define a setting specific to a given run using this mechanism.  This may be useful when tuning the installation to a given system, for example experimenting MPI mapping and "bind to core" type features.
`CS_RENUMBER`          | Deactivating mesh renumbering in the Solver is possible by setting `CS_RENUMBER=off`.
`CATALYST_ROOT_DIR`    | Indicate where the ParaView Catalyst libraries are installed; the associated library path is added to `LD_LIBRARY_PATH` by the low-level Solver launch script, but does not otherwise interfere with the user's normal environment
`CATALYST_LD_ADD_PATH` | Indicate where additional libraries needed by ParaView Catalyst may be found; this path is added to `LD_LIBRARY_PATH` by the low-level Solver launch script, but does not otherwise interfere with the user's normal environment
`CATALYST_PLUGIN_PATH` | Indicate where additional libraries needed by ParaView Catalyst may be found; the low-level Solver launch script sets `PV_PLUGIN_PATH` to this value, which does not otherwise interfere with the user's normal environment

Environment variables specific to the Preprocessor {#sec_env_var_pcs}
--------------------------------------------------

Note that in general, if a given behavior is modifiable through an environment
variable rather than by a command-line option, it has little interest for a
non-developer, or is expected to be needed only in exceptional cases.

Variable                                 | Role
-----------------------------------------|----------------------------------
`CS_PREPROCESS_MEM_LOG`                  | Allows defining a file name in which memory allocation, reallocation, and freeing is logged.
`CS_PREPROCESS_MIN_EDGE_LEN`             | Under the indicated length ( *10<sup>-15</sup>* by default), an edge is considered to be degenerate and its vertices will be merged after the transformation to descending connectivity. Degenerate edges and faces will thus be removed. Hence, the post-processed element does not change, but the Solver may handle a prism where the preprocessor input contained a hexahedron with two identical vertex couples (and thus a face of zero surface). If the Preprocessor does not print any information relative to this type of correction, it means that it has not been necessary. To completely deactivate this automatic correction, a negative value may be assigned to this environment variable.
`CS_PREPROCESS_MEM_IGNORE_IDEAS_COO_SYS` | If this variable is defined and is a strictly positive integer, coordinate systems in [I-deas universal](@ref sec_fmtdesc_unv) format files will be ignored. The behavior of the Preprocessor will thus be the same as that of versions 1.0 and 1.1. Note that in any case, non Cartesian coordinate systems are not handled yet.


Running the solver directly {#sec_ug_cs_solver_direct}
===========================

In the standard cases, the compilation of code_saturne main program (`cs_solver`)
and its execution are entirely controlled by the launch script (using the `code_saturne run`
or `code_saturne submit` command). The command line options for (`cs_solver`) are passed
through user modifiable variables at the beginning of the `cs_user_scripts.py` file,
if present (if needed, it may be copied from a case's `DATA/REFERENCE` subdirectory to its
`DATA` directory; it may also be copied from the `${install_prefix}/share/code_saturne`
directory). This way, the user only has to fill these variables and does not need
to search deeply in the script for the solver command line. For more advanced
usage, the main options are described below.

Solver executable program {#sec_ug_opt_cs_solver_exe}
-------------------------

The executable of the main solver is named `cs_solver`. This program may usually
be found under `${install_prefix}/libexec/code_saturne` on most Linux systems
(or `${install_prefix}/lib/` on some systems such as OpenSUSE).

Whenever user-defined functions are used, a modified version of this program
is generated, including the user-defined functions. This is usually handled
by the run script, but can be done direclty using the `code_saturne compile`
command, and generates a `cs_solver` executable program directly in the
current directory.

The list of Command-line options for the `code_saturne compile` command
may be obtained by calling `code_saturne compile --help`. The most useful
options are `-t` or `--test` (test compilation only, without generating the program)
and `-s <source_path>` or `--source <source_path>`.

The `code_saturne compile` command always compiles all C, Fortran, and C++ files
from the specified directory (or current directory if not specified)
and links them with the installed library so as to produce a new `cs_solver`
with the compiled files. Unless forced, when no such files are present, no program
is generated, unless the `-f` or `--force` option is used.

Additional options allow adding additional compile or link flags.

Solver command line options {#sec_ug_opt_cs_solver}
---------------------------

It may be practical (especially when debugging) to run the `cs_solver`
directly. A complete and up-to-date list of command-line options
may be obtained by running `./cs_solver --help` (or
`${install_prefix}/libexec/code_saturne --help`).

The various command-line options are detailed here:

* `--app-name`

  Specifies the application name. This is useful only in the case of code coupling,
  where the application name is used to distinguish between different code instances
  launched together.

* `-wdir`, `--wdir`

  Specifies the working directory that should be used (used for code coupling, where
  each computational domain is assigned a specific sub-directory).

* `--mpi`

  Specifies that the calculation is running with MPI communications. The number of processes
  used will be determined automatically by the solver. With most MPI implementations, the
  code will detect the presence of an MPI environment automatically, and this option
  is redundant. It is only kept for the rare case in which the MPI environment might
  not be detected.

* `--sig-defaults`

  Use default runtime behavior when signals are received rather than code_saturne's
  own signal handlers.

* `--preprocess`

  Triggers the preprocessing-only mode.
  The code may run without any additional setup (parameter file or user function),
  though the user-defined setup may modify its default behavior.
  Only the initial operations such as mesh joining and modification are
  executed.

* `--quality`

  Triggers the mesh verification mode.

  This mode is similar to the preprocessing-only mode, and includes the preprocessing
  stages, adding elementary tests:
  - Quality criteria of the mesh are calculated (non-orthogonality angles,
    internal faces off-set, ... and corresponding visualizable post-processing output
    is generated.
  - A few additional mesh consistency tests are run.

* `--benchmark`

  Triggers the benchmark mode, for a timing of elementary operations on the machine.

  A secondary `--mpitrace` can be added. It is to be activated when the benchmark mode
  is used in association with an MPI trace utility. It restricts the elementary
  operations to those implying MPI communications and does only one of each
  elementary operation, to avoid overfilling the MPI trace report.
  This command can also be placed in the `domain.solver_args` variable
  in the `cs_user_scripts.py` when using this mode with the high-level run script is
  desired.

* `--logp`

  Activates the output for the processors of rank 1 to *n-1* in a calculation in
  parallel on *N* processors to files `run_solver_rxxxx.log` where `xxxx` is the
  rank id.

* `--trace`

  Activates the tracing of the output to the standard output.

* `--system-info`

  Print system information and exit.

* `--version`

  Print version number and exit.

* `-h`, `--help`

  Displays a summary of the different command line options.
