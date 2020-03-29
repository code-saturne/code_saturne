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

\page cs_run_computation Running a calculation

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
  unchecking "Use unmodified checkpoint mesh in case of restart",
  or setting `domain.preprocess_on_restart = True` in `cs_user_scripts.py`).
  This may be useful in 2 cases:
  - when the preprocessed mesh has not been saved from a previous run
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
on. Indeed, when running the same pre-processing stages such as mesh joining
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
graph-based partitioning may best be run in a serial pre-processing stage.
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

* __execute__
  - Execute `run_solver` mini-script:
    sets module or environement variables, runs MPI or serial command;

* __finalize__
  - Cleanup only when previous stages successful: remove executable and
    possibly otther user-specified files.

By default, all steps are executed. If some stages are specified, by adding
`--stage`, `--initialize`, `--execute`, and/or `--finalize` to the
`code_saturne run` options, only steps between the first and last one
specified are executed.

### Job submission on cluster {#sec_prg_exec_stages_hpc}

The *initialize* step itself can be split into two substeps, so that when
submitting a run to a batch system (using the GUI or the `code_saturne submit`
command), *stage* is run interactively, and the following steps are executed:

* run the __stage__ step immediately;
* switch to `RESU/<run-id>` directory;
* generate `runcase` script based on `run.cfg` parameters, containing
  batch directives and a `code_sature run --id <run_id>` command removing the
  (already done) `--stage` step;
* submit the `runcase` script

This allows errors in data paths (such as incorrect restart paths) or user-defined
sources to be detected immediately, before the job is sumbitted, and
also modifying the base setup without impacting an already sumbitted job.
On some HPC systems, the compilers may also be available only on the front-end nodes,
so this also avoids possible issues related to trying to compile user-defined
functions and re-link the code from a a compute node.

Environment variables {#sec_envcs}
---------------------

Setting a few environment variables specific to code_saturne allows modifying
its default behavior. The environment variables used by code_saturne
are described here:

Variable             | Role
---------------------|------------------------------------------------------------
`CS_SCRATCHDIR`      | Allows defining the execution directory (see [temporary directory] (@ref sec_prg_temporarydirectory)),overriding the default path or settings from the global or user `code_saturne.cfg`.
`CS_MPIEXEC_OPTIONS` | This variable allows defining extra arguments to be passed to the MPI execution command by the run scripts.  If this option is defined, it will have priority over the value defined in the preferences file (or by computed defaults), so if necessary, it is possible to define a setting specific to a given run using this mechanism.  This may be useful when tuning the installation to a given system, for example experimenting MPI mapping and "bind to core" type features.
