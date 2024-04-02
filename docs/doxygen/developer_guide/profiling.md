<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

<!--

Colors
------

blueedf          rgb(0,91,187)
blueedf_light    rgb(0,91,187)
blueedf_dark     rgb(9,53,122)

greenedf         rgb(0,158,47)
greenedf_light   rgb(127,175,16)
greenbed_dark    rgb(48,119,16)

orangeedf        rgb(255,160,47)
orangeedf_light  rgb(255,122,0)
orangeedf_dark   rgb(254,88,21)

antiquewhite  rgb(0.98, 0.92, 0.84)
whitesmoke    rgb(0.96, 0.96, 0.96)

-->

\page cs_dg_profiling profiling

[TOC]

Introduction
============

For solvers such as code_saturne, which can run on large resources for
long durations, improving performance is always essential, to reduce
both user wait times and IT costs (of which a large part is nowadays
energy cost).

Performance gains are usually a combination of progress in computing power and in algorithms. For similar software over the last few decades, both factors have been important, with algorithm progress being a major driver.

So although the first step in improving program performance is on the algorithmic (and theory) side, detailed analysis of the behavior of those algorithms on actual hardware is important.

To assist developers and users in performance optimization, code_saturne includes many timers, and tries to log synthetic performance information.

This allows comparing the performance of numerical options and checking that no "unexpected" performance bottlenecks are present.

For more detailed analysis, the use of profiling tools is recommended.

Preliminary knowledge
---------------------

To be able to understand the performance behavior of code_saturne, the user should have at least introductory knowlege of several hardware and programming model related aspects.

- code_saturne is mostly memory bound, so understanding of the [roofline model](https://en.wikipedia.org/wiki/Roofline_model) is important.
  * This implies some understading of memory [models](https://en.wikipedia.org/wiki/Memory_model_(programming)) and [hierarchies](https://en.wikipedia.org/wiki/Memory_hierarchy).
  * Merging operators so as to increase arithmetic intensity (reuse of variables in memory for multiple computations) is essential to obtaining proper performance on modern hardware.

- Regarding parallelism, for both MPI (first level) and OpenMP (second level), understanding implicit and explicit barriers, and their relation to load balance is also important.

Built-in performance measurement
================================

Performance of each time step
-----------------------------

When running, the code_saturne solver generates a <span style="color:rgb(48,119,16)"><b>`timer_stats.csv`</b></span> file, which traces the elapsed time for each major operation type (mesh modification, post-processing, gradient reconstructions, linear solvers, and such). This information may be easily plotted using a spreadsheet or a visualization tool such as ParaView.

The code also generates a <span style="color:rgb(48,119,16)"><b>`performance.log`</b></span> file, which summarizes timings for various operations, in a manner independent of the number of time steps actually run (so this file is complete only after a successful run).

Profiling tools
===============

To obtain more detailed performance information, use of a **profiling** tool is needed.

Use <span style="color:rgb(48,119,16)">`--enable-profile`</span>
to configure builds for profiling.

Common profiling tools
----------------------

Several types of tools may be available. We list a few commonly available tools, though the list is far from exhaustive:

- Various profilers are commonly found on Linux machines.
  - *gprof* is obsolete, and does not work well with non-static builds. Please forget about it outside of purely historical interest.
  - The *gprofng* tool is a new profiler, available as part of recent Linux distributions base *binutils* package.
    * As this tool is not broadly available on older systems, we do not describe it here yet, but will recommend it and provide links to additional resources as its availability progresses.

- The Valgrind tool suite includes several tool which are very useful for profiling. Note that as usual when running under Valgrind, there is an overhead relative to actual performance, and the obtained timing results may be simulated as much as measured, but the information obtained is very similar to that obtained with less ubiquitous tools.

  - [Massif](https://valgrind.org/info/tools.html#massif) is a heap profiler,
    allowing measuring the evolution of memory over time.
    Using this tool simply requires defining:
```
    valgrind --tool=massif
```
    as a tool in the advanced run options in the GUI, or running
```
    code_saturne run --tool-args="valgrind --tool=massif"
```
    on the command line.
    The files produced are text files, but can also be visualized with the [Massif-vizualizer](https://apps.kde.org/fr/massif-visualizer/) application where available.

  - [Callgrind](https://valgrind.org/docs/manual/cl-manual.html) functions as a "general purpose" profiler, and can be used in a similar manner, using
```
    valgrind --tool=callgrind
```

 Combined with the [kcachegrind](https://kcachegrind.github.io/html/Home.html) visualization tool, it is extremely easy to use on a Linux workstation. It allows easy visualization of call trees and hot spots, as illustrated below:

  \image html dg/kcachegrind.png "GDB in terminal mode" width=80%

Vendor profiling tools
------------------------

Other advanced profiling tools may be provided by various vendors, for example:
- [VTune](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler-documentation.html) for Intel systems.
- [NVIDIA Nsight Systems](https://developer.nvidia.com/nsight-systems) for NVIDIA GPUs.

### Using VTune with code_saturne

If Intel's <span style="color: rgb(48,119,16)">VTune</span> is available, the following procedure may be used:

- First, initialize the execution directory, using one of the following methods:
  - From the GUI, in the advanced run/job submission options, check "initialize only", then submit the computation.
  - Outside the GUI, run `code_saturne submit --initialize`
  In either case, the code will prepare the execution directory,
  and preprocess the mesh if needed, but not remove the executable and temporary script.

- Once the stage has finished, `cd` to the execution directory, and edit the <span style="color: rgb(48,119,16)">`run_solver`</span> script script:
  - Add commands necessary to load the profiler's environment.
    For example, on the EDF Cronos cluster, this requires adding
```{.sh}
module load intel-Basekit
```
    in the section where other modules are added.
  - Search for the line actually launching the solver, near the end of the script;
    * Before that line, add:
```{.sh}
# Select a different profiling directory at each run
export DIR_PROF=profiling.$SLURM_NTASKS.$SLURM_JOB_ID

# Avoid using /tmp, which may be too small
export TMPDIR=$SCRATCH
```
    - Before the `cs_solver` command, insert the profiling commands.
      For example, replace
```{.sh}
mpiexec <options> ./cs_solver --mpi "$@"
```
      with:
```{.sh}
mpiexec <options> -q -collect hotspots -r $DIR_PROF -data-limit=4000 -- ./cs_solver --mpi "$@"
```
      (`mpiexec` might be replaced by another command, such as `srun` depending on the system, but the logic remains the same).

- If using a batch system (usually true on a cluster), you will need to submit the `run_solver` script rather than running it directly.
  * If you are familiar with batch commands, pass the required commands to the submission command.
  * Otherwise, you can copy the batch headers at the beginning of the `runcase` file to the same position (starting at line 2) in the `run_solver` file.
    - With the SLURM resource manager, this means:
      * Copy all lines of `runcase` starting with `#SLURM` to `run_solver`.
      * Run
```{.sh}
sbatch run_solver
```
    - With other systems, the syntax will be slightly different but the principle remains the same.

- Once the code has finished running, simply run
```{.sh}
vtune-gui <profiling_directory_name>
```
This may require loading an environment module (as in the `run_solver` file for the VTune installation (`intel-Basekit` in the previous example).

VTune allows many exploration views, for example:

\image html dg/vtune_screen.png "VTune summary page" width=80%

or

\image html dg/vtune_screen_hotspots.png "VTune hotspots view" width=80%
