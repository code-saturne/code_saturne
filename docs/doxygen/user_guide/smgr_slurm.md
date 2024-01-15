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

\page cs_ug_smgr_slurm How to submit cases on cluster with SLURM

[TOC]

SLURM command lines
===================

This page explains how to submit cases on a cluster using the SLURM resource
manager with the studymanager tool.

In order to activate the submission of cases on cluster with SLURM, it is
necessary to specify at least one of the two following options:
- `--slurm-batch-size=N`: maximum number of cases per batch in SLURM batch mode
  (50 by default).
- `--slurm-batch-wtime=M`: maximum computation time in hours per batch in SLURM
  batch mode (12 hours by default).

For instance, the following command
  ```
  $ code_saturne smgr -f smgr.xml -r --with-tags=coarse --slurm-batch-size=20 --slurm-batch-wtime=8
  ```
will submit batch of cases with a maximum of 20 cases per batch and a maximum
total computation time of 8 hours. The number of cases per batch could then be
inferior to 20 if the total computation time exceeds 8h.

In order to compute the total computation time per batch, it is necessary to
specify an expected computation time per case (HH:MM) in the smgr xml file.

```{.xml}
    <study label="MyStudy1" status="on">
        <case label="Grid1" run_id="Grid1" status="on" compute="on" post="on" expected_time="2:15"/>
        <case label="Grid2" run_id="Grid2" status="on" compute="on" post="on" expected_time="5:00"/>
    </study>
```
In this case, the run _Grid1_ has an expected computation time of 2 hours and 15
minutes and the one of _Grid2_ is 5 hours.

SLURM batch options
-------------------

SLURM batch files are automatically generated in the folder _slurm_files_ in
**destination**. The following file is an example of a SLURM file used with the
SLURM batch mode.

```{.sh}
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=0:12:00
#SBATCH --output=vnv_6
#SBATCH --error=vnv_6
#SBATCH --job-name=saturne_vnv_6
```

- ntasks : number of cores required for the computation.
  * All cases are automatically sorted by number of required cores so that the
    number of tasks per batch is the same.
- time : sum of the computation times of the cases in the batch

Batch cases which require 5 or more cores will be executed in exclusive mode
(i.e. no other submission will run on the node).

Additional SLURM batch parameters can be also specified at run time using the 
`--slurm-batch-arg` option. This option only takes into account one argument at
a time. For example, to add the "exclusive" and send an e-mail notification use
the following command-line option:
`--slurm-batch-arg=--exclusive --slurm-batch-arg=--mail-user=name.last@email.com`

\warning
For EDF users, the `wckey` argument should be defined. It can be done by either
using `--slurm-batch-arg=--wckey=<key>` during run time, or by setting an
environment variable with the following command: `export SBATCH_WCKEY=<key>`.

All ouput and error files are also in the folder _slurm_files_ in destination.

Dependency between cases
------------------------

Job-dependencies are defined automatically such that blocks of dependency level
`M` will wait until all blocks of level `M-1` are successfully finished.

Three methods are available to define a dependency between cases:
- Set a restart in the data settings of a case using the graphical user
  interface of code_saturne;
- Add a restart in parametric options of a SMGR parameter file using
  `-r` or `--restart` argument
```{.xml}
  <study label='STUDY' status='on'>
      <case label='CASE1' status='on' compute="on" post="on">
          <parametric args="-r run_id"/>
      </case>
  </study>
```
- Add a `<depends>` node to a case in SMGR parameter file :
```{.xml}
  <study label='STUDY' status='on'>
      <case label='CASE1' status='on' compute="on" post="on">
          <depends args="STUDY/CASE/RESU/run_id"/>
      </case>
  </study>
```

Dependencies defined using a `<depends>` node have priority over those deduced
from parametric arguments. They both have priority over restarts in code_saturne
data settings.

In the rare cases where dependencies are not related to restarts, the `<depends>`
approach allow fine-grained control. In other cases, dependencies are deduced from
restart definitions, so no additional user settings are needed.

Postprocessing and analysis
---------------------------

The state analysis is automatically added with the slurm batch mode. This final
batch will depend on all previous submissions. It can also include
postprocessing and comparison steps if these options are activated.

Example
-------

In the following example, a list of 8 cases is launched in SLURM batch mode:

```
$ code_saturne smgr -f smgr.xml -rp --slurm-batch-size=2 --slurm-batch-wtime=5
```

\image html smgr_dependency.png "Exemple of cases allocation per batch"

Here are some explanations on cases allocation per batch :
- In level 0, batch 1 is limited by the number of cases with the same number of
  cores. Batch 2 is limited to one case (CASE1/run3) as the maximum time in this
  batch would have exceeded 5 hours including the next case (CASE2/run1). Batch
  3 only includes the last case at level 0. 
- In level 1, batch 4 is limited by the maximum number of cases per batch (2).
  Batch 5 only includes the last case at level 1. Batches 4 and 5 depends on all
  batches from level 0.
- In level 2, batch 6 only includes the last case of the level. It depends on all
  batches from level 0 and 1.
- In level 3, batch 7 includes postprocessing and state analysis. The final
  batch depends on all cases from all steps.
