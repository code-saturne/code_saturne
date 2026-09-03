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

\page ug_bpg_vnv_test_case Best Practice Guide for creating a V&V test case

[TOC]

# BPG - Creation of Verification and Validation Test Cases

## Study Architecture

A Verification or Validation (V&V) test case should follow the standard
architecture of a code_saturne study (see \subpage cs_ug_case_structure). In
addition, it must include an XML configuration file required by the
studymanager (SMGR) tool (see \subpage cs_ug_studymanager).

All stages of the workflow, including pre-processing, simulation execution,
post-processing, and report generation, should be managed through SMGR.

When executed, SMGR automatically creates a RUN_STUDY directory alongside the
study directory and handles the execution of all simulations, either on a local
workstation or on a computing cluster. It also performs the required
post-processing tasks, generates figures, and produces the final report
according to the information specified in the SMGR XML file.

Consequently, a REPORT directory containing the LaTeX source files of the report
is mandatory. The report should provide a comprehensive description of the test
case and an analysis of the obtained results. For validation cases, the results
are typically compared against experimental measurements, reference numerical
solutions and, when relevant, against results obtained with previous versions of
the code. For verification cases, the results are generally compared with
analytical solutions.

The following sections provide more detailed information about the required
files and directory structure.

### Studymanager xml file

The studymanager configuration file must be named `smgr.xml`. Additional XML
files may exist within the study directory, but only `smgr.xml` is considered
part of the V&V workflow.

This file defines the reference data (experimental measurements or reference
results), the list of simulation cases to be executed, the post-processing
scripts to be run, and the figures to be generated.

The inclusion of study metadata is mandatory (see example below).

```{.xml}
<study_keywords>
  Laminar, 2D, Unsteady, Incompressible
</study_keywords>
```

Simulation results displayed in figures must be clearly labeled with the version
of code_saturne used to generate them (e.g. v9.0).

### List of Cases

The total number of cases included in a study should be kept as small as
possible. Variants of a given test case can be generated directly by
studymanager without duplicating the entire case definition.

For example, notebook variables, numerical parameters, or mesh files can be
modified through the notebook, parametric, and kw_args nodes, respectively (see
example below).

```{.xml}
<study label='STUDY' status='on'>
  <case label='CASE1' status='on' compute="on" post="on">
    <notebook args="u_inlet_1=0.1 u_inlet_2=0.2"/>
    <parametric args="-m grid2.med --iter-dt 0.005"/>
    <kw_args args="--my-gradient=lsq --my-restart-100-iter"/>
  </case>
</study>
```

### Data settings

To minimize maintenance costs and improve code reusability, common user source
files should be grouped in the `GENERIC_SRC` directory, while shared input data
should be stored in the `GENERIC_DATA` directory.

```{text}
STUDY/
├── CASE1/
│   ├── DATA/
│   │   ├── code_saturne
│   │   ├── run.cfg
│   │   ├── setup.xml
│   │   ├── catalyst.py -> ../../GENERIC_DATA/catalyst.py
│   │   └── input.dat -> ../../GENERIC_DATA/input.dat
│   └── SRC/
│       ├── cs_user_cs_user_extra_operations.cpp
│       ├── cs_user_initialization.cpp -> ../../GENERIC_SRC/cs_user_initialization.cpp
│       └── cs_user_parameters.cpp -> ../../GENERIC_SRC/cs_user_parameters.cpp
├── CASE2/
│   ├── DATA/
│   │   ├── code_saturne
│   │   ├── run.cfg
│   │   ├── setup.xml
│   │   ├── catalyst.py -> ../../GENERIC_DATA/catalyst.py
│   │   └── input.dat -> ../../GENERIC_DATA/input.dat
│   └── SRC/
│       ├── cs_user_cs_user_extra_operations.cpp
│       ├── cs_user_initialization.cpp -> ../../GENERIC_SRC/cs_user_initialization.cpp
│       └── cs_user_parameters.cpp -> ../../GENERIC_SRC/cs_user_parameters.cpp
├── GENERIC_DATA/
│   ├── catalyst.py
│   └── input.dat
├── GENERIC_SRC/
│   ├── cs_user_initialization.cpp
│   └── cs_user_parameters.cpp
├── MESH/
├── POST/
└── REPORT/
```

The use of the Graphical User Interface (GUI) should be maximized whenever
possible to define simulation settings. This improves readability,
maintainability, and reproducibility of the test case.

All user source files must include the appropriate copyright headers and comply
with the code formatting and coding style guidelines (see coding style
guidelines).

Particular attention should be paid to the selection of numerical parameters.
Whenever possible, the default code_saturne settings should be used.

### POST Folder

The `POST` folder should contain all post-processing scripts, reference data,
and previously generated code_saturne results required for result analysis and
comparison.

```{text}
STUDY/
├── CASE1/
├── CASE2/
├── GENERIC_SRC/
├── MESH/
├── POST/
│   ├── EXP/
│   │   ├── experiment.dat
│   ├── v9.0
│   │   ├── CASE1/
│   │   │   └── run1/
│   │   │       └── results.dat
│   │   └── CASE2/
│   │       └── run1/
│   │           └── results.dat
│   └── post.py
└── REPORT/
```

### REPORT Folder

The `REPORT` folder must contain:
- a LaTeX report file describing the test case,
- a `Makefile` used to build the documentation,
- a `FIG` directory containing all permanent figures and supporting files
  referenced in the report (see example below).

```{text}
STUDY/
├── CASE1/
├── CASE2/
├── GENERIC_SRC/
├── MESH/
├── POST/
└── REPORT/
    ├── FIG/
    │   ├── geom.png
    │   └── mesh.png
    ├── report.tex
    └── Makefile
```

The report should provide a comprehensive and self-contained description of the
test case, including:
- the geometry and physical setup,
- mesh generation and discretization strategies,
- simulation settings and input data,
- validation and analysis of the results through appropriate figures,
  comparisons and discussions.

Any numerical parameter differing from the default configuration shall be
explicitly justified and documented to ensure traceability and reproducibility
of the simulation setup.
