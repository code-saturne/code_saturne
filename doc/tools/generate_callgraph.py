#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

# -------------------------------------------------------------------------------

"""
code_saturne Call Graph Generator
=================================

Generates an interactive, standalone HTML call graph and Graphviz DOT/SVG/Mermaid
architectural representations of the main code_saturne solver execution stages:
  1. Setup
  2. Mesh Preprocessing
  3. Initialization
  4. Time Loop (Navier-Stokes, Turbulence, Radiation, Scalars, Lagrangian, etc.)
  5. Finalization

Identifies major solver milestones (iterative solves, velocity prediction,
pressure correction, mass flux updates, Lagrangian steps, couplings) and all
user-defined function call sites (cs_user_*), clearly indicating whether each
function is called:
  - Once only [1x]
  - Once per time step [1x / Step]
  - Multiple times per time step [Nx / Step]

Outputs:
  - Standalone Interactive HTML with settable detail levels (Summary vs Detailed),
    zoom/pan, search, filter by frequency, and inspector sidebar.
  - Graphviz DOT files (summary and detailed).
  - Markdown with Mermaid diagram syntax.
"""

import os
import sys
import json
import argparse
import subprocess
from enum import Enum
from typing import Dict, List, Optional, Set, Any


# =============================================================================
# Enums and Data Schema
# =============================================================================

class Stage(str, Enum):
    SETUP = "SETUP"
    PREPROCESS = "PREPROCESS"
    INITIALIZATION = "INITIALIZATION"
    TIME_LOOP = "TIME_LOOP"
    FINALIZATION = "FINALIZATION"


class CallFrequency(str, Enum):
    ONCE = "ONCE"                       # Called once during run [1x]
    PER_TIME_STEP = "PER_TIME_STEP"     # Called once per time step [1x / Step]
    MULTIPLE_PER_STEP = "MULTIPLE"      # Called multiple times per step [Nx / Step]


class NodeCategory(str, Enum):
    CORE = "CORE"                       # Architectural driver / core function
    SOLVER_STAGE = "SOLVER_STAGE"       # Major numerical/algorithmic milestone
    PHYSICS_MODULE = "PHYSICS_MODULE"   # Specialized physics model
    USER_HOOK = "USER_HOOK"             # User-defined function (cs_user_*)


class DetailLevel(str, Enum):
    SUMMARY = "SUMMARY"                 # Shown in high-level executive view
    DETAILED = "DETAILED"               # Shown in detailed view
    BOTH = "BOTH"                       # Shown in both views


class CallNode:
    """Represents a node in the code_saturne execution call graph."""
    def __init__(
        self,
        node_id: str,
        name: str,
        label: str,
        stage: Stage,
        cluster: str,
        category: NodeCategory,
        frequency: CallFrequency,
        frequency_detail: str,
        file_ref: str,
        description: str,
        detail_level: DetailLevel = DetailLevel.BOTH,
        callees: Optional[List[str]] = None,
        is_iterative_solve: bool = False,
    ):
        self.node_id = node_id
        self.name = name
        self.label = label
        self.stage = stage
        self.cluster = cluster
        self.category = category
        self.frequency = frequency
        self.frequency_detail = frequency_detail
        self.file_ref = file_ref
        self.description = description
        self.detail_level = detail_level
        self.callees = callees or []
        self.is_iterative_solve = is_iterative_solve

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.node_id,
            "name": self.name,
            "label": self.label,
            "stage": self.stage.value,
            "cluster": self.cluster,
            "category": self.category.value,
            "frequency": self.frequency.value,
            "frequency_detail": self.frequency_detail,
            "file_ref": self.file_ref,
            "description": self.description,
            "detail_level": self.detail_level.value,
            "callees": self.callees,
            "is_iterative_solve": self.is_iterative_solve,
        }


# =============================================================================
# Graph Definition
# =============================================================================

def build_saturne_call_graph() -> Dict[str, CallNode]:
    """Constructs the comprehensive data model of code_saturne operations."""
    nodes: Dict[str, CallNode] = {}

    def add_node(node: CallNode):
        nodes[node.node_id] = node

    # -------------------------------------------------------------------------
    # TOP-LEVEL ENTRY POINT
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="main_entry",
        name="main / _run",
        label="code_saturne solver entry (_run)",
        stage=Stage.SETUP,
        cluster="setup_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Single execution at process startup",
        file_ref="src/apps/cs_solver.cxx:220",
        description="Main entry point for code_saturne solver. Initializes environment, "
                    "timers, halo buffers, coordinates setup, preprocessing, time stepping, "
                    "and clean finalization.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "cs_user_parallel_io", "cs_user_time_table", "cs_user_zones",
            "cs_setup", "cs_preprocess_mesh", "cs_time_stepping", "finalization_core"
        ]
    ))

    # -------------------------------------------------------------------------
    # STAGE 1: SETUP
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="cs_user_parallel_io",
        name="cs_user_parallel_io",
        label="cs_user_parallel_io",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at startup if rank count > 1",
        file_ref="src/apps/cs_solver.cxx:267",
        description="Allows configuring parallel IO heuristics, buffer sizes, "
                    "and MPI-IO algorithm parameters.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_time_table",
        name="cs_user_time_table",
        label="cs_user_time_table",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at startup during setup",
        file_ref="src/apps/cs_solver.cxx:274",
        description="Defines custom user time tables and interpolation curves for "
                    "time-dependent parameters.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_zones",
        name="cs_user_zones",
        label="cs_user_zones",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during domain setup",
        file_ref="src/apps/cs_solver.cxx:314",
        description="Defines custom volume and boundary geometric zones (head-loss zones, "
                    "mass source zones, boundary patches).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_setup",
        name="cs_setup",
        label="cs_setup (Calculation Setup)",
        stage=Stage.SETUP,
        cluster="setup_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Executed once during calculation preparation",
        file_ref="src/base/cs_setup.cpp:3856",
        description="Central calculation setup orchestrator. Reads user GUI XML parameters, "
                    "initializes flow, turbulence, and thermal models, instantiates variable "
                    "and property fields, and invokes early user configuration hooks.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "_init_setup", "_init_user", "cs_user_model", "cs_user_parameters",
            "cs_user_linear_solvers", "cs_user_time_moments", "cs_user_finalize_setup"
        ]
    ))

    add_node(CallNode(
        node_id="_init_setup",
        name="_init_setup",
        label="_init_setup (Module Defaults)",
        stage=Stage.SETUP,
        cluster="setup_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once inside cs_setup",
        file_ref="src/base/cs_setup.cpp:2730",
        description="Initializes internal module defaults before user parameters are read.",
        detail_level=DetailLevel.DETAILED
    ))

    add_node(CallNode(
        node_id="_init_user",
        name="_init_user",
        label="_init_user (Read GUI & Models)",
        stage=Stage.SETUP,
        cluster="setup_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once inside cs_setup",
        file_ref="src/base/cs_setup.cpp:2784",
        description="Parses GUI setup XML data, initializes time-step options, selects physical "
                    "models, and invokes cs_user_model.",
        detail_level=DetailLevel.DETAILED,
        callees=["cs_user_model"]
    ))

    add_node(CallNode(
        node_id="cs_user_model",
        name="cs_user_model",
        label="cs_user_model",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at earliest stage before field allocation",
        file_ref="src/base/cs_setup.cpp:2828",
        description="Selects physical models (Navier-Stokes, compressible, atmospheric, "
                    "turbulence RANS/LES, combustion, radiation, CDO, etc.) and user variables.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_parameters",
        name="cs_user_parameters",
        label="cs_user_parameters",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once after model fields are defined",
        file_ref="src/base/cs_setup.cpp:2981",
        description="User configuration of numerical parameters, relaxation factors, "
                    "time-stepping schemes, convergence tolerances, and physical constants.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_linear_solvers",
        name="cs_user_linear_solvers",
        label="cs_user_linear_solvers",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during solver setup",
        file_ref="src/apps/cs_solver.cxx:349",
        description="Configures native iterative linear solvers (Conjugate Gradient, BiCGStab, "
                    "GMRES, Multigrid AMG) or external solver libraries (PETSc, HYPRE, MUMPS, AmgX).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_time_moments",
        name="cs_user_time_moments",
        label="cs_user_time_moments",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during setup",
        file_ref="src/base/cs_setup.cpp:3000",
        description="Defines temporal statistics and moments (mean velocity, Reynolds stresses, "
                    "RMS fluctuations, temperature correlations).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_finalize_setup",
        name="cs_user_finalize_setup",
        label="cs_user_finalize_setup",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at the end of setup",
        file_ref="src/base/cs_setup.cpp:3020",
        description="Finalizes setup of equations, boundary conditions definitions, and CDO modules.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_matrix_tuning",
        name="cs_user_matrix_tuning",
        label="cs_user_matrix_tuning",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once before matrix initialization",
        file_ref="src/apps/cs_solver.cxx:475",
        description="Fine-tunes sparse matrix storage formats (MSR, CSR, native) and thread/GPU options.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_couplings",
        name="cs_user_syrthes/saturne_coupling",
        label="cs_user_syrthes / saturne_coupling",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during MPI coupling initialization",
        file_ref="src/apps/cs_solver.cxx:325",
        description="Defines coupled simulation interfaces with SYRTHES (thermal) or other "
                    "code_saturne instances.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_postprocess_setup",
        name="cs_user_postprocess_writers/meshes/probes",
        label="cs_user_postprocess_writers / meshes / probes",
        stage=Stage.SETUP,
        cluster="setup_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during post-processing setup",
        file_ref="src/apps/cs_solver.cxx:387",
        description="Defines output formats (EnSight, MED, CGNS), post-processing sub-meshes, "
                    "and monitoring point/line probes.",
        detail_level=DetailLevel.BOTH
    ))

    # -------------------------------------------------------------------------
    # STAGE 2: MESH PREPROCESSING
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="cs_preprocess_mesh",
        name="cs_preprocess_mesh",
        label="cs_preprocess_mesh (Mesh Engine)",
        stage=Stage.PREPROCESS,
        cluster="preprocess_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during mesh preparation phase",
        file_ref="src/base/cs_preprocess.cpp:150",
        description="Orchestrates mesh file reading, non-conformal face joining, periodic "
                    "boundary transformations, partitioning across MPI ranks, renumbering, "
                    "and mesh quality evaluation.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "cs_user_mesh_input", "cs_user_join", "cs_user_periodicity",
            "cs_user_mesh_modify", "cs_user_partition", "cs_mesh_quality"
        ]
    ))

    add_node(CallNode(
        node_id="cs_user_mesh_input",
        name="cs_user_mesh_input",
        label="cs_user_mesh_input",
        stage=Stage.PREPROCESS,
        cluster="preprocess_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at mesh loading time",
        file_ref="src/base/cs_preprocess.cpp:168",
        description="Specifies input mesh files, formats, and associated initial coordinate "
                    "transformations.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_join",
        name="cs_user_join",
        label="cs_user_join",
        stage=Stage.PREPROCESS,
        cluster="preprocess_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during mesh joining stage",
        file_ref="src/base/cs_preprocess.cpp:184",
        description="Defines tolerance and criteria for joining non-conformal mesh interfaces.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_periodicity",
        name="cs_user_periodicity",
        label="cs_user_periodicity",
        stage=Stage.PREPROCESS,
        cluster="preprocess_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during periodicity definition",
        file_ref="src/base/cs_preprocess.cpp:187",
        description="Defines periodic boundary transformations (translation, rotation).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_mesh_modify",
        name="cs_user_mesh_modify",
        label="cs_user_mesh_modify / smoothe",
        stage=Stage.PREPROCESS,
        cluster="preprocess_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during geometric preprocessing",
        file_ref="src/base/cs_preprocess.cpp:275",
        description="Modifies mesh coordinates, inserts internal boundaries, or triggers mesh smoothing.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_partition",
        name="cs_user_partition / numbering",
        label="cs_user_partition / numbering",
        stage=Stage.PREPROCESS,
        cluster="preprocess_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once prior to domain decomposition",
        file_ref="src/base/cs_preprocess.cpp:237",
        description="Configures parallel domain decomposition (SCOTCH, METIS, Hilbert) and renumbering.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_mesh_quality",
        name="cs_mesh_quality",
        label="cs_mesh_quality & bad cells postprocess",
        stage=Stage.PREPROCESS,
        cluster="preprocess_core",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during mesh validation",
        file_ref="src/mesh/cs_mesh_quality.cpp:210",
        description="Evaluates geometric quality criteria (orthogonality, non-conjunction, "
                    "expansion ratio) and detects degenerate cells.",
        detail_level=DetailLevel.BOTH
    ))

    # -------------------------------------------------------------------------
    # STAGE 3: INITIALIZATION
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="init_stage",
        name="Initialization Stage",
        label="Field & Model Initialization",
        stage=Stage.INITIALIZATION,
        cluster="init_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Executed once before entering time loop",
        file_ref="src/base/cs_time_stepping.cpp:180",
        description="Allocates variable fields, reads restart checkpoints, initializes fluid "
                    "properties, sets up porosity, wall models, and invokes user initialization.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "cs_initialize_fields", "cs_restart_read", "cs_porous_model_init",
            "cs_user_initialization", "cs_user_extra_operations_initialize",
            "cs_user_porosity", "cs_lagr_solve_initialize"
        ]
    ))

    add_node(CallNode(
        node_id="cs_initialize_fields",
        name="cs_initialize_fields_stage_0/1",
        label="cs_initialize_fields_stage_0/1",
        stage=Stage.INITIALIZATION,
        cluster="init_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during field setup",
        file_ref="src/base/cs_initialize_fields.cpp:120",
        description="Allocates and maps fields, sets default values for velocity, pressure, "
                    "turbulent quantities, temperature, and passive scalars.",
        detail_level=DetailLevel.DETAILED,
        callees=["cs_user_initialization"]
    ))

    add_node(CallNode(
        node_id="cs_restart_read",
        name="cs_time_stepping_read_checkpoint",
        label="Read Checkpoint (Restart)",
        stage=Stage.INITIALIZATION,
        cluster="init_core",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once if restart file is present",
        file_ref="src/base/cs_time_stepping.cpp:330",
        description="Reads previous solution state, time step counters, and auxiliary field "
                    "data from checkpoint files.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_porous_model_init",
        name="cs_porous_model_init",
        label="Porosity Model Initialization",
        stage=Stage.INITIALIZATION,
        cluster="init_physics",
        category=NodeCategory.PHYSICS_MODULE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once if porous media activated",
        file_ref="src/base/cs_time_stepping.cpp:355",
        description="Computes cell porosity and permeability tensors from scan, immersed boundaries, "
                    "or user definition.",
        detail_level=DetailLevel.DETAILED,
        callees=["cs_user_porosity"]
    ))

    add_node(CallNode(
        node_id="cs_user_porosity",
        name="cs_user_porosity",
        label="cs_user_porosity",
        stage=Stage.INITIALIZATION,
        cluster="init_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once during initialization if porous model is active",
        file_ref="src/base/cs_time_stepping.cpp:380",
        description="User definition of cell porosity factor (fluid volume fraction) and "
                    "directional resistance.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_initialization",
        name="cs_user_initialization",
        label="cs_user_initialization",
        stage=Stage.INITIALIZATION,
        cluster="init_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once before time loop (or upon restart to override)",
        file_ref="src/base/cs_initialize_fields.cpp:526",
        description="User routine for setting initial field distributions (velocity profiles, "
                    "temperatures, turbulent intensities, species concentrations).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_extra_operations_initialize",
        name="cs_user_extra_operations_initialize",
        label="cs_user_extra_operations_initialize",
        stage=Stage.INITIALIZATION,
        cluster="init_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at initial calculation state",
        file_ref="src/base/cs_initialize_fields.cpp:831",
        description="Pre-computation initialization of custom user tracking structures, "
                    "files, and supplementary memory allocations.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_lagr_solve_initialize",
        name="cs_lagr_solve_initialize",
        label="Lagrangian Initialization (cs_lagr_solve_initialize)",
        stage=Stage.INITIALIZATION,
        cluster="init_physics",
        category=NodeCategory.PHYSICS_MODULE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once if Lagrangian particle model is active",
        file_ref="src/base/cs_time_stepping.cpp:495",
        description="Initializes particle tracking arrays, boundary models, statistics, and "
                    "injects initial particle population.",
        detail_level=DetailLevel.BOTH
    ))

    # -------------------------------------------------------------------------
    # STAGE 4: TIME LOOP
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="cs_time_stepping",
        name="cs_time_stepping",
        label="Main Time Loop (cs_time_stepping)",
        stage=Stage.TIME_LOOP,
        cluster="timeloop_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Iterates over all time steps: do ... while(nt_cur < nt_max)",
        file_ref="src/base/cs_time_stepping.cpp:540",
        description="Governs time progression, adaptive time step evaluation, execution of "
                    "cs_solve_all(), Lagrangian tracking, temporal moments, checkpointing, and postprocessing.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "time_step_control", "cs_solve_all", "cs_lagr_solve_time_step",
            "cs_time_moment_update_all", "cs_user_extra_operations",
            "cs_time_stepping_write_checkpoint", "cs_post_default_write_variables"
        ]
    ))

    add_node(CallNode(
        node_id="time_step_control",
        name="Time Step Control & AMR",
        label="Time Step Control & AMR",
        stage=Stage.TIME_LOOP,
        cluster="timeloop_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at beginning of each time step",
        file_ref="src/base/cs_time_stepping.cpp:560",
        description="Evaluates CFL/Fourier criteria, updates dt, checks runaway convergence "
                    "and performs adaptive mesh refinement (AMR) if triggered.",
        detail_level=DetailLevel.DETAILED
    ))

    add_node(CallNode(
        node_id="cs_solve_all",
        name="cs_solve_all",
        label="Flow & Transport Equations (cs_solve_all)",
        stage=Stage.TIME_LOOP,
        cluster="solve_all",
        category=NodeCategory.CORE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once per time step",
        file_ref="src/base/cs_solve_all.cpp:890",
        description="Top-level solver orchestrator for incompressible/dilatable Navier-Stokes, "
                    "physical properties, velocity-pressure coupling, turbulence, radiative "
                    "transfer, and scalar transport equations.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "cs_physical_properties_update", "cs_head_losses_compute",
            "_solve_most", "_solve_turbulence", "cs_rad_transfer_solve",
            "cs_solve_transported_variables"
        ]
    ))

    add_node(CallNode(
        node_id="cs_physical_properties_update",
        name="cs_physical_properties_update",
        label="Physical Properties Update",
        stage=Stage.TIME_LOOP,
        cluster="solve_all",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at start of time step and upon property variation",
        file_ref="src/base/cs_physical_properties_default.cpp:1050",
        description="Updates density, laminar and turbulent viscosity, specific heat capacity, "
                    "thermal conductivity, and invokes user property laws.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_physical_properties", "cs_user_physical_properties_turb_viscosity"]
    ))

    add_node(CallNode(
        node_id="cs_user_physical_properties",
        name="cs_user_physical_properties",
        label="cs_user_physical_properties",
        stage=Stage.TIME_LOOP,
        cluster="solve_all_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at each time step to define variable physical properties",
        file_ref="src/base/cs_physical_properties_default.cpp:1080",
        description="User definition of temperature/pressure-dependent fluid properties "
                    "(density, viscosity, specific heat Cp, diffusivity).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_physical_properties_turb_viscosity",
        name="cs_user_physical_properties_turb_viscosity",
        label="cs_user_physical_properties_turb_viscosity",
        stage=Stage.TIME_LOOP,
        cluster="solve_all_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once per time step",
        file_ref="src/base/cs_physical_properties_default.cpp:1205",
        description="Allows custom user modification or clipping of turbulent eddy viscosity (mu_t).",
        detail_level=DetailLevel.DETAILED
    ))

    add_node(CallNode(
        node_id="cs_head_losses_compute",
        name="cs_head_losses_compute",
        label="Head Losses Computation",
        stage=Stage.TIME_LOOP,
        cluster="solve_all",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at each time step if head-loss zones are defined",
        file_ref="src/base/cs_head_losses.cpp:140",
        description="Computes anisotropic quadratic head loss tensor coefficients in volume zones.",
        detail_level=DetailLevel.DETAILED,
        callees=["cs_user_head_losses"]
    ))

    add_node(CallNode(
        node_id="cs_user_head_losses",
        name="cs_user_head_losses",
        label="cs_user_head_losses",
        stage=Stage.TIME_LOOP,
        cluster="solve_all_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at each time step for active head loss zones",
        file_ref="src/base/cs_head_losses.cpp:149",
        description="Defines directional head loss coefficients (K11, K22, K33, K12, K23, K13) "
                    "as functions of local velocity or flow conditions.",
        detail_level=DetailLevel.BOTH
    ))

    # --- Navier-Stokes Sub-iterations (_solve_most) ---
    add_node(CallNode(
        node_id="_solve_most",
        name="_solve_most",
        label="Velocity-Pressure Coupling Loop (_solve_most)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.CORE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Executes nterup inner sub-iterations (iterns = 1 .. nterup)",
        file_ref="src/base/cs_solve_all.cpp:445",
        description="Coordinates velocity-pressure sub-iterations, boundary conditions update, "
                    "Navier-Stokes equation assembly, and wall distance evaluation.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_boundary_conditions_set_coeffs", "cs_solve_navier_stokes"]
    ))

    add_node(CallNode(
        node_id="cs_boundary_conditions_set_coeffs",
        name="cs_boundary_conditions_set_coeffs",
        label="Boundary Conditions Setup (cs_boundary_conditions_set_coeffs)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called at each Navier-Stokes sub-iteration iterns",
        file_ref="src/base/cs_boundary_conditions_set_coeffs.cpp:840",
        description="Computes Dirichlet, Neumann, and Robin coefficients for all variables, "
                    "evaluating wall friction laws, roughness, and calling user boundary hooks.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_boundary_conditions"]
    ))

    add_node(CallNode(
        node_id="cs_user_boundary_conditions",
        name="cs_user_boundary_conditions",
        label="cs_user_boundary_conditions",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called at each time step and sub-iteration iterns",
        file_ref="src/base/cs_boundary_conditions_set_coeffs.cpp:843",
        description="User assignment of boundary types (inlet, outlet, wall, symmetry) and "
                    "boundary values (velocity profiles, pressures, temperatures, fluxes).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_solve_navier_stokes",
        name="cs_solve_navier_stokes",
        label="Navier-Stokes Solver (cs_solve_navier_stokes)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called each sub-iteration iterns for velocity-pressure resolution",
        file_ref="src/base/cs_solve_navier_stokes.cpp:3793",
        description="Performs the fractional step projection method: velocity prediction, "
                    "pressure Poisson equation solve, mass flux update on faces, and velocity correction.",
        detail_level=DetailLevel.BOTH,
        callees=[
            "_cs_mass_flux_prediction", "cs_velocity_prediction_1",
            "cs_pressure_correction", "cs_mass_flux", "cs_velocity_prediction_2"
        ]
    ))

    add_node(CallNode(
        node_id="_cs_mass_flux_prediction",
        name="_cs_mass_flux_prediction",
        label="Mass Flux Prediction (Low Mach)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called in Low Mach / dilatable algorithm",
        file_ref="src/base/cs_solve_navier_stokes.cpp:165",
        description="Predicts face mass flux from previous step for variable-density flows.",
        detail_level=DetailLevel.DETAILED
    ))

    add_node(CallNode(
        node_id="cs_velocity_prediction_1",
        name="cs_velocity_prediction (iappel = 1)",
        label="Velocity Prediction Step (iappel = 1)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called at sub-iteration iterns during prediction stage",
        file_ref="src/base/cs_solve_navier_stokes.cpp:1840",
        description="Assembles convective, diffusive, and explicit source terms for momentum; "
                    "solves the vector momentum equation for predicted velocity (u*).",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_source_terms", "cs_equation_iterative_solve_vector"]
    ))

    add_node(CallNode(
        node_id="cs_pressure_correction",
        name="cs_pressure_correction",
        label="Pressure Correction / Poisson Solve",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called at each sub-iteration iterns",
        file_ref="src/base/cs_pressure_correction.cpp:150",
        description="Applies Rhie-Chow interpolation to predicted velocity, computes divergence, "
                    "assembles and solves the Poisson pressure increment equation via SLES/AMG.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_sles_solve_pressure"],
        is_iterative_solve=True
    ))

    add_node(CallNode(
        node_id="cs_sles_solve_pressure",
        name="cs_sles_solve (Pressure Poisson)",
        label="Iterative Solve: Pressure Poisson",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Iterative linear solver execution",
        file_ref="src/alge/cs_sles.cpp:410",
        description="Solves the elliptic pressure Poisson system using Multigrid (AMG) or Conjugate Gradient.",
        detail_level=DetailLevel.DETAILED,
        is_iterative_solve=True
    ))

    add_node(CallNode(
        node_id="cs_mass_flux",
        name="cs_mass_flux",
        label="Mass Flux Update (cs_mass_flux)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called after pressure correction to enforce continuity",
        file_ref="src/base/cs_solve_navier_stokes.cpp:4167",
        description="Updates interior and boundary face mass fluxes using the corrected pressure "
                    "gradient, guaranteeing discrete mass conservation.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_velocity_prediction_2",
        name="cs_velocity_prediction (iappel = 2)",
        label="Velocity Correction Step (iappel = 2)",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called after pressure correction",
        file_ref="src/base/cs_solve_navier_stokes.cpp:4680",
        description="Corrects predicted cell velocities using the pressure increment gradient "
                    "and computes Navier-Stokes residual estimators.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_equation_iterative_solve_vector",
        name="cs_equation_iterative_solve_vector",
        label="Iterative Solve: Momentum",
        stage=Stage.TIME_LOOP,
        cluster="navier_stokes",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Vector iterative linear solve for velocity components",
        file_ref="src/base/cs_solve_equation.cpp:2400",
        description="Iterative linear solution of the convection-diffusion momentum matrix (BiCGStab / Jacobi).",
        detail_level=DetailLevel.DETAILED,
        is_iterative_solve=True
    ))

    # --- Turbulence ---
    add_node(CallNode(
        node_id="_solve_turbulence",
        name="_solve_turbulence",
        label="Turbulence Equations (_solve_turbulence)",
        stage=Stage.TIME_LOOP,
        cluster="turbulence",
        category=NodeCategory.PHYSICS_MODULE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once per time step after Navier-Stokes convergence",
        file_ref="src/base/cs_solve_all.cpp:320",
        description="Solves turbulent transport equations for RANS (k-epsilon, k-omega SST, "
                    "v2-f, Reynolds Stress Model Rij/DRSM) or LES subgrid models.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_source_terms", "cs_turbulence_iterative_solve"]
    ))

    add_node(CallNode(
        node_id="cs_turbulence_iterative_solve",
        name="cs_turbulence_iterative_solve",
        label="Iterative Solve: Turbulence (k, eps, Rij)",
        stage=Stage.TIME_LOOP,
        cluster="turbulence",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Solves linear systems for each turbulent scalar or stress component",
        file_ref="src/turb/cs_turbulence_ke.cpp:1480",
        description="Iterative linear solution of coupled or segregated turbulence equations.",
        detail_level=DetailLevel.DETAILED,
        is_iterative_solve=True
    ))

    # --- Radiative Transfer ---
    add_node(CallNode(
        node_id="cs_rad_transfer_solve",
        name="cs_rad_transfer_solve",
        label="Radiative Transfer Solve (cs_rad_transfer_solve)",
        stage=Stage.TIME_LOOP,
        cluster="radiation",
        category=NodeCategory.PHYSICS_MODULE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called at specified radiative time step frequency",
        file_ref="src/rayt/cs_rad_transfer_solve.cpp:200",
        description="Solves radiative transfer equation using P-1 or Discrete Ordinates Method (DOM), "
                    "computing incident radiation and radiative divergence source terms.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_rad_transfer_bcs", "cs_user_rad_transfer_absorption"]
    ))

    add_node(CallNode(
        node_id="cs_user_rad_transfer_bcs",
        name="cs_user_radiative_transfer_bcs",
        label="cs_user_radiative_transfer_bcs",
        stage=Stage.TIME_LOOP,
        cluster="radiation_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called per radiative solve step",
        file_ref="src/base/cs_prototypes.h:1280",
        description="Specifies wall emissivity, outside temperature, and gray/reflecting boundary conditions.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_rad_transfer_absorption",
        name="cs_user_rad_transfer_absorption",
        label="cs_user_rad_transfer_absorption",
        stage=Stage.TIME_LOOP,
        cluster="radiation_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called per radiative solve step",
        file_ref="src/base/cs_prototypes.h:1180",
        description="Defines the fluid's absorption coefficient Ck (transparent medium, gas mixtures, soot).",
        detail_level=DetailLevel.BOTH
    ))

    # --- Transported Scalars ---
    add_node(CallNode(
        node_id="cs_solve_transported_variables",
        name="cs_solve_transported_variables",
        label="Scalar & Energy Equations (cs_solve_transported_variables)",
        stage=Stage.TIME_LOOP,
        cluster="scalars",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called per time step for each active thermal and user scalar",
        file_ref="src/base/cs_solve_transported_variables.cpp:110",
        description="Assembles and solves convection-diffusion equations for temperature, "
                    "enthalpy, and user-defined transported passive or reactive scalars.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_source_terms", "cs_solve_equation_scalar"]
    ))

    add_node(CallNode(
        node_id="cs_solve_equation_scalar",
        name="cs_solve_equation (Scalar)",
        label="Iterative Solve: Transported Scalars",
        stage=Stage.TIME_LOOP,
        cluster="scalars",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Executed for each scalar equation",
        file_ref="src/base/cs_solve_equation.cpp:1310",
        description="Assembles system matrix and solves linear equation for scalar fields.",
        detail_level=DetailLevel.DETAILED,
        is_iterative_solve=True
    ))

    # --- Common User Source Terms Hook ---
    add_node(CallNode(
        node_id="cs_user_source_terms",
        name="cs_user_source_terms",
        label="cs_user_source_terms",
        stage=Stage.TIME_LOOP,
        cluster="timeloop_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.MULTIPLE_PER_STEP,
        frequency_detail="Called multiple times per time step (for momentum, turbulence, and each scalar)",
        file_ref="src/base/cs_solve_equation.cpp:1315",
        description="Defines explicit (st_exp) and implicit (st_imp) volumetric source terms for "
                    "momentum, thermal energy, turbulence variables, and user scalars: "
                    "rho*vol*dF/dt = st_imp*F + st_exp.",
        detail_level=DetailLevel.BOTH
    ))

    # --- Lagrangian Tracking ---
    add_node(CallNode(
        node_id="cs_lagr_solve_time_step",
        name="cs_lagr_solve_time_step",
        label="Lagrangian Particle Tracking (cs_lagr_solve_time_step)",
        stage=Stage.TIME_LOOP,
        cluster="lagrangian",
        category=NodeCategory.PHYSICS_MODULE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once per time step if Lagrangian module is active",
        file_ref="src/base/cs_time_stepping.cpp:710",
        description="Integrates particle equations of motion (drag, gravity, Saffman lift, "
                    "Brownian motion), handles particle-wall collisions, deposition, and two-way coupling.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_lagr_bcs", "cs_user_lagr_particle"]
    ))

    add_node(CallNode(
        node_id="cs_user_lagr_bcs",
        name="cs_user_lagr_boundary_conditions",
        label="cs_user_lagr_boundary_conditions",
        stage=Stage.TIME_LOOP,
        cluster="lagrangian_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called per time step for boundary interaction",
        file_ref="src/user/cs_user_lagr_boundary_conditions.cpp:40",
        description="Specifies particle boundary behavior (rebound coefficient, absorption, deposition, fouling).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_user_lagr_particle",
        name="cs_user_lagr_particle / volume_conditions",
        label="cs_user_lagr_particle / volume_conditions",
        stage=Stage.TIME_LOOP,
        cluster="lagrangian_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called during particle tracking step",
        file_ref="src/user/cs_user_lagr_particle.cpp:50",
        description="Defines continuous particle injection, diameter distributions, and custom particle forces.",
        detail_level=DetailLevel.BOTH
    ))

    # --- Post-Solve / Diagnostics / Post-processing ---
    add_node(CallNode(
        node_id="cs_time_moment_update_all",
        name="cs_time_moment_update_all",
        label="Temporal Moments Accumulation",
        stage=Stage.TIME_LOOP,
        cluster="post_solve",
        category=NodeCategory.SOLVER_STAGE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once per time step",
        file_ref="src/base/cs_time_stepping.cpp:728",
        description="Accumulates temporal means, variances, and cross-correlations across the time domain.",
        detail_level=DetailLevel.DETAILED
    ))

    add_node(CallNode(
        node_id="cs_user_extra_operations",
        name="cs_user_extra_operations",
        label="cs_user_extra_operations",
        stage=Stage.TIME_LOOP,
        cluster="post_solve_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called once at the end of each time step",
        file_ref="src/base/cs_time_stepping.cpp:758",
        description="General-purpose user routine called at the end of each time step for "
                    "custom data extraction, online integrals, monitoring, and simulation steering.",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_time_stepping_write_checkpoint",
        name="cs_time_stepping_write_checkpoint",
        label="Write Checkpoint (Restart)",
        stage=Stage.TIME_LOOP,
        cluster="post_solve",
        category=NodeCategory.CORE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Triggered based on checkpoint time frequency or user signal",
        file_ref="src/base/cs_time_stepping.cpp:805",
        description="Writes binary restart checkpoint files (main variables, auxiliary fields, "
                    "mesh topology if moving mesh).",
        detail_level=DetailLevel.BOTH
    ))

    add_node(CallNode(
        node_id="cs_post_default_write_variables",
        name="cs_post_default_write_variables",
        label="Visualization & Probes Output",
        stage=Stage.TIME_LOOP,
        cluster="post_solve",
        category=NodeCategory.CORE,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Evaluated every time step, outputs according to configured writer frequency",
        file_ref="src/base/cs_time_stepping.cpp:830",
        description="Writes field visualization outputs (EnSight, MED, CGNS) and writes time-history "
                    "probe monitoring files.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_postprocess_values"]
    ))

    add_node(CallNode(
        node_id="cs_user_postprocess_values",
        name="cs_user_postprocess_values",
        label="cs_user_postprocess_values / activate",
        stage=Stage.TIME_LOOP,
        cluster="post_solve_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.PER_TIME_STEP,
        frequency_detail="Called when a postprocessing writer is activated",
        file_ref="src/base/cs_post.cpp:8176",
        description="Allows calculating and writing custom derived user fields onto visualization meshes.",
        detail_level=DetailLevel.BOTH
    ))

    # -------------------------------------------------------------------------
    # STAGE 5: FINALIZATION
    # -------------------------------------------------------------------------
    add_node(CallNode(
        node_id="finalization_core",
        name="Finalization Stage",
        label="Clean Finalization & Deallocation",
        stage=Stage.FINALIZATION,
        cluster="finalize_core",
        category=NodeCategory.CORE,
        frequency=CallFrequency.ONCE,
        frequency_detail="Executed once after exiting time loop",
        file_ref="src/apps/cs_solver.cxx:540",
        description="Calls user finalization, frees solver memory, matrices, SLES structures, "
                    "couplings, probes, and logs final performance statistics.",
        detail_level=DetailLevel.BOTH,
        callees=["cs_user_extra_operations_finalize"]
    ))

    add_node(CallNode(
        node_id="cs_user_extra_operations_finalize",
        name="cs_user_extra_operations_finalize",
        label="cs_user_extra_operations_finalize",
        stage=Stage.FINALIZATION,
        cluster="finalize_user",
        category=NodeCategory.USER_HOOK,
        frequency=CallFrequency.ONCE,
        frequency_detail="Called once at the very end of calculation",
        file_ref="src/apps/cs_solver.cxx:542",
        description="Final user hook for closing custom files, printing execution balance sheets, "
                    "and releasing user-allocated memory.",
        detail_level=DetailLevel.BOTH
    ))

    return nodes


# =============================================================================
# Exporter: Interactive HTML
# =============================================================================

def generate_interactive_html(nodes: Dict[str, CallNode], output_path: str):
    """Generates a rich, self-contained interactive HTML visualization."""
    nodes_data = [node.to_dict() for node in nodes.values()]
    nodes_json = json.dumps(nodes_data, indent=2)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>code_saturne Execution Flow & Call Graph</title>
  <style>
    :root {{
      --bg: #0f172a;
      --card-bg: #1e293b;
      --card-border: #334155;
      --text: #f8fafc;
      --text-muted: #94a3b8;
      --primary: #38bdf8;
      --primary-dark: #0284c7;
      
      --freq-once: #38bdf8;
      --freq-once-bg: rgba(56, 189, 248, 0.15);
      --freq-step: #4ade80;
      --freq-step-bg: rgba(74, 222, 128, 0.15);
      --freq-multi: #fb923c;
      --freq-multi-bg: rgba(251, 146, 60, 0.15);

      --user-color: #fbbf24;
      --user-bg: rgba(251, 191, 36, 0.15);
      --core-color: #94a3b8;
      --solve-color: #c084fc;
      --physics-color: #38bdf8;

      --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    }}

    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}

    body {{
      font-family: var(--font-sans);
      background-color: var(--bg);
      color: var(--text);
      line-height: 1.5;
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }}

    /* Top Navigation Header */
    header {{
      background-color: var(--card-bg);
      border-bottom: 1px solid var(--card-border);
      padding: 0.75rem 1.5rem;
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 1rem;
      z-index: 10;
    }}

    .header-title h1 {{
      font-size: 1.25rem;
      font-weight: 700;
      color: var(--primary);
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }}

    .header-title span {{
      font-size: 0.8rem;
      color: var(--text-muted);
      font-weight: 400;
    }}

    .controls {{
      display: flex;
      align-items: center;
      gap: 1rem;
      flex-wrap: wrap;
    }}

    .control-group {{
      display: flex;
      align-items: center;
      background: var(--bg);
      border: 1px solid var(--card-border);
      border-radius: 6px;
      padding: 2px;
    }}

    .segmented-btn {{
      background: transparent;
      border: none;
      color: var(--text-muted);
      padding: 0.35rem 0.75rem;
      font-size: 0.8rem;
      font-weight: 500;
      cursor: pointer;
      border-radius: 4px;
      transition: all 0.2s ease;
    }}

    .segmented-btn.active {{
      background: var(--primary-dark);
      color: #fff;
    }}

    .search-box {{
      position: relative;
    }}

    .search-box input {{
      background: var(--bg);
      border: 1px solid var(--card-border);
      color: var(--text);
      padding: 0.35rem 0.75rem 0.35rem 2rem;
      font-size: 0.825rem;
      border-radius: 6px;
      width: 220px;
      outline: none;
      transition: border-color 0.2s;
    }}

    .search-box input:focus {{
      border-color: var(--primary);
    }}

    .search-icon {{
      position: absolute;
      left: 0.6rem;
      top: 50%;
      transform: translateY(-50%);
      color: var(--text-muted);
      font-size: 0.8rem;
      pointer-events: none;
    }}

    /* Legend Bar */
    .legend-bar {{
      background: rgba(30, 41, 59, 0.75);
      border-bottom: 1px solid var(--card-border);
      padding: 0.4rem 1.5rem;
      display: flex;
      align-items: center;
      gap: 1.5rem;
      font-size: 0.75rem;
      overflow-x: auto;
    }}

    .legend-item {{
      display: flex;
      align-items: center;
      gap: 0.4rem;
      cursor: pointer;
      user-select: none;
    }}

    .legend-item.dimmed {{
      opacity: 0.35;
    }}

    .badge {{
      display: inline-block;
      padding: 0.15rem 0.45rem;
      font-size: 0.68rem;
      font-weight: 700;
      border-radius: 4px;
      letter-spacing: 0.02em;
    }}

    .badge-once {{
      color: var(--freq-once);
      background: var(--freq-once-bg);
      border: 1px solid var(--freq-once);
    }}

    .badge-step {{
      color: var(--freq-step);
      background: var(--freq-step-bg);
      border: 1px solid var(--freq-step);
    }}

    .badge-multi {{
      color: var(--freq-multi);
      background: var(--freq-multi-bg);
      border: 1px solid var(--freq-multi);
    }}

    .badge-user {{
      color: var(--user-color);
      background: var(--user-bg);
      border: 1px solid var(--user-color);
    }}

    /* Main Content Area */
    .main-wrapper {{
      display: flex;
      flex: 1;
      height: calc(100vh - 95px);
      position: relative;
      overflow: hidden;
    }}

    /* Graph Canvas */
    .canvas-container {{
      flex: 1;
      position: relative;
      overflow: auto;
      background: radial-gradient(circle at 50% 50%, #1e293b 0%, #0f172a 100%);
      cursor: grab;
      display: flex;
      justify-content: center;
      padding: 2.5rem;
    }}

    .canvas-container:active {{
      cursor: grabbing;
    }}

    .graph-pipeline {{
      display: flex;
      flex-direction: column;
      gap: 2rem;
      width: 100%;
      max-width: 1100px;
      margin: 0 auto;
      transition: transform 0.2s ease;
      transform-origin: top center;
    }}

    /* Stage Cluster Box */
    .stage-cluster {{
      background: var(--card-bg);
      border: 1px solid var(--card-border);
      border-radius: 8px;
      padding: 1.25rem;
      box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
      position: relative;
    }}

    .stage-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      border-bottom: 1px solid var(--card-border);
      padding-bottom: 0.5rem;
      margin-bottom: 1rem;
    }}

    .stage-title {{
      font-size: 0.95rem;
      font-weight: 700;
      color: var(--primary);
      text-transform: uppercase;
      letter-spacing: 0.05em;
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }}

    .stage-badge {{
      font-size: 0.7rem;
      color: var(--text-muted);
      font-weight: 500;
    }}

    .nodes-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
      gap: 0.85rem;
    }}

    /* Node Card */
    .node-card {{
      background: rgba(15, 23, 42, 0.65);
      border: 1px solid var(--card-border);
      border-radius: 6px;
      padding: 0.75rem 0.85rem;
      cursor: pointer;
      transition: all 0.2s ease;
      display: flex;
      flex-direction: column;
      gap: 0.4rem;
      position: relative;
    }}

    .node-card:hover {{
      border-color: var(--primary);
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(56, 189, 248, 0.15);
    }}

    .node-card.selected {{
      border-color: var(--primary);
      background: rgba(56, 189, 248, 0.08);
      box-shadow: 0 0 0 2px var(--primary-dark);
    }}

    .node-card.user-hook {{
      border-left: 3px solid var(--user-color);
    }}

    .node-card.solver-stage {{
      border-left: 3px solid var(--solve-color);
    }}

    .node-card.iterative-solve {{
      border-left: 3px solid #ec4899;
    }}

    .node-card.dimmed {{
      opacity: 0.2;
      pointer-events: none;
    }}

    .node-header {{
      display: flex;
      justify-content: space-between;
      align-items: flex-start;
      gap: 0.5rem;
    }}

    .node-name {{
      font-family: var(--font-mono);
      font-size: 0.85rem;
      font-weight: 600;
      color: #fff;
      word-break: break-all;
    }}

    .node-label {{
      font-size: 0.75rem;
      color: var(--text-muted);
    }}

    .node-footer {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-top: auto;
      padding-top: 0.4rem;
    }}

    .node-file {{
      font-family: var(--font-mono);
      font-size: 0.65rem;
      color: var(--text-muted);
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      max-width: 170px;
    }}

    /* Inspector Sidebar */
    .inspector {{
      width: 380px;
      background: var(--card-bg);
      border-left: 1px solid var(--card-border);
      padding: 1.5rem;
      overflow-y: auto;
      display: flex;
      flex-direction: column;
      gap: 1.25rem;
      transition: transform 0.3s ease;
      z-index: 5;
    }}

    .inspector.collapsed {{
      transform: translateX(100%);
      display: none;
    }}

    .inspector-header {{
      display: flex;
      justify-content: space-between;
      align-items: flex-start;
      border-bottom: 1px solid var(--card-border);
      padding-bottom: 0.75rem;
    }}

    .inspector-title h2 {{
      font-size: 1.05rem;
      color: var(--primary);
      font-family: var(--font-mono);
      word-break: break-all;
    }}

    .inspector-title p {{
      font-size: 0.8rem;
      color: var(--text-muted);
      margin-top: 2px;
    }}

    .close-btn {{
      background: transparent;
      border: none;
      color: var(--text-muted);
      font-size: 1.25rem;
      cursor: pointer;
      line-height: 1;
    }}

    .close-btn:hover {{
      color: #fff;
    }}

    .inspector-section {{
      display: flex;
      flex-direction: column;
      gap: 0.4rem;
    }}

    .section-label {{
      font-size: 0.7rem;
      text-transform: uppercase;
      font-weight: 700;
      color: var(--text-muted);
      letter-spacing: 0.05em;
    }}

    .inspector-desc {{
      font-size: 0.85rem;
      color: var(--text);
      line-height: 1.6;
    }}

    .info-table {{
      width: 100%;
      font-size: 0.8rem;
      border-collapse: collapse;
    }}

    .info-table td {{
      padding: 0.35rem 0;
      border-bottom: 1px solid rgba(255, 255, 255, 0.05);
    }}

    .info-table td.prop {{
      color: var(--text-muted);
      width: 35%;
    }}

    .info-table td.val {{
      font-family: var(--font-mono);
      font-size: 0.75rem;
    }}

    .callees-list {{
      display: flex;
      flex-direction: column;
      gap: 0.3rem;
    }}

    .callee-pill {{
      background: rgba(15, 23, 42, 0.5);
      border: 1px solid var(--card-border);
      padding: 0.3rem 0.5rem;
      border-radius: 4px;
      font-family: var(--font-mono);
      font-size: 0.75rem;
      cursor: pointer;
    }}

    .callee-pill:hover {{
      border-color: var(--primary);
      color: var(--primary);
    }}

    /* Zoom Controls Floating Box */
    .zoom-toolbar {{
      position: absolute;
      bottom: 1.5rem;
      left: 1.5rem;
      background: var(--card-bg);
      border: 1px solid var(--card-border);
      border-radius: 6px;
      display: flex;
      flex-direction: column;
      padding: 2px;
      box-shadow: 0 4px 10px rgba(0, 0, 0, 0.3);
      z-index: 5;
    }}

    .zoom-btn {{
      background: transparent;
      border: none;
      color: var(--text);
      padding: 0.4rem 0.6rem;
      cursor: pointer;
      font-size: 0.9rem;
      font-weight: bold;
    }}

    .zoom-btn:hover {{
      color: var(--primary);
    }}

    /* Arrow Connector Styling between Stages */
    .stage-connector {{
      display: flex;
      justify-content: center;
      align-items: center;
      color: var(--text-muted);
      font-size: 1.25rem;
      margin: -1rem 0;
      z-index: 1;
    }}
  </style>
</head>
<body>

  <!-- Header -->
  <header>
    <div class="header-title">
      <h1>code_saturne User-Defined Functions Call Graph</h1>
      <span>Interactive Solver Workflow & User Call Site Reference</span>
    </div>

    <div class="controls">
      <!-- Search Input -->
      <div class="search-box">
        <span class="search-icon">🔍</span>
        <input type="text" id="searchInput" placeholder="Search functions (e.g. source_terms)...">
      </div>

      <!-- Detail Switcher -->
      <div class="control-group">
        <button class="segmented-btn" id="btnSummary" onclick="setDetailLevel('SUMMARY')">Executive Overview</button>
        <button class="segmented-btn active" id="btnDetailed" onclick="setDetailLevel('DETAILED')">Detailed Pipeline</button>
      </div>

      <!-- Quick Toggles -->
      <div class="control-group">
        <button class="segmented-btn" id="btnToggleUser" onclick="toggleUserOnly()">Show [USER] Only</button>
      </div>
    </div>
  </header>

  <!-- Legend Bar -->
  <div class="legend-bar">
    <span style="font-weight: 600; color: var(--text-muted);">EXECUTION FREQUENCY:</span>
    <div class="legend-item" id="legOnce" onclick="toggleFrequencyFilter('ONCE')">
      <span class="badge badge-once">[1x]</span>
      <span>Once Only (Setup / Init / Finish)</span>
    </div>
    <div class="legend-item" id="legStep" onclick="toggleFrequencyFilter('PER_TIME_STEP')">
      <span class="badge badge-step">[1x / Step]</span>
      <span>Once per Time Step</span>
    </div>
    <div class="legend-item" id="legMulti" onclick="toggleFrequencyFilter('MULTIPLE')">
      <span class="badge badge-multi">[Nx / Step]</span>
      <span>Multiple per Step (Sub-iterations / Equations)</span>
    </div>
    <span style="border-left: 1px solid var(--card-border); height: 16px; margin: 0 0.5rem;"></span>
    <span style="font-weight: 600; color: var(--text-muted);">NODE TYPES:</span>
    <div class="legend-item">
      <span class="badge badge-user">[USER]</span>
      <span>User Call Site (cs_user_*)</span>
    </div>
    <div class="legend-item">
      <span style="width: 12px; height: 12px; border-left: 3px solid var(--solve-color); display: inline-block;"></span>
      <span>Major Solver Milestone</span>
    </div>
    <div class="legend-item">
      <span style="width: 12px; height: 12px; border-left: 3px solid #ec4899; display: inline-block;"></span>
      <span>Iterative Linear Solve</span>
    </div>
  </div>

  <!-- Main View -->
  <div class="main-wrapper">
    <!-- Graph Canvas Container -->
    <div class="canvas-container" id="canvasContainer">
      <div class="graph-pipeline" id="graphPipeline">
        <!-- Stages will be populated dynamically -->
      </div>

      <!-- Zoom Toolbar -->
      <div class="zoom-toolbar">
        <button class="zoom-btn" onclick="zoomIn()" title="Zoom In">+</button>
        <button class="zoom-btn" onclick="zoomReset()" title="Reset Zoom">1:1</button>
        <button class="zoom-btn" onclick="zoomOut()" title="Zoom Out">-</button>
      </div>
    </div>

    <!-- Inspector Sidebar -->
    <div class="inspector collapsed" id="inspector">
      <div class="inspector-header">
        <div class="inspector-title">
          <h2 id="insName">cs_user_source_terms</h2>
          <p id="insLabel">User source terms hook</p>
        </div>
        <button class="close-btn" onclick="closeInspector()">&times;</button>
      </div>

      <div class="inspector-section">
        <div class="section-label">Classification</div>
        <table class="info-table">
          <tr><td class="prop">Stage</td><td class="val" id="insStage">TIME_LOOP</td></tr>
          <tr><td class="prop">Category</td><td class="val" id="insCategory">USER_HOOK</td></tr>
          <tr><td class="prop">Frequency</td><td class="val" id="insFreq">[Nx / Step]</td></tr>
          <tr><td class="prop">Source File</td><td class="val" id="insFile">src/...</td></tr>
        </table>
      </div>

      <div class="inspector-section">
        <div class="section-label">Frequency Rationale</div>
        <p class="inspector-desc" id="insFreqDetail">Detailed explanation...</p>
      </div>

      <div class="inspector-section">
        <div class="section-label">Function Purpose & Theory</div>
        <p class="inspector-desc" id="insDesc">Description...</p>
      </div>

      <div class="inspector-section" id="calleesSection">
        <div class="section-label">Direct Calls / Downstream Sub-steps</div>
        <div class="callees-list" id="insCallees">
          <!-- Callee pills inserted here -->
        </div>
      </div>
    </div>
  </div>

  <!-- JavaScript Engine -->
  <script>
    const allNodes = {nodes_json};
    let currentDetail = 'DETAILED';
    let filterUserOnly = false;
    let activeFrequencies = new Set(['ONCE', 'PER_TIME_STEP', 'MULTIPLE']);
    let currentSearchTerm = '';
    let currentZoom = 1.0;
    let selectedNodeId = null;

    const stageDefinitions = [
      {{ id: 'SETUP', title: '1. Setup Stage (cs_setup)', badge: 'Process startup & model configuration' }},
      {{ id: 'PREPROCESS', title: '2. Mesh Preprocessing (cs_preprocess_mesh)', badge: 'Mesh input, joining, partitioning, numbering' }},
      {{ id: 'INITIALIZATION', title: '3. Field & Physics Initialization', badge: 'Restart read, initial conditions, porosity' }},
      {{ id: 'TIME_LOOP', title: '4. Main Time Loop (cs_time_stepping & cs_solve_all)', badge: 'Time progression, velocity prediction, Poisson solve, scalars, Lagrangian' }},
      {{ id: 'FINALIZATION', title: '5. Finalization Stage', badge: 'Memory release, couplings & MPI teardown' }}
    ];

    function renderPipeline() {{
      const container = document.getElementById('graphPipeline');
      container.innerHTML = '';

      stageDefinitions.forEach((stageDef, index) => {{
        // Filter nodes belonging to this stage and current filters
        const stageNodes = allNodes.filter(n => {{
          if (n.stage !== stageDef.id) return false;
          if (currentDetail === 'SUMMARY' && n.detail_level === 'DETAILED') return false;
          if (filterUserOnly && n.category !== 'USER_HOOK') return false;
          if (!activeFrequencies.has(n.frequency)) return false;
          if (currentSearchTerm) {{
            const term = currentSearchTerm.toLowerCase();
            const matches = n.name.toLowerCase().includes(term) ||
                            n.label.toLowerCase().includes(term) ||
                            n.description.toLowerCase().includes(term);
            if (!matches) return false;
          }}
          return true;
        }});

        // If no nodes match this stage under current filters, still show container or hide
        if (stageNodes.length === 0 && (filterUserOnly || currentSearchTerm)) {{
          return;
        }}

        const stageCluster = document.createElement('div');
        stageCluster.className = 'stage-cluster';

        stageCluster.innerHTML = `
          <div class="stage-header">
            <div class="stage-title">
              <span>${{stageDef.title}}</span>
            </div>
            <div class="stage-badge">${{stageDef.badge}} (${{stageNodes.length}} functions)</div>
          </div>
          <div class="nodes-grid" id="grid-${{stageDef.id}}"></div>
        `;

        const grid = stageCluster.querySelector(`#grid-${{stageDef.id}}`);

        stageNodes.forEach(node => {{
          const card = document.createElement('div');
          card.className = 'node-card';
          if (node.id === selectedNodeId) card.classList.add('selected');
          if (node.category === 'USER_HOOK') card.classList.add('user-hook');
          if (node.category === 'SOLVER_STAGE') card.classList.add('solver-stage');
          if (node.is_iterative_solve) card.classList.add('iterative-solve');

          let freqBadge = '';
          if (node.frequency === 'ONCE') freqBadge = '<span class="badge badge-once">[1x]</span>';
          else if (node.frequency === 'PER_TIME_STEP') freqBadge = '<span class="badge badge-step">[1x / Step]</span>';
          else if (node.frequency === 'MULTIPLE') freqBadge = '<span class="badge badge-multi">[Nx / Step]</span>';

          let userBadge = node.category === 'USER_HOOK' ? '<span class="badge badge-user">[USER]</span>' : '';

          card.innerHTML = `
            <div class="node-header">
              <span class="node-name">${{node.name}}</span>
              <div style="display: flex; gap: 4px;">
                ${{userBadge}}
                ${{freqBadge}}
              </div>
            </div>
            <div class="node-label">${{node.label}}</div>
            <div class="node-footer">
              <span class="node-file" title="${{node.file_ref}}">${{node.file_ref}}</span>
            </div>
          `;

          card.onclick = () => selectNode(node);
          grid.appendChild(card);
        }});

        container.appendChild(stageCluster);

        // Add visual stage connector arrow between stages
        if (index < stageDefinitions.length - 1) {{
          const connector = document.createElement('div');
          connector.className = 'stage-connector';
          connector.innerHTML = '⬇';
          container.appendChild(connector);
        }}
      }});
    }}

    function selectNode(node) {{
      selectedNodeId = node.id;
      const inspector = document.getElementById('inspector');
      inspector.classList.remove('collapsed');

      document.getElementById('insName').textContent = node.name;
      document.getElementById('insLabel').textContent = node.label;
      document.getElementById('insStage').textContent = node.stage;
      document.getElementById('insCategory').textContent = node.category;

      let freqText = node.frequency === 'ONCE' ? '[1x] Once Only' :
                     node.frequency === 'PER_TIME_STEP' ? '[1x / Step] Per Time Step' :
                     '[Nx / Step] Multiple times per Step';
      document.getElementById('insFreq').textContent = freqText;
      document.getElementById('insFile').textContent = node.file_ref;
      document.getElementById('insFreqDetail').textContent = node.frequency_detail;
      document.getElementById('insDesc').textContent = node.description;

      const calleesContainer = document.getElementById('insCallees');
      calleesContainer.innerHTML = '';
      if (node.callees && node.callees.length > 0) {{
        document.getElementById('calleesSection').style.display = 'flex';
        node.callees.forEach(calleeId => {{
          const calleeNode = allNodes.find(n => n.id === calleeId);
          const pill = document.createElement('div');
          pill.className = 'callee-pill';
          pill.textContent = calleeNode ? `↳ ${{calleeNode.name}}` : `↳ ${{calleeId}}`;
          if (calleeNode) {{
            pill.onclick = () => selectNode(calleeNode);
          }}
          calleesContainer.appendChild(pill);
        }});
      }} else {{
        document.getElementById('calleesSection').style.display = 'none';
      }}

      // Highlight selected card visually
      document.querySelectorAll('.node-card').forEach(c => c.classList.remove('selected'));
      renderPipeline();
    }}

    function closeInspector() {{
      document.getElementById('inspector').classList.add('collapsed');
      selectedNodeId = null;
      renderPipeline();
    }}

    function setDetailLevel(level) {{
      currentDetail = level;
      document.getElementById('btnSummary').classList.toggle('active', level === 'SUMMARY');
      document.getElementById('btnDetailed').classList.toggle('active', level === 'DETAILED');
      renderPipeline();
    }}

    function toggleUserOnly() {{
      filterUserOnly = !filterUserOnly;
      document.getElementById('btnToggleUser').classList.toggle('active', filterUserOnly);
      renderPipeline();
    }}

    function toggleFrequencyFilter(freq) {{
      const legElem = freq === 'ONCE' ? document.getElementById('legOnce') :
                      freq === 'PER_TIME_STEP' ? document.getElementById('legStep') :
                      document.getElementById('legMulti');

      if (activeFrequencies.has(freq)) {{
        if (activeFrequencies.size > 1) {{
          activeFrequencies.delete(freq);
          legElem.classList.add('dimmed');
        }}
      }} else {{
        activeFrequencies.add(freq);
        legElem.classList.remove('dimmed');
      }}
      renderPipeline();
    }}

    document.getElementById('searchInput').addEventListener('input', (e) => {{
      currentSearchTerm = e.target.value.trim();
      renderPipeline();
    }});

    function zoomIn() {{
      currentZoom = Math.min(currentZoom + 0.15, 1.8);
      applyZoom();
    }}

    function zoomOut() {{
      currentZoom = Math.max(currentZoom - 0.15, 0.5);
      applyZoom();
    }}

    function zoomReset() {{
      currentZoom = 1.0;
      applyZoom();
    }}

    function applyZoom() {{
      document.getElementById('graphPipeline').style.transform = `scale(${{currentZoom}})`;
    }}

    // Initial render
    renderPipeline();
  </script>
</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"[OK] Generated Interactive HTML: {output_path}")


# =============================================================================
# Exporter: Graphviz DOT (Calling-Order Aligned Vertical Layout)
# =============================================================================

def generate_graphviz_dot(
    nodes: Dict[str, CallNode],
    output_path: str,
    detail_level: DetailLevel = DetailLevel.BOTH
):
    """Generates a structured Graphviz DOT file with calling-order aligned ranks."""
    dot_lines = [
        "digraph CodeSaturneCallGraph {",
        "  rankdir=TB;",
        "  compound=true;",
        '  ratio="auto";',
        '  size="8.5,14";',
        "  nodesep=0.4;",
        "  ranksep=0.45;",
        "  ordering=out;",
        '  fontname="Helvetica,Arial,sans-serif";',
        '  node [fontname="Helvetica,Arial,sans-serif", fontsize=10, shape=box, style="rounded,filled"];',
        '  edge [fontname="Helvetica,Arial,sans-serif", fontsize=8, color="#64748b"];',
        '  bgcolor="#f8fafc";',
        ""
    ]

    filtered_nodes = {}
    for nid, node in nodes.items():
        if detail_level == DetailLevel.SUMMARY and node.detail_level == DetailLevel.DETAILED:
            continue
        filtered_nodes[nid] = node

    stage_clusters = {
        Stage.SETUP: ("cluster_setup", "1. Setup Stage (cs_setup)", "#e0f2fe", "#0284c7"),
        Stage.PREPROCESS: ("cluster_preprocess", "2. Preprocessing Stage (cs_preprocess_mesh)", "#f1f5f9", "#475569"),
        Stage.INITIALIZATION: ("cluster_init", "3. Initialization Stage", "#ecfdf5", "#059669"),
        Stage.TIME_LOOP: ("cluster_timeloop", "4. Time Loop (cs_time_stepping & cs_solve_all)", "#fffbeb", "#d97706"),
        Stage.FINALIZATION: ("cluster_finalize", "5. Finalization Stage", "#fdf2f8", "#db2777"),
    }

    for stage, (cid, clabel, cbg, cborder) in stage_clusters.items():
        dot_lines.append(f"  subgraph {cid} {{")
        dot_lines.append(f'    label = "{clabel}";')
        dot_lines.append(f'    color = "{cborder}";')
        dot_lines.append(f'    bgcolor = "{cbg}";')
        dot_lines.append("    style = \"filled,rounded\";")
        dot_lines.append("    margin = 12;")
        dot_lines.append("")

        for nid, node in filtered_nodes.items():
            if node.stage != stage:
                continue

            fillcolor = "#ffffff"
            border_color = "#334155"
            shape = "box"
            penwidth = 1.0

            freq_badge = ""
            if node.frequency == CallFrequency.ONCE:
                freq_badge = "[1x]"
            elif node.frequency == CallFrequency.PER_TIME_STEP:
                freq_badge = "[1x / Step]"
            elif node.frequency == CallFrequency.MULTIPLE_PER_STEP:
                freq_badge = "[Nx / Step]"

            if node.category == NodeCategory.USER_HOOK:
                fillcolor = "#fef3c7"
                border_color = "#d97706"
                penwidth = 2.0
                label = f"<{user_dot_label(node, freq_badge)}>"
            elif node.is_iterative_solve:
                fillcolor = "#fce7f3"
                border_color = "#db2777"
                penwidth = 1.8
                label = f"<{iterative_dot_label(node, freq_badge)}>"
            else:
                label = f"<{standard_dot_label(node, freq_badge)}>"

            dot_lines.append(
                f'    "{node.node_id}" [label={label}, fillcolor="{fillcolor}", '
                f'color="{border_color}", shape={shape}, penwidth={penwidth}];'
            )

        dot_lines.append("  }")
        dot_lines.append("")

    # Primary Vertical Spine (Chronological calling order down the left column)
    dot_lines.append("  // Primary vertical spine (Chronological calling order)")
    spine_pairs = [
        ("main_entry", "cs_setup"),
        ("cs_setup", "cs_preprocess_mesh"),
        ("cs_preprocess_mesh", "init_stage"),
        ("init_stage", "cs_time_stepping"),
        ("cs_time_stepping", "finalization_core"),
    ]
    for src, dst in spine_pairs:
        if src in filtered_nodes and dst in filtered_nodes:
            dot_lines.append(f'  "{src}" -> "{dst}" [weight=10, penwidth=2.2, color="#0284c7"];')

    if detail_level == DetailLevel.DETAILED:
        timeloop_spine = [
            ("cs_time_stepping", "time_step_control"),
            ("time_step_control", "cs_solve_all"),
            ("cs_solve_all", "cs_physical_properties_update"),
            ("cs_physical_properties_update", "cs_head_losses_compute"),
            ("cs_head_losses_compute", "_solve_most"),
            ("_solve_most", "cs_boundary_conditions_set_coeffs"),
            ("cs_boundary_conditions_set_coeffs", "cs_solve_navier_stokes"),
            ("cs_solve_navier_stokes", "cs_velocity_prediction_1"),
            ("cs_velocity_prediction_1", "cs_pressure_correction"),
            ("cs_pressure_correction", "cs_mass_flux"),
            ("cs_mass_flux", "cs_velocity_prediction_2"),
            ("cs_velocity_prediction_2", "_solve_turbulence"),
            ("_solve_turbulence", "cs_rad_transfer_solve"),
            ("cs_rad_transfer_solve", "cs_solve_transported_variables"),
            ("cs_solve_transported_variables", "cs_lagr_solve_time_step"),
            ("cs_lagr_solve_time_step", "cs_time_moment_update_all"),
            ("cs_time_moment_update_all", "cs_time_stepping_write_checkpoint"),
            ("cs_time_stepping_write_checkpoint", "cs_post_default_write_variables"),
        ]
        for src, dst in timeloop_spine:
            if src in filtered_nodes and dst in filtered_nodes:
                dot_lines.append(f'  "{src}" -> "{dst}" [weight=8, penwidth=1.8, color="#0284c7"];')

    # Rank Alignment: Pair Core Callers and User Hooks side-by-side relative to calling order
    dot_lines.append("")
    dot_lines.append("  // Aligned ranks: Core caller on left, User hook on right")
    caller_hook_pairs = [
        ("main_entry", "cs_user_parallel_io"),
        ("cs_setup", "cs_user_model"),
        ("cs_setup", "cs_user_parameters"),
        ("cs_setup", "cs_user_linear_solvers"),
        ("cs_setup", "cs_user_finalize_setup"),
        ("cs_preprocess_mesh", "cs_user_mesh_input"),
        ("cs_preprocess_mesh", "cs_user_join"),
        ("cs_preprocess_mesh", "cs_user_periodicity"),
        ("cs_preprocess_mesh", "cs_user_partition"),
        ("init_stage", "cs_user_initialization"),
        ("init_stage", "cs_user_porosity"),
        ("init_stage", "cs_user_extra_operations_initialize"),
        ("cs_physical_properties_update", "cs_user_physical_properties"),
        ("cs_head_losses_compute", "cs_user_head_losses"),
        ("cs_boundary_conditions_set_coeffs", "cs_user_boundary_conditions"),
        ("cs_velocity_prediction_1", "cs_user_source_terms"),
        ("cs_rad_transfer_solve", "cs_user_rad_transfer_bcs"),
        ("cs_lagr_solve_time_step", "cs_user_lagr_bcs"),
        ("cs_time_stepping", "cs_user_extra_operations"),
        ("cs_post_default_write_variables", "cs_user_postprocess_values"),
        ("finalization_core", "cs_user_extra_operations_finalize"),
    ]

    for caller_id, hook_id in caller_hook_pairs:
        if caller_id in filtered_nodes and hook_id in filtered_nodes:
            dot_lines.append(f'  {{ rank=same; "{caller_id}"; "{hook_id}"; }}')
            dot_lines.append(f'  "{caller_id}" -> "{hook_id}" [color="#d97706", style=dashed, penwidth=1.6, label="calls", fontcolor="#b45309"];')

    dot_lines.append("}")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(dot_lines))
    print(f"[OK] Generated Graphviz DOT (Calling-Order Aligned): {output_path}")

    try:
        svg_path = os.path.splitext(output_path)[0] + ".svg"
        res = subprocess.run(["dot", "-Tsvg", output_path, "-o", svg_path], capture_output=True, text=True)
        if res.returncode == 0:
            print(f"[OK] Rendered SVG via Graphviz: {svg_path}")
    except Exception:
        pass


def user_dot_label(node: CallNode, freq_badge: str) -> str:
    return (
        f'<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1">'
        f'<TR><TD ALIGN="CENTER"><B><FONT COLOR="#b45309">[USER]</FONT> {node.name}</B></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="8" COLOR="#78350f">{node.label}</FONT></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="8" COLOR="#b45309">{freq_badge}</FONT></B></TD></TR>'
        f'</TABLE>'
    )


def iterative_dot_label(node: CallNode, freq_badge: str) -> str:
    return (
        f'<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1">'
        f'<TR><TD ALIGN="CENTER"><B><FONT COLOR="#be185d">[SOLVE]</FONT> {node.name}</B></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="8" COLOR="#831843">{node.label}</FONT></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="8" COLOR="#be185d">{freq_badge}</FONT></B></TD></TR>'
        f'</TABLE>'
    )


def standard_dot_label(node: CallNode, freq_badge: str) -> str:
    return (
        f'<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1">'
        f'<TR><TD ALIGN="CENTER"><B>{node.name}</B></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="8" COLOR="#64748b">{node.label}</FONT></TD></TR>'
        f'<TR><TD ALIGN="CENTER"><B><FONT POINT-SIZE="8" COLOR="#0284c7">{freq_badge}</FONT></B></TD></TR>'
        f'</TABLE>'
    )


# =============================================================================
# Exporter: Calling-Order Aligned Standalone Vertical SVG Generator
# =============================================================================

def generate_standalone_svg(
    nodes: Dict[str, CallNode],
    output_path: str,
    detail_level: DetailLevel = DetailLevel.SUMMARY
):
    """Generates a portrait-oriented standalone SVG with left and right columns
    strictly aligned according to the actual calling sequence.
    """
    canvas_width = 920
    margin_x = 40
    content_width = canvas_width - 2 * margin_x

    card_w = 360
    card_h = 56
    row_gap = 18
    stage_header_h = 42
    stage_padding = 18

    col1_x = margin_x + stage_padding
    col2_x = col1_x + card_w + 84  # 84px gap for horizontal calling connector

    # Define stage rows aligned chronologically: (core_node_id, [user_node_ids], call_label)
    if detail_level == DetailLevel.SUMMARY:
        stage_definitions = [
            (
                Stage.SETUP,
                "1. SETUP STAGE (cs_setup)",
                "#0284c7",
                "#e0f2fe",
                [
                    ("main_entry", ["cs_user_parallel_io"], "calls"),
                    ("cs_setup", ["cs_user_time_table", "cs_user_zones"], "calls"),
                    ("cs_setup", ["cs_user_model"], "calls"),
                    ("cs_setup", ["cs_user_parameters"], "calls"),
                    ("cs_setup", ["cs_user_linear_solvers"], "calls"),
                    ("cs_setup", ["cs_user_time_moments", "cs_user_couplings"], "calls"),
                    ("cs_setup", ["cs_user_finalize_setup", "cs_user_matrix_tuning"], "calls"),
                    ("cs_setup", ["cs_user_postprocess_setup"], "calls"),
                ]
            ),
            (
                Stage.PREPROCESS,
                "2. MESH PREPROCESSING (cs_preprocess_mesh)",
                "#475569",
                "#f1f5f9",
                [
                    ("cs_preprocess_mesh", ["cs_user_mesh_input"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_join"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_periodicity"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_mesh_modify"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_partition"], "calls"),
                    ("cs_mesh_quality", [], ""),
                ]
            ),
            (
                Stage.INITIALIZATION,
                "3. INITIALIZATION STAGE",
                "#059669",
                "#ecfdf5",
                [
                    ("cs_restart_read", [], ""),
                    ("init_stage", ["cs_user_initialization"], "calls"),
                    ("init_stage", ["cs_user_porosity"], "calls"),
                    ("init_stage", ["cs_user_extra_operations_initialize"], "calls"),
                    ("cs_lagr_solve_initialize", [], ""),
                ]
            ),
            (
                Stage.TIME_LOOP,
                "4. MAIN TIME LOOP (cs_time_stepping & cs_solve_all)",
                "#d97706",
                "#fffbeb",
                [
                    ("cs_time_stepping", [], ""),
                    ("cs_physical_properties_update", ["cs_user_physical_properties"], "calls"),
                    ("cs_head_losses_compute", ["cs_user_head_losses"], "calls"),
                    ("cs_boundary_conditions_set_coeffs", ["cs_user_boundary_conditions"], "calls"),
                    ("cs_velocity_prediction_1", ["cs_user_source_terms"], "calls"),
                    ("cs_pressure_correction", [], ""),
                    ("cs_mass_flux", [], ""),
                    ("cs_velocity_prediction_2", [], ""),
                    ("_solve_turbulence", ["cs_user_source_terms"], "calls"),
                    ("cs_rad_transfer_solve", ["cs_user_rad_transfer_bcs", "cs_user_rad_transfer_absorption"], "calls"),
                    ("cs_solve_transported_variables", ["cs_user_source_terms"], "calls"),
                    ("cs_lagr_solve_time_step", ["cs_user_lagr_bcs", "cs_user_lagr_particle"], "calls"),
                    ("cs_time_moment_update_all", ["cs_user_extra_operations"], "calls"),
                    ("cs_post_default_write_variables", ["cs_user_postprocess_values"], "calls"),
                    ("cs_time_stepping_write_checkpoint", [], ""),
                ]
            ),
            (
                Stage.FINALIZATION,
                "5. FINALIZATION STAGE",
                "#db2777",
                "#fdf2f8",
                [
                    ("finalization_core", ["cs_user_extra_operations_finalize"], "calls"),
                ]
            )
        ]
    else:  # DETAILED
        stage_definitions = [
            (
                Stage.SETUP,
                "1. SETUP STAGE (cs_setup)",
                "#0284c7",
                "#e0f2fe",
                [
                    ("main_entry", ["cs_user_parallel_io"], "calls"),
                    ("_init_setup", ["cs_user_time_table", "cs_user_zones"], "calls"),
                    ("_init_user", ["cs_user_model"], "calls"),
                    ("cs_setup", ["cs_user_parameters"], "calls"),
                    ("cs_setup", ["cs_user_linear_solvers"], "calls"),
                    ("cs_setup", ["cs_user_time_moments", "cs_user_couplings"], "calls"),
                    ("cs_setup", ["cs_user_finalize_setup", "cs_user_matrix_tuning"], "calls"),
                    ("cs_setup", ["cs_user_postprocess_setup"], "calls"),
                ]
            ),
            (
                Stage.PREPROCESS,
                "2. MESH PREPROCESSING (cs_preprocess_mesh)",
                "#475569",
                "#f1f5f9",
                [
                    ("cs_preprocess_mesh", ["cs_user_mesh_input"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_join"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_periodicity"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_mesh_modify"], "calls"),
                    ("cs_preprocess_mesh", ["cs_user_partition"], "calls"),
                    ("cs_mesh_quality", [], ""),
                ]
            ),
            (
                Stage.INITIALIZATION,
                "3. INITIALIZATION STAGE",
                "#059669",
                "#ecfdf5",
                [
                    ("cs_restart_read", [], ""),
                    ("cs_initialize_fields", ["cs_user_initialization"], "calls"),
                    ("cs_porous_model_init", ["cs_user_porosity"], "calls"),
                    ("init_stage", ["cs_user_extra_operations_initialize"], "calls"),
                    ("cs_lagr_solve_initialize", [], ""),
                ]
            ),
            (
                Stage.TIME_LOOP,
                "4. MAIN TIME LOOP (cs_time_stepping & cs_solve_all)",
                "#d97706",
                "#fffbeb",
                [
                    ("time_step_control", [], ""),
                    ("cs_physical_properties_update", ["cs_user_physical_properties", "cs_user_physical_properties_turb_viscosity"], "calls"),
                    ("cs_head_losses_compute", ["cs_user_head_losses"], "calls"),
                    ("cs_boundary_conditions_set_coeffs", ["cs_user_boundary_conditions"], "calls"),
                    ("cs_velocity_prediction_1", ["cs_user_source_terms"], "calls"),
                    ("cs_equation_iterative_solve_vector", [], ""),
                    ("cs_pressure_correction", [], ""),
                    ("cs_sles_solve_pressure", [], ""),
                    ("cs_mass_flux", [], ""),
                    ("cs_velocity_prediction_2", [], ""),
                    ("_solve_turbulence", ["cs_user_source_terms"], "calls"),
                    ("cs_turbulence_iterative_solve", [], ""),
                    ("cs_rad_transfer_solve", ["cs_user_rad_transfer_bcs", "cs_user_rad_transfer_absorption"], "calls"),
                    ("cs_solve_transported_variables", ["cs_user_source_terms"], "calls"),
                    ("cs_solve_equation_scalar", [], ""),
                    ("cs_lagr_solve_time_step", ["cs_user_lagr_bcs", "cs_user_lagr_particle"], "calls"),
                    ("cs_time_moment_update_all", ["cs_user_extra_operations"], "calls"),
                    ("cs_post_default_write_variables", ["cs_user_postprocess_values"], "calls"),
                    ("cs_time_stepping_write_checkpoint", [], ""),
                ]
            ),
            (
                Stage.FINALIZATION,
                "5. FINALIZATION STAGE",
                "#db2777",
                "#fdf2f8",
                [
                    ("finalization_core", ["cs_user_extra_operations_finalize"], "calls"),
                ]
            )
        ]

    # Pre-calculate Y coordinates and box heights
    current_y = 115
    rendered_stages = []

    for stage_enum, title, border_col, bg_col, rows in stage_definitions:
        stage_y = current_y
        row_y = stage_y + stage_header_h + stage_padding
        rendered_rows = []

        for core_id, user_ids, call_label in rows:
            num_users = max(len(user_ids), 1)
            row_height = num_users * card_h + (num_users - 1) * 10
            rendered_rows.append({
                "core_id": core_id,
                "user_ids": user_ids,
                "call_label": call_label,
                "y": row_y,
                "height": row_height
            })
            row_y += row_height + row_gap

        cluster_h = (row_y - stage_y) - row_gap + stage_padding
        rendered_stages.append({
            "stage": stage_enum,
            "title": title,
            "border_col": border_col,
            "bg_col": bg_col,
            "y": stage_y,
            "height": cluster_h,
            "rows": rendered_rows
        })
        current_y += cluster_h + 38

    total_canvas_height = current_y + 20

    # Build SVG XML
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {canvas_width} {total_canvas_height}" width="{canvas_width}" height="{total_canvas_height}">',
        '  <defs>',
        '    <marker id="spine-arrow" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">',
        '      <path d="M 0 1 L 10 5 L 0 9 z" fill="#0284c7" />',
        '    </marker>',
        '    <marker id="user-call-arrow" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">',
        '      <path d="M 0 1 L 10 5 L 0 9 z" fill="#d97706" />',
        '    </marker>',
        '  </defs>',
        f'  <rect width="{canvas_width}" height="{total_canvas_height}" fill="#f8fafc" />',
        '',
        '  <!-- Title Header -->',
        f'  <text x="{canvas_width / 2}" y="42" text-anchor="middle" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="22" font-weight="700" fill="#0f172a">code_saturne Execution Call Graph</text>',
        f'  <text x="{canvas_width / 2}" y="66" text-anchor="middle" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="12" fill="#64748b">Calling-Order Aligned Architecture &amp; User Call Sites ({detail_level.value})</text>',
        '',
        '  <!-- Legend Pills -->',
        f'  <g transform="translate({margin_x}, 82)">',
        '    <rect x="0" y="0" width="70" height="18" rx="4" fill="#e0f2fe" stroke="#0284c7" stroke-width="1"/>',
        '    <text x="35" y="13" text-anchor="middle" font-family="monospace" font-size="10" font-weight="700" fill="#0284c7">[1x] Once</text>',
        '    <rect x="80" y="0" width="105" height="18" rx="4" fill="#f0fdf4" stroke="#16a34a" stroke-width="1"/>',
        '    <text x="132" y="13" text-anchor="middle" font-family="monospace" font-size="10" font-weight="700" fill="#16a34a">[1x/Step] Per Step</text>',
        '    <rect x="195" y="0" width="115" height="18" rx="4" fill="#fff7ed" stroke="#ea580c" stroke-width="1"/>',
        '    <text x="252" y="13" text-anchor="middle" font-family="monospace" font-size="10" font-weight="700" fill="#ea580c">[Nx/Step] Multiple</text>',
        '    <rect x="320" y="0" width="90" height="18" rx="4" fill="#fef3c7" stroke="#d97706" stroke-width="1.5"/>',
        '    <text x="365" y="13" text-anchor="middle" font-family="monospace" font-size="10" font-weight="700" fill="#b45309">[USER] Hook</text>',
        '  </g>',
        ''
    ]

    for s_idx, stage in enumerate(rendered_stages):
        sy = stage["y"]
        sh = stage["height"]

        # Stage background container
        svg.append(f'  <!-- Stage: {stage["title"]} -->')
        svg.append(f'  <rect x="{margin_x}" y="{sy}" width="{content_width}" height="{sh}" rx="8" fill="{stage["bg_col"]}" stroke="{stage["border_col"]}" stroke-width="1.5" />')
        svg.append(f'  <text x="{margin_x + 16}" y="{sy + 24}" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="13" font-weight="700" fill="{stage["border_col"]}">{stage["title"]}</text>')
        svg.append(f'  <text x="{col1_x}" y="{sy + 42}" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="10" font-weight="700" fill="#475569">CORE SOLVER PIPELINE (CHRONOLOGICAL)</text>')
        svg.append(f'  <text x="{col2_x}" y="{sy + 42}" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="10" font-weight="700" fill="#b45309">ALIGNED USER-DEFINED CALL SITES (cs_user_*)</text>')

        # Render rows
        for r_idx, row in enumerate(stage["rows"]):
            ry = row["y"]
            core_id = row["core_id"]
            core_node = nodes.get(core_id)

            if core_node:
                svg.append(render_svg_node(core_node, col1_x, ry, card_w, card_h))

            # Render vertical spine arrow to next core row within this stage
            if r_idx < len(stage["rows"]) - 1:
                next_ry = stage["rows"][r_idx + 1]["y"]
                spine_x = col1_x + card_w // 2
                svg.append(f'  <line x1="{spine_x}" y1="{ry + card_h}" x2="{spine_x}" y2="{next_ry - 2}" stroke="#0284c7" stroke-width="2.0" marker-end="url(#spine-arrow)" />')

            # Render User hooks and calling arrows
            for u_idx, uid in enumerate(row["user_ids"]):
                u_node = nodes.get(uid)
                if u_node:
                    uy = ry + u_idx * (card_h + 10)
                    svg.append(render_svg_node(u_node, col2_x, uy, card_w, card_h))

                    # Horizontal calling connector arrow between core caller and user hook
                    arrow_y = uy + card_h // 2
                    svg.append(f'  <line x1="{col1_x + card_w}" y1="{arrow_y}" x2="{col2_x - 6}" y2="{arrow_y}" stroke="#d97706" stroke-width="1.8" stroke-dasharray="4,2" marker-end="url(#user-call-arrow)" />')
                    svg.append(f'  <rect x="{col1_x + card_w + 16}" y="{arrow_y - 8}" width="42" height="16" rx="3" fill="#fef3c7" stroke="#d97706" stroke-width="0.8"/>')
                    svg.append(f'  <text x="{col1_x + card_w + 37}" y="{arrow_y + 4}" text-anchor="middle" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="8" font-weight="700" fill="#b45309">calls ➔</text>')

        # Downward spine arrow to next stage
        if s_idx < len(rendered_stages) - 1:
            arrow_x = col1_x + card_w // 2
            arrow_y1 = sy + sh + 2
            arrow_y2 = sy + sh + 32
            svg.append(f'  <line x1="{arrow_x}" y1="{arrow_y1}" x2="{arrow_x}" y2="{arrow_y2}" stroke="#0284c7" stroke-width="2.5" marker-end="url(#spine-arrow)" />')

    svg.append('</svg>')

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(svg))
    print(f"[OK] Generated Calling-Order Aligned Standalone SVG: {output_path}")


def render_svg_node(node: CallNode, x: int, y: int, w: int, h: int) -> str:
    """Renders a single node card in SVG."""
    is_user = node.category == NodeCategory.USER_HOOK
    is_solve = node.is_iterative_solve

    bg_color = "#ffffff"
    border_color = "#cbd5e1"
    border_w = 1.0

    if is_user:
        bg_color = "#fffbeb"
        border_color = "#d97706"
        border_w = 1.8
    elif is_solve:
        bg_color = "#fdf2f8"
        border_color = "#db2777"
        border_w = 1.5

    # Frequency badge
    freq_text = "[1x]"
    freq_bg = "#e0f2fe"
    freq_fg = "#0284c7"
    if node.frequency == CallFrequency.PER_TIME_STEP:
        freq_text = "[1x/Step]"
        freq_bg = "#f0fdf4"
        freq_fg = "#16a34a"
    elif node.frequency == CallFrequency.MULTIPLE_PER_STEP:
        freq_text = "[Nx/Step]"
        freq_bg = "#fff7ed"
        freq_fg = "#ea580c"

    short_name = node.name if len(node.name) < 32 else node.name[:29] + "..."
    short_label = node.label if len(node.label) < 38 else node.label[:35] + "..."

    badge_x = x + w - 70
    badge_y = y + 8

    lines = [
        f'  <g id="node-{node.node_id}">',
        f'    <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{bg_color}" stroke="{border_color}" stroke-width="{border_w}" />',
        f'    <text x="{x + 10}" y="{y + 20}" font-family="monospace" font-size="11" font-weight="700" fill="#0f172a">{short_name}</text>',
        f'    <text x="{x + 10}" y="{y + 36}" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif" font-size="9" fill="#64748b">{short_label}</text>',
        f'    <text x="{x + 10}" y="{y + 48}" font-family="monospace" font-size="7.5" fill="#94a3b8">{node.file_ref}</text>',
        f'    <rect x="{badge_x}" y="{badge_y}" width="60" height="15" rx="3" fill="{freq_bg}" stroke="{freq_fg}" stroke-width="0.7"/>',
        f'    <text x="{badge_x + 30}" y="{badge_y + 11}" text-anchor="middle" font-family="monospace" font-size="8" font-weight="700" fill="{freq_fg}">{freq_text}</text>',
        '  </g>'
    ]
    return "\n".join(lines)


# =============================================================================
# Exporter: Markdown / Mermaid
# =============================================================================

def generate_markdown_mermaid(nodes: Dict[str, CallNode], output_path: str):
    """Generates documentation in Markdown format with structured Mermaid diagrams."""
    lines = [
        "# code_saturne Architectural Call Graph & Execution Workflow",
        "",
        "This document details the main execution stages of **code_saturne**, "
        "highlighting major numerical steps and all user-defined function call sites.",
        "",
        "### Call Frequency Nomenclature",
        "- 🔵 **`[1x]` (Once Only)**: Setup, mesh preprocessing, initial conditions, finalization.",
        "- 🟢 **`[1x / Step]` (Once per Time Step)**: Physical properties update, Lagrangian tracking, "
        "temporal statistics, post-processing writers.",
        "- 🟠 **`[Nx / Step]` (Multiple per Time Step)**: Inner velocity-pressure iterations ($iterns \\le nterup$), "
        "boundary condition evaluations, per-equation volumetric source terms (`cs_user_source_terms`).",
        "",
        "---",
        "",
        "## 1. High-Level Execution Pipeline",
        "",
        "```mermaid",
        "graph TD",
        "    classDef user fill:#fef3c7,stroke:#d97706,stroke-width:2px;",
        "    classDef core fill:#e2e8f0,stroke:#475569,stroke-width:1px;",
        "    classDef solve fill:#fce7f3,stroke:#db2777,stroke-width:2px;",
        "",
        "    subgraph S1 [1. Setup Stage]",
        "        MAIN[main / _run] --> SETUP[cs_setup]",
        "        SETUP --> U_MOD[cs_user_model <br/> 🔵 <b>[1x]</b>]:::user",
        "        SETUP --> U_PAR[cs_user_parameters <br/> 🔵 <b>[1x]</b>]:::user",
        "        SETUP --> U_SLES[cs_user_linear_solvers <br/> 🔵 <b>[1x]</b>]:::user",
        "    end",
        "",
        "    subgraph S2 [2. Preprocessing Stage]",
        "        MAIN --> PREP[cs_preprocess_mesh]",
        "        PREP --> U_MESH[cs_user_mesh_input <br/> 🔵 <b>[1x]</b>]:::user",
        "        PREP --> U_JOIN[cs_user_join <br/> 🔵 <b>[1x]</b>]:::user",
        "        PREP --> QUAL[cs_mesh_quality <br/> 🔵 <b>[1x]</b>]",
        "    end",
        "",
        "    subgraph S3 [3. Initialization Stage]",
        "        MAIN --> INIT[Initialization Stage]",
        "        INIT --> RESTART[Read Checkpoint <br/> 🔵 <b>[1x]</b>]",
        "        INIT --> U_INIT[cs_user_initialization <br/> 🔵 <b>[1x]</b>]:::user",
        "        INIT --> U_PORO[cs_user_porosity <br/> 🔵 <b>[1x]</b>]:::user",
        "    end",
        "",
        "    subgraph S4 [4. Time Loop: cs_time_stepping]",
        "        INIT --> TLOOP[Time Loop: do ... while nt_cur < nt_max]",
        "        TLOOP --> SOLVE_ALL[Flow & Scalars: cs_solve_all]",
        "        SOLVE_ALL --> PROP[cs_user_physical_properties <br/> 🟢 <b>[1x / Step]</b>]:::user",
        "        SOLVE_ALL --> NS[_solve_most / Navier-Stokes]",
        "        NS --> BC[cs_user_boundary_conditions <br/> 🟠 <b>[Nx / Step]</b>]:::user",
        "        NS --> VPRED[Velocity Prediction: iappel=1 <br/> 🟠 <b>[Nx / Step]</b>]",
        "        VPRED --> ST_U[cs_user_source_terms <br/> 🟠 <b>[Nx / Step]</b>]:::user",
        "        NS --> PCORR[Pressure Correction / Poisson <br/> 🟠 <b>[Nx / Step]</b>]:::solve",
        "        NS --> MFLUX[Mass Flux Update: cs_mass_flux <br/> 🟠 <b>[Nx / Step]</b>]",
        "        NS --> VCORR[Velocity Correction: iappel=2 <br/> 🟠 <b>[Nx / Step]</b>]",
        "        SOLVE_ALL --> TURB[Turbulence: _solve_turbulence <br/> 🟢 <b>[1x / Step]</b>]",
        "        TURB --> ST_TURB[cs_user_source_terms <br/> 🟠 <b>[Nx / Step]</b>]:::user",
        "        SOLVE_ALL --> SCAL[Transported Scalars <br/> 🟢 <b>[1x / Step]</b>]",
        "        SCAL --> ST_SCAL[cs_user_source_terms <br/> 🟠 <b>[Nx / Step]</b>]:::user",
        "        TLOOP --> LAGR[Lagrangian: cs_lagr_solve_time_step <br/> 🟢 <b>[1x / Step]</b>]",
        "        TLOOP --> U_EXTRA[cs_user_extra_operations <br/> 🟢 <b>[1x / Step]</b>]:::user",
        "        TLOOP --> POST[Post-processing & Checkpoint <br/> 🟢 <b>[1x / Step]</b>]",
        "    end",
        "",
        "    subgraph S5 [5. Finalization Stage]",
        "        TLOOP --> FIN[Finalization & Teardown]",
        "        FIN --> U_FIN[cs_user_extra_operations_finalize <br/> 🔵 <b>[1x]</b>]:::user",
        "    end",
        "```",
        "",
        "---",
        "",
        "## 2. Inventory of User Call Sites (`cs_user_*`)",
        "",
        "| Function Name | Stage | Frequency | Purpose | Source File |",
        "| :--- | :--- | :--- | :--- | :--- |"
    ]

    for node in nodes.values():
        if node.category == NodeCategory.USER_HOOK:
            freq_str = (
                "🔵 `[1x]`" if node.frequency == CallFrequency.ONCE else
                "🟢 `[1x / Step]`" if node.frequency == CallFrequency.PER_TIME_STEP else
                "🟠 `[Nx / Step]`"
            )
            lines.append(
                f"| `{node.name}` | `{node.stage.value}` | {freq_str} | {node.description} | `{node.file_ref}` |"
            )

    lines.append("")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[OK] Generated Markdown & Mermaid documentation: {output_path}")


# =============================================================================
# CLI Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate code_saturne Call Graph in HTML, Graphviz DOT, and Markdown."
    )
    parser.add_argument(
        "--output-html",
        default="cs_callgraph.html",
        help="Path for interactive HTML output (default: cs_callgraph.html)"
    )
    parser.add_argument(
        "--output-dot-summary",
        default="doc/cs_callgraph_summary.dot",
        help="Path for summary Graphviz DOT file (default: doc/cs_callgraph_summary.dot)"
    )
    parser.add_argument(
        "--output-dot-detailed",
        default="doc/cs_callgraph_detailed.dot",
        help="Path for detailed Graphviz DOT file (default: doc/cs_callgraph_detailed.dot)"
    )
    parser.add_argument(
        "--output-svg-summary",
        default="doc/cs_callgraph_summary.svg",
        help="Path for summary vertical SVG file (default: doc/cs_callgraph_summary.svg)"
    )
    parser.add_argument(
        "--output-svg-detailed",
        default="doc/cs_callgraph_detailed.svg",
        help="Path for detailed vertical SVG file (default: doc/cs_callgraph_detailed.svg)"
    )
    parser.add_argument(
        "--output-md",
        default="doc/callgraph.md",
        help="Path for Markdown / Mermaid documentation (default: doc/callgraph.md)"
    )

    args = parser.parse_args()

    # Ensure output directories exist
    all_paths = [
        args.output_html, args.output_dot_summary, args.output_dot_detailed,
        args.output_svg_summary, args.output_svg_detailed, args.output_md
    ]
    for path in all_paths:
        parent = os.path.dirname(path)
        if parent and not os.path.exists(parent):
            os.makedirs(parent, exist_ok=True)

    print("Building code_saturne execution model...")
    nodes = build_saturne_call_graph()
    print(f"Loaded {len(nodes)} call graph nodes.")

    # 1. Interactive HTML
    generate_interactive_html(nodes, args.output_html)

    # 2. Graphviz DOT files (Vertical Layout)
    generate_graphviz_dot(nodes, args.output_dot_summary, detail_level=DetailLevel.SUMMARY)
    generate_graphviz_dot(nodes, args.output_dot_detailed, detail_level=DetailLevel.DETAILED)

    # 3. Standalone Portrait-Oriented Vector SVG
    generate_standalone_svg(nodes, args.output_svg_summary, detail_level=DetailLevel.SUMMARY)
    generate_standalone_svg(nodes, args.output_svg_detailed, detail_level=DetailLevel.DETAILED)

    # 4. Markdown / Mermaid
    generate_markdown_mermaid(nodes, args.output_md)

    print("\nAll outputs generated successfully!")


if __name__ == "__main__":
    main()
