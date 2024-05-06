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

\page cs_ug_cdo_navsto Using CDO schemes for solving the Navier-Stokes system

[TOC]

<!--
    References used in this page
-->

[Mila20]: https://hal.science/tel-03080530v2
[MiBoE22]: https://hal.science/hal-03215118v1

Introduction {#sec_ug_cdo_navsto_intro}
============

The Navier--Stokes equations can be solved using **CDO face-based**
discretizations. Other space discretizations are not available. Up to
now, there is **no turbulence modelling** when using CDO schemes. This
is a work in progress.

The rationale to set up the computation in the case of the Navier-Stokes
equations is the very near of the one explained for user-defined equations
[here](@ref cs_ug_cdo_hho_base). One only focuses on the specifities arising
from the Navier-Stokes case.


Settings done in cs_user_model()
===================

The activation of the NavSto module is done thanks to the function \ref cs_navsto_system_activate
For instance, here are two examples to activate and set the main parameters for the
NavSto module

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_activate

or

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_activate2

The first parameter is the structure managing the domain boundaries. An example
of settings for the domain boundaries is described [here](@ref ug_cdo_sec_domain_boundaries).

The second parameter specifies the type of model to
consider among the following choice:
    - \ref CS_NAVSTO_MODEL_STOKES,
    - \ref CS_NAVSTO_MODEL_OSEEN,
    - \ref CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,

The *base* model can be updated thanks to the third parameter which is a flag
built from the following elemental bit (\ref cs_navsto_param_model_bit_t):
    - \ref CS_NAVSTO_MODEL_STEADY to specify the model is steady-state (by
      default, this is not the case).
    - \ref CS_NAVSTO_MODEL_GRAVITY_EFFECTS
    - \ref CS_NAVSTO_MODEL_BOUSSINESQ

The fourth parameter specifies which type of velocity-pressure algorithm will be
used (\ref cs_navsto_param_coupling_t). The choice is done among:
    - \ref CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY (cf. \cite MiBoE22 and
      \cite Mila20),
    - \ref CS_NAVSTO_COUPLING_MONOLITHIC
    - \ref CS_NAVSTO_COUPLING_PROJECTION (work in progress),

The last parameter specifies predefined post-processing operations. As the third
parameter, this a flag built from the following elemental bit
(\ref cs_navsto_param_post_bit_t):
    - \ref CS_NAVSTO_POST_VELOCITY_DIVERGENCE
    - \ref CS_NAVSTO_POST_KINETIC_ENERGY
    - \ref CS_NAVSTO_POST_VORTICITY
    - \ref CS_NAVSTO_POST_VELOCITY_GRADIENT
    - \ref CS_NAVSTO_POST_STREAM_FUNCTION (This adds an equation named
      *streamfunction_eq* and its variable field named *stream_function*. This
      equation relies on a CDO vertex-based equation.)
    - \ref CS_NAVSTO_POST_HELICITY
    - \ref CS_NAVSTO_POST_ENSTROPHY
    - \ref CS_NAVSTO_POST_MASS_DENSITY
    - \ref CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE
    - \ref CS_NAVSTO_POST_PRESSURE_GRADIENT




Predefined equations associated to the Navier--Stokes equations are
    - the *momentum* equation is automatically added when activated with a
    **monolithic** or **artificial compressibility** velocity-pressure coupling



In all cases, a vector-valued field named _velocity_ and a scalar-valued field
named _"pressure"_ are created. Moreover, several properties are added:
    - the property _mass_density_ (a \ref cs_property_t structure);
    - the property _laminar viscosity_ (a \ref cs_property_t structure);
    - the properties _turbulent_viscosity_ and the _total_viscosity_ are added
      if the model of turbulence is different from the laminar one (cf. \ref cs_turb_model_t);
    - the property _graddiv_coef_ (a \ref cs_property_t structure) when the
**artificial compressibility** is set;

along with the advection field _mass_flux_ (a \ref cs_adv_field_t structure)


Settings done in cs_user_finalize_setup()
===================

Boundary conditions
-------------------

In the case of the NavSto module, this is done as follows (One does not access
to the equation directly since it depends on the way the velocity/pressure
coupling is done).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_bc1

where the function `_vel_def` is defined as follows

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_vel_function

Definition of source terms
---------------------

Since the equation on which this source term applies, depends on the
choice of the velocity-pressure algorithm, the way to proceed varies
slightly of the way used on a user-defined equation.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st1

where the function `_src_def` is defined as follows

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_st_function


Settings done in cs_user_parameters() {#cs_ug_cdo_navsto_set_param}
===================

The rationale is similar to the one detailed in [this section](@ref cs_ug_cdo_hho_base_user_param).

Here are listed some specificities related to the NavSto module.

Set the strategy to solve the Navier-Stokes system when a monolithic coupling is used
-------------------

When a *monolithic* velocity-pressure is set, the linear system to solve is a
_saddle-point_ problem. This class of linear systems needs specific choices of
preconditioner/solver. The default settings is not always the optimal choice in
terms of efficiency. Some examples of settings involving differents strategies
of resolution are presented hereafter.

### Augmented Lagrangian Uzawa algorithm (ALU)

Here is another example settings a strategy to solve a saddle-point problem
arising from a Navier-Stokes equation with a monolithic velocity-pressure
coupling. One assumes that the external libraries have been installed and have
been configured with code_saturne (see the installation guide for more
details).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_alu

Since the velocity block is _augmented_ by the term
\f$\underline{\mathsf{grad}}\,\mathsf{div}\f$ (with a scaling coefficient given
by `CS_EQKEY_SADDLE_AUGMENT_SCALING`), the resulting linear system is hard to
solve for an iterative method (or one needs a very specific preconditionner for
\f$H(\mathsf{div})\f$).  Our strategy is to consider the sparse direct solver
MUMPS \cite MUMPS01. Here is an example of an advanced usage of the MUMPS
settings.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_mumps

The structure storing all the parameters related to the way to solve a **SLES**
(= _Sparse Linear Equation Solver_) is \ref cs_param_sles_t One can easily
retrieve this structure when one has a pointer to a \ref cs_equation_param_t
calling

```c
  cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEqName");

  cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(eqp);
```

This will be useful for more advanced settings.
The function \ref cs_param_sles_mumps allows one to set the two
following main parameters:
1. The factorization is performed using a *single* (parameter =
   `"true"`) ou a *double* precision (parameter = `"false"`). If the
   usage of MUMPS is not as a preconditioner, then one advises to
   consider a *double*-precision factorization.
2. Which type of factorization to perform ? There are the following choices:
   + \ref CS_PARAM_MUMPS_FACTO_LU : the default choice.
   + \ref CS_PARAM_MUMPS_FACTO_LDLT_SYM : Only for symmetric matrices
   + \ref CS_PARAM_MUMPS_FACTO_LDLT_SPD : Cholesky factorization only for
     Symmetric Positive Definite matrices

In addition, one can set more advanced options calling \ref cs_param_sles_mumps
to reach a better performance but this is case dependent. Since the velocity
block is vector-valued, one can benefit from a block analysis (the third
parameter is set to 3). According to your MUMPS installation, we can have
access to external libraries such as (PT-)Scotch ou (Par)Metis to speed-up the
_analysis_ step.

Here is a second example.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_alu_mumps



### Golub-Kahan Bidiagonalization algorithm (GKB)

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_gkb_kcycle

The linear system may be augmented to improve the convergence rate of the
algorithm (but the system is harder to solve). Here is another example:

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_gkb_mumps



### Block preconditioners with a Krylov solver

The two in-house Krylov solvers available with block preconditioning are
- **Minres** (Minimal residual) algorithm for symmetric indefinite systems such
  as systems encountered in the case of Stokes systems
- **GCR** (Generalized conjugate residual) algorithm for general indefinite
  systems. This is a flexible solver (preconditioner may vary between two
  iterations).

These two algorithms are optimized to handle saddle-point problems in
code_saturne since the (1,2) and (2,1) which are transposed is stored only
once. Moreover, this block is stored in an unassembled way.

Other block preconditioners with a Krylov solver can be set using the external
librairy [PETSc](https://petsc.org)

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_minres



### Uzawa algorithm with a CG acceleration

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_uzacg


To go further
=============

The detailed description of CDO schemes, the mathmatical analysis and numerical
results on benchmarks are available in the following PhD thesis:

* [**PhD thesis**: "Compatible Discrete Operator schemes for the unsteady incompessible Navier-Stokes equations"][Mila20] \cite Mila20


Additional publications :
* [**Article**: "Artificial compressibility methods for the incompressible Navier-Stokes equations using lowest-order face-based schemes on polytopal meshes"][MiBoE22] \cite MiBoE22
