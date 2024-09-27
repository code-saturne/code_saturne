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

\page cs_ug_cdo_sles How to set linear solvers with CDO/HHO schemes

[TOC]

# Introduction

In code_saturne, a **SLES** means a **S** parse **L** inear **E** quation
**S** olver.  When using CDO/HHO schemes, the settings related to a
SLES is stored in the structure \ref cs_param_sles_t

To access more advanced settings, it is possible to retrieve this
structure when one has a pointer to a \ref cs_equation_param_t
Simply write:
```c
  cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEqName");

  cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(eqp);
```

To get a full access to the low-level structures handling the linear
algebra in code_saturne, please consider to edit the function \ref
cs_user_linear_solvers For a Finite Volume scheme, this is the main
way to proceed.

# Solve a linear system {#cs_ug_cdo_sles_set}

Several examples are detailed hereafter to modify the default settings
related to a linear solver.

## A first example

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_user_simple

- More details about the available solvers are gathered [in this section](@ref cs_ug_cdo_sles_solver)
- More details about the available preconditioners are gathered [in this section](@ref cs_ug_cdo_sles_precond)

## A second example

For a symmetric positive definite (SPD) system, the solver of choice is the
conjugate gradient (\ref CS_PARAM_SOLVER_CG) or its flexible variant (\ref
CS_PARAM_SOLVER_FCG). In this case, a good preconditioner is a multilevel
algorithm such as an algebraic multigrid (AMG). Different multigrid techniques
are available when using code_saturne depending on the installation
configuration (in-house implementations or algorithms available from PETSc
and/or HYPRE libraries). Please refer to [this section](@ref cs_ug_cdo_sles_amg),
for more details about AMG techniques.

Here is a second example:
\snippet cs_user_parameters-cdo-condif.c param_cdo_sles_settings1

If the external library [PETSc](https://petsc.org) is available, one uses its
algebraic multigrid, otherwise one uses the in-house solver. Here a flexible
conjugate gradient (`fcg`) with a diagonal preconditionning (`jacobi`).


## Set the family of linear solvers

By default, *in-house* solvers are considered. With some choices of
solver or preconditioner, one has to use an external library. The
external libraries which can be linked to code_saturne are:
- [HYPRE](https://hypre.readthedocs.io)
- [MUMPS](https://mumps-solver.org), a robust sparse direct solver \cite MUMPS01
- [PETSc](https://petsc.org)

\ref CS_EQKEY_SOLVER_FAMILY allows one to specify the family of linear solvers
to consider. This can be useful if an automatic switch is not done or when
there are several possibilities. For instance,
- using HYPRE solver either directly from the HYPRE library or through
the PETSc library;
- using the MUMPS solver either directly from the MUMPS library or
through the PETSc library.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_solver_family


## Set the linear solver {#cs_ug_cdo_sles_solver}

Available linear solvers are listed in \ref cs_param_solver_type_t and
can be set using the key \ref CS_EQKEY_SOLVER and the key values
gathered in the following table:

key value | description | associated type | available family
:--- | :--- | :--- | :---:
`"amg"` | Algebraic MultiGrid (AMG) iterative solver. This is a good choice to achieve a scalable solver related to symmetric positive definite system when used with the default settings. Usually, a better choice is to use it as a preconditioner of a Krylov solver (for instance a conjugate gradient). For convection-dominated systems, one needs to adapt the smoothers, the coarsening algorithm or the coarse solver. | \ref CS_PARAM_SOLVER_AMG |  saturne, PETSc, HYPRE
`"bicgs"`, `"bicgstab"` | stabilized BiCG algorithm (for non-symmetric linear systems). It may lead to breakdown. | \ref CS_PARAM_SOLVER_BICGS | saturne, PETSc
`"bicgstab2"` | variant of the stabilized BiCG algorithm (for non-symmetric linear systems). It may lead to breakdown. It should be more robust than the BiCGstab algorithm | \ref CS_PARAM_SOLVER_BICGS2 | saturne, PETSc
`"cg"`| Conjuguate gradient algorithm for symmetric positive definite (SPD) linear systems | \ref CS_PARAM_SOLVER_CG  | saturne, PETSc, HYPRE
`"cr3"`| 3-layer conjuguate residual algorithm for non-symmetric linear systems | \ref CS_PARAM_SOLVER_CR3  | saturne
`"fcg"` | Flexible version of the conjuguate gradient algorithm. This is useful when the preconditioner changes between two successive iterations | \ref CS_PARAM_SOLVER_FCG | saturne, PETSc
`"fgmres"` | Flexible Generalized Minimal Residual: robust and flexible iterative solver. More efficient solver can be found on easy-to-solve systems. Additional settings related to the number of directions stored before restarting the algorithm (cf. the key \ref CS_EQKEY_SOLVER_RESTART) | \ref CS_PARAM_SOLVER_FGMRES | PETSc
`"gauss_seidel"`, `"gs"` | Gauss-Seidel algorithm (Jacobi on interfaces between MPI ranks in case of parallel computing). No preconditioner. | \ref CS_PARAM_SOLVER_GAUSS_SEIDEL | saturne, PETSc, HYPRE
`"gcr"` | Generalized Conjugate Residual: **robust** and flexible iterative solver. More efficient solver can be found on easy-to-solve systems. This is close to a FGMRES algorithm. Additional settings related to the number of directions stored before restarting the algorithm (cf. the key \ref CS_EQKEY_SOLVER_RESTART). This is the **default solver for a user-defined equation**. | \ref CS_PARAM_SOLVER_GCR | saturne, PETSc
`"gmres"` | Generalized Conjugate Residual: **robust** iterative solver. More efficient solver can be found on easy-to-solve systems. Additional settings related to the number of directions stored before restarting the algorithm (cf. the key \ref CS_EQKEY_SOLVER_RESTART) | \ref CS_PARAM_SOLVER_GMRES | saturne, PETSc, HYPRE
`"jacobi"`, `"diag"`, `"diagonal"` | Jacobi algorithm = reciprocal of the diagonal values (**simplest iterative solver**). This is an efficient solver if the linear system is dominated by the diagonal. No preconditioner. | \ref CS_PARAM_SOLVER_JACOBI | saturne, PETSc, HYPRE
`"mumps"` | sparse direct solver from the MUMPS library (cf. \cite MUMPS01). This is the most robust solver but it may require a huge memory (especially for 3D problems). By default, it performs a LU factorization on general matrices using double precision. More options are available using the function \ref cs_param_sles_mumps and the function \ref cs_param_sles_mumps_advanced. More options are available through the user-defined function \ref cs_user_sles_mumps_hook | \ref CS_PARAM_SOLVER_MUMPS | mumps, PETSc (according to the configuration)
`"sym_gauss_seidel"`, `"sgs"` | Symmetric Gauss-Seidel algorithm. Jacobi on interfaces between MPI ranks in case of parallel computing). No preconditioner. | \ref CS_PARAM_SOLVER_SYM_GAUSS_SEIDEL | saturne, PETSc
`"user"` | User-defined iterative solver (rely on the function \ref cs_user_sles_it_solver) | \ref CS_PARAM_SOLVER_USER_DEFINED | saturne


## Set the preconditioner {#cs_ug_cdo_sles_precond}

Available preconditioners are listed in \ref cs_param_precond_type_t
and can be set using the key \ref CS_EQKEY_PRECOND and the key values
gathered in the following table:

key value | description | type | family
:--- | :--- | :--- | :---:
`"amg"` | Algebraic MultiGrid (AMG) iterative solver. See [this section](@ref cs_ug_cdo_sles_amg) for more details. | \ref CS_PARAM_PRECOND_AMG |  saturne, PETSc, HYPRE
`"block_jacobi"`, `"bjacobi"` | Block Jacobi with an ILU(0) factorization in each block. By default, there is one block per MPI rank. | \ref CS_PARAM_PRECOND_BJACOB_ILU0 | PETSc
`"bjacobi_sgs"`, `"bjacobi_ssor"` | Block Jacobi with a SSOR algorithm in each block. By default, there is one block per MPI rank. For a sequential run, this is the SSOr algorithm. | \ref CS_PARAM_PRECOND_BJACOB_SGS | PETSc
`"diag"`, `"jacobi"` | Diagonal preconditioning | \ref CS_PARAM_PRECOND_DIAG | saturne, PETSc, HYPRE
`"ilu0"` | Incomplete LU factorization (zero fill-in meaning that the same sparsity level is used). | \ref CS_PARAM_PRECOND_ILU0 | PETSc, HYPRE
`"icc0"` | Incomplete Choleski factorization (zero fill-in meaning that the same sparsity level is used). For SPD matrices. | \ref CS_PARAM_PRECOND_ICC0 | PETSc, HYPRE
`"mumps"` | the sparse direct solver MUMPS used as preconditioner. Please refer to [the MUMPS section](@ref cs_ug_cdo_sles_mumps) for more details. | \ref CS_PARAM_PRECOND_MUMPS | mumps
`"none"` | No preconditioner. | \ref CS_PARAM_PRECOND_NONE | saturne, PETSc, HYPRE
`"poly1"`, `"poly_1"` | 1st order Neumann polynomial preconditioning | \ref CS_PARAM_PRECOND_POLY1 | saturne
`"poly2"`, `"poly_2"` | 2nd order Neumann polynomial preconditioning | \ref CS_PARAM_PRECOND_POLY2 | saturne
`"ssor"` | Symmetric Successive OverRelaxation (SSOR) algorithm. In case of parallel computations, each MPI rank performs a SSOR on the local matrix. | \ref CS_PARAM_PRECOND_SSOR | PETSc




## Algebraic multigrids {#cs_ug_cdo_sles_amg}

Algebraic multigrids (AMG) can be used either as a solver or as a
preconditioner. According to this choice, the settings of the main components
should be different. The main options related to an AMG algorithm are:
- The choice of the cycle (`V`, `W`, `K` for instance)
- The choice of the down and/or up smoothers (type and number of sweeps)
- The choice of the coarse solver
- The choice of the coarsening algorithm

AMG algorithm is a good choice to achieve a scalable solver. The default
settings are related to a symmetric positive definite (SPD) system. Usually, a
better efficiency is reached when the AMG is used as a preconditioner of a
Krylov solver rather than as a solver. For convection-dominated systems, one
needs to adapt the main ingredients.

Several multigrid algorithms are available in code_saturne. The default one is
the in-house `V-cycle`. Use the key \ref CS_EQKEY_AMG_TYPE to change the type of
multigrid. Possible values for this key are gathered in the next table.

key value | description
:--- | :---
`"v_cycle"` | This is the default choice. This corresponds to an in-house V-cycle detailed in \cite MeFoH09 The main ingredients can set using the functions \ref cs_param_sles_amg_inhouse and \ref cs_param_sles_amg_inhouse_advanced Please refer to [this section](@ref cs_ug_cdo_sles_amg_inhouse) for more details.
`"k_cycle"` or `"kamg"` | Switch to a K-cycle strategy. This type of cycle has been detailed in \cite Notay05 Be aware that a lot of work can be done in the coarsest levels yielding more communications between MPI ranks during a parallel computation. The main ingredients can set using the functions \ref cs_param_sles_amg_inhouse and \ref cs_param_sles_amg_inhouse_advanced Please refer to [this section](@ref cs_ug_cdo_sles_amg_inhouse) for more details.
`"boomer"`, `"bamg"` or "boomer_v" | Switch to use the V-cycle of the BoomerAMG algorithm from the [HYPRE library](https://hypre.readthedocs.io)  Please refer to [this section](@ref cs_ug_cdo_sles_amg_boomer) for more details on the way to set the main ingredients.
`"boomer_w"` or `"bamg_w"` | Switch to a W-cycle and the BoomerAMG algorithm from the [HYPRE library](https://hypre.readthedocs.io). Please refer to [this section](@ref cs_ug_cdo_sles_amg_boomer) for more details on the way to set the main ingredients.
`"gamg"` or `"gamg_v"` | Switch to a V-cycle and the GAMG algorithm from the [PETSc](https://petsc.org) library. GAMG means "Geometric Algebraic multigrid. Please refer to [this section](@ref cs_ug_cdo_sles_amg_gamg) for more details on the way to set the main ingredients.
`"gamg_w"` | Switch to a W-cycle and the GAMG algorithm from the [PETSc](https://petsc.org) library. GAMG means "Geometric Algebraic multigrid. Please refer to [this section](@ref cs_ug_cdo_sles_amg_gamg) for more details on the way to set the main ingredients.


### In-house multigrids {#cs_ug_cdo_sles_amg_inhouse}

The in-house multigrid algorithms can be easily tuned thanks to the function
\ref cs_param_sles_amg_inhouse More advanced settings can also be added using
the function \ref cs_param_sles_amg_inhouse_advanced

\snippet cs_user_parameters-linear_solvers.c sles_kamg_momentum

In the case of an advanced settings, it is possible to keep the default
settings associated to a parameter using the value \ref CS_CDO_KEEP_DEFAULT

It is also possible to access directly the structures and set the
parameters in the function \ref cs_user_linear_solvers . Here is an
example of a such usage:

\snippet cs_user_parameters-linear_solvers.c sles_mgp_1

or

\snippet cs_user_parameters-linear_solvers.c sles_mgp_2

In order to get the best efficiency in case of HPC computations, parallel grid
merging can be optimized as follows:

\snippet cs_user_parameters-linear_solvers.c sles_mg_parall


### BoomerAMG {#cs_ug_cdo_sles_amg_boomer}

The main ingredients of the BoomerAMG algorithm are defined thanks to the
functions \ref cs_param_sles_boomeramg and \ref cs_param_sles_boomeramg_advanced
for additional parameters.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_boomer

We strongly invite users to read the HYPRE documentation and especially the
part dedicated to BoomerAMG to understand the different choices
[BoomerAMG documentation](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html)

### GAMG {#cs_ug_cdo_sles_amg_gamg}

The main ingredients of the GAMG algorithm are defined thanks to the
functions \ref cs_param_sles_gamg and \ref cs_param_sles_gamg_advanced
for additional parameters.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_gamg

We strongly invite users to read the PETSc documentation dedicated to algebraic multigrid to understand the different choices
[GAMG documentation](https://petsc.org/release/manual/ksp/#algebraic-multigrid-amg-preconditioners)

#### Alternative way (low-level API)

HYPRE setting functions can be called thanks to a lower-level API (relying on a
setup hook function). Here is an example of a such usage

In the function \ref cs_user_linear_solvers,

\snippet cs_user_parameters-linear_solvers.c sles_hypre_3

where the function `_hypre_p_setup_hook` is defined as follows:

\snippet cs_user_parameters-linear_solvers.c sles_hypre_hook_1


### GAMG {#cs_ug_cdo_sles_amg_gamg}

Up to now, the advanced settings of GAMG is possible inside \ref cs_user_linear_solvers
Here is an example of a such usage

\snippet cs_user_parameters-linear_solvers.c sles_petsc_gamg_1


## Sparse direct solver: MUMPS {#cs_ug_cdo_sles_mumps}

For the most difficult to solve linear systems or in order to check a
numerical scheme, it may be useful to solve the linear systems with a
direct solver. In code_saturne, this is possible using the MUMPS
library.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_user_mumps

### Additional options

Additional options may be set to improve the performance of the sparse
direct solver. Please refer to the MUMPS documentation to know more
about the different options.


\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_user_mumps_advanced

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

Further options are available through the function
\ref cs_param_sles_mumps_advanced in order to reach a better performance
but this is case dependent.

There are three main stages in the resolution:
1. The *analysis* stage : renumbering to improve the fill-in of the factorization
2. the *factorization* stage
3. The *solve* stage (apply the factorization)

#### Renumbering

During the *analysis* stage, an important step is the
*renumbering*. The quality of the analysis and in particular of the
renumbering has a direct impact on the efficiency of the factorization
and its fill-in (memory consumption).

Several algorithms are available according to your
MUMPS installation. Some of these algorithms rely on external
libraries such as (PT-)Scotch ou (Par)Metis

The type of algorithm which can be set from code_saturne is stored
in \ref  cs_param_mumps_analysis_algo_t

| parameter | description | dependency |
| ---       | ---         | ---        |
| \ref CS_PARAM_MUMPS_ANALYSIS_AMD  | *AMD* sequential algorithm. Well-suited for 2D cases. | no |
| \ref CS_PARAM_MUMPS_ANALYSIS_QAMD | *QAMD* sequential algorithmn. Well-suited for 2D cases. | no |
| \ref CS_PARAM_MUMPS_ANALYSIS_PORD | Sequential algorithmn. For 3D cases if there is no dependency in your MUMPS installation | no |
| \ref CS_PARAM_MUMPS_ANALYSIS_SCOTCH | Sequential algorithm for 3D cases. | Scotch or PT-Scotch |
| \ref CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH | Parallel algorithm for 3D cases. | PT-Scotch |
| \ref CS_PARAM_MUMPS_ANALYSIS_METIS | Sequential algorithm for 3D cases. Often the best choice | METIS or PARMETIS |
| \ref CS_PARAM_MUMPS_ANALYSIS_PARMETIS | Parallel algorithm for 3D cases. | PARMETIS |
| \ref CS_PARAM_MUMPS_ANALYSIS_AUTO | The choice is performed by MUMPS according to the available algorithms and the type of matrices |  |

#### Iterative refinement

When the system is really ill-conditioned, it can be useful to add one
or several steps of *iterative refinement* to improvement the
solution. By default, no iteration is performed. A fixed number of
iterations is always performed.

#### Block Low Rank (BLR) algorithms

There is no BLR activated if the value of the parameter related to the
BLR algorithm is set to 0. When the value is positive, BLR is
activated for the factorization and solve stages using the `UCFS`
variant and the value is used to indicate the level of
compression. When the value is negative, the `UFSC` algorithm is
used instead.

If the option \ref CS_PARAM_MUMPS_MEMORY_CONSTRAINED is set, other
compressions are activated (cleaning of the workspace, `ICNTL(37)=1`
and `ICNTL(40)=1`) in addition to the BLR compression.

#### Memory increase

In some cases (for instance with 2D cases), the estimation of the
memory consumption is not sufficient and an error message is
returned. To circumvent this issue, one needs to set a correct
`mem_coef` parameter. This is directly given to `ICNTL(14)` in MUMPS.

### Interfaces with MUMPS structures

To get the finest control on the tuning of the MUMPS solver, one has
to use the function \ref cs_user_sles_mumps_hook

This function is called before the main stages:
1. Before the analysis step
2. Before the factorization step

\snippet cs_user_parameters-cdo-linear_solvers.c mumps_user_hook

According to the configuration *single* or *double* precision, the
structure storing the MUMPS solver differs. In case of a *single*
precision factorization, the cast is done to a pointer to a
`SMUMPS_STRUC_C` structure


# Solve a saddle-point linear system {#cs_ug_cdo_sles_saddle}

Saddle-point systems arise from the monolithic velocity-pressure
coupling of the Navier-Stokes equations or from the discretization of
CDO cell-based schemes (steady diffusion equation equation for
instance).  Some examples of settings involving differents strategies
of resolution are presented hereafter.

One considers the following type of saddle-point systems \f$\mathsf{Mx = b}\f$ with
\f[
 \mathsf{M} = \begin{bmatrix}
    A & B^T \\
    B & 0
 \end{bmatrix}
 \qquad
 \text{and}
 \qquad
 \mathsf{b} =
 \begin{bmatrix}
    f \\
    g
 \end{bmatrix}
\f]
where \f$A\f$ is a square non-singular matrix (symmetric or not
according to the numerical and/or modelling choices).

The key `CS_EQKEY_SADDLE_SOLVER` enables to set the solver (or the strategy) to
solve the saddle-point problem associated to an equation. More details are
gathered in the following table.

key value | description | type | family
:--- | :--- | :--- | :---:
`"none"` | No solver. No saddle-point system to solve. | \ref CS_PARAM_SADDLE_SOLVER_NONE | saturne
`"alu"` | Uzawa algorithm with an Augmented Lagrangian acceleration. Segregated technique. | \ref CS_PARAM_SADDLE_SOLVER_ALU | saturne
`"fgmres"` | Flexible GMRES Krylov solver on the full system. Possibility to consider block preconditioning. The choice to approximate the Schur complement has an impact on the effeciency of the algorithm. Monolithic technique. | \ref CS_PARAM_SADDLE_SOLVER_FGMRES | PETSc
`"gcr"` | GCR (Generalized Conjugate Residual) Krylov solver on the full system along with block preconditioning. The choice to approximate the Schur complement has an impact on the effeciency of the algorithm. Monolithic technique. | \ref CS_PARAM_SADDLE_SOLVER_GCR | saturne
`"gkb"` | Golub-Kahan bidiagonalization.  An augmentation of the (1,1) block is possible. Use this technique on **symmetric** saddle-point systems. Segregated technique. | \ref CS_PARAM_SADDLE_SOLVER_GKB | saturne
`"minres"` | MINRES Krylov solver on the full system along with block preconditioning. Close to the choice `"gcr"` but only for **symmetric** saddle-point systems. Monolithic technique. | \ref CS_PARAM_SADDLE_SOLVER_MINRES | saturne
`"mumps"` | MUMPS sparse direct solver on the full system. Monolithic technique. | \ref CS_PARAM_SADDLE_SOLVER_MUMPS | MUMPS
`"notay"` | Notay's algebraic transformation. Monolithic technique relying on a FMGRES on the modified system. Block preconditioning can be specified on the transformed blocks. (Experimental) | \ref CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM | PETSc
`"uzawa_cg"` | Uzawa algorithm accelerated by a CG technique. Only for **symmetric** saddle-point systems. The choice to approximate the Schur complement has an impact on the effeciency of the algorithm. Segregated technique. | \ref CS_PARAM_SADDLE_SOLVER_UZAWA_CG | saturne

In segregated technique, the (1,1)-block solver refers to the main solver. It
can be set following the rationale described in \ref cs_ug_cdo_sles_set In
monolothic technique, the (1,1)-block solver refers to the approximation of
this block for preconditioning.  Moreover, for some strategies, additional SLES
such as the Schur complement system can also be specified using the same
rationale described in \ref cs_ug_cdo_sles_set.



## Augmented Lagrangian Uzawa algorithm (ALU)

Here is a first example to set a saddle-point solver. One assumes that
the external library MUMPS has been installed and has been configured
with code_saturne (see the installation guide for more details
[here](@ref cs_dg_build_system)).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_alu

Since the velocity block is _augmented_ by the term
\f$\underline{\mathsf{grad}}\,\mathsf{div}\f$ (with a scaling coefficient given
by `CS_EQKEY_SADDLE_AUGMENT_SCALING`), the resulting linear system is hard to
solve for an iterative method (or one needs a very specific preconditionner for
\f$H(\mathsf{div})\f$).  Our strategy is to consider the sparse direct solver
MUMPS \cite MUMPS01. Here is an example of an advanced usage of the MUMPS
settings.

Additional options may be set to improve the performance of the sparse
direct solver.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_mumps

Since the velocity block is vector-valued, one can benefit from a
block analysis (the third parameter is set to 3).

Here is a second example.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_alu_mumps


## Block preconditioners with a Krylov solver

### In-house solvers

The two in-house Krylov solvers available with block preconditioning are
- **Minres** (Minimal residual) algorithm for symmetric indefinite systems such
  as systems encountered in the case of Stokes systems
- **GCR** (Generalized conjugate residual) algorithm for general indefinite
  systems. This is a flexible solver (preconditioner may vary between two
  iterations).

These two algorithms are optimized to handle saddle-point problems in
code_saturne since the (1,2) and (2,1) which are transposed is stored only
once. Moreover, this block is stored in an unassembled way.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_minres

### PETSc solvers

Other block preconditioners with a Krylov solver can be set using the external
librairy [PETSc](https://petsc.org). The two main differences (in addition to
the use of the PETSc library) are:
1. The Krylov solver is a flexible GMRES algorithm
2. The Schur complement approximation can differ since one hinges on the
   choices available in PETSc. Please refer to the PETSc documentation for more
   details.

One uses the `FIELDSPLIT` functionnality to set the block preconditioning strategy.



## Golub-Kahan Bidiagonalization algorithm (GKB)

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_gkb_kcycle

The linear system may be augmented to improve the convergence rate of the
algorithm (but the system is harder to solve). Here is another example:

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_gkb_mumps

## Notay's algebraic transformation (Experimental)

Here is an example in the case of a saddle-point system stemming from the
Navier-Stokes equations.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_sles_notay

## Uzawa algorithm with a CG acceleration

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_navsto_uzacg



## Schur complement approximation

The Schur complement related to the saddle-point system \f$\mathsf{Mx = b}\f$
detailed above is
\f$\mathsf{S} = \mathsf{-B\cdot A^{-1} \cdot B}\f$
This is a square matrix. The cost to build exactly this matrix
is prohibitive. Thus, one needs an approximation denoted by
\f$\mathsf{\widetilde{S}}\f$

An approximation of the Schur complement is needed by the following
strategies:
- Block preconditionner with a Krylov method as in `"fgmres"` (only
  with PETSc), `"gcr"` or `"minres"`
- Uzawa algorithm with a CG acceleration: `"uzawa_cg"`

The efficiency of the strategy to solve the saddle-point problem is
related to the quality of the approximation. Several types of
approximation of the Schur complement are available and a synthesis of
the different options is gathered in the following table.

key value | description | related type | family
:--- | :--- | :--- | :---:
`"none"` | No approximation. The Schur complement approximation is not used. | \ref CS_PARAM_SADDLE_SCHUR_NONE | all
`"diag_inv"` | \f$\mathsf{\widetilde{S}} = \mathsf{-B\cdot diag(A)^{-1} \cdot B}\f$ | \ref CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE | saturne, petsc
`"identity"` | The identity matrix is used. | \ref CS_PARAM_SADDLE_SCHUR_IDENTITY | saturne, petsc
`"lumped_inv"` | \f$\mathsf{\widetilde{S}} = \mathsf{-B\cdot L \cdot B}\f$ where \f$\mathsf{A\,L = 1}\f$, \f$\mathsf{1}\f$ is the vector fill with the value 1. | \ref CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE | saturne
`"mass_scaled"` | The pressure mass matrix is used with a scaling coefficient. This scaling coefficient is related to the viscosity for instance when dealing with the Navier-Stokes system. | \ref CS_PARAM_SADDLE_SCHUR_MASS_SCALED | saturne
`"mass_scaled_diag_inv"` | This is a combination of the `"mass_scaled"` and `"diag_inv"` approximation. | \ref CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE | saturne
`"mass_scaled_lumped_inv"` | This is a combination of the `"mass_scaled"` and `"lumped_inv"` approximation. | \ref CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE | saturne

To set the way to approximate the Schur complement use
```c
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEquationName");

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, key_value);
```

For the Schur approximation implying the `"diag_inv"` or the
`"lumped_inv"` approximation, one needs to build and then solve a
linear system associated to the Schur approximation. In the case of
the `"lumped_inv"`, an additionnal linear system is solved to define
the vector \f$\mathsf{L}\f$ from the resolution of \f$\mathsf{A\,L =
1}\f$. Please consider the examples below for more details.

### Example 1: "mass_scaled"

In the first example, one assumes that the saddle-point solver is
associated to the structure \cs_equation_param_t nammed `mom_eqp`. For
instance, one uses a `gcr` strategy for solving the saddle-point
problem. This is one of the simpliest settings for the Schur
complement approximation. The `"mass_scaled"` approximation is
well-suited for solving the Stokes problem since the pressure mass
matrix is a rather good approximation (in terms of spectral behavior)
of the Schur complement.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_schur_mass_scaled


### Example 2: "mass_scaled_diag_inv"

In the second example, one assumes that the saddle-point solver is
associated to the structure \cs_equation_param_t nammed `mom_eqp`. One
now considers a more elaborated approximation which needs to solve a
system implying the Schur complement approximation.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_schur_mass_scaled_diag_inv


### Example 3: "lumped_inv"

In the third and last example, one assumes that the saddle-point solver is
associated to the structure \cs_equation_param_t nammed `mom_eqp`. One
now considers a more elaborated approximation which needs to solve a
system implying the Schur complement approximation and an additional
system to build the Schur complement approximation (the system related
to \f$\mathsf{A\,L = 1}\f$).

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_schur_lumped_inv


# Additional settings

## Immediate exit

In some cases, the initial residual is very small leading to an
*immediate exit*.  If this happens and it should not be the case, one
can specify a stronger criterion to trigger the immediate exit. In the
function \ref cs_user_linear_solvers, add the following lines for
instance:

\snippet cs_user_parameters-cdo-linear_solvers.c linear_solver_immediate_exit

## Allow no operation

An immediate exit is possible (i.e. there is no solve) when the norm of the
right-hand side is equal to zero or very close to zero (\ref cs_sles_set_epzero).
To allow this kind of behavior for the SLES associated to an equation, one has
to set the key `CS_EQKEY_SOLVER_NO_OP` to `true`.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_allow_no_op
