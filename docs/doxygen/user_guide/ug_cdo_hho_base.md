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

\page cs_ug_cdo_hho_base Using CDO/HHO schemes

[TOC]

<!--
    References used in this page
-->

[Bonel14]: https://hal.archives-ouvertes.fr/tel-01116527
[Cant16]: https://hal.science/tel-01419312v1
[Mila20]: https://hal.science/tel-03080530v2
[MiBoE22]: https://hal.science/hal-03215118v1


Introduction {#sec_ug_cdo_hho_intro}
============

**Compatible Discrete Operator** (CDO) and **Hybrid High-Order** (HHO) schemes
are available in code_saturne as alternative space discretizations to the
legacy Finite Volumes (FV) schemes. These are new generation of space
discretizations dedicated to polyhedral meshes.  The capabilities of these
schemes in terms of physical modelling are limited compared to the legacy FV
schemes.

Overview of physical models
------------

With **CDO schemes**, the resolution of the following models/equations are
available:
    - **User-defined scalar-valued or vector-valued equations** with
      - unsteady term
      - advection term (the advection-field can be user-defined or one
        retrieves the mass flux from the Navier--Stokes equations or the Darcy
        velocity from the Richards equation for instance)
      - diffusion term with a homogeneous or heterogeneous +
        isotropic/orthotropic or anisotropic diffusivity tensor
      - reaction term (also called *implicit source term* in the FV parts of
        this documentation)
      - source terms (also called *explicit source term* in the FV parts of
        this documentation)
    - **Steady/unsteady Stokes or Navier--Stokes equations** in a laminar flow
      regime : NavSto module
      - Several velocity/pressure coupling algorithms are available: monolithic,
        artifical compressibility or prediction/correction algorithm
      - Please refer to the [dedicated page](@ref cs_ug_cdo_navsto) for more details.
    - **Thermal systems** with the possibility to model phase-change
      (liquid/solid)
    - **Groundwater flows**: GWF module
      -  This module computes the Darcy velocity from the hydraulic state of
         the soil seen as porous media and can advect tracers (possibly
         radioactive tracers) in that soil.
      -  Please refer to [the dedicated page](@ref cs_ug_cdo_gwf)
         for more details
    - **Solidification process**:
        - This module relies on the NavSto module and the Thermal module. Phase
          transitions between solid and liquid states (solidification or
          melting) are taken into account.
        - Segregation phenomena can be considered with binary alloys
        - Please refer to [the dedicated page](@ref cs_ug_cdo_solidification)
          for more details
    - **Wall distance** computations.
    - **Maxwell module**. This is an on-going work. Up to now, electro-static
      and magneto-static problems are available.


To get started
------------

A simple example relying on the resolution of the Laplace equation (steady
isotropic diffusion equation with Dirichlet boundary conditions) in a cube with
CDO schemes is available [here](@ref cs_cdo_laplacian). This is a a good
starting point for beginers.


Case settings for CDO/HHO schemes
=============

To set-up a CDO computation, one has to update the cs_user_parameters.c file and
edit successively the following functions
    - \ref cs_user_model,
    - \ref cs_user_parameters
    - \ref cs_user_finalize_setup

In addition, volume or boundary zones useful to set properties, equations,
etc. are defined either thanks to the GUI or by editing the function \ref
cs_user_zones for more advanced definitions from the file
`cs_user_zones.c`. This is similar to what is done with the legacy FV schemes.

More precisely, in each function above-mentioned the following settings can be
done:
    - \ref cs_user_model
      + Activation of CDO/HHO schemes,
      + Definition of the domain boundaries
      - Add new user-defined equations
      - Add new user-defined advection fields
      - Add new user-defined properties
      - Activate predefined modules (GWF, solidification, Navier--Stokes,
        thermal system or wall-distance computation)
    - \ref cs_user_parameters
      +  Definition of the time stepping strategy (if different from the one
         pre-defined in th GUI)
      +  Definition of the logging/restart frequency
      + Update of the numerical settings of each equation thanks to \ref
        cs_equation_param_set (including the way to solve the linear or the
        non-linear equations).
    - \ref cs_user_finalize_setup
      + Complete the definition of advection fields
      + Complete the definition of properties.
      + Complete the definition of an equation
        - Definition of boundary conditions
        - Definition of initial conditions
        - Defintion of source terms,
        - Activate predefined terms of a user-defined equation
          - The unsteady term by adding a link between an equation to a property
          - The advection term by adding a link to an advection field
          - The diffusion term by adding a link to a property
          - The reaction term by adding a link to a property


Settings done in cs_user_model()
===================

Activation of the CDO/HHO part
---------------------

The very first step is to activate the CDO module in the function
\ref cs_user_model There are two ways to switch on CDO/HHO schemes:
    - \ref CS_DOMAIN_CDO_MODE_ONLY for a usage of CDO or HHO in stand-lone
      (i.e. without the legacy Finite Volume approach)
    - \ref CS_DOMAIN_CDO_MODE_WITH_FV for a usage which can share some
    equations/models solved with CDO/HHO schemes and some other equations/models
    solved with the legacy Finite Volume approach.

CDO/HHO schemes can be activated within this function as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_activation

Domain boundaries {#ug_cdo_sec_domain_boundaries}
-------------------

Domain boundaries are useful for the (Navier-)Stokes equations or for the
computation of the wall distance. Several types of domain boundaries can be
defined. They are gathered in \ref cs_boundary_type_t (a flag built from a set
of predefined bits)

The definition of the domain boundaries for CDO/HHO schemes is performed in two
steps: (1) set the default boundary and then add other boundaries which do not
fit the default one. The two possible default domain boundaries are
\ref CS_BOUNDARY_WALL or \ref CS_BOUNDARY_SYMMETRY Here is a first example.

\snippet cs_user_parameters-cdo-condif.c param_cdo_domain_boundary

Here is a second example.

\snippet cs_user_parameters-cdo-navsto.c param_cdo_navsto_boundary


Activate predefined equations
-----------------

Several predefined models are available (and thus equation or set of
equations). When activating a predefined module, this implies:
    - adding equations and variable fields;
    - adding properties (and sometimes an advection field);
    - adapting the default numerical settings of the added equations to get more
      relevant settings.

### Navier-Stokes (NavSto module)

Please refer to the dedicated documentation available [here](@ref cs_ug_cdo_navsto).

### Groundwater flows (GWF)

Please refer to the dedicated documentation available [here](@ref cs_ug_cdo_gwf).

### Solidification/melting module

Please refer to the dedicated documentation available [here](@ref cs_ug_cdo_solidification).

### Thermal system

The activation of the thermal module yields the following actions:
- Add an equation called "thermal_equation" (the name of the associated variable
  depends on the choice of the variable, by default, the variable is
  "temperature")
- Add the properties "mass_density", "thermal_capacity" and "thermal_conductivity"

For the thermal equation, the default boundary condition is a no flux
(i.e. a homogeneous Neumann boundary condition).  Here is the simplest
example to activate the thermal module.

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_activate_thermal_solver

The first parameter is a flag to describe the thermal model to consider. This
flag can be built from the following tags (\ref cs_thermal_model_type_bit_t)

- \ref CS_THERMAL_MODEL_STEADY
- \ref CS_THERMAL_MODEL_NAVSTO_ADVECTION

To specify the choice of the variable used in the thermal model (by default, the
temperature in Kelvin). This can be modified by adding the tag

- \ref CS_THERMAL_MODEL_IN_CELSIUS

and to change the main variable use the one of the following tag

- \ref CS_THERMAL_MODEL_USE_ENTHALPY
- \ref CS_THERMAL_MODEL_USE_TOTAL_ENERGY

The two previous options are not fully implemented.

The second parameter is a flag related to the numerical aspects. There is no tag
available up to now.

The third parameter is a flag related to the activation of automatic
post-processings.

- \ref CS_THERMAL_POST_ENTHALPY

### Wall-distance computation

It is possible to activate the computation of the distance to the wall using CDO
schemes (CDO vertex-based schemes are used) as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_wall_distance



Add a user-defined equation
-------------------------

User-defined equation with CDO/HHO schemes are added thanks to a call to the
function \ref cs_equation_add_user in \ref cs_user_model Here are several
examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_equation


Add a property
---------------

User-defined properties are added thanks to a call to the function
\ref cs_property_add in \ref cs_user_model Here are several examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties

The function \ref cs_property_add returns a pointer to a \ref cs_property_t
structure which can be used to set advanced parameters. If the pointer to a \ref
cs_property_t structure is not available, this can be easily recover thanks to a
call to the function \ref cs_property_by_name

To enable the computation of the **Fourier number** related to a given property in
an unsteady simulation proceed as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_properties_opt


Add an advection field
---------------------

The definition of an advection field allows one to handle flows with a frozen
velocity field or the transport of scalar quantities without solving the
Navier-Stokes system. The add of a new user-defined advection field with CDO/HHO
schemes is specified as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field

When an advection field is defined, it is possible to retrieve it and then set
a post-processing operation. For instance, to activate the post-processing of
the **CFL number** proceed as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_user_adv_field_post



Settings done in cs_user_finalize_setup()
===================

Boundary conditions
-------------------

### User-defined equations

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_bcs

### Thermal system

The mechanism is the same as for user-defined equations. The name of the
equation is defined in the macro \ref CS_THERMAL_EQNAME

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_bc


Initial conditions
------------------

### Thermal system

The mechanism is detailed for the thermal system but the same mechanism can be
used for user-defined equations. In our example, the name of the equation is
defined in the macro \ref CS_THERMAL_EQNAME

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_ic

where the function `_initial_temperature` has a predefined prototype (cf. the
definition of the function pointer \ref cs_analytic_func_t)

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_initial_temperature_function


### GWF module

Here is another example extracted from the file `cs_user_parameters-cdo-gwf.c`
(in src/user_examples) but this is readily applicable to any equation.

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_set_ic

where the function `get_tracer_ic` has a predefined prototype (cf. the
definition of the function pointer \ref cs_analytic_func_t)

\snippet cs_user_parameters-cdo-gwf.c param_cdo_gwf_get_tracer_ic


Definition of properties
------------------

### User-defined properties

When a property has been added, the second step is to define this
property. According to the type of property (isotropic, orthotropic or
anisotropic) definitions differ.  Here are two examples:

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_property

### Thermal system

The mechanism is the same as for user-defined properties. The named are
predefined and associated to the following macros:

- \ref CS_THERMAL_LAMBDA_NAME for the thermal conductivity
- \ref CS_THERMAL_CP_NAME for the thermal capacity
- \ref CS_PROPERTY_MASS_DENSITY for the mass density. This property is also used
  when the NavSto module is activated.

\snippet cs_user_parameters-cdo-thermal-solver.c param_cdo_define_thermal_properties


Definition of an advection field
-----------------

When an advection field has been added, the second step is to define this
advection field. Here are is an example of definition using an anlytic function
and the activation of optional features:

\snippet cs_user_parameters-cdo-condif.c param_cdo_setup_advfield


Definition of source terms {#cs_ug_cdo_hho_base_source_term}
---------------------

### User-defined equation

A first simple example.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_simple_source_terms

The second example shows an advanced usage relying on an analytic function and
an input structure. Moreover, a more accurate quadrature rule is specified.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_source_terms

The user-defined function to compute the source term is specified as follows

\snippet cs_user_parameters-cdo-condif.c param_cdo_condif_analytic_st

and the function for the memory management of a \ref cs_xdef_t structure is

\snippet cs_user_parameters-cdo-condif.c param_cdo_condif_free_input




Add diffusion, advection, etc. to a user-defined equation
---------------------

Add terms to an equation like a diffusion term, an advection term,
unsteady term, reaction terms.

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_terms

In some cases, one can also add less common terms such as a
\f$\mathsf{\underline{grad}\cdot div}\f$ or
\f$\mathsf{\underline{curl}\cdot\underline{curl} }\f$ (only available
with **CDO edge_based** schemes).

\snippet cs_user_parameters-cdo-condif.c param_cdo_add_terms_2







Settings done in cs_user_parameters() {#cs_ug_cdo_hho_base_user_param}
===================

Logging options
----------------

The management of the level and frequency of details written by the code can be
specified for CDO/HHO schemes as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_domain_output


Time stepping strategy
----------------

The management of the time step with CDO/HHO schemes can be specified as
follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_time_step

This can be completed with numerical options to specify the time scheme. See
the next section)

Discretization parameters associated to an equation {#cs_ug_cdo_hho_base_set_eqparam}
-------------------------

To modify the numerical settings, the mechanism relies on a `(key, value)`
rationale.  The function to use is \ref cs_equation_param_set For instance, with
an equation called `"MyEquation"`, and a key called `CS_EQKEY_PARAM_TO_SET` with
the value `"value_to_set"`.

```c
cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEquation");

cs_equation_param_set(eqp, CS_EQKEY_PARAM_TO_SET, "value_to_set");
```

If one wants to modify the settings for all equations, this is possible using
the function \ref cs_equation_set_default_param

```c
cs_equation_set_default_param(CS_EQKEY_PARAM_TO_SET, "value_to_set");
```

All available keys are listed and described with \ref cs_equation_key_t
One gives some examples for some of them.

### Set the space discretization scheme

The key is \ref CS_EQKEY_SPACE_SCHEME with the possible value gathered in the
following table

key_value | description
:--- | :---
`"cdo_vb"` or `"cdovb"` | Switch to a **CDO vertex-based discretization** (degrees of freedom are located at the mesh vertices). One value per vertex in the case of scalar-valued equation and three values per vertex in the case of a vector-valued equation.
`"cdo_vcb"` or `"cdovcb"` | Switch to a **CDO vertex+cell-based discretization** (degrees of freedom are located at the mesh vertices and at the mesh cells). One value per vertex and per cell in the case of scalar-valued equation and three values per vertex and pêr cell in the case of a vector-valued equation. Thanks to a static condensation operation, the algebraic system is reduced to only vertex unknows.
`"cdo_fb"` or `"cdofb"` | Switch to a **CDO face-based discretization** (degrees of freedom are located at interior and boundary faces and at mesh cell). One value per face and mesh cell in the case of scalar-valued equation and three values per face and mesh cell in the case of a vector-valued equation. Thanks to a static condensation operation, the algebraic system is reduced to only face unknows.
`"cdo_cb"` or `"cdocb"` | Switch to a **CDO cell-based discretization** (degrees of freedom are located at mesh cells for the potential and at faces for the flux). Only scalar-valued equation are possible. One value per cell for the potential. One value per face for the flux unknown (the normal component of the flux).
`"cdo_eb"` or `"cdoeb"` | Switch to **CDO edge-based discretization** (degrees of freedom are located at mesh edges, one scalar per edge corresponding to the circulation). Only vector-valued equation are handled with this discretization.
`"hho_p0"` | Switch to a **HHO(k=0)** discretization relying on \f$P_0\f$ polynomial approximation. (degrees of freedom are located at interior and boundary faces and at mesh cells). One value per face and per mesh cell in the case of scalar-valued equation and three values per face and mesh cell in the case of a vector-valued equation.  Thanks to a static condensation operation, the algebraic system is reduced to only face unknows.
`"hho_p1"` | Switch to a **HHO(k=1)** discretization relying on \f$P_1\f$ polynomial approximation. (degrees of freedom are located at interior and boundary faces and at mesh cells). Three values per face and four values per cell in the case of scalar-valued equation and nine values per face and 12 values per cell in the case of a vector-valued equation.  Thanks to a static condensation operation, the algebraic system is reduced to only face unknows.
`"hho_p2"` | Switch to a **HHO(k=2)** discretization relying on \f$P_2\f$ polynomial approximation. (degrees of freedom are located at interior and boundary faces and at mesh cells). Six values per face and ten values per cell in the case of scalar-valued equation and 18 values per face and 30 values per cell in the case of a vector-valued equation. Thanks to a static condensation operation, the algebraic system is reduced to only face unknows.

An example of usage:

\snippet cs_user_parameters-cdo-condif.c param_cdo_numerics

More details can be found in \cite Bonel14 for CDO-Vb, CDO-Fb, CDO-Cb and
CDO-Eb. CDO-Fb are also detailed in \cite Mila20. CDO-VCb are detailed in
\cite Cant16


### Set the advection scheme

\snippet cs_user_parameters-cdo-condif.c param_cdo_conv_numerics

The available advection schemes are listed in the description of the key
\ref CS_EQKEY_ADV_SCHEME

\snippet cs_user_parameters-cdo-condif.c param_cdo_conv_schemes

### Set the time scheme

When the equation to solve is unsteady, one has to specify a time
discretization scheme. Available time schemes are listed in the description of
the key \ref CS_EQKEY_TIME_SCHEME
By default, a first order implicit Euler scheme is used. To modify this default
settings, please proceed as follows:

\snippet cs_user_parameters-cdo-condif.c param_cdo_time_schemes

The mass matrix associated to the unsteady term is also a parameter. One can
use either a mass matrix like in FV scheme using a `"voronoi"` algorithm
(default) or the `"wbs"` algorithm like in Finite Element (FE) schemes.

\snippet cs_user_parameters-cdo-condif.c param_cdo_time_hodge

### Modify the numerical scheme for the diffusion term

Several algorithms are available to build the diffusion term. They are all
listed in the description of the key \ref CS_EQKEY_HODGE_DIFF_ALGO Please refer
to \cite Bonel14 for more details In the case of the `cost` (or `ocs`), one can
specify the value of the scaling coefficient in front of the stabilization
part. This is done using the key \ref CS_EQKEY_HODGE_DIFF_COEF

\snippet cs_user_parameters-cdo-condif.c param_cdo_diff_numerics


Linear algebra settings {#cs_ug_cdo_hho_base_linalg}
-------------------

Several examples are detailed hereafter to specify a solver different from the
default one integrated to code_saturne.


### A first example

For a symmetric positive definite (SPD) system, the solver of choice is the
conjugate gradient (\ref CS_PARAM_SOLVER_CG) or its flexible variant (\ref
CS_PARAM_SOLVER_FCG). In this case, a good preconditioner is a multilevel
algorithm such as an algebraic multigrid (AMG). Different multigrid techniques
are available when using code_saturne depending on the installation
configuration (in-house implementations or algorithms available from PETSc
and/or HYPRE libraries).

Here is a first example:
\snippet cs_user_parameters-cdo-condif.c param_cdo_sles_settings1

If the external library [PETSc](https://petsc.org) is available, one uses its
algebraic multigrid, otherwise one uses the in-house solver. Here a flexible
conjugate gradient (`fcg`) with a diagonal preconditionning (`jacobi`).


### Set the family of solvers

Some of the choices are only available with external libraries. The external
libraries which can be linked to code_saturne are:
- [HYPRE library](https://hypre.readthedocs.io)
- [MUMPS](https://mumps-solver.org), a robust sparse direct solver \cite MUMPS01
- [PETSc](https://petsc.org)

\ref CS_EQKEY_SOLVER_FAMILY allows one to specify the family of linear solvers
to consider. This can be useful if an automatic switch is not done or when
there are several possibilities (for instance when using HYPRE solver either
directly from the library or through the PETSc library or using the MUMPS
solver either directly or through the PETSc library).

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_solver_family


### Set the solver {#cs_ug_cdo_hho_base_solver}

Available linear solvers are listed in \ref cs_param_solver_type_t and can be
set using the key \ref CS_EQKEY_SOLVER

key value | description | type | family
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


### Set the preconditioner {#cs_ug_cdo_hho_base_precond}

Available preconditioners are listed in \ref cs_param_precond_type_t and can be
set using the key \ref CS_EQKEY_PRECOND

key value | description | type | family
:--- | :--- | :--- | :---:
`"amg"` | Algebraic MultiGrid (AMG) iterative solver. See [this section](@ref cs_ug_cdo_hho_base_amg) for more details. | \ref CS_PARAM_PRECOND_AMG |  saturne, PETSc, HYPRE
`"block_jacobi"`, `"bjacobi"` | Block Jacobi with an ILU(0) factorization in each block. By default, there is one block per MPI rank. | \ref CS_PARAM_PRECOND_BJACOB_ILU0 | PETSc
`"bjacobi_sgs"`, `"bjacobi_ssor"` | Block Jacobi with a SSOR algorithm in each block. By default, there is one block per MPI rank. For a sequential run, this is the SSOr algorithm. | \ref CS_PARAM_PRECOND_BJACOB_SGS | PETSc
`"diag"`, `"jacobi"` | Diagonal preconditioning | \ref CS_PARAM_PRECOND_DIAG | saturne, PETSc, HYPRE
`"ilu0"` | Incomplete LU factorization (zero fill-in meaning that the same sparsity level is used). | \ref CS_PARAM_PRECOND_ILU0 | PETSc, HYPRE
`"icc0"` | Incomplete Choleski factorization (zero fill-in meaning that the same sparsity level is used). For SPD matrices. | \ref CS_PARAM_PRECOND_ICC0 | PETSc, HYPRE
`"mumps"` | the sparse direct solver MUMPS used as preconditioner. Please refer to [the MUMPS section](@ref cs_ug_cdo_hho_base_mumps) for more details. | \ref CS_PARAM_PRECOND_MUMPS | mumps
`"none"` | No preconditioner. | \ref CS_PARAM_PRECOND_NONE | saturne, PETSc, HYPRE
`"poly1"` | 1st order Neumann polynomial preconditioning | \ref CS_PARAM_PRECOND_POLY1 | saturne
`"poly2"` | 2nd order Neumann polynomial preconditioning | \ref CS_PARAM_PRECOND_POLY2 | saturne
`"ssor"` | Symmetric Successive OverRelaxation (SSOR) algorithm. In case of parallel computations, each MPI rank performs a SSOR on the local matrix. | \ref CS_PARAM_PRECOND_SSOR | PETSc




### Algebraic multigrids {#cs_ug_cdo_hho_base_amg}

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
`"v_cycle"` | This is the default choice. This corresponds to an in-house V-cycle detailed in \cite MeFoH09 The main ingredients can set using the functions \ref cs_param_sles_amg_inhouse and \ref cs_param_sles_amg_inhouse_advanced Please refer to [this section](@ref cs_ug_cdo_hho_base_amg_inhouse) for more details.
`"k_cycle"` or `"kamg"` | Switch to a K-cycle strategy. This type of cycle has been detailed in \cite Notay05 Be aware that a lot of work can be done in the coarsest levels yielding more communications between MPI ranks during a parallel computation. The main ingredients can set using the functions \ref cs_param_sles_amg_inhouse and \ref cs_param_sles_amg_inhouse_advanced Please refer to [this section](@ref cs_ug_cdo_hho_base_amg_inhouse) for more details.
`"boomer"`, `"bamg"` or "boomer_v" | Switch to use the V-cycle of the BoomerAMG algorithm from the [HYPRE library](https://hypre.readthedocs.io)  Please refer to [this section](@ref cs_ug_cdo_hho_base_amg_boomer) for more details on the way to set the main ingredients.
`"boomer_w"` or `"bamg_w"` | Switch to a W-cycle and the BoomerAMG algorithm from the [HYPRE library](https://hypre.readthedocs.io). Please refer to [this section](@ref cs_ug_cdo_hho_base_amg_boomer) for more details on the way to set the main ingredients.
`"gamg"` or `"gamg_v"` | Switch to a V-cycle and the GAMG algorithm from the [PETSc](https://petsc.org) library. iGAMG means "Geometric Algebraic multigrid. Please refer to [this section](@ref cs_ug_cdo_hho_base_amggamg) for more details on the way to set the main ingredients.
`"gamg_w"` | Switch to a W-cycle and the GAMG algorithm from the [PETSc](https://petsc.org) library. GAMG means "Geometric Algebraic multigrid. Please refer to [this section](@ref cs_ug_cdo_hho_base_amg_gamg) for more details on the way to set the main ingredients.


#### In-house multigrids {#cs_ug_cdo_hho_base_amg_inhouse}

The in-house multigrid algorithms can be easily tuned thanks to the function
\ref cs_param_sles_amg_inhouse More advanced settings can also be added using
the function \ref cs_param_sles_amg_inhouse_advanced

\snippet cs_user_parameters-linear_solvers.c sles_kamg_momentum

In the case of an advanced settings, it is possible to keep the default
settings associated to a parameter using the value \ref CS_CDO_KEEP_DEFAULT

It is also possible to access directly the structures and set the
parameters. Here is an example of a such usage:

\snippet cs_user_parameters-linear_solvers.c sles_mgp_1

or

\snippet cs_user_parameters-linear_solvers.c sles_mgp_2

In order to get the best efficiency in case of HPC computations, parallel grid
merging can be optimized as follows:

\snippet cs_user_parameters-linear_solvers.c sles_mg_parall


#### BoomerAMG {#cs_ug_cdo_hho_base_amg_boomer}

The main ingredients of the BoomerAMG algorithm are defined thanks to the
functions \ref cs_param_sles_boomeramg and \ref cs_param_sles_boomeramg_advanced
for additional parameters.

\snippet cs_user_parameters-cdo-linear_solvers.c cdo_sles_boomer

We strongly invite users to read the HYPRE documentation and especially the
part dedicated to BoomerAMG to understand the different choices
[BoomerAMG documentation(https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html)

HYPRE setting functions can be called thanks to a lower-level API (relying on a
setup hook function). Here is an example of a such usage

\snippet cs_user_parameters-linear_solvers.c sles_hypre_hook_1

#### GAMG {#cs_ug_cdo_hho_base_amg_gamg}

Up to now, the advanced settings of GAMG is possible inside \ref cs_user_linear_solvers
Here is an example of a such usage

\snippet cs_user_parameters-linear_solvers.c sles_petsc_gamg_1


### Sparse direct solver: MUMPS {#cs_ug_cdo_hho_base_mumps}

For the hardest to solve linear systems or in order to check a numerical
scheme, it may be useful to solve the linear systems with a direct solver. In
code_saturne, this is possible using the MUMPS library.



To go further
=============

The detailed description of CDO schemes, the mathmatical analysis and numerical
results on benchmarks are available in the following PhD thesis:

* [**PhD thesis**: "Compatible Discrete Operator schemes on polyhedral meshes for elliptic and Stokes equations"][Bonel14] \cite Bonel14
* [**PhD thesis**: "Approximation of scalar and vector transport problems on polyhedral meshes"][Cant16] \cite Cant16
* [**PhD thesis**: "Compatible Discrete Operator schemes for the unsteady incompessible Navier-Stokes equations"][Mila20] \cite Mila20


Additional publications :
* [**Article**: "Artificial compressibility methods for the incompressible Navier-Stokes equations using lowest-order face-based schemes on polytopal meshes"][MiBoE22] \cite MiBoE22
