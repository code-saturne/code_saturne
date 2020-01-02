/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
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
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /*--------------------------------------------------------------------------*/

  /*! [turbulence_model_choice] */

  /* Example: Chose a turbulence model
   *   CS_TURB_NONE: no turbulence model (laminar flow)
   *   CS_TURB_MIXING_LENGTH: mixing length model
   *   CS_TURB_K_EPSILON: standard k-epsilon model
   *   CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
   *     Linear Production (LP) correction
   *   CS_TURB_K_EPSILON_LS: Launder-Sharma low Re k-epsilon model
   *   CS_TURB_K_EPSILON_QUAD: Baglietto et al. low Re
   *   CS_TURB_RIJ_EPSILON_LRR: Rij-epsilon (LRR)
   *   CS_TURB_RIJ_EPSILON_SSG: Rij-epsilon (SSG)
   *   CS_TURB_RIJ_EPSILON_EBRSM: Rij-epsilon (EBRSM)
   *   CS_TURB_LES_SMAGO_CONST: LES (constant Smagorinsky model)
   *   CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic Smagorisky model)
   *   CS_TURB_LES_WALE: LES (WALE)
   *   CS_TURB_V2F_PHI: v2f phi-model
   *   CS_TURB_V2F_BL_V2K: v2f BL-v2-k
   *   CS_TURB_K_OMEGA: k-omega SST
   *   CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */

  cs_turb_model_t *turb_model = cs_get_glob_turb_model();
  turb_model->iturb = CS_TURB_K_EPSILON_LIN_PROD;

  /*! [turbulence_model_choice] */

  /*--------------------------------------------------------------------------*/

  /* Advanced choice of Wall function */
  {
    cs_wall_functions_t *wf = cs_get_glob_wall_functions();
     wf->iwallf = CS_WALL_F_2SCALES_VDRIEST;
  }


  /*--------------------------------------------------------------------------*/

  /*! [Rij_coupled_solver_choice] */
  /* Example: Coupled solver for Rij components (when iturb=30, 31 or 32)
   *   0: switch off
   *   1: switch on (default)
   */

  cs_turb_rans_model_t *rans_model = cs_get_glob_turb_rans_model();
  rans_model->irijco = 1;

  /* Advanced re-initialization for EBRSM or k-omega models
     - 0: switch off
     - 1: switch on (default)
     */

  rans_model->reinit_turb = 1;

  /* Turbulent diffusion model for second moment closure (CS_TURB_RIJ_*)
     - 0: scalar diffusivity (Shir model)
     - 1: tensorial diffusivity (Daly and Harlow model, default model)
     */
  rans_model->idirsm = 0;

  /*! [Rij_coupled_solver_choice] */

  /*--------------------------------------------------------------------------*/

  /* Reference velocity for turbulence initialization (m2/s)
    (useful only with turbulence) */

  cs_turb_ref_values_t *turb_ref_values = cs_get_glob_turb_ref_values();
  turb_ref_values->uref = 1.;

  /* Reference length scale in meters for initialization
     of epsilon (and specific clipping of turbulence, but
     this is not the default option)
     Assign a value of the order of the largest dimension of the
     physical domain in which the flow may develop.
     If a negative value is set here, or no value set and the GUI not
     used, the cubic root of the domain will be used.
     (useful only for turbulence). */

    turb_ref_values->almax = 0.5;

  /*--------------------------------------------------------------------------*/

  /*! [thermal_model_choice] */

  /* Example: Choose a thermal model
   *
   * 0: none
   * 1: temperature
   * 2: enthalpy
   * 3: total energy (only for compressible module)
   *
   *  For temperature, the temperature scale may be set later using itpscl
   *  (1 for Kelvin, 2 for Celsius).
   *
   *  Warning: When using specific physics, this value is
   *           set automatically by the physics model.
   *
   */

  cs_thermal_model_t *thermal_model = cs_get_glob_thermal_model();
  thermal_model->itherm = 1;

  /*! [thermal_model_choice] */

  /*--------------------------------------------------------------------------*/

  /* Volume of Fluid model with mass transfer Merkle model (cavitating flow)
   * to take into account vaporization / condensation */

  /*! [enable_cavit] */
  cs_vof_parameters_t *vof_param = cs_get_glob_vof_parameters();
  vof_param->vof_model = CS_VOF_ENABLED | CS_VOF_MERKLE_MASS_TRANSFER;
  /*! [enable_cavit] */

  /*--------------------------------------------------------------------------*/

  /*! [ALE_activation] */

  /* Example: activate ALE (Arbitrary Lagrangian Eulerian) method
   *   CS_ALE_NONE: switch off
   *   CS_ALE_LEGACY: legacy solver
   *   CS_ALE_CDO: CDO solver
   */

  cs_glob_ale = CS_ALE_LEGACY;

  /*! [ALE_activation] */

  /*--------------------------------------------------------------------------*/

  /*! [scalars_addition] */

  /* Example: add 2 scalar variables ("species" in the GUI nomenclature).
   *
   * Note that at this (very early) stage of the data setup, fields are
   * not defined yet. Associated fields will be defined later (after
   * model-defined fields) in the same order as that used here, and
   * after user-defined variables defined throught the GUI, if used.
   *
   * Currently, only 1d (scalar) fields are handled.
   *
   * parameters for cs_parameters_add_variable():
   *   name             <-- name of variable and associated field
   *   dim              <-- variable dimension
   */

  cs_parameters_add_variable("species_1", 1);
  cs_parameters_add_variable("tracer", 1);

  /*! [scalars_addition] */

  /*--------------------------------------------------------------------------*/

  /*! [scalars_variance_addition] */

  /* Example: add the variance of a user variable.
   *
   * parameters for cs_parameters_add_variable_variance():
   *   name          <-- name of variance and associated field
   *   variable_name <-- name of associated variable
   */

  cs_parameters_add_variable_variance("variance_1",
                                      "species_1");

  /*! [scalars_variance_addition] */

  /*--------------------------------------------------------------------------*/

  /*! [user_property_addition] */

  /* Example: add a user property defined on boundary faces.
   *
   * parameters for cs_parameters_add_property():
   *   name        <-- name of property and associated field
   *   dim         <-- property dimension
   *   location_id <-- id of associated mesh location, which must be one of:
   *                     CS_MESH_LOCATION_CELLS
   *                     CS_MESH_LOCATION_INTERIOR_FACES
   *                     CS_MESH_LOCATION_BOUNDARY_FACES
   *                     CS_MESH_LOCATION_VERTICES
   */

  cs_parameters_add_property("user_b_property_1",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);

  /*! [user_property_addition] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t *domain)
{

  /*! [time_stepping_options] */

  /* Time stepping  (0 : uniform and constant
                     1 : variable in time, uniform in space
                     2 : variable in time and space
                    -1 : steady algorithm) */


  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  time_opt->idtvar = 0;

  /*! [time_stepping_options] */

  /*! [ref_time_step] */

  /* Reference time step dt_ref
     The example given below is probably not adapted to your case. */

  cs_real_t dt_ref = 0.005;
  domain->time_step->dt_ref = dt_ref;

  /*! [ref_time_step] */

  /*! [duration] */

  /* Duration
     nt_max is absolute number of the last time step required
     if we have already run 10 time steps and want to
     run 10 more, nt_max must be set to 10 + 10 = 20 */

  domain->time_step->nt_max = (int) (40./dt_ref);

  /*! [duration] */

  /* Example: set options for Stokes solving */
  /*-----------------------------------------*/

  /* Members of the structure cs_glob_stokes_model:
   *  ivisse: take viscous term of transposed velocity
   *          gradient into account in momentum equation
   *                   - 1: true (default)
   *                   - 0: false
   *  arak: Arakawa multiplicator for the Rhie and Chow
   *        filter (1 by default)
   *  ipucou: pseudo coupled pressure-velocity solver
   *                   - 0: false (default)
   *                   - 1: true
   *  irevmc: reconstruct the velocity field after correction step
   *                   - 0: with the pressure increment gradient (default)
   *                   - 1: with an RT0 like formula using the mass fluxes
   *  iccvfg: calculation with a fixed velocity field
   *                   - 0: false (default)
   *                   - 1: true
   *  idilat: algorithm to take into account the density
   *          variation in time
   *                   - 0: Boussinesq approximation
   *                   - 1: dilatable steady algorithm (default)
   *                   - 2: dilatable unsteady algorithm
   *                   - 3: low-Mach algorithm
   *                   - 4: algorithm for fire
   *                   - 0: boussinesq algorithm with constant
   *                   density (not yet available)
   *  iphydr: improve hydrostatic pressure algorithm
   *                   - 0: no treatment (default)
   *                   - 1: impose the equilibrium of the hydrostaic
   *                     part of the pressure with any external force,
   *                     even head losses
   *                   - 2: compute an hydrostatic pressure due to
   *                     buoyancy forces before the prediction step
   *  igprij: improve static pressure algorithm
   *                   - 0: no treatment (default)
   *                   - 1: take -div(rho R) in the static pressure
   *                     treatment IF iphydr=1
   *  igpust: improve static pressure algorithm
   *                   - 0: no treatment
   *                   - 1: take user momemtum source terms in the
   *                     static pressure treatment IF iphydr=1 (default)
   *  fluid_solid: Has a solid zone where dynamics must be killed?
   *                   - false (default)
   *                   - true
   *
   */

  /*! [param_stokes_model] */
  {
    cs_stokes_model_t *stokes = cs_get_glob_stokes_model();
    stokes->arak = 0.;
  }
  /*! [param_stokes_model] */

  /* Example: change Reference fluid properties options */
  /*----------------------------------------------------*/

  /* Members of the structure cs_fluid_properties_t
   */

  /*! [param_fluid_properties] */
  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

    /*      ro0        : density in kg/m3
            viscl0     : dynamic viscosity in kg/(m s)
            cp0        : specific heat in J/(Kelvin kg)
            t0         : reference temperature in Kelvin
            p0         : total reference pressure in Pascal
                         the calculation is based on a
                         reduced pressure P*=Ptot-ro0*g.(x-xref)
                         (except in compressible case)
            xyzp0(3)   : coordinates of the reference point for
                         the total pressure (where it is equal to p0)

          In general, it is not necessary to furnish a reference point xyz0.
            If there are outlets, the code will take the center of the
            reference outlet face.
            On the other hand, if we plan to explicitly fix Dirichlet conditions
            for pressure, it is better to indicate to which reference the
            values relate (for a better resolution of reduced pressure).


          Other properties are given by default in all cases.

          Nonetheless, we may note that:

            In the standard case (no gas combustion, coal, electric arcs,
                                  compressibility):
            ---------------------
              ro0, viscl0 and cp0
                  are useful and represent either the fluid properties if they
                  are constant, either simple mean values for the initialization
                  if properties are variable and defined in cs_user_physical_properties.
              t0  is not useful
              p0  is useful but is not used in an equation of state. p0
                  is a reference value for the incompressible solver
                  which will serve to set the (possible) domain outlet pressure.
                  We may also take it as 0 or as a physical value in Pascals.

            With the electric module:
            ------------------------
              ro0, viscl0 and cp0
                  are useful but simply represent mean initial values;
                  the density, molecular dynamic viscosity, and specific
                  heat are necessarily defined as fields (whether they are
                  physically variable or not): see cs_user_physical_properties
                  for the Joule effect
                  module and the electric arcs dp_ELE data file.
              t0  is useful an must be in Kelvin (> 0) but represents a simple
                  initialization value.
              p0  is useful bu is not used in the equation of state. p0
                  is a reference value for the incompressible solver which
                  will be used to calibrate the (possible) outlet pressure
                  of the domain. We may take it as zero or as a physical
                  value in Pascals.

            With gas combustion:
            --------------------
              ro0 is not useful (it is automatically recalculated by the
                  law of ideal gases from t0 and p0).
              viscl0 is indispensable: it is the molecular dynamic viscosity,
                  assumed constant for the fluid.
              cp0 is indispensable: it is the heat capacity, assumed constant,
                  (modelization of source terms involving a local Nusselt in
                  the Lagrangian module, reference value allowing the
                  calculation of a radiative
                  (temperature, exchange coefficient) couple).
              t0  is indispensible and must be in Kelvin (> 0).
              p0  is indispensable and must be in Pascal (> 0).

            With pulverized coal:
            ---------------------
              ro0 is not useful (it is automatically recalculated by the
                  law of ideal gases from t0 and p0).
              viscl0 is indispensable: it is the molecular dynamic viscosity,
                  assumed constant for the fluid (its effect is expected to
                  be small compared to turbulent effects).
              cp0 is indispensable: it is the heat capacity, assumed constant,
                  (modelization of source terms involving a local Nusselt in
                  the coal or Lagrangian module, reference value allowing the
                  calculation of a radiative
                  (temperature, exchange coefficient) couple).
              t0  is indispensable and must be in Kelvin (> 0).
              p0  is indispensable and must be in Pascal (> 0).

            With compressibility:
            ---------------------
              ro0 is not useful, stricto sensu; nonetheless, as experience
                  shows that users often use this variable, it is required
                  to assign to it a strictly positive value (for example,
                  an initial value).
              viscl0 is useful and represents the molecular dynamic viscosity,
                  when it is constant, or a value which will be used during
                  initializations (or in inlet turbulence conditions,
                  depending on the user choice.
              cp0 is indispensable: it is the heat capacity, assumed constant
                  in the thermodynamics available by default
              t0  is indispensable and must be in Kelvin (> 0).
              p0  is indispensable and must be in Pascal (> 0).
                  With the thermodynamic law available by default,
                  t0 and p0 are used for the initialization of the density.
              xyzp0 is not useful because the pressure variable directly
                  represents the total pressure.
    */

    fp->ro0    = 1.17862;
    fp->viscl0 = 1.83337e-5;
    fp->cp0    = 1017.24;

    fp->t0 = 20. + 273.15;
    fp->p0 = 1.01325e5;

    /* We only specify XYZ0 if we explicitely fix Dirichlet conditions
       for the pressure. */

    fp->xyzp0[0] = 0.;
    fp->xyzp0[1] = 0.;
    fp->xyzp0[2] = 0.;


    /* irovar, ivivar, icp: constant or variable density,
                              viscosity/diffusivity, and specific heat

         When a specific physics module is active
           (coal, combustion, electric arcs, compressible: see usppmo)
           we MUST NOT set variables 'irovar', 'ivivar', and 'icp' here, as
           they are defined automatically.
         Nonetheless, for the compressible case, ivivar may be modified
           in the uscfx2 user subroutine.

         When no specific physics module is active, we may specify if the
           density, specific heat, and the molecular viscosity
           are constant (irovar=0, ivivar=0, icp=-1), which is the default
           or variable (irovar=1, ivivar=1, icp=0)

         For those properties we choose as variable, the corresponding law
           must be defined in cs_user_physical_properties
           (incs_user_physical_properties.f90);
           if they are constant, they take values ro0, viscl0, and cp0.

    */

    fp->irovar = 1;
    fp->ivivar = 1;
    fp->icp = -1;
  }
  /*! [param_fluid_properties] */

  /* Example: Change physical constants
   *------------------------------------*/

  /*! [param_physical_constants] */
  {
    cs_physical_constants_t *pc = cs_get_glob_physical_constants();

    pc->gravity[0] = 0.;
    pc->gravity[1] = 0.;
    pc->gravity[2] = -9.81; /* gravity  (m/s2) in the z direction */

  }

  /*! [param_physical_constants] */

  /* Example: Change options relative to the PISO-like sub-iterations
   * over prediction-correction.
   * - nterup: number of sub-iterations (default 1)
   * - epsup: relative precision (default 10e-5)
   *------------------------------------*/

  /*! [param_piso] */
  {
    cs_piso_t *piso = cs_get_glob_piso();
    piso->nterup = 3;
  }
  /*! [param_piso] */

  /* Example: activate the porous model
   * - 0 No porosity taken into account (Default)
   * - 1 Porosity taken into account
   *------------------------------------*/

  /*! [param_porous_model] */
  {
    cs_glob_porous_model = 1;
  }
  /*! [param_porous_model] */

  /* Example: log verbosity */
  /*------------------------*/

  /*! [param_log_verbosity] */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_var_cal_opt_t vcopt;
        int key_cal_opt_id = cs_field_key_id("var_cal_opt");

        cs_field_get_key_struct(f, key_cal_opt_id, &vcopt);
        vcopt.iwarni = 2;
        cs_field_set_key_struct(f, key_cal_opt_id, &vcopt);
      }
    }
  }
  /*! [param_log_verbosity] */

  /*Convective scheme

    blencv = 0 for upwind (order 1 in space, "stable but diffusive")
           = 1 for centered/second order (order 2 in space)
      we may use intermediate real values.
      Here we choose:
        for the velocity and user scalars:
          an upwind-centered scheme with 100% centering (blencv=1)
        for other variables
          the default code value (upwind standard, centered in LES)

    Specifically, for user scalars
      if we suspect an excessive level of numerical diffusion on
        a variable ivar representing a user scalar
        iscal (with ivar=isca(iscal)), it may be useful to set
        blencv = 1.0d0 to use a second-order scheme in space for
        convection. For temperature or enthalpy in particular, we
        may thus choose in this case:

        call field_get_key_struct_var_cal_opt(ivarfl(isca(iscalt)), vcopt)
        vcopt%blencv = 1.0d0
        call field_set_key_struct_var_cal_opt(ivarfl(isca(iscalt)), vcopt)

      For non-user scalars relative to specific physics
        implicitly defined by the model,
        the corresponding information is set automatically elsewhere:
        we do not modify blencv here. */

  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_var_cal_opt_t vcopt;
        int key_cal_opt_id = cs_field_key_id("var_cal_opt");

        cs_field_get_key_struct(f, key_cal_opt_id, &vcopt);
        vcopt.blencv = 1.;
        cs_field_set_key_struct(f, key_cal_opt_id, &vcopt);
      }
    }
  }

  /* Linear solver parameters (for each unknown)
     epsilo: relative precision for the solution of the linear system. */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_var_cal_opt_t vcopt;
        int key_cal_opt_id = cs_field_key_id("var_cal_opt");

        cs_field_get_key_struct(f, key_cal_opt_id, &vcopt);
        vcopt.epsilo = 1.e-6;
        cs_field_set_key_struct(f, key_cal_opt_id, &vcopt);
      }
    }
  }

  /* Dynamic reconstruction sweeps to handle non-orthogonlaities
     This parameter computes automatically a dynamic relax factor,
     and can be activated for any variable.
      - iswdyn = 1: means that the last increment is relaxed
      - iswdyn = 2: means that the last two increments are used to
                         relax
     NB: when iswdyn is greater than 1, then the number of
         non-orthogonality sweeps is increased to 20.*/
  {

    cs_var_cal_opt_t vcopt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");

    cs_field_get_key_struct(CS_F_(p), key_cal_opt_id, &vcopt);
    vcopt.iswdyn = 2;
    cs_field_set_key_struct(CS_F_(p), key_cal_opt_id, &vcopt);
  }

  /* Stabilization in turbulent regime

    For difficult cases, a stabilization may be obtained by not
    reconstructing the convective and diffusive flux for variables
    of the turbulence model, that is for k-epsilon models:
    */
  {

    cs_var_cal_opt_t vcopt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");

    cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &vcopt);
    vcopt.ircflu = 0;
    cs_field_set_key_struct(CS_F_(k), key_cal_opt_id, &vcopt);

    cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &vcopt);
    vcopt.ircflu = 0;
    cs_field_set_key_struct(CS_F_(eps), key_cal_opt_id, &vcopt);
  }


  /* Example: choose a convective scheme and
   * a limiter for a given variable (user and non-user) */
  /*----------------------------------------------*/

  /*! [param_var_limiter_choice] */
  {
    /* retrieve scalar field by its name */
    cs_field_t *sca1 = cs_field_by_name("scalar1");

    /* ischcv is the type of convective scheme:
       0: second order linear upwind
       1: centered
       2: pure upwind gradient in SOLU */

    /* isstpc:
      0: swich on the slope test
      1: swich off the slope test (default)
      2: continuous limiter ensuring boundedness (beta limiter)
      3: NVD/TVD Scheme */

    cs_var_cal_opt_t vcopt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");

    cs_field_get_key_struct(sca1, key_cal_opt_id, &vcopt);
    vcopt.ischcv = 0;
    vcopt.isstpc = 3;
    cs_field_set_key_struct(sca1, key_cal_opt_id, &vcopt);

    /* Min/Max limiter or NVD/TVD limiters
     * then "limiter_choice" keyword must be set:
     *   0: Gamma
     *   1: SMART
     *   2: CUBISTA
     *   3: SUPERBEE
     *   4: MUSCL
     *   5: MINMOD
     *   6: CLAM
     *   7: STOIC
     *   8: OSHER
     *   9: WASEB
     *   --- VOF scheme ---
     *   10: M-HRIC
     *   11: M-CICSAM       */

    int key_lim_id = cs_field_key_id("limiter_choice");
    cs_field_set_key_int(sca1, key_lim_id, CS_NVD_SUPERBEE);

    /* Get the Key for the Sup and Inf for the convective scheme */
    int kccmin = cs_field_key_id("min_scalar");
    int kccmax = cs_field_key_id("max_scalar");

    /* Set the Value for the Sup and Inf of the studied scalar
     * for the Gamma beta limiter for the temperature */
    cs_field_set_key_double(CS_F_(t), kccmin, 0.);
    cs_field_set_key_double(CS_F_(t), kccmax, 1.);

  }
  /*! [param_var_limiter_choice] */

  /* Example: Minimum and maximum admissible values for each USER scalar
   * Results are clipped at the end of each time step.
   *
   * If min > max, we do not clip
   *
   * For a scalar jj representing the variance of another, we may
   * abstain from defining these values
   * (a default clipping is set in place).
   * This is the purpose of the test on iscavr(jj) in the example below.
   *
   * For non-user scalars relative to specific physics (coal, combustion,
   * electric arcs: see usppmo) implicitly defined according to the
   * model, the information is automatically set elsewhere: we
   * do not set min or max values here. */

  {

    /* We define the min and max bounds */
    cs_field_set_key_double(CS_F_(t),
                            cs_field_key_id("min_scalar_clipping"),
                            0.);
    cs_field_set_key_double(CS_F_(t),
                            cs_field_key_id("max_scalar_clipping"),
                            1.);

  }


  /*-----------------------------------------------------------------------*/

  /*! [param_var_blend_st] */
  {
    /* retrieve scalar field by its name */
    cs_field_t *sca1 = cs_field_by_name("scalar1");

    /* blend_st (can be between 0 and 1):
      0: full upwind (default)
      1: scheme without upwind */

    cs_var_cal_opt_t vcopt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");

    cs_field_get_key_struct(sca1, key_cal_opt_id, &vcopt);
    vcopt.blend_st = 0.1;
    cs_field_set_key_struct(sca1, key_cal_opt_id, &vcopt);

  }
  /*! [param_var_blend_st] */


  /* Example: declare a scalar as buoyant so that it is
   * included in the velocity pressure PISO loop  */
  /*----------------------------------------------*/

  /*! [param_var_is_buoyant] */
  {
    /* retrieve scalar field by its name */
    cs_field_t *sca1 = cs_field_by_name("scalar1");

    int key_is_buoyant = cs_field_key_id("is_buoyant");

    cs_field_set_key_int(sca1, key_is_buoyant, 1);

  }
  /*! [param_var_is_buoyant] */


  /* Example: Scalar with a drift (key work "drift_scalar_model">0) or without drift
     ((key work "drift_scalar_model"=0, default option) for each USER scalar.
       - to specify that a scalar have a drift and need the drift computation:
         drift |= CS_DRIFT_SCALAR_ADD_DRIFT_FLUX

     Then, for each scalar with a drift, add a flag to specify if
     specific terms have to be taken into account:
       - thermophoresis terms:
         drift |= CS_DRIFT_SCALAR_THERMOPHORESIS
       - turbophoresis terms:
         drift |= CS_DRIFT_SCALAR_TURBOPHORESIS
       - centrifugal force terms:
         drift |= CS_DRIFT_SCALAR_CENTRIFUGALFORCE

      To kill boundary fluxes:
         drift |= CS_DRIFT_SCALAR_ZERO_BNDY_FLUX
   */

  /*! [param_var_drift] */
  {
    /* retrieve scalar field by its name */
    cs_field_t *sca1 = cs_field_by_name("scalar1");

    /* Key id for drift scalar */
    int key_drift = cs_field_key_id("drift_scalar_model");

    int drift = CS_DRIFT_SCALAR_ON;
    drift |= CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

    if (false)
      drift |= CS_DRIFT_SCALAR_THERMOPHORESIS;

    if (false)
      drift |= CS_DRIFT_SCALAR_TURBOPHORESIS;

    if (false)
      drift |= CS_DRIFT_SCALAR_CENTRIFUGALFORCE;

    if (false)
      drift |= CS_DRIFT_SCALAR_ZERO_BNDY_FLUX;

    /* Set the key word "drift_scalar_model" into the field structure */
    cs_field_set_key_int(sca1, key_drift, drift);

  }
  /*! [param_var_drift] */

  /* Example: activate mesh robustness options */
  /*-------------------------------------------*/

  /*! [mesh_tag_bad_cells_correction] */

  cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_WARPED_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_REGULARISATION;
  cs_glob_mesh_quantities_flag |= CS_CELL_FACE_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_CELL_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_FACE_DISTANCE_CLIP;
  cs_glob_mesh_quantities_flag |= CS_FACE_RECONSTRUCTION_CLIP;
  cs_glob_mesh_quantities_flag |= CS_CELL_VOLUME_RATIO_CORRECTION;

  /*! [mesh_tag_bad_cells_correction] */


  /* Postprocessing-related fields
     ============================= */

  /* Example: enforce existence of 'yplus', 'tplus' and 'tstar' fields, so that
              yplus may be saved, or a local Nusselt number may be computed using
              the post_boundary_nusselt subroutine.
              When postprocessing of these quantities is activated, those fields
              are present, but if we need to compute them in the
              cs_user_extra_operations user subroutine without postprocessing them,
              forcing the definition of these fields to save the values computed
              for the boundary layer is necessary. */

  /*! [param_force_yplus] */
  {
    cs_field_t *f;

    f = cs_field_by_name_try("yplus");
    if (f != NULL)
      cs_parameters_add_property("yplus",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    f = cs_field_by_name_try("tplus");
    if (f != NULL)
      cs_parameters_add_property("tplus",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

    f = cs_field_by_name_try("tstar");
    if (f != NULL)
      cs_parameters_add_property("tstar",
                                 1,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);
  }
  /*! [param_force_yplus] */

  /*--------------------------------------------------------------------------*/

  /* Example: add field to post-process the predicted-velocity divergence
   * and the pressure gradient in the momentum equation.
   */

  cs_parameters_add_property("predicted_vel_divergence",
                             1,
                             CS_MESH_LOCATION_CELLS);

  cs_parameters_add_property("pressure_gradient",
                             3,
                             CS_MESH_LOCATION_CELLS);
  /*--------------------------------------------------------------------------*/

  /* Example: add field to post-process the Reynolds stress production tensor
   * (for DRSM models only)
   */

  cs_parameters_add_property("rij_production",
                             6,
                             CS_MESH_LOCATION_CELLS);


  /* Example: force presence of boundary temperature field */
  /*-------------------------------------------------------*/

  /*! [param_force_b_temperature] */
  {
    cs_parameters_add_boundary_temperature();
  }

  /*! [param_force_b_temperature] */

  /* Example: add boundary values for all scalars */
  /*----------------------------------------------*/

  /*! [param_var_boundary_vals_1] */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE)
        cs_parameters_add_boundary_values(f);

    }
  }
  /*! [param_var_boundary_vals_1] */

  /* Example: add handling (storage) of previous values for a field */
  /*-------------------------------------------------------------------*/

  /*! [user_field_n_time_vals] */
  {
    /* add previous values handling for boundary temperature */
    cs_field_set_n_time_vals(CS_F_(t_b), 2);
  }
  /*! [user_field_n_time_vals] */

  /* Example: post-process clippings for slope test */
  /*------------------------------------------------*/

  /*! [param_post_slop_test] */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE)
        cs_field_set_key_int(f, cs_field_key_id("slope_test_upwind_id"), 0);
    }
  }
  /*! [param_post_slop_test] */

  /* Example: post-process clippings for Rij tensor */
  /*------------------------------------------------*/

  /*! [param_var_rij_clipping] */

  cs_field_set_key_int(CS_F_(rij), cs_field_key_id("clipping_id"), 1);
  cs_field_set_key_int(CS_F_(eps), cs_field_key_id("clipping_id"), 1);

  /*! [param_var_rij_clipping] */

  /* Example: post-process the Q-criterion on the whole domain mesh */
  /*----------------------------------------------------------------*/

  /*! [param_var_q_criterion] */

  cs_glob_post_util_flag[CS_POST_UTIL_Q_CRITERION] = 1;

  /*! [param_var_q_criterion] */

  /* Example: homogeneous mixture physical properties */
  /*--------------------------------------------------*/

  /*! [phprop] */

  cs_vof_parameters_t *vof_param = cs_get_glob_vof_parameters();

  /* Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
     liquid phase */

  vof_param->rho1 = 1.e3;
  vof_param->mu1 = 1.e-3;

  /* Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
     gas phase */

  vof_param->rho2 = 1.;
  vof_param->mu2 = 1.e-5;

  /*! [phprop] */

  /* Example: retrieve cavitation parameters structure */
  /*---------------------------------------------------*/

  /*! [cavit_param] */

  cs_cavitation_parameters_t *cavit_param =
    cs_get_glob_cavitation_parameters();

  /*! [cavit_param] */

  /* Example: Model parameters of the vaporization term (Merkle model) */
  /*-------------------------------------------------------------------*/

  /* Reference saturation pressure in kg/(m s2) */

  /*! [presat] */
  cavit_param->presat = 2.e3;
  /*! [presat] */

  /* Reference length, in meters, and velocity scales, in m/s, of the flow */

  /*! [scales_inf] */
  cavit_param->linf = 0.1;
  cavit_param->uinf = 1.;
  /*! [scales_inf] */

  /* Example: Interaction with turbulence */
  /*--------------------------------------*/

  /* Eddy-viscosity correction (Reboud et al. correction)
     0: deactivated
     1: activated */

  /*! [reboud_activ] */
  cavit_param->icvevm = 1;
  /*! [reboud_activ] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void)
{
  /* Example: define coupling between one volume zone and the rest of the
     mesh; this will automatically transform the selection boundaries
     to actual mesh boundaries.
     --------------------------------------------------------------------*/

  /*! [param_internal_coupling_add_volume] */

  cs_internal_coupling_add_volume(NULL,
                                  "x<.5"); /* Solid volume criterion */

  /*! [param_internal_coupling_add_volume] */

  /* Example: define coupling along an existing mesh boundary.
     ---------------------------------------------------------*/

  /*! [param_internal_coupling_add] */

  cs_internal_coupling_add(NULL,
                           "solid_volume_criterion",
                           "interface_criterion");

  /*! [param_internal_coupling_add] */

  /* Example: couple field whose name is "scalar1"
     ---------------------------------------------*/

  /*! [param_internal_coupling] */

  int f_id = cs_field_id_by_name("scalar1");

  cs_internal_coupling_add_entity(f_id);  /* Field to be coupled */

  /*! [param_internal_coupling] */
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Define or modify log user parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t     *domain)
{

  /* Frequency of log output */
  /*! [setup_log] */
  cs_glob_log_frequency = 1;
  /*! [setup_log] */

  /* Change a property's label
     (here for specific heat, first checking if it is variable) */
  /*! [setup_label] */
  if (CS_F_(cp) != NULL)
    cs_field_set_key_str(CS_F_(cp), cs_field_key_id("label"), "Cp");
  /*! [setup_label] */

  /* Probes for variables and properties
   * (example for velocity) */
  /*! [setup_post] */
  cs_field_set_key_int_bits(CS_F_(vel),
                            cs_field_key_id("post_vis"),
                            CS_POST_MONITOR);
  /*! [setup_post] */

  /* Probes for Radiative Transfer
   * (Luminance and radiative density flux vector) */

  /*! [setup_post_lum] */
   cs_field_t *f = cs_field_by_name_try("luminance");
   if (f != NULL)
     cs_field_set_key_int_bits(f,
                               cs_field_key_id("post_vis"),
                               CS_POST_MONITOR);

   f = cs_field_by_name_try("radiative_flux");
   if (f != NULL)
     cs_field_set_key_int_bits(f,
                               cs_field_key_id("post_vis"),
                               CS_POST_MONITOR);
  /*! [setup_post_lum] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
