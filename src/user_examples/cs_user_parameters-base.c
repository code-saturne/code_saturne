/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
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
 * See \ref parameters for examples.
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

  /* Activate Atmospheric flow model
   *  CS_ATMO_OFF:                module not activated
   *  CS_ATMO_CONSTANT_DENSITY:   standard modelling
   *  CS_ATMO_DRY:                dry atmosphere
   *  CS_ATMO_HUMID:              humid atmosphere (experimental)
   * */

  /*! [atmo_user_model_1] */
  cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_DRY;

  cs_atmo_set_meteo_file_name("meteo");

  /*--------------------------------------------------------------------------*/

  /* Atmospheric module options
   */

  /*! [atmo_module] */

  /*  Microphysics parameterization options */

  /* Option for nucleation for humid atmosphere
   *  0: without nucleation
   *  1: Pruppacher and Klett 1997
   *  2: Cohard et al. 1998,1999
   *  3: Abdul-Razzak et al. 1998,2000
   *  logarithmic standard deviation of the log-normal law of the droplet spectrum
   */
  cs_glob_atmo_option->nucleation_model = 3;

  /* Option for liquid water content distribution models
   *  1: all or nothing
   *  2: Gaussian distribution
   */
  cs_glob_atmo_option->distribution_model = 1;

  /*  Option for subgrid models
   *   0: the simplest parameterization (for numerical verifications)
   *   1: Bechtold et al. 1995 (Luc Musson-Genon)
   *   2: Bouzereau et al. 2004
   *   3: Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria and
   *                 Deardorff 1977
   */
  cs_glob_atmo_option->subgrid_model = 0;

  /* Sedimentation flag */
  cs_glob_atmo_option->sedimentation_model = 1;

  /* Deposition flag */
  cs_glob_atmo_option->deposition_model = 1;

  /* Read the meteo file (1) or impose directly the input values to compute it
   * in code_saturne (2) */
  cs_glob_atmo_option->meteo_profile = 2;

  /* Inverse LMO length (m^-1) */
  cs_glob_atmo_option->meteo_dlmo = 0.;
  /* Large scale roughness (m) */
  cs_glob_atmo_option->meteo_z0 = 0.1;
  /* Elevation for reference velocity (m) */
  cs_glob_atmo_option->meteo_zref = 10.;
  /* Friction velocity (m/s) */
  cs_glob_atmo_option->meteo_ustar0 = 1.;
  /* Velocity direction (degrees from North) */
  cs_glob_atmo_option->meteo_angle = 270.;
  /* Temperature at 2m (K) */
  cs_glob_atmo_option->meteo_t0 = 288.;
  /* Pressure at sea level (Pa) */
  cs_glob_atmo_option->meteo_psea = 1.01325e5;

  /* Option to compute ground elevation in the domain */
  cs_glob_atmo_option->compute_z_ground = true;

  /* Automatic open boundary conditions
   *   1: meteo mass flow rate is imposed with a constant large scale
   *      pressure gradient
   *   2: same plus velocity profile imposed at ingoing faces
   */
  cs_glob_atmo_option->open_bcs_treatment = 1;

  /* Time of the simulation (for radiative model or chemistry)
   * syear:  starting year
   * squant: starting quantile
   * shour:  starting hour (UTC)
   * smin:   starting minute
   * ssec:   starting second
   */
  cs_glob_atmo_option->syear = 2020;
  cs_glob_atmo_option->squant = 1;
  cs_glob_atmo_option->shour = 1;
  cs_glob_atmo_option->smin = 0;
  cs_glob_atmo_option->ssec = 0.;

  /* Geographic position
   *  longitude: longitude of the domain origin
   *  latitude: latitude of the domain origin
   */
  cs_glob_atmo_option->longitude = 0.;
  cs_glob_atmo_option->latitude = 45.0;

  /* Chemistry:
   *   model: choice of chemistry resolution scheme
   *     0: no atmospheric chemistry
   *     1: quasi steady equilibrium NOx scheme with 4 species and 5 reactions
   *     2: scheme with 20 species and 34 reactions
   *     3: scheme CB05 with 52 species and 155 reactions
   *     4: user defined scheme
   *      for model = 4, a SPACK file must be provided using
   *        cs_atmo_chemistry_set_spack_file_name("species.spack.dat")
   *      the following sources generated by SPACK should be included in the SRC folder
   *        kinetic.f90, fexchem.f90, jacdchemdc.f90, rates.f90, dratedc.f90
   *        dimensions.f90, LU_decompose.f90, LU_solve.f90
   */
   cs_glob_atmo_chemistry->model = 0;

   /* Default file for the chemistry profile is "chemistry" */
   cs_atmo_set_chem_conc_file_name("chem_01_01_2000");

   /* Chemistry with photolysis: inclusion (true) or not (false) of photolysis reactions
    * warning: photolysis is not compatible with space-variable time step
    */
   cs_glob_atmo_chemistry->chemistry_with_photolysis = true;

   /* Aerosol chemistry
    * -----------------*/

   /*   aerosol_model: flag to activate aerosol chemistry
    *   CS_ATMO_AEROSOL_OFF: aerosol chemistry deactivated, default
    *   CS_ATMO_AEROSOL_SSH: model automatically set to 4
    *     External library SSH-aerosol is used
    *     Corresponding SPACK files must be provided
    *     The SSH namelist file can be specified using
    *       call atmo_chemistry_set_aerosol_file_name("namelist_coag.ssh")
    *       if no namelist file is specified, "namelist.ssh" is used
    */
   cs_glob_atmo_chemistry->aerosol_model = CS_ATMO_AEROSOL_SSH;

   /* Default file for the aerosol profile is "aerosols" */
   cs_atmo_set_aero_conc_file_name("aero_01_01_2001");

   /* Frozen gaseous chemistry
    *   false: gaseous chemistry is activated (default)
    *   true: gaseous chemistry is frozen
    */
   cs_glob_atmo_chemistry->frozen_gas_chem = false;

  /*! [atmo_module] */


  /*! [atmo_user_model_1] */

  /*--------------------------------------------------------------------------*/

  /* Activate compressible model */

  /*! [compressible_user_model_1] */
  /* -1: not active, 0: activated, 1: barotropic version,
     2: homogeneous two phase model, 3: by pressure increment */
  cs_glob_physical_model_flag[CS_COMPRESSIBLE] = -1;
  /*! [compressible_user_model_1] */

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

  /* Advanced choice of Wall function:
   *  Wall function for the dynamic, iwallf can be:
   *   CS_WALL_F_DISABLED: no wall function
   *   CS_WALL_F_1SCALE_POWER: deprecated 1 velocity scale power law
   *   CS_WALL_F_1SCALE_LOG: 1 velocity scale log law
   *   CS_WALL_F_2SCALES_LOG: 2 velocity scales log law (default for many models)
   *   CS_WALL_F_SCALABLE_2SCALES_LOG: Scalable log wll function with two
   *                                   velocity scales
   *   CS_WALL_F_2SCALES_VDRIEST: 2 velocity scales with Van Driest damping
   *   CS_WALL_F_2SCALES_SMOOTH_ROUGH: 2 velocity scales valid for smooth and
   *                                   rough regimes
   *   CS_WALL_F_2SCALES_CONTINUOUS: 2 velcoity scales (all y+), default for
   *                                 low Reynolds turbulence models.
   */
  {
    cs_wall_functions_t *wf = cs_get_glob_wall_functions();
     wf->iwallf = CS_WALL_F_2SCALES_VDRIEST;
  }

  /*--------------------------------------------------------------------------*/

  /* Advanced choice of scalars wall function:
   *   CS_WALL_F_S_ARPACI_LARSEN: 3 layer Arpaci and Larsen wall function
   *                              (default for most of the models)
   *   CS_WALL_F_S_VDRIEST: Van Driest wall function
   *   CS_WALL_F_S_LOUIS: Atmospheric Louis wall function
   *   CS_WALL_F_S_MONIN_OBUKHOV: Monin Obukhov atmospheric wall function
   *   CS_WALL_F_S_SMOOTH_ROUGH: smooth rough wall function
   * */
  {
    cs_wall_functions_t *wf = cs_get_glob_wall_functions();
     wf->iwalfs = CS_WALL_F_S_MONIN_OBUKHOV;
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
  rans_model->idirsm = 1;

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
  thermal_model->itherm = CS_THERMAL_MODEL_TEMPERATURE;

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

  /*! [gax_mix_activation] */

  cs_gas_mix_type_t gas_mix_type          = CS_GAS_MIX_AIR_HELIUM;
  cs_glob_physical_model_flag[CS_GAS_MIX] = gas_mix_type;

  /*! [gax_mix_activation] */

  /*! [wall_condensation] */

  /* Activated wall condensation model for natural convection
   *   CS_WALL_COND_MODEL_COPAIN : Legacy implementation of COPAIN correlation
   * (default) CS_WALL_COND_MODEL_COPAIN_BD : Update of COPAIN correlation from
   * Benteboula and Dabbene CS_WALL_COND_MODEL_UCHIDA : Correlation of Uchida
   *   CS_WALL_COND_MODEL_DEHBI  : Correlation of Dehbi
   */

  cs_wall_cond_t *wall_cond = cs_get_glob_wall_cond();
  wall_cond->icondb         = 0; // activate wall codnensation
  wall_cond->natural_conv_model
    = CS_WALL_COND_MODEL_DEHBI; // choose correlation

  /*! [wall_condensation] */

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
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t *domain)
{
  /*! [time_stepping_options] */

  /* Time step type */

  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  time_opt->idtvar = CS_TIME_STEP_CONSTANT;

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

  /*! [param_vp_arak] */
  {
    cs_velocity_pressure_param_t *vp_param
      = cs_get_glob_velocity_pressure_param();
    /* Switch off Rhie & Chow filter */
    vp_param->arak = 0.;

    /* Activate RT0 reconstruction for the velocity from
     * mass fluxes */
    vp_param->irevmc = 1;
  }
  /*! [param_vp_arak] */

  /* Example: change Reference fluid properties options */
  /*----------------------------------------------------*/

  /* Members of the structure cs_fluid_properties_t
   */

  /*! [param_fluid_properties] */
  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

    /*
      ro0        : density in kg/m3
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

      In the standard case (no combustion, electric arcs, compressibility):
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

  /* Example: Change options relative to the inner iterations
   * over prediction-correction.
   * - nterup: number of sub-iterations (default 1)
   * - epsup: relative precision (default 10e-5)
   *------------------------------------*/

  /*! [param_vp_netrup] */
  {
    cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
    vp_param->nterup = 3;
  }
  /*! [param_vp_netrup] */

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
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->verbosity = 2;
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
        blencv = 1.0to use a second-order scheme in space for
        convection. For temperature or enthalpy in particular, we
        may thus choose in this case:

        cs_field_t *f = cs_thermal_model_field();
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->blencv = 1.;

      For non-user scalars relative to specific physics
        implicitly defined by the model,
        the corresponding information is set automatically elsewhere:
        we do not modify blencv here. */

  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->blencv = 1.;
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
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->epsilo = 1.e-6;
      }
    }
  }

  /* Dynamic reconstruction sweeps to handle non-orthogonlaities
     This parameter computes automatically a dynamic relax factor,
     and can be activated for any variable.
      - iswdyn = 0: no relaxation
      - iswdyn = 1: means that the last increment is relaxed
      - iswdyn = 2: (default) means that the last two increments are used
                    to relax.
  */
  {
    cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(p));
    eqp->iswdyn = 2;
  }

  /* Stabilization in turbulent regime

    For difficult cases, a stabilization may be obtained by not
    reconstructing the convective and diffusive flux for variables
    of the turbulence model, that is for k-epsilon models:
    */
  {
    cs_equation_param_t *eqp;

    eqp = cs_field_get_equation_param(CS_F_(k));
    eqp->ircflu = 0;

    eqp = cs_field_get_equation_param(CS_F_(eps));
    eqp->ircflu = 0;
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
       2: pure upwind gradient in SOLU
       3: blending SOLU and centered
       4: NVD/TVD Scheme */

    /* isstpc:
      0: slope test enabled
      1: slope test disabled (default)
      2: continuous limiter ensuring boundedness (beta limiter) enabled */

    cs_equation_param_t *eqp;

    eqp = cs_field_get_equation_param(sca1);
    eqp->ischcv = 1;
    eqp->isstpc = 0;

    /* Min/Max limiter or NVD/TVD limiters
     * then "limiter_choice" keyword must be set:
     *   CS_NVD_Gamma
     *   CS_NVD_SMART
     *   CS_NVD_CUBISTA
     *   CS_NVD_SUPERBEE
     *   CS_NVD_MUSCL
     *   CS_NVD_MINMOD
     *   CS_NVD_CLAM
     *   CS_NVD_STOIC
     *   CS_NVD_OSHER
     *   CS_NVD_WASEB
     *   --- VOF scheme ---
     *   CS_NVD_VOF_HRIC
     *   CS_NVD_VOF_CICSAM
     *   CS_NVD_VOF_STACS        */

    int key_lim_id = cs_field_key_id("limiter_choice");
    cs_field_set_key_int(sca1, key_lim_id, CS_NVD_VOF_CICSAM);

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

    cs_equation_param_t *eqp = cs_field_get_equation_param(sca1);
    eqp->blend_st = 0.1;

  }
  /*! [param_var_blend_st] */

  /* Example: declare a scalar as buoyant so that it is
   * included in the velocity pressure inner loop  */
  /*----------------------------------------------*/

  /*! [param_var_is_buoyant] */
  {
    /* retrieve scalar field by its name */
    cs_field_t *sca1 = cs_field_by_name("scalar1");

    int key_is_buoyant = cs_field_key_id("is_buoyant");

    cs_field_set_key_int(sca1, key_is_buoyant, 1);
  }
  /*! [param_var_is_buoyant] */

  /* Example: Scalar with a drift (key work "drift_scalar_model">0)
              or without drift
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
              cs_user_extra_operations user subroutine without postprocessing
              them, forcing the definition of these fields to save the
              values computed for the boundary layer is necessary. */

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

  /* Example: Change the number of checkpoint files which are saved. */
  /*-----------------------------------------------------------------*/

  /* By default, only one checkpoint file (last, end of run) is kept.
   * If it is wished to keep more (debugging purposes for example),
   * it can be chaned.
   * Herebelow is an example showing how to keep the last two files:
   */

  /*! [change_nsave_checkpoint_files] */
  cs_restart_set_n_max_checkpoints(2);
  /*! [change_nsave_checkpoint_files] */

  /* Example: compute porosity from a scan of points
   * ------------------------------------------------*/

  cs_porosity_from_scan_set_file_name("chbre_chbre33.pts");

  /* Apply a transformation to the scanned points */
  /* Translation part */
  cs_glob_porosity_from_scan_opt->transformation_matrix[0][3] = 4.;
  cs_glob_porosity_from_scan_opt->transformation_matrix[1][3] = 1.98;
  cs_glob_porosity_from_scan_opt->transformation_matrix[2][3] = 37.5477;
  /* Rotation part arround z axis */
  cs_real_t angle = 15. /180. * cs_math_pi;
  cs_glob_porosity_from_scan_opt->transformation_matrix[0][0] =  cos(angle);
  cs_glob_porosity_from_scan_opt->transformation_matrix[0][1] =  sin(angle);
  cs_glob_porosity_from_scan_opt->transformation_matrix[1][0] = -sin(angle);
  cs_glob_porosity_from_scan_opt->transformation_matrix[1][1] =  cos(angle);
  cs_glob_porosity_from_scan_opt->transformation_matrix[2][2] = 1.;

  /* Add some sources to fill fluid space */
  {
    cs_real_3_t source = {4.295, 1.15326, 0.5};
    /* If a transformation matrix has been applied
     * chose if if has to be applied to the source */
    bool transform = true ;
    cs_porosity_from_scan_add_source(source, transform);
  }
  {
    cs_real_3_t source = {4.295, 3.2, 0.5};
    /* If a transformation matrix has been applied
     * chose if if has to be applied to the source */
    bool transform = true ;
    cs_porosity_from_scan_add_source(source, transform);
  }

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

  const cs_zone_t *sz = cs_volume_zone_by_name("solid");
  cs_internal_coupling_add_volume_zone(sz);

  /* Activate fluid-solid mode to kill dynamics in the solid */
  cs_velocity_pressure_model_t *vp_model
    = cs_get_glob_velocity_pressure_model();
  vp_model->fluid_solid = true;

  /*! [param_internal_coupling_add_volume] */

  /* Example: define coupling along an existing mesh boundary.
     ---------------------------------------------------------*/

  /*! [param_internal_coupling_add] */

  cs_internal_coupling_add("solid_volume_criterion",
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
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t     *domain)
{
  CS_UNUSED(domain);

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
