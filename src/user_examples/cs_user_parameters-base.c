/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
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
  /*! [time_stepping_options] */

  /* Time step type */

  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  time_opt->idtvar = CS_TIME_STEP_CONSTANT;

  /*! [time_stepping_options] */

  /*! [activate_user_model] */

  /* Activate Atmospheric flow model
   *  CS_ATMO_OFF:                module not activated
   *  CS_ATMO_CONSTANT_DENSITY:   standard modelling
   *  CS_ATMO_DRY:                dry atmosphere
   *  CS_ATMO_HUMID:              humid atmosphere (experimental)
   * */
  cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_DRY;

  /* Activate compressible model
   *  -1: not active, 0: activated, 1: barotropic version,
   *   2: homogeneous two phase model, 3: by pressure increment */
  cs_glob_physical_model_flag[CS_COMPRESSIBLE] = 3;

  /* Activate Eddy Break Up pre-mixed flame combustion model
   * -1: not active
   *  0: adiabatic conditions at constant richness
   *  1: permeatic conditions at constant richness
   *  2: adiabatic conditions at variable richness
   *  3: permeatic conditions at variable richness */
  cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = -1;

  /* Activate 3-pointcombustion model
   * -1: not active, 0: adiabatic, 1: permeatic */
  cs_glob_physical_model_flag[CS_COMBUSTION_3PT] = -1;

  /* Activate Libby-Williams pre-mixed flame combustion model
   * -1: not active
   *  0: two peak model with: adiabiatic conditions
   *  1: two peak model with: permeatic conditions
   *  2: three peak model: with adiabiatic conditions
   *  3: three peak model: with permeatic conditions
   *  4: four peak model with: adiabiatic conditions
   *  5: four peak model with: permeatic conditions*/
  cs_glob_physical_model_flag[CS_COMBUSTION_LW] = -1;

  /* Activate pulverized coal combustion model
   * -1: not active
   *  0: active
   *  1: with drying */
  cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = 1;

  /* Activate the drift (for combustion)
   * 0 (no activation),
   * 1 (transported particle velocity)
   * 2 (limit drop particle velocity) */
  cs_glob_combustion_model->idrift = 1;

  /* Cooling towers model
   * -1: not active
   *  1: Poppe's model
   *  2: Merkel's model */
  cs_glob_physical_model_flag[CS_COOLING_TOWERS] = 1;

  /* Electric arcs model
   * -1: not active
   *  2: electric potential and vector potential */
  cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] = 2;

  /* Joule effcets model
   * -1: not active
   *  1: real potential
   *  2: complex potential
   *  3: real potential + Transfo
   *  4: complex potential + Transfo */
  cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 2;

  /* Radiative transfer module
   * CS_RAD_TRANSFER_NONE: not active
   * CS_RAD_TRANSFER_DOM: discrete ordinates model
   * CS_RAD_TRANSFER_P1: P1 model */
  cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_DOM;

  /*  Activate gas mix model
   *  CS_GAS_MIX_OFF                      gas mix model off
   *  CS_GAS_MIX_AIR_HELIUM               air/helium, helium deduced
   *  CS_GAS_MIX_AIR_HYDROGEN             air/hydrogen, hydrogen deduced
   *  CS_GAS_MIX_AIR_STEAM                air/steam, steam deduced
   *  CS_GAS_MIX_AIR_HELIUM_STEAM         helium/steam, steam deduced
   *  CS_GAS_MIX_AIR_HYDROGEN_STEAM       hydrogen/steam, steam deduced
   *  CS_GAS_MIX_HELIUM_AIR               helium/air, O2 from air deduced
   *  CS_GAS_MIX_USER                     user defined
   * */
  cs_glob_physical_model_flag[CS_GAS_MIX] = CS_GAS_MIX_AIR_HELIUM;

  /*! [activate_user_model] */

  /*! [atmo_user_model_1] */

  cs_atmo_set_meteo_file_name("meteo");

  /*--------------------------------------------------------------------------*/

  /* Atmospheric module options
   */

  /*! [atmo_module] */

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  /*  Microphysics parameterization options */

  /* Option for nucleation for humid atmosphere
   *  0: without nucleation
   *  1: Pruppacher and Klett 1997
   *  2: Cohard et al. 1998,1999
   *  3: Abdul-Razzak et al. 1998,2000
   *  logarithmic standard deviation of the log-normal law of the droplet
   *  spectrum
   */
  at_opt->nucleation_model = 3;

  /* Option for liquid water content distribution models
   *  1: all or nothing
   *  2: Gaussian distribution
   */
  at_opt->distribution_model = 1;

  /*  Option for subgrid models
   *   0: the simplest parameterization (for numerical verifications)
   *   1: Bechtold et al. 1995 (Luc Musson-Genon)
   *   2: Bouzereau et al. 2004
   *   3: Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria and
   *                 Deardorff 1977
   */
  at_opt->subgrid_model = 0;

  /* Sedimentation flag */
  at_opt->sedimentation_model = 1;

  /* Deposition flag */
  at_opt->deposition_model = 1;

  /* Read the meteo file (1) or impose directly the input values to compute it
   * in code_saturne (2) */
  at_opt->meteo_profile = 2;

  /* Advanced choice of universal functions among for stable
   *  - CS_ATMO_UNIV_FN_CHENG (default)
   *  - CS_ATMO_UNIV_FN_HOGSTROM
   *  - CS_ATMO_UNIV_FN_BUSINGER
   *  - CS_ATMO_UNIV_FN_HARTOGENSIS
   *  for unstable
   *  - CS_ATMO_UNIV_FN_HOGSTROM (default)
   *  - CS_ATMO_UNIV_FN_BUSINGER
   *
   * */

  /* Hartogensis for stable */
  at_opt->meteo_phim_s = CS_ATMO_UNIV_FN_HARTOGENSIS;
  at_opt->meteo_phih_s = CS_ATMO_UNIV_FN_HARTOGENSIS;

  /* Inverse LMO length (m^-1) */
  at_opt->meteo_dlmo = 0.;
  /* Large scale roughness (m) */
  at_opt->meteo_z0 = 0.1;
  /* Elevation for reference velocity (m) */
  at_opt->meteo_zref = 10.;
  /* Friction velocity (m/s) */
  at_opt->meteo_ustar0 = 1.;
  /* Velocity direction (degrees from North) */
  at_opt->meteo_angle = 270.;
  /* Temperature at 2m (K) */
  at_opt->meteo_t0 = 288.;
  /* Pressure at sea level (Pa) */
  at_opt->meteo_psea = 1.01325e5;

  /* Option to compute ground elevation in the domain */
  at_opt->compute_z_ground = true;

  /* Automatic open boundary conditions
   *   1: meteo mass flow rate is imposed with a constant large scale
   *      pressure gradient
   *   2: same plus velocity profile imposed at ingoing faces
   */
  at_opt->open_bcs_treatment = 1;

  /* Time of the simulation (for radiative model or chemistry)
   * syear:  starting year
   * squant: starting quantile
   * shour:  starting hour (UTC)
   * smin:   starting minute
   * ssec:   starting second
   */
  at_opt->syear = 2020;
  at_opt->squant = 1;
  at_opt->shour = 1;
  at_opt->smin = 0;
  at_opt->ssec = 0.;

  /* Geographic position
   *  longitude: longitude of the domain origin
   *  latitude: latitude of the domain origin
   */
  at_opt->longitude = 0.;
  at_opt->latitude = 45.0;

  /* Chemistry:
   *   model: choice of chemistry resolution scheme
   *     0: no atmospheric chemistry
   *     1: quasi steady equilibrium NOx scheme with 4 species and 5 reactions
   *     2: scheme with 20 species and 34 reactions
   *     3: scheme CB05 with 52 species and 155 reactions
   *     4: user defined scheme
   *      for model = 4, a SPACK file must be provided using
   *        cs_atmo_chemistry_set_spack_file_name("species.spack.dat")

   *      the following sources generated by SPACK should be included in the
   *      SRC folder
   *        kinetic.f90, fexchem.f90, jacdchemdc.f90, rates.f90, dratedc.f90
   *        dimensions.f90, LU_decompose.f90, LU_solve.f90
   */
  cs_glob_atmo_chemistry->model = 0;

  /* Default file for the chemistry profile is "chemistry" */
  cs_atmo_set_chem_conc_file_name("chem_01_01_2000");

  /* Chemistry with photolysis: inclusion (true) or not (false) of photolysis
   * reactions
   * warning: photolysis is not compatible with space-variable time step
   */
  cs_glob_atmo_chemistry->chemistry_with_photolysis = true;

  /* chemistry_sep_mode: split (=1) or semi-coupled (=2, pu-sun)
   * resolution of chemistry.
   * Split (=1) mandatory for aerosols.
   * Semi-coupled (=2) by default. */
  cs_glob_atmo_chemistry->chemistry_sep_mode = 1;

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

  /* Soil Atmosphere model
   * ---------------------*/

  /*! [atmo_soil_set] */
  at_opt->soil_model = 1; /* Switch on soil model */

  /* Set the number of predefined categories (+1 which is the default one)
   * among:
   *  - CS_ATMO_SOIL_5_CAT
   *  - CS_ATMO_SOIL_7_CAT
   * */
  at_opt->soil_cat= CS_ATMO_SOIL_5_CAT; /* Switch on soil model */

  /* Specify the boundary zone which is modeled */
  at_opt->soil_zone_id = cs_boundary_zone_by_name("Sol")->id;
  /*! [atmo_soil_set] */

  /* 1-D radiative transfer
   * ---------------------*/

  /* Activate 1-D radiative transfer model */
  at_opt->radiative_model_1d = 1;

  /* Specify the number of verticals and the number of levels.
   * The mesh levels can be specified in cs_user_parameters function.
   * */
  at_opt->rad_1d_nvert = 1;
  at_opt->rad_1d_nlevels = 50;
  at_opt->rad_1d_nlevels_max = at_opt->rad_1d_nlevels;

  /* Complete 1-D mesh to ztop in case of radiative transfer */
  if (at_opt->rad_1d_nvert > 0) {
    cs_real_t zvmax = 1975.;/* top of the domain */
    cs_real_t ztop = 11000.;/* top of the troposphere */
    for (cs_real_t zzmax = (((int) zvmax)/1000)*1000.;
        zzmax <= (ztop -1000);
        zzmax += 1000.) {
      (at_opt->rad_1d_nlevels_max)++;
    }

  }

  /*! [atmo_module] */

  /*! [atmo_user_model_1] */

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

  /*! [turb_rans_model_settings] */
  /* Coupled solver for Rij components (CS_TURB_RIJ_*)
   *   0: switch off
   *   1: switch on (default)
   */

  cs_turb_rans_model_t *rans_model = cs_get_glob_turb_rans_model();
  rans_model->irijco = 0;

  /* Advanced re-initialization for EBRSM or k-omega models
     - 0: switch off
     - 1: switch on (default)
     */

  rans_model->reinit_turb = 1;

  /* Turbulent diffusion model for second moment closure (CS_TURB_RIJ_*)
     - 0: scalar diffusivity (Shir model, default)
     - 1: tensorial diffusivity (Daly and Harlow model)
     */
  rans_model->idirsm = 1;

  /* Rotation/curvature correction for eddy-viscosity turbulence models
     (k-epsilon, k-omega) */

  rans_model->irccor = 1;

  /*! [turb_rans_model_settings] */

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
   * CS_THERMAL_MODEL_NONE        : none
   * CS_THERMAL_MODEL_TEMPERATURE : temperature
   * CS_THERMAL_MODEL_ENTHALPY    : enthalpy
   * CS_THERMAL_MODEL_TOTAL_ENERGY: total energy (only for compressible module)
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

  /* Number of iterations for fluid initialization.
     cs_glob_ale_n_ini_f is not an absolute iteration number, so for a
     restarted calculation, it corresponds to the number of iterations
     for fuid initialization relative to the first restart iteration.
     In general it should be set to 0 in that case. */

  cs_glob_ale_n_ini_f = 75;

  /* Maximum number of iterations in case of implicit Fluid Structure Coupling
     with structural calculations (internal and/or external
     (i.e. using code_aster).
     For an explicit FSI scheme, set cs_glob_mobile_structures_i_max = 1; */

  cs_glob_mobile_structures_i_max = 15;

  /* Relative precision of sub-cycling Fluid Structure Coupling algorithm */

  cs_glob_mobile_structures_i_eps = 1.e-5;

  /*! [ALE_activation] */

  /*--------------------------------------------------------------------------*/

  /*! [wall_condensation] */

  /* Activated wall condensation model for natural convection
   * - CS_WALL_COND_MODEL_COPAIN:     Legacy implementation of COPAIN
   *                                  correlation (default)
   * - CS_WALL_COND_MODEL_COPAIN_BD:  Update of COPAIN correlation from
   *                                  Benteboula and Dabbene
   * - CS_WALL_COND_MODEL_UCHIDA:     Uchida correlation
   * - CS_WALL_COND_MODEL_DEHBI:      Dehbi correlation
   */

  cs_wall_condensation_t *wall_cond = cs_get_glob_wall_condensation();
  wall_cond->icondb         = 0; /* activate wall condensation */
  wall_cond->natural_conv_model
    = CS_WALL_COND_MODEL_DEHBI;  /* choose correlation */

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

  /*--------------------------------------------------------------------------*/

  /* Physical constants / body forces */

  /*! [user_model_gravity] */
  {
    cs_physical_constants_t *pc = cs_get_glob_physical_constants();

    pc->gravity[0] = 0.;
    pc->gravity[1] = 0.;
    pc->gravity[2] = -9.81; /* gravity  (m/s2) in the z direction */
  }
  /*! [user_model_gravity] */

  /*! [user_model_coriolis] */
  {
    cs_physical_constants_t *pc = cs_get_glob_physical_constants();

    cs_real_t omega_r[3] = {0., 0., 1.};    /* rotation vector, in rad/s */
    cs_real_t invariant[3] = {0., 0., 0.};  /* invariant point */

    cs_rotation_define(omega_r[0], omega_r[1], omega_r[2],
                       invariant[0], invariant[1], invariant[2]);

    pc->icorio = 1;  /* take Coriolis source terms into account */
  }
  /*! [user_model_coriolis] */

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

  /* Add some sources from a file
   * The file contains lines with x,y,z (in meters) format */
  cs_ibm_add_sources_by_file_name("sources.csv");

  /* Example: setup options for radiative transfer
   * --------------------------------------------- */

  /*! [cs_user_radiative_transfer_parameters] */

  /* Local pointer to global parameters, for conciseness */

  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  /* indicate whether the radiation variables should be
     initialized (=0) or read from a restart file (=1) */

  rt_params->restart = (cs_restart_present()) ? 1 : 0;

  /* Update period of the radiation module */

  cs_time_control_init_by_time_step
    (&( rt_params->time_control),
     - 1,     /* nt_start */
     -1,      /* nt_end */
     5,       /* interval */
     true,    /* at start */
     false);  /* at end */

  /* Quadrature Sn (n(n+2) directions)

     1: S4 (24 directions)
     2: S6 (48 directions)
     3: S8 (80 directions)

     Quadrature Tn (8n^2 directions)

     4: T2 (32 directions)
     5: T4 (128 directions)
     6: Tn (8*ndirec^2 directions)
  */

  rt_params->i_quadrature = 4;

  /* Number of directions, only for Tn quadrature */
  rt_params->ndirec = 3;

  /* Method used to calculate the radiative source term:
     - 0: semi-analytic calculation (required with transparent media)
     - 1: conservative calculation
     - 2: semi-analytic calculation corrected
          in order to be globally conservative
     (If the medium is transparent, the choice has no effect) */

  rt_params->idiver = 2;

  /* Verbosity level in the log concerning the calculation of
     the wall temperatures (0, 1 or 2) */

  rt_params->iimpar = 1;

  /* Verbosity mode for the radiance (0, 1 or 2) */

  rt_params->verbosity = 1;

  /* Compute the absorption coefficient through a model (if different from 0),
     or use a constant absorption coefficient (if 0).
     Useful ONLY when gas or coal combustion is activated
     - imodak = 1: ADF model with 8 wave length intervals
     - imodak = 2: Magnussen et al. and Kent and Honnery models */

  rt_params->imodak = 2;

  /* Compute the absorption coefficient via ADF model
     Useful ONLY when coal combustion is activated
     imoadf = 0: switch off the ADF model
     imoadf = 1: switch on the ADF model (with 8 bands ADF08)
     imoadf = 2: switch on the ADF model (with 50 bands ADF50) */

  rt_params->imoadf = 1;

  /* Compute the absorption coefficient through FSCK model (if 1)
     Useful ONLY when coal combustion is activated
     imfsck = 1: activated
     imfsck = 0: not activated */

  rt_params->imfsck = 1;

  /* Activate  3D radiative models for  atmospheric flows
       atmo_model |=  CS_RAD_ATMO_3D_DIRECT_SOLAR: direct solar
       atmo_model |=  CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND: direct solar
       atmo_model |=  CS_RAD_ATMO_3D_DIFFUSE_SOLAR: diffuse solar
       atmo_model |=  CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND: diffuse solar
       atmo_model |=  CS_RAD_ATMO_3D_INFRARED: Infrared
  */

  rt_params->atmo_model |= CS_RAD_ATMO_3D_DIRECT_SOLAR;
  rt_params->atmo_model |= CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND;
  rt_params->atmo_model |= CS_RAD_ATMO_3D_DIFFUSE_SOLAR;
  rt_params->atmo_model |= CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND;
  rt_params->atmo_model |= CS_RAD_ATMO_3D_INFRARED;

  /*! [cs_user_radiative_transfer_parameters] */
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

  /* Calculation (restart) with frozen velocity field (1 yes, 0 no) */

  /*! [param_iccvfg] */
  {
    cs_time_scheme_t *t_sc = cs_get_glob_time_scheme();

    t_sc->iccvfg = 1;
  }
  /*! [param_iccvfg] */

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

    /* Change the hydrostatic pressure algorithm */
    vp_param->iphydr = 1;
  }
  /*! [param_vp_arak] */

  /* Example: change Reference fluid properties options */
  /*----------------------------------------------------*/

  /* Members of the structure cs_fluid_properties_t */

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

    /* Account the thermodynamical pressure variation in time
       (only if cs_glob_velocity_pressure_param->idilat = 3)

      By default:
      ----------
      - The thermodynamic pressure (pther) is initialized with p0 = p_atmos.
      - The maximum thermodynamic pressure (pthermax) is initialized with -1
        (no maximum by default, this term is used to model a venting effect when
        a positive value is given by the user).
      - A global leak can be set through a leakage surface sleak with a head
       loss kleak of 2.9 (Idelcick) */

    fp->ipthrm = 1;

    fp->pthermax= -1.;

    fp->sleak = 0.;
    fp->kleak = 2.9;
  }
  /*! [param_fluid_properties] */

  /*! [param_diffusivity_id] */
  {
    const int kivisl = cs_field_key_id("diffusivity_id");

    /* For thermal scalar */
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
      cs_field_set_key_int(cs_thermal_model_field(), kivisl, -1);
    else
      cs_field_set_key_int(cs_field_by_name("temperature"), kivisl, -1);

    /* For user-defined scalars */
    const int n_fields = cs_field_n_fields();
    const int k_scal = cs_field_key_id("scalar_id");
    const int kscavr = cs_field_key_id("first_moment_id");

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & (CS_FIELD_VARIABLE | CS_FIELD_USER)) {
        int s_num = cs_field_get_key_int(f, k_scal);
        int iscavr = cs_field_get_key_int(f, kscavr);
        if (s_num > 0 && iscavr <= 0)
          cs_field_set_key_int(f, kivisl, 0);
      }
    }
  }
  /*! [param_diffusivity_id] */

  /*! [param_density_id] */
  {
    const int kromsl = cs_field_key_id("density_id");

    /* Example user-defined scalar */
    cs_field_t *f = cs_field_by_name("scalar1");
    cs_field_set_key_int(f, kromsl, 0);
  }
  /*! [param_density_id] */

  /* Set temperature scale */

  /*! [param_itpscl] */
  {
    cs_thermal_model_t *thm = cs_get_glob_thermal_model();

    thm->itpscl = CS_TEMPERATURE_SCALE_CELSIUS;
  }
  /*! [param_itpscl] */

  /* If a user-defined scalar behaves like a temperature (relative to Cp):
     we set the "is_temperature" keyword to 1. */

  /*! [param_kscacp] */
  {
    const int kscacp = cs_field_key_id("is_temperature");

    cs_field_set_key_int(cs_field_by_name("t_1"), kscacp, 0);
  }
  /*! [param_kscacp] */

  /* Physical properties for compressible flows */

  /*! [param_compressible] */
  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

    /* Reference molecular thermal conductivity
       visls0 = lambda0 (molecular thermal conductivity, W/(m K)) */

    const int kvisl0 = cs_field_key_id("diffusivity_ref");

    cs_field_set_key_double(CS_F_(t), kvisl0, 3.e-2);

    /* If the molecular thermal conductivity is variable, its values
       must be provided in 'cs_user_physical_properties'. */

    /* Reference volumetric molecular viscosity
       viscv0 = kappa0  (volumetric molecular viscosity, kg/(m s)) */

    fp->viscv0 = 0;

    /* If the volumetric molecular viscosity is variable, its values
       must be provided in 'cs_user_physical_properties' */

    /* Molar mass of the gas (kg/mol)
       For example with dry air, xmasml is around 28.8d-3 kg/mol */

    fp->xmasmr = 0.028966;

    /* Hydrostatic equilibrium at boundaries
       Specify if the hydrostatic equilibrium must be accounted for
       (yes = 1 , no = 0) */

    cs_cf_model_t *cf_model = cs_get_glob_cf_model();

    cf_model->icfgrp = 1;
  }
  /*! [param_compressible] */

  /* Example: Change options relative to the inner iterations
   * over prediction-correction.
   * - nterup: number of sub-iterations (default 1)
   * - epsup: relative precision (default 10e-5)
   *------------------------------------*/

  /*! [param_vp_netrup] */
  {
    cs_velocity_pressure_param_t *vp_param
      = cs_get_glob_velocity_pressure_param();
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

  /* Linear solver parameters (for each unknown)
     epsilo: relative precision for the solution of the linear system. */

  /*! [param_linear_solver_epsilo] */
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
  /*! [param_linear_solver_epsilo] */

  /* Dynamic reconstruction sweeps to handle non-orthogonalities
     This parameter automatically computes a dynamic relax factor,
     and can be activated for any variable.
      - iswdyn = 0: no relaxation
      - iswdyn = 1: means that the last increment is relaxed
      - iswdyn = 2: (default) means that the last two increments are used
                    to relax.
  */

  /*! [param_iswydn] */
  {
    cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(p));
    eqp->iswdyn = 2;
  }
  /*! [param_iswydn] */

  /* Stabilization in turbulent regime

    For difficult cases, a stabilization may be obtained by not
    reconstructing the convective and diffusive flux for variables
    of the turbulence model, that is for k-epsilon models: */

  /*! [param_ircflu] */
  {
    cs_equation_param_t *eqp;

    eqp = cs_field_get_equation_param(CS_F_(k));
    eqp->ircflu = 0;

    eqp = cs_field_get_equation_param(CS_F_(eps));
    eqp->ircflu = 0;
  }
  /*! [param_ircflu] */

  /*
   * Turbulent flux model u'T' for the scalar T
   *   Algebraic Model
   *      0  SGDH
   *     10 GGDH
   *     11 EB-GGDH (Elliptic Blending)
   *     20 AFM
   *     21 EB-AFM (Elliptic Blending)
   *   Model with transport equations
   *     30 DFM
   *     31 EB-DFM (Elliptic Blending)
   */

  /*! [param_iturt] */
  {
    const int kturt = cs_field_key_id("turbulent_flux_model");

    const int turb_flux_model = 10;

    /* GGDH for thermal scalar */
    cs_field_t *f_t = cs_thermal_model_field();
    if (f_t != NULL)
      cs_field_set_key_int(f_t, kturt, turb_flux_model);

    /* GGDH for all user-defined scalars */
    const int n_fields = cs_field_n_fields();
    const int k_scal = cs_field_key_id("scalar_id");
    const int kscavr = cs_field_key_id("first_moment_id");

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & (CS_FIELD_VARIABLE | CS_FIELD_USER)) {
        int s_num = cs_field_get_key_int(f, k_scal);
        int iscavr = cs_field_get_key_int(f, kscavr);
        if (s_num > 0 && iscavr <= 0)
          cs_field_set_key_int(f, kturt, turb_flux_model);
      }
    }
  }
  /*! [param_iturt] */

  /* Example: choose a convective scheme and
   * a limiter for a given variable (user and non-user) */
  /*----------------------------------------------*/

  /* Convective scheme

    blencv = 0 for upwind (order 1 in space, "stable but diffusive")
           = 1 for centered/second order (order 2 in space)
    We may use intermediate real values. */

  /*! [param_convective_scheme] */
  {
    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);

      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->blencv = 0.98;
      }

    }
  }
  /*! [param_convective_scheme] */

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

    cs_equation_param_t *eqp = cs_field_get_equation_param(sca1);

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

  /* Example: to change the turbulent Schmidt number */
  /*-------------------------------------------------*/
  {
    cs_field_set_key_double(CS_F_(t),
                            cs_field_key_id("turbulent_schmidt"),
                            0.7);
  }

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

  /*--------------------------------------------------------------------------*/

  /* Example: add field to post-process subgrid TKE and epsilon
   * (for LES models only)
   */

  cs_parameters_add_property("k_sgs",
                             1,
                             CS_MESH_LOCATION_CELLS);

  cs_parameters_add_property("epsilon_sgs",
                             1,
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

  cs_function_define_q_criterion();

  /*! [param_var_q_criterion] */

  /* Example: post-process error estimators */
  /*----------------------------------------*/

  /*! [error_indicators] */

  /* We recommend running a calculation restart on a few time steps
     with the activation of the most interesting error indicators. */

  {
    const int field_type = CS_FIELD_INTENSIVE | CS_FIELD_POSTPROCESS;

    const char *name[] = {"est_error_cor_2",
                          "est_error_tot_2"};

    for (int i = 0; i < 2; i++) {
      cs_field_t *f = cs_field_create(name[i],
                                      field_type,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      false);
      cs_field_set_key_int(f, cs_field_key_id("log"), 1);
      cs_field_set_key_int(f, cs_field_key_id("post_vis"), 1);
    }
  }

  /*! [error_indicators] */

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

  /* Example: deactivate output of auxiliary checkpoint */
  /*----------------------------------------------------*/

  /* By default, this file is read, but it may be useful to deactivate
     its use when restarting after a preprocessing stage possibly leading
     to a different number of faces (such as simply joining meshes on
     a different architecture or optimization level or with different options).

     Writing of auxiliary restart files may also be deactivated using
  */

  /*! [deactivate_aux_checkpoint_write] */
  cs_glob_restart_auxiliary->read_auxiliary = 0;
  /*! [deactivate_aux_checkpoint_write] */

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

  /*-----------------------------------------------------------------*/

  /* Cooling tower:
   * Definition of cooling tower model and exchange zones
   * Air and liquid properties
   * */
  /*! [cs_user_cooling_towers] */
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  fp->ro0    = 1.17862;
  fp->viscl0 = 1.765e-5;

  cs_air_fluid_props_t  *air_prop = cs_glob_air_props;
  air_prop->cp_a    = 1006.0;
  air_prop->cp_v    = 1831.0;
  air_prop->cp_l    = 4179.0;
  air_prop->hv0     = 2501600.0;
  air_prop->rho_l   = 997.85615;
  air_prop->lambda_l = 0.02493;
  air_prop->humidity0 = 0.0;
  air_prop->droplet_diam = 0.005;
  /*! [cs_user_cooling_towers] */
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
  CS_NO_WARN_IF_UNUSED(domain);

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

  /* Output and probes for Radiative Transfer
   * (Luminance and radiative density flux vector) */

  /*! [setup_post_lum] */
  {
    cs_field_t *f = cs_field_by_name_try("rad_energy");
    if (f != NULL)
      cs_field_set_key_int_bits(f,
                                cs_field_key_id("post_vis"),
                                CS_POST_ON_LOCATION | CS_POST_MONITOR);

    f = cs_field_by_name_try("radiative_flux");
    if (f != NULL)
      cs_field_set_key_int_bits(f,
                                cs_field_key_id("post_vis"),
                                CS_POST_ON_LOCATION | CS_POST_MONITOR);
  }
  /*! [setup_post_lum] */

  /* Example: define 1-D radiative transfer mesh for
   * the atmospheric module */
  /*-----------------------------------------------------------------*/
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  at_opt->rad_1d_z[0 ] = 0.;
  at_opt->rad_1d_z[1 ] = 5.;
  at_opt->rad_1d_z[2 ] = 20.5;
  at_opt->rad_1d_z[3 ] = 42.0;
  at_opt->rad_1d_z[4 ] = 65.0;
  at_opt->rad_1d_z[5 ] = 89.5;
  at_opt->rad_1d_z[6 ] = 115.0;
  at_opt->rad_1d_z[7 ] = 142.0;
  at_opt->rad_1d_z[8 ] = 170.5;
  at_opt->rad_1d_z[9 ] = 199.5;
  at_opt->rad_1d_z[10] = 230.0;
  at_opt->rad_1d_z[11] = 262.0;
  at_opt->rad_1d_z[12] = 294.5;
  at_opt->rad_1d_z[13] = 328.5;
  at_opt->rad_1d_z[14] = 363.5;
  at_opt->rad_1d_z[15] = 399.0;
  at_opt->rad_1d_z[16] = 435.5;
  at_opt->rad_1d_z[17] = 473.5;
  at_opt->rad_1d_z[18] = 512.0;
  at_opt->rad_1d_z[19] = 551.0;
  at_opt->rad_1d_z[20] = 591.5;
  at_opt->rad_1d_z[21] = 632.5;
  at_opt->rad_1d_z[22] = 674.0;
  at_opt->rad_1d_z[23] = 716.0;
  at_opt->rad_1d_z[24] = 759.0;
  at_opt->rad_1d_z[25] = 802.5;
  at_opt->rad_1d_z[26] = 846.5;
  at_opt->rad_1d_z[27] = 891.5;
  at_opt->rad_1d_z[28] = 936.5;
  at_opt->rad_1d_z[29] = 982.0;
  at_opt->rad_1d_z[30] = 1028.0;
  at_opt->rad_1d_z[31] = 1074.5;
  at_opt->rad_1d_z[32] = 1122.0;
  at_opt->rad_1d_z[33] = 1169.5;
  at_opt->rad_1d_z[34] = 1217.0;
  at_opt->rad_1d_z[35] = 1265.5;
  at_opt->rad_1d_z[36] = 1314.5;
  at_opt->rad_1d_z[37] = 1363.5;
  at_opt->rad_1d_z[38] = 1413.0;
  at_opt->rad_1d_z[39] = 1462.5;
  at_opt->rad_1d_z[40] = 1512.5;
  at_opt->rad_1d_z[41] = 1563.0;
  at_opt->rad_1d_z[42] = 1613.5;
  at_opt->rad_1d_z[43] = 1664.5;
  at_opt->rad_1d_z[44] = 1715.5;
  at_opt->rad_1d_z[45] = 1767.0;
  at_opt->rad_1d_z[46] = 1818.5;
  at_opt->rad_1d_z[47] = 1870.0;
  at_opt->rad_1d_z[48] = 1922.5;
  at_opt->rad_1d_z[49] = 1975.0;

  /* Complete 1-D mesh to ztop in case of radiative transfer */
  if (at_opt->rad_1d_nvert > 0) {
    int i = at_opt->rad_1d_nlevels;
    cs_real_t zvmax = 1975.;/* top of the domain */
    cs_real_t ztop = 11000.;/* top of the troposphere */
    for (cs_real_t zzmax = (((int) zvmax)/1000)*1000.;
         zzmax <= (ztop -1000);
         i++) {
      zzmax += 1000.;
      at_opt->rad_1d_z[i] = zzmax;
    }
  }

  /* Initialize position of each vertical */
  for (int i = 0; i < at_opt->rad_1d_nvert; i++) {
    at_opt->rad_1d_xy[0 * at_opt->rad_1d_nvert + i] = 50.; /* X coord */
    at_opt->rad_1d_xy[1 * at_opt->rad_1d_nvert + i] = 50.; /* Y coord */
    at_opt->rad_1d_xy[2 * at_opt->rad_1d_nvert + i] = 1.; /* kmin in case of
                                                             non-flat terrain */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
