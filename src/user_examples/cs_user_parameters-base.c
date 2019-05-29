/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_internal_coupling.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_stokes_model.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

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
  /* If the GUI is used, user fields should preferably be defined with the GUI,
     so that associated numerical options, boundary conditions, initializations,
     and such may also be defined using the GUI. */

  if (cs_gui_file_is_loaded())
    return;

  /*--------------------------------------------------------------------------*/

  /* Example: Chose a turbulence model
   *   CS_TURB_NONE: no turbulence model (laminar flow)
   *   CS_TURB_MIXING_LENGTH: mixing length model
   *   CS_TURB_K_EPSILON: standard k-epsilon model
   *   CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
   *     Linear Production (LP) correction
   *   CS_TURB_K_EPSILON_LS: Launder-Sharma low Re k-epsilon model
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

  /*--------------------------------------------------------------------------*/

  /* Example: Coupled solver for Rij components (when iturb=30, 31 or 32)
   *   0: switch off
   *   1: switch on (default)
   */

  cs_turb_rans_model_t *rans_model = cs_get_glob_turb_rans_model();
  rans_model->irijco = 1;

  /*--------------------------------------------------------------------------*/

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

  /*--------------------------------------------------------------------------*/

  /* Example: activate ALE (Arbitrary Lagrangian Eulerian) method
   *   CS_ALE_NONE: switch off
   *   CS_ALE_LEGACY: legacy solver
   *   CS_ALE_CDO: CDO solver
   */

  cs_glob_ale = CS_ALE_LEGACY;

  /*--------------------------------------------------------------------------*/

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

  /*--------------------------------------------------------------------------*/

  /* Example: add the variance of a user variable.
   *
   * parameters for cs_parameters_add_variable_variance():
   *   name          <-- name of variance and associated field
   *   variable_name <-- name of associated variable
   */

  cs_parameters_add_variable_variance("variance_1",
                                      "species_1");

  /*--------------------------------------------------------------------------*/

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

  /* Time stepping  (0 : uniform and constant
                     1 : variable in time, uniform in space
                     2 : variable in time and space
                    -1 : steady algorithm) */


  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  time_opt->idtvar = 0;


  /* Reference time step dt_ref
     The example given below is probably not adapted to your case. */

  cs_real_t dt_ref = 0.005; // was 0.05
  domain->time_step->dt_ref = dt_ref;

  /* Duration
     nt_max is absolute number of the last time step required
     if we have already run 10 time steps and want to
     run 10 more, nt_max must be set to 10 + 10 = 20 */

  domain->time_step->nt_max = (int) (40./dt_ref);


  /* Example: set options for Stokes solving */
  /*-----------------------------------------*/

  /* Members of the strucutre cs_glob_stokes_model:
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
   */

  /*! [param_stokes_model] */
  {
    cs_stokes_model_t *stokes = cs_get_glob_stokes_model();
    stokes->arak = 0.;
  }
  /*! [param_stokes_model] */

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
   cs_field_t *f = cs_field_by_name_try('luminance');
   if (f != NULL)
     cs_field_set_key_int_bits(f,
                               cs_field_key_id("post_vis"),
                               CS_POST_MONITOR);

   f = cs_field_by_name_try('radiative_flux');
   if (f != NULL)
     cs_field_set_key_int_bits(f,
                               cs_field_key_id("post_vis"),
                               CS_POST_MONITOR);
  /*! [setup_post_lum] */



}

END_C_DECLS
