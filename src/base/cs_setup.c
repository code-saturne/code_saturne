/*============================================================================
 * Setup computation based on provided user data and functions.
 *============================================================================*/

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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "cs_ale.h"
#include "cs_at_data_assim.h"
#include "cs_atmo.h"
#include "cs_atmo_variables.h"
#include "cs_cf_thermo.h"
#include "cs_coal_read_data.h"
#include "cs_ctwr.h"
#include "cs_ctwr_variables.h"
#include "cs_domain_setup.h"
#include "cs_elec_model.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_function_default.h"
#include "cs_gas_mix.h"
#include "cs_gui.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_gui_output.h"
#include "cs_gui_radiative_transfer.h"
#include "cs_gui_specific_physics.h"
#include "cs_gwf.h"
#include "cs_ibm.h"
#include "cs_internal_coupling.h"
#include "cs_lagr.h"
#include "cs_lagr_options.h"
#include "cs_mesh_location.h"
#include "cs_mobile_structures.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties.h"
#include "cs_porous_model.h"
#include "cs_porosity_from_scan.h"
#include "cs_post.h"
#include "cs_pressure_correction.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_options.h"
#include "cs_restart.h"
#include "cs_runaway_check.h"
#include "cs_thermal_model.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_wall_condensation.h"
#include "cs_wall_distance.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_steady_laminar_flamelet_read_base(void);

void
cs_f_usipes(int *nmodpp);

void
cs_f_impini(void);

void
cs_f_indsui(void);

void
cs_f_colecd(void);

void
cs_f_usppmo(void);

void
cs_f_fldvar(int *nmodpp);

void
cs_f_atini1(void);

void
cs_f_solcat(int iappel);

void
cs_f_fldprp(void);

void
cs_f_ppprop(void);

void
cs_f_usipsu(int *nmodpp);

void
cs_f_varpos(void);

void
cs_f_iniini(void);

void
cs_f_ppinii(void);

void
cs_f_ppini1(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief First initialization stages
 */
/*----------------------------------------------------------------------------*/

static void
_init_setup(void)
{
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n\n"
       "===============================================================\n\n\n"
       "                   CALCULATION PREPARATION\n"
       "                   =======================\n\n\n"
       "===============================================================\n\n\n"));

  /* File for some specific physical models */

  cs_atmo_set_meteo_file_name("meteo");

  /* Handle some reference and physical values */

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  fp->pther = fp->p0;

  /* Other mappings */

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    cs_set_glob_turb_model(); /* set global pointer to turbulence model */
  }

  cs_f_iniini();
  cs_f_ppinii();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Disable logging and visualization for a given field.
 *
 * \param (in, out]  f  pointer to field
 */
/*----------------------------------------------------------------------------*/

static void
_hide_field(cs_field_t *f)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property field at cells.
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_property_field
(
  const char *f_name,       /*!< field name */
  const char *f_label,      /*!< field label, or NULL */
  int         dim,          /*!< field dimension */
  bool        has_previous  /*!< do we maintain time step values */
)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  if (cs_field_by_name_try(f_name) != NULL)
    cs_parameters_error(CS_ABORT_IMMEDIATE,
                        _("initial data setup"),
                        _("Field %s has already been assigned.\n"),
                        f_name);

  cs_physical_property_define_from_field(f_name,
                                         CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                         CS_MESH_LOCATION_CELLS,
                                         dim,
                                         has_previous);

  int f_id = cs_physical_property_field_id_by_name(f_name);
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  if (f_label != NULL)
    cs_field_set_key_str(f, keylbl, f_label);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an associted property field at boundary.
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_property_field_boundary
(
  const char *f_name,       /*!< field name */
  const char *f_label,      /*!< field label, or NULL */
  int         dim,          /*!< field dimension */
  bool        has_previous  /*!< do we maintain time step values */
)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  if (cs_field_by_name_try(f_name) != NULL)
    cs_parameters_error(CS_ABORT_IMMEDIATE,
                        _("initial data setup"),
                        _("Field %s has already been assigned.\n"),
                        f_name);

  cs_field_t *f = cs_field_create(f_name,
                                  CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                  CS_MESH_LOCATION_BOUNDARY_FACES,
                                  dim,
                                  has_previous);

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  if (f_label != NULL)
    cs_field_set_key_str(f, keylbl, f_label);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a simple variable field.
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_variable_field
(
  const char *f_name,
  const char *f_label,
  int         dim
)
{
  cs_field_t *f
    = cs_field_by_id(cs_variable_field_create(f_name,
                                              f_label,
                                              CS_MESH_LOCATION_CELLS,
                                              dim));

  cs_add_variable_field_indexes(f->id);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a model scalar field.
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_model_scalar_field
(
  const char *f_name,
  const char *f_label,
  int         dim
)
{
  cs_field_t *f
    = cs_field_by_id(cs_variable_field_create(f_name,
                                              f_label,
                                              CS_MESH_LOCATION_CELLS,
                                              dim));

  cs_add_model_field_indexes(f->id);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Variables definition initialization based on active model.
 */
/*----------------------------------------------------------------------------*/

static void
_init_variable_fields(void)
{
  /* Set turbulence model type to allow for simpler tests */
  cs_turb_model_t *turb_model = cs_get_glob_turb_model();
  turb_model->itytur = turb_model->iturb / 10;

  /* Check models compatibility (this test should be improved,
     as allowed ranges vary from one model to another). */

  int nmodpp = 0;
  for (int ipp = CS_PHYSICAL_MODEL_FLAG +1;
       ipp < CS_N_PHYSICAL_MODEL_TYPES;
       ipp++) {
    if (cs_glob_physical_model_flag[ipp] != -1) {
      nmodpp += 1;
      if (   cs_glob_physical_model_flag[ipp] < -1
          || cs_glob_physical_model_flag[ipp] > 5) {
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("initial data verification"),
           _("Incorrect physical model selection.\n\n"
             "Out-of bounds value (%d) for cs_glob_physical model flag[%d].\n"),
           cs_glob_physical_model_flag[ipp], ipp);
      }
    }
  }

  int nmodpp_compatibility = nmodpp;

  /* Compressible module and gas mix are compatible */
  if (   cs_glob_physical_model_flag[CS_GAS_MIX] != -1
      && cs_glob_physical_model_flag[CS_COMPRESSIBLE] != -1) {
    nmodpp_compatibility -= 1;
  }

  /* Atmo in humid atmosphere et Couling tower (iaeros) coupling */
  if (   cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2
      && cs_glob_physical_model_flag[CS_COOLING_TOWERS] != -1) {
    nmodpp_compatibility -= 1;
  }

  if (nmodpp_compatibility > 1) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("initial data verification"),
       _("%d incompatible physical models are enabled simultaneously.\n"),
       nmodpp_compatibility);
  }

  cs_parameters_error_barrier();

  /* Set global indicator */

  cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] = (nmodpp > 0) ? 1 : 0;

  /* In case ideal gas mix specific physics was enabled by the user
     together with the compressible module, the equation of state
     indicator is reset to the approprate value automatically (ieos=3)
     and the user is warned. */
  if (   cs_glob_physical_model_flag[CS_GAS_MIX] >= 0
      && cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
      && cs_glob_cf_model->ieos != CS_EOS_GAS_MIX) {
    cs_cf_model_t *cf_model = cs_get_glob_cf_model();
    cf_model->ieos = CS_EOS_GAS_MIX;
    cs_parameters_error
      (CS_WARNING,
       _("initial data verification"),
       _("Equation of state incompatible with selected physical model.\n"
         "\n"
         "The compressible and gas mix models are  enabled but the selected\n"
         "equation of state is not ideal gas mix.\n"
         "\n"
         "cs_glob_cf_model->ieos is forced to CS_EOS_GAS_MIX.\n"));
  }

  /* Enable VoF model if free surface or mass transfer modeling enabled */
  int vof_mask = CS_VOF_FREE_SURFACE | CS_VOF_MERKLE_MASS_TRANSFER;
  if ((cs_glob_vof_parameters->vof_model & vof_mask) != 0) {
    cs_vof_parameters_t *vof_parameters = cs_get_glob_vof_parameters();
    vof_parameters->vof_model |= CS_VOF_ENABLED;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create variable fields based on active model.
 */
/*----------------------------------------------------------------------------*/

static void
_create_variable_fields(void)
{
  const int keycpl = cs_field_key_id("coupled");

  int dyn_field_id_start = cs_field_n_fields();

  /*! Velocity */

  {
    cs_field_t *f = _add_variable_field("velocity", "Velocity", 3);
    cs_field_set_key_int(f, keycpl, 1);
    cs_field_pointer_map(CS_ENUMF_(vel), f);
  }

  /* Pressure */

  {
    cs_field_t *f = _add_variable_field("pressure", "Pressure", 1);
    cs_field_pointer_map(CS_ENUMF_(p), f);

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);

    // elliptic equation
    eqp->iconv = 0;

    // compressible algorithm
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
        || (cs_glob_velocity_pressure_model->idilat == 2
            && cs_glob_cf_model->ieos != CS_EOS_NONE))
      eqp->istat = 1;
    else
      eqp->istat = 0;

    // VoF algorithm: activate the weighting for the pressure
    if (cs_glob_vof_parameters->vof_model > 0)
      eqp->iwgrec = 1;
  }

  /* void fraction (VoF algorithm) */

  if (cs_glob_vof_parameters->vof_model > 0) {

    cs_field_t *f = _add_variable_field("void_fraction", "Void Fraction", 1);
    cs_field_pointer_map(CS_ENUMF_(void_f),  f);

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    eqp->idiff = 0;  // pure convection equation

    // NVD/TVD scheme
    eqp->ischcv = 4;
    // (CICSAM limiter)
    cs_field_set_key_int(f, cs_field_key_id("limiter_choice"), 11);
    // Beta Limiter
    eqp->isstpc = 2;

    // Bounds for the beta limiter
    cs_field_set_key_double(f, cs_field_key_id("min_scalar"), 0.);
    cs_field_set_key_double(f, cs_field_key_id("max_scalar"), 1.);

  }

  /* Turbulence */

  const int itytur = cs_glob_turb_model->itytur;
  const cs_turb_model_type_t iturb = cs_glob_turb_model->iturb;

  if (itytur == 2) {
    cs_field_pointer_map(CS_ENUMF_(k),
                         _add_variable_field("k", "Turb Kinetic Energy", 1));
    cs_field_pointer_map(CS_ENUMF_(eps),
                         _add_variable_field("epsilon", "Turb Dissipation", 1));
  }
  else if (itytur == 3) {
    cs_field_pointer_map(CS_ENUMF_(rij),
                         _add_variable_field("rij", "Rij", 6));
    cs_field_set_key_int(CS_F_(rij), keycpl, 1);

    cs_field_pointer_map(CS_ENUMF_(eps),
                         _add_variable_field("epsilon", "Turb Dissipation", 1));

    if (iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_field_pointer_map(CS_ENUMF_(alp_bl),
                           _add_variable_field("alpha", "Alphap", 1));

      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(alp_bl));
      // Elliptic equation (no convection, no time term)
      eqp->istat = 0;
      eqp->iconv = 0;
      // For alpha, we always have a diagonal term, so do not shift the diagonal
      eqp->idircl = 0;
    }
  }
  else if (itytur == 5) {
    cs_field_pointer_map(CS_ENUMF_(k),
                         _add_variable_field("k", "Turb Kinetic Energy", 1));
    cs_field_pointer_map(CS_ENUMF_(eps),
                         _add_variable_field("epsilon", "Turb Dissipation", 1));
    cs_field_pointer_map(CS_ENUMF_(phi),
                         _add_variable_field("phi", "Phi", 1));
    if (iturb == CS_TURB_V2F_PHI) {
      cs_field_pointer_map(CS_ENUMF_(f_bar),
                           _add_variable_field("f_bar", "f_bar", 1));
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(f_bar));
      // Elliptic equation (no convection, no time term)
      eqp->istat = 0;
      eqp->iconv = 0;
      // For f_bar, we always have a diagonal term, so do not shift the diagonal
      eqp->idircl = 0;
    }
    else if (iturb == CS_TURB_V2F_BL_V2K) {
      cs_field_pointer_map(CS_ENUMF_(alp_bl),
                           _add_variable_field("alpha", "Alpha", 1));
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(alp_bl));
      // Elliptic equation (no convection, no time term)
      eqp->istat = 0;
      eqp->iconv = 0;
      // For alpha, we always have a diagonal term, so do not shift the diagonal
      eqp->idircl = 0;
    }
  }
  else if (iturb == CS_TURB_K_OMEGA) {
    cs_field_pointer_map(CS_ENUMF_(k),
                         _add_variable_field("k", "Turb Kinetic Energy", 1));
    cs_field_pointer_map(CS_ENUMF_(omg),
                         _add_variable_field("omega", "Omega", 1));
  }
  else if (iturb == CS_TURB_SPALART_ALLMARAS) {
    cs_field_pointer_map(CS_ENUMF_(nusa),
                         _add_variable_field("nu_tilda", "NuTilda", 1));
  }

  /* Mesh velocity with ALE */

  if (cs_glob_ale != CS_ALE_NONE) {
    cs_field_t *f;

    // field defined on vertices if CDO-Vb scheme is used
    if (cs_glob_ale == CS_ALE_CDO) {
      int f_id = cs_variable_cdo_field_create("mesh_velocity",
                                              "Mesh Velocity",
                                              CS_MESH_LOCATION_VERTICES,
                                              3,
                                              1);
      f = cs_field_by_id(f_id);
      // TODO remove this once iuma is not referenced in Fortran anymore
      cs_add_variable_field_indexes(f->id);
    }
    else
      f = _add_variable_field("mesh_velocity", "Mesh Velocity", 3);

    cs_field_pointer_map(CS_ENUMF_(mesh_u), f);

    cs_field_set_key_int(f, keycpl, 1);

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    eqp->istat = 0;
    eqp->iconv = 0;
    eqp->idifft = 0;
    eqp->relaxv = 1;
  }

  int dyn_field_id_end = cs_field_n_fields();

  /* Lock time step multiplier for all "dynamic" model fields.
     (used only for user or model scalars or similar fields). */

  const int keycdt = cs_field_key_id("time_step_factor");

  for (int i = dyn_field_id_start; i < dyn_field_id_end; i++) {
    cs_field_t *f = cs_field_by_id(i);
    cs_field_lock_key(f, keycdt);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create variable fields based on active model.
 */
/*----------------------------------------------------------------------------*/

static void
_create_property_fields(void)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  cs_field_t *f;

  /* Main properties
     -------------- */

  // Base properties, always present

  {
    f = _add_property_field("density", "Density", 1, false);
    cs_field_pointer_map(CS_ENUMF_(rho), f);

    // Postprocessed and in the log file by default, hidden later in
    // cs_parameters_*_complete if constant.
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    cs_field_set_key_int(f, keylog, 1);

    f = _add_property_field_boundary("boundary_density", "Boundary Density",
                                     1, false);
    cs_field_pointer_map(CS_ENUMF_(rho_b), f);
  }

  {
    f = _add_property_field("molecular_viscosity", "Laminar Viscosity",
                            1, false);
    cs_field_pointer_map(CS_ENUMF_(mu), f);

    f = _add_property_field("turbulent_viscosity", "Turb Viscosity",
                            1, false);
    cs_field_pointer_map(CS_ENUMF_(mu_t), f);
    if (turb_model->iturb == CS_TURB_NONE)
      _hide_field(f);
  }

  /* Hybrid RANS/LES function f_d is stored for Post Processing in hybrid_blend.
     If the hybrid spatial scheme is activated for the velocity (ischcv=3)
     create field hybrid_blend which contains the local blending factor. */

  {
    cs_equation_param_t *eqp_u = cs_field_get_equation_param(CS_F_(vel));
    if (eqp_u->ischcv == 3 || turb_model->hybrid_turb > 0) {
      _add_property_field("hybrid_blend", "Hybrid blending function", 1, false);
    }

    if (turb_model->hybrid_turb == 3) {
      _add_property_field("hybrid_sas_source_term",
                          "SAS hybrid source term",
                          1, false);
    }
    else if (turb_model->hybrid_turb == 4) {
      _add_property_field("k_tot",   "Energy total",     1, false);
      _add_property_field("k_mod",   "Modelised Energy", 1, false);
      _add_property_field("k_res",   "Resolved Energy",  1, false);
      _add_property_field("eps_mod", "Mean Dissipation", 1, false);
      if (turb_model->iturb == CS_TURB_K_OMEGA) {
        _add_property_field("omg_mod",  "Mean Specific Dissipation", 1, false);
        _add_property_field("f1_kwsst", "Function F1 of k-omg SST",  1, false);
      }
      _add_property_field("htles_psi", "Psi HTLES",          1, false);
      _add_property_field("htles_r",   "Energy ratio",       1, false);
      _add_property_field("htles_t",   "Time scale HTLES",   1, false);
      _add_property_field("htles_icc", "ICC coefficient",    1, false);
      _add_property_field("htles_fs",  "Shielding function", 1, false);
      _add_property_field("Delta_max", "Delta max",          1, false);

      // Time averaged with exponential filtering, TODO use standard time moment
      _add_property_field("vel_mag_mean","Mean velocity mag.",1, false);
      _add_property_field("velocity_mean", "Vel Tavg", 3, false);

      // Diagonal part of time moment of uiuj
      _add_property_field("ui2_mean", "Vel Tavg", 3, false);
    }
  }

  if (turb_model->iturb == CS_TURB_K_OMEGA) {
    // Square of the norm of the deviatoric part of the deformation rate
    // tensor (\f$S^2=2S_{ij}^D S_{ij}^D\f$).
    f = _add_property_field("s2", "S2", 1, false);
    _hide_field(f);

    // Divergence of the velocity. More precisely, trace of the velocity gradient
    // (and not a finite volume divergence term). Defined only for k-omega SST
    // (because in this case it may be calculated at the same time as \f$S^2\f$)
    f = _add_property_field("vel_gradient_trace", "Vel. Gradient Trace",
                            1, false);
    _hide_field(f);
  }

  {
    int idtvar = cs_glob_time_step_options->idtvar;

    f = _add_property_field("courant_number", "CFL", 1, false);
    if (idtvar < 0)
      _hide_field(f);

    if (cs_glob_vof_parameters->vof_model > 0) {
      f = _add_property_field("volume_courant_number", "CourantNbVol", 1, false);
      if (idtvar < 0)
        _hide_field(f);
    }

    f = _add_property_field("fourier_number", "Fourier Number", 1, false);
    if (idtvar < 0)
      _hide_field(f);
  }

  // Total pressure is stored in a property field.
  // if the compressible module is not enabled (otherwise Ptot=P*).
  // only used if the gravity is set.

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
    f = _add_property_field("total_pressure", "Total Pressure", 1, false);

    // Save total pressure in auxiliary restart file
    int k_restart_id = cs_field_key_id("restart_file");
    cs_field_set_key_int(f, k_restart_id, CS_RESTART_AUXILIARY);
  }

  //! Cs^2 for dynamic LES model
  if (turb_model->iturb == CS_TURB_LES_SMAGO_DYN) {
    _add_property_field("smagorinsky_constant^2", "Csdyn2", 1, false);
  }

  /* Additions for specific models
     ----------------------------- */

  cs_f_ppprop();

  // Compressible model
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
    cs_cf_set_thermo_options();
    cs_cf_add_property_fields();
  }

  // Electric models
  if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 1
      || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_elec_add_property_fields();

  // Atmospheric modules
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0)
    cs_atmo_add_property_fields();

  // Cooling towers model
  if (cs_glob_physical_model_flag[CS_COOLING_TOWERS] >= 0)
    cs_ctwr_add_property_fields();

  // Add the mixture molar mass fraction field
  if (cs_glob_physical_model_flag[CS_GAS_MIX] >= 0)
    cs_gas_mix_add_property_fields();

  if (cs_glob_vof_parameters->vof_model & CS_VOF_FREE_SURFACE) {
    cs_vof_parameters_t *vof_parameters = cs_get_glob_vof_parameters();
    vof_parameters->idrift = 2;
    _add_property_field("drift_velocity", "Drift Velocity", 3, false);
  }

  // Auxiliary property fields dedicated to ALE model
  if (cs_glob_ale != CS_ALE_NONE)
    cs_ale_add_property_fields();

  /* Other properties and fields
     --------------------------- */

  cs_parameters_define_auxiliary_fields();

  // User-defined properties
  cs_parameters_create_added_properties();

  /* Ensure some field pointers are mapped */

  cs_field_pointer_map_base();
  cs_field_pointer_map_boundary();

  /* Map some Fortran field ids */

  cs_f_fldprp();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create additional fields based on user options.
 */
/*----------------------------------------------------------------------------*/

static void
_additional_fields_stage_1(void)
{
  /* Initialization
     -------------- */

  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  const int pflag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  cs_turb_rans_model_t *turb_rans_model = cs_get_glob_turb_rans_model();
  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();

  /* Determine itycor now that irccor is known (iturb/itytur known much earlier)
     type of rotation/curvature correction for turbulent viscosity models. */

  if (turb_rans_model->irccor == 1) {
    if (turb_model->itytur == 2 || turb_model->itytur == 5)
      turb_rans_model->itycor = 1;
    else if (   turb_model->iturb == CS_TURB_K_OMEGA
             || turb_model->iturb == CS_TURB_SPALART_ALLMARAS)
      turb_rans_model->itycor = 2;
  }

  /* Additional physical properties
     ------------------------------ */

  if (cs_glob_vof_parameters->vof_model != 0) {
    // variable density
    fluid_props->irovar = 1;
    fluid_props->ivivar = 1;
  }

  // CP when variable
  if (fluid_props->icp >= 0) {
    cs_field_t *f_cp = _add_property_field("specific_heat",
                                           "Specific Heat",
                                           1,
                                           false);

    cs_field_set_key_int(f_cp, keyvis, 1);
    cs_field_set_key_int(f_cp, keylog, 1);

    fluid_props->icp = f_cp->id;
  }

  // ALE mesh viscosity
  if (cs_glob_ale != CS_ALE_NONE) {
    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(CS_F_(mesh_u));
    const int idftnp = eqp->idften;

    if (idftnp & CS_ISOTROPIC_DIFFUSION)
      _add_property_field("mesh_viscosity", "Mesh Visc", 1, false);

    else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION)
      _add_property_field("mesh_viscosity", "Mesh Visc", 6, false);
  }

  if (turb_rans_model->irccor == 1 && cs_glob_time_step_options->idtvar >= 0) {
    // Strain rate tensor at the previous time step
    cs_field_t *f = _add_property_field("strain_rate_tensor",
                                        "Strain Rate Tensor",
                                        6,
                                        false);
    cs_field_set_key_int(f, keyvis, 1);
  }

  /* TODO: migrate varpos to C */

  cs_f_varpos();

  /* Porosity
     -------- */

  if (cs_glob_porous_model >= 1) {

    const char porosity_name[] = "porosity";
    cs_field_t *f = NULL;

    if (   cs_glob_porosity_from_scan_opt->compute_porosity_from_scan
        || cs_glob_porosity_ibm_opt->porosity_mode > 0) {

      // TODO move it to _create_variable_fields() ?
      f = _add_variable_field(porosity_name, NULL, 1);

      // Pure convection equation (no time term)
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      eqp->iconv = 1;
      eqp->blencv= 0.; // Pure upwind
      eqp->istat = 0;
      eqp->nswrsm = 1;
      eqp->idiff  = 0;
      eqp->idifft = 0;
      eqp->relaxv = 1.; // No relaxation, even for steady algorithm.

      // Activate the drift for all scalars with key "drift" > 0
      int iscdri =   CS_DRIFT_SCALAR_ON
                   | CS_DRIFT_SCALAR_ADD_DRIFT_FLUX
                   | CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;

      int keydri = cs_field_key_id("drift_scalar_model");
      cs_field_set_key_int(f, keydri, iscdri);

    }
    else {
      f = cs_field_create(porosity_name,
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_int(f, keyvis, pflag);
    }

    f = cs_field_create("cell_f_vol",
                        CS_FIELD_EXTENSIVE,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        false);

    if (cs_glob_porous_model == 2) {
      f = cs_field_create("tensorial_porosity",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_CELLS,
                          6,
                          false);
    }
    else if (cs_glob_porous_model == 3) {
      const int key_restart_file = cs_field_key_id("restart_file");

      f = cs_field_create("poro_div_duq",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_CELLS,
                          3,
                          false);

      f = cs_field_create("i_poro_duq_0",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          1,
                          false);
      f = cs_field_create("i_poro_duq_1",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          1,
                          false);

      f = cs_field_create("b_poro_duq",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1,
                          false);

      f = cs_field_create("i_f_face_normal",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          3,
                          false);

      f = cs_field_create("i_f_face_surf",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          1,
                          false);

      f = cs_field_create("i_f_face_cog",
                          0,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          3,
                          false);

      f = cs_field_create("b_f_face_normal",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          3,
                          false);

      f = cs_field_create("b_f_face_surf",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1,
                          false);

      f = cs_field_create("b_f_face_cog",
                          0,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          3,
                          false);

      f = cs_field_create("i_f_face_factor",
                          CS_FIELD_INTENSIVE,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          2,  // 2 values per face
                          false);

      f = cs_field_create("b_f_face_factor",
                          CS_FIELD_INTENSIVE,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1,
                          false);

      // Interior faces weighting factor with new cell cog
      f = cs_field_create("i_f_weight",
                          CS_FIELD_INTENSIVE,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          2,  // 2 values per face
                          false);

      // Solid surface normal immersed in the cells
      f = cs_field_create("c_w_face_normal",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_CELLS,
                          3,
                          false);
      cs_field_set_key_int(f, key_restart_file, CS_RESTART_IBM);

      // Center of gravity of solid face immersed in the cells
      f = cs_field_create("c_w_face_cog",
                          0,
                          CS_MESH_LOCATION_CELLS,
                          3,
                          false);

      // Solid surface of cells
      f = cs_field_create("c_w_face_surf",
                          CS_FIELD_EXTENSIVE,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);

      // Distance between the centers of the cell and the solid face
      f = cs_field_create("c_w_dist_inv",
                          0,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);

      // Cell fluid center coordinates
      f = cs_field_create("cell_f_cen",
                          0,
                          CS_MESH_LOCATION_CELLS,
                          3,
                          false);

      // Cell solid center coordinates
      f = cs_field_create("cell_s_cen",
                          0,
                          CS_MESH_LOCATION_CELLS,
                          3,
                          false);

      // Porosity at internal faces
      f = cs_field_create("i_face_porosity",
                          CS_FIELD_INTENSIVE,
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          1,
                          false);

      // Porosity at boundary faces
      f = cs_field_create("b_face_porosity",
                          CS_FIELD_INTENSIVE,
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1,
                          false);
    }

  }

  /* Local time step and postprocessing fields
     ----------------------------------------- */

  /* Local time step */

  {
    cs_field_t *f = cs_field_create("dt",
                                    CS_FIELD_INTENSIVE,
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    false);

    cs_field_set_key_str(f, keylbl, "Local Time Step");

    if (cs_glob_time_step_options->idtvar == 2) {
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_int(f, keyvis, pflag);
    }
  }

  /* Transient velocity/pressure coupling, postprocessing field
     (variant used for computation is a tensorial field, not this one) */

  int ncpdct = cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_HEAD_LOSS);

  const cs_velocity_pressure_param_t *vp_param
    = cs_glob_velocity_pressure_param;

  if (vp_param->ipucou != 0 || ncpdct > 0 || cs_glob_porous_model == 2) {

    cs_field_t *f = cs_field_create("dttens",
                                    CS_FIELD_INTENSIVE,
                                    CS_MESH_LOCATION_CELLS,
                                    6,
                                    false);
    if (vp_param->ipucou != 0 || ncpdct > 0)
      cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    cs_field_set_key_int(f, keylog, 1);
    if (cs_glob_porous_model == 2) {
      int kwgrec = cs_field_key_id("gradient_weighting_id");
      cs_field_set_key_int(CS_F_(p), kwgrec, f->id);
    }

    /* Tensorial diffusivity */

    if (cs_glob_porous_model == 2) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));
      eqp->idften = CS_ANISOTROPIC_LEFT_DIFFUSION;
    }

    /* Diagonal cell tensor for the pressure solving when needed */

    if (vp_param->ipucou == 1 || ncpdct > 0 || cs_glob_porous_model == 2) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(p));
      eqp->idften = CS_ANISOTROPIC_LEFT_DIFFUSION;
    }

  }

  /* Map to field pointers
     --------------------- */

  cs_field_pointer_map_base();
  cs_field_pointer_map_boundary();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create additional fields based on user options.
 */
/*----------------------------------------------------------------------------*/

static void
_additional_fields_stage_2(void)
{
  const int n_fields = cs_field_n_fields();
  cs_turb_les_model_t *turb_les_param = cs_get_glob_turb_les_model();

  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

  /* Keys not stored globally */
  const int kturt = cs_field_key_id("turbulent_flux_model");
  const int kfturt = cs_field_key_id("turbulent_flux_id");
  const int kfturt_alpha = cs_field_key_id("alpha_turbulent_flux_id");
  const int keycpl = cs_field_key_id("coupled");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int ks = cs_field_key_id_try("scalar_id");
  const int keyvis = cs_field_key_id("post_vis");
  const int kwgrec = cs_field_key_id("gradient_weighting_id");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kscacp = cs_field_key_id("is_temperature");
  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int key_sgs_sca_coef = cs_field_key_id("sgs_scalar_flux_coef_id");
  const int key_clipping_id = cs_field_key_id("clipping_id");
  const int key_turb_schmidt = cs_field_key_id("turbulent_schmidt_id");
  const int kromsl = cs_field_key_id_try("density_id");

  /* Key id for drift scalar */
  const int keydri = cs_field_key_id("drift_scalar_model");

  /* Time interpolation ? */
  const int key_t_ext_id = cs_field_key_id("time_extrapolated");

  /* Restart file key */
  const int key_restart_id = cs_field_key_id("restart_file");

  /* Initialization */
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  /* Additional variable fields
     -------------------------- */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = (ks > -1) ? cs_field_get_key_int(f, ks) -1 : -1;
    if (scalar_id < 0)
      continue;

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    if (eqp != NULL) {
      int post_flag = cs_field_get_key_int(f, keyvis);
      int log_flag = cs_field_get_key_int(f, keylog);
      int turb_flux_model = cs_field_get_key_int(f, kturt);
      int turb_flux_model_type = turb_flux_model / 10.;

      if (turb_flux_model_type > 0) {
        char f_tf_name[128];
        snprintf(f_tf_name, 127, "%s_turbulent_flux", f->name);
        f_tf_name[127] = '\0';
        cs_field_t *f_turb_flux = NULL;

        if (turb_flux_model_type == 3) {
          f_turb_flux = _add_variable_field(f_tf_name,
                                            f_tf_name,
                                            3);
          cs_field_set_key_int(f_turb_flux, keycpl, 1);
          cs_field_set_key_int(f_turb_flux, key_clipping_id, 1);

          /* Tensorial diffusivity */
          cs_equation_param_t *eqp_turb_flux
            = cs_field_get_equation_param(f_turb_flux);
          eqp_turb_flux->idften = CS_ANISOTROPIC_RIGHT_DIFFUSION;

          /* If variance is required */
          cs_real_t grav = cs_math_3_norm(cs_glob_physical_constants->gravity);
          if (   (   cs_glob_fluid_properties->irovar > 0
                  || cs_glob_velocity_pressure_model->idilat == 0)
              && grav > cs_math_epzero) {
            char f_var_name[128];
            snprintf(f_var_name, 127, "%s_variance", f->name);
            f_var_name[127] = '\0';

            cs_field_t *f_var = _add_model_scalar_field(f_var_name,
                                                        f_var_name,
                                                        1);
            cs_field_set_key_int(f_var, kscavr, f_id);
          }
        }
        else {
          f_turb_flux = cs_field_create(f_tf_name,
                                        CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_CELLS,
                                        3,
                                        true);
          cs_field_set_key_int(f_turb_flux, keyvis, post_flag);
          cs_field_set_key_int(f_turb_flux, keylog, log_flag);
        }

        cs_field_set_key_int(f, kfturt, f_turb_flux->id);

        /* Elliptic Blending (AFM or DFM) */
        if (   turb_flux_model == 11
            || turb_flux_model == 21
            || turb_flux_model == 31) {
          char f_alp_name[128];
          snprintf(f_alp_name, 127, "%s_alpha", f->name);
          f_alp_name[127] = '\0';

          cs_field_t *f_alpha = _add_variable_field(f_alp_name,
                                                    f_alp_name,
                                                    1);

          cs_equation_param_t *eqp_alp = cs_field_get_equation_param(f_alpha);
          eqp_alp->iconv = 0;
          eqp_alp->istat = 0;

          cs_field_set_key_int(f, kfturt_alpha, f_alpha->id);
        }
      }
    }
  }

  /* Hydrostatic pressure used to update pressure BCs */
  if (cs_glob_velocity_pressure_param->icalhy == 1) {
    cs_field_t *f = _add_variable_field("hydrostatic_pressure",
                                        "Hydrostatic Pressure",
                                        1);
    cs_field_set_key_int(f, keyvis, 0);

    /* Elliptic equation (no convection, no time term) */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    eqp->iconv = 0;
    eqp->istat = 0;
    eqp->nswrsm = 2;
    eqp->idifft = 0;
    eqp->relaxv = 1.; /* No relaxation even for steady algorithm */
  }

  /* Head losses weighting field in case of Lagrangian deposition and
   * reentrainment model (general case in varpos, but the Lagrangian
   * options are not known yet at the call site, so we have a similar
   * code block here for this special case. */
  cs_field_t *f_dttens = cs_field_by_name_try("dttens");
  if (   cs_glob_lagr_reentrained_model->iflow > 0
      && f_dttens == NULL) {
    cs_field_t *f = cs_field_create("dttens",
                                    field_type,
                                    CS_MESH_LOCATION_CELLS,
                                    6,
                                    false);
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_int(CS_F_(p), kwgrec, f->id);
  }

  /* Additional property fields
     --------------------------

     Add a scalar diffusivity when defined as variable.
     The kivisl key should be equal to -1 for constant diffusivity,
     and f_id for a variable diffusivity defined by field f_id
     Assuming the first field created is not a diffusivity property
     (we define variables first), f_id > 0, so we use 0 to indicate
     the diffusivity is variable but its field has not been created yet. */

  cs_field_t *thm_field = cs_thermal_model_field();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = (ks > -1) ? cs_field_get_key_int(f, ks) -1 : -1;
    if (scalar_id < 0)
      continue;
    int ifcvsl = cs_field_get_key_int(f, kivisl);
    int var_f_id = cs_field_get_key_int(f, kscavr);
    if (ifcvsl == 0 && var_f_id < 0) {
      /* Build name and label, using a general rule, with
       * a fixed name for temperature or enthalpy */
      const char th_name[] = "thermal";
      const char th_label[] = "Th";
      const char *f_name = f->name;
      const char *f_label = cs_field_get_label(f);
      if (f == thm_field) {
        f_name = th_name;
        f_label = th_label;
      }
      char s_name[128];
      char s_label[128];
      int iscacp = cs_field_get_key_int(f, kscacp);
      if (iscacp > 0) {
        snprintf(s_name, 127, "%s_conductivity", f_name);
        snprintf(s_label, 127, "%s Cond", f_label);
      }
      else {
        snprintf(s_name, 127, "%s_diffusivity", f_name);
        snprintf(s_label, 127, "%s Diff", f_label);
      }
      /* Special case for electric arcs: real and imaginary electric
       * conductivity is the same (and ipotr < ipoti) */
      if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 2
          || cs_glob_physical_model_flag[CS_JOULE_EFFECT] == 4) {
        cs_field_t *f_elec_port_r = cs_field_by_name_try("elec_port_r");
        cs_field_t *f_elec_port_i = cs_field_by_name_try("elec_port_i");
        if (f == f_elec_port_r) {
          snprintf(s_name, 127, "elec_sigma");
          snprintf(s_label, 127, "Sigma");
        }
        else if(f == f_elec_port_i) {
          int potr_ifcvsl = cs_field_get_key_int(f_elec_port_r, kivisl);
          cs_field_set_key_int(f, kivisl, potr_ifcvsl);
          continue; /* go to next scalar in loop, avoid creating a property */
        }
      }
      s_name[127] = '\0';
      s_label[127] = '\0';

      /* Now create matching property */
      cs_field_t *f_s = _add_property_field(s_name,
                                            s_label,
                                            1,
                                            false);
      cs_field_set_key_int(f, kivisl, f_s->id);
    }
  }

  /* For variances, the diffusivity is that of the associated scalar,
   * and must not be initialized first. */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    const int iscavr = cs_field_get_key_int(f, kscavr);
    if (iscavr < 0)
      continue;
    int ifcvsl = cs_field_get_key_int(cs_field_by_id(iscavr), kivisl);
    if (cs_field_is_key_set(f, kivisl) == true)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("initial data verification"),
         _("The field %d, represents the variance\n"
           "of fluctuations of the field %d.\n"
           "The diffusivity must not be set, it will\n"
           "automatically set equal to that of the associated scalar.\n"),
         f_id, iscavr);
    else
      cs_field_set_key_int(f, kivisl, ifcvsl);
  }

  /* Add a scalar turbulent diffusivity field */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int ifcdep = cs_field_get_key_int(f, key_turb_diff);
    int var_f_id = cs_field_get_key_int(f, kscavr);
    if (ifcdep >= 0 && var_f_id < 0) {
      /* Build name and label, using a general rule, with
       * a fixed name for temperature or enthalpy */
      char s_name[128];
      char s_label[128];
      snprintf(s_name, 127, "%s_turb_diffusivity", f->name);
      snprintf(s_label, 127, "%s Turb Diff", cs_field_get_label(f));
      s_name[127] = '\0';
      s_label[127] = '\0';

      if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0) {
        int ifm_ifcdep = cs_field_get_key_int(CS_F_(fm), key_turb_diff);
        if (f != CS_F_(fm))
          cs_field_set_key_int(f, key_turb_diff, ifm_ifcdep);
        continue;
      }

      /* Now create matching property */
      cs_field_t *f_s = _add_property_field(s_name,
                                            s_label,
                                            1,
                                            false);
      cs_field_set_key_int(f, key_turb_diff, f_s->id);
    }
  }

  /* For variances, the turbulent diffusivity is that of the associated scalar,
   * and must not be initialized first */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    const int iscavr = cs_field_get_key_int(f, kscavr);
    if (iscavr < 0)
      continue;
    int ifcdep = cs_field_get_key_int(cs_field_by_id(iscavr),
                                      key_turb_diff);
    if (cs_field_is_key_set(f, key_turb_diff) == true)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("initial data verification"),
         _("The field %d, represents the variance\n"
           "of fluctuations of the field %d.\n"
           "The diffusivity must not be set, it will\n"
           "automatically set equal to that of the associated scalar.\n"),
         f_id, iscavr);
    else
      cs_field_set_key_int(f, key_turb_diff, ifcdep);
  }

  if (cs_glob_turb_model->iturb == CS_TURB_LES_SMAGO_DYN) {
    /* Add a subgrid-scale scalar flux coefficient field */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
        continue;
      int ifcdep = cs_field_get_key_int(f, key_sgs_sca_coef);
      int var_f_id = cs_field_get_key_int(f, kscavr);
      if (ifcdep >= 0 && var_f_id < 0) {
        /* Build name and label using a general rule */
        char s_name[128];
        char s_label[128];
        snprintf(s_name, 127, "%s_sgs_flux_coef", f->name);
        snprintf(s_label, 127, "%s SGS Flux Coef", cs_field_get_label(f));
        s_name[127] = '\0';
        s_label[127] = '\0';

        if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0) {
          int ifm_ifcdep = cs_field_get_key_int(CS_F_(fm), key_sgs_sca_coef);
          if (f != CS_F_(fm))
            cs_field_set_key_int(f, key_sgs_sca_coef, ifm_ifcdep);
          continue;
        }
        /* Now create matching property */
        cs_field_t *f_s = _add_property_field(s_name,
                                              s_label,
                                              1,
                                              false);
        cs_field_set_key_int(f, key_sgs_sca_coef, f_s->id);
      }
    }

    /* For variances, the subgrid scale flux is that of the associated scalar,
     * and must not be initialized first */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
        continue;
      const int iscavr = cs_field_get_key_int(f, kscavr);
      if (iscavr > -1) {
        int ifcdep = cs_field_get_key_int(cs_field_by_id(iscavr),
                                          key_sgs_sca_coef);
        if (cs_field_is_key_set(f, key_sgs_sca_coef) == true)
          cs_parameters_error
            (CS_ABORT_DELAYED,
             _("initial data verification"),
             _("The field %d, represents the variance\n"
               "of fluctuations of the field %d.\n"
               "The diffusivity must not be set, it will\n"
               "automatically set equal to that of the associated scalar.\n"),
             f_id, iscavr);
        else
          cs_field_set_key_int(f, key_sgs_sca_coef, ifcdep);
      }
    }
  } /* End for CS_TURB_LES_SMAGO_DYN) */

  /* Add a scalar density when defined as variable and different from the bulk.
   * WARNING: it must be consitent with continuity equation, this is used
   * for fluid solid computation with passive scalars with different
   * density in the solid.
   * The kromsl key should be equal to -1 for constant density
   * and f_id for a variable density defined by field f_id
   * Assuming the first field created is not a density property
   * (we define variables first), f_id > 0, so we use 0 to indicate */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = (ks > -1) ? cs_field_get_key_int(f, ks) -1 : -1;
    if (scalar_id < 0)
      continue;
    int ifcdep = cs_field_get_key_int(f, kromsl);
    int var_f_id = cs_field_get_key_int(f, kscavr);
    if (ifcdep == 0 && var_f_id < 0) {
      char f_name[128];
      char f_label[128];
      snprintf(f_name, 127, "%s_density", f->name);
      snprintf(f_label, 127, "%s Rho", cs_field_get_label(f));
      f_name[127] = '\0';
      f_label[127] = '\0';

      /* Now create matching property */
      cs_field_t *f_s = _add_property_field(f_name,
                                            f_label,
                                            1,
                                            false);
      cs_field_set_key_int(f, kromsl, f_s->id);
    }
  }

  /* For variances, the density is that associated to the scalar,
   * and must not be initialized first */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    const int iscavr = cs_field_get_key_int(f, kscavr);
    if (iscavr < 0)
      continue;
    int ifcdep = cs_field_get_key_int(cs_field_by_id(iscavr), kromsl);
    if (cs_field_is_key_set(f, kromsl) == true)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("initial data verification"),
         _("The field %d, represents the variance\n"
           "of fluctuations of the field %d.\n"
           "The diffusivity must not be set, it will\n"
           "automatically set equal to that of the associated scalar.\n"),
         f_id, iscavr);
    else
      cs_field_set_key_int(f, kromsl, ifcdep);
  }

  /* Add a scalar turbulent Schmidt field. */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    const int iscavr = cs_field_get_key_int(f, kscavr);
    if (iscavr < 0)
      continue;
    int ifcdep = cs_field_get_key_int(f, key_turb_schmidt);
    int var_f_id = cs_field_get_key_int(f, kscavr);
    if (ifcdep == 0 && var_f_id < 0) {
      char f_name[128];
      char f_label[128];
      snprintf(f_name, 127, "%s_turb_schmidt", f->name);
      snprintf(f_label, 127, "%s ScT", cs_field_get_label(f));
      f_name[127] = '\0';
      f_label[127] = '\0';

      /* Now create matching property */
      cs_field_t *f_s = _add_property_field(f_name,
                                            f_label,
                                            1,
                                            false);
      cs_field_set_key_int(f, key_turb_schmidt, f_s->id);
    }
  }

  /* For variances, the Schmidt is that associated to the scalar,
   * and must not be initialized first */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    const int iscavr = cs_field_get_key_int(f, kscavr);
    if (iscavr < 0)
      continue;
    int ifcdep = cs_field_get_key_int(cs_field_by_id(iscavr),
                                      key_turb_schmidt);
    if (cs_field_is_key_set(f, key_turb_schmidt) == true)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("initial data verification"),
         _("The field %d, represents the variance\n"
           "of fluctuations of the field %d.\n"
           "The diffusivity must not be set, it will\n"
           "automatically set equal to that of the associated scalar.\n"),
         f_id, iscavr);
    else
      cs_field_set_key_int(f, key_turb_schmidt, ifcdep);
  }

  /* Boundary roughness (may be already created by the atmospheric module) */
  if (   cs_glob_wall_functions->iwallf == 5
      || cs_glob_wall_functions->iwallf == 6) {
    cs_field_find_or_create("boundary_rougness",
                            CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                            CS_MESH_LOCATION_BOUNDARY_FACES,
                            1,
                            false);
  }

  /* Van Driest damping */
  if (cs_glob_turb_les_model->idries == -1) {
    if (cs_glob_turb_model->iturb == 40)
      turb_les_param->idries = 1;
    else if (   cs_glob_turb_model->iturb == CS_TURB_LES_SMAGO_DYN
             || cs_glob_turb_model->iturb == CS_TURB_LES_WALE)
      turb_les_param->idries = 0;
  }

  /* Wall distance for some turbulence models
   * and for Lagrangian multilayer deposition for DRSM models, needed for inlets */
  cs_wall_distance_options_t *wdo = cs_get_glob_wall_distance_options();
  if (   cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_QUAD
      || cs_glob_turb_model->itytur == 3
      || (   cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR
          && cs_glob_turb_rans_model->irijec == 1)
      || (cs_glob_turb_model->itytur == 4 && cs_glob_turb_les_model->idries == 1)
      || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA
      || cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS
      || cs_glob_lagr_reentrained_model->iflow == 1
      || cs_glob_turb_model->hybrid_turb == 4)
    wdo->need_compute = 1;

  if (wdo->need_compute == 1) {
    cs_field_t *f_wd = _add_variable_field("wall_distance",
                                           "Wall distance",
                                           1);
    cs_field_set_key_int(f_wd, key_restart_id, CS_RESTART_AUXILIARY);

    /* Elliptic equation (no convection, no time term) */
    cs_equation_param_t *eqp_wd = cs_field_get_equation_param(f_wd);
    eqp_wd->iconv = 0;
    eqp_wd->istat = 0;
    eqp_wd->nswrsm = 2;
    eqp_wd->idifft = 0;
    eqp_wd->relaxv = 1.; // No relaxation, even for steady algorithm.

    /* Working field to store value of the solved variable at the previous time step
     * if needed (ALE) */
    int model = cs_turbomachinery_get_model();
    if (cs_glob_ale != 0 || model > 0) {
      cs_field_t *f_wd_aux_pre = _add_property_field("wall_distance_aux_pre",
                                                     NULL,
                                                     1,
                                                     false);
      cs_field_set_key_int(f_wd_aux_pre, keyvis, 0);
      cs_field_set_key_int(f_wd_aux_pre, keylog, 0);
    }

    /* Dimensionless wall distance "y+"
     * non-dimensional distance \f$y^+\f$ between a given volume and the
     * closest wall, when it is necesary (LES with van Driest wall damping */
    if (   cs_glob_turb_model->itytur == 4
        && cs_glob_turb_les_model->idries == 1) {
      cs_field_t *f_yp = _add_variable_field("wall_yplus",
                                             "Wall Y+",
                                             1);
      cs_field_set_key_int(f_yp, keyvis, CS_POST_ON_LOCATION);
      cs_field_set_key_int(f_yp, keylog, 1);

      /* Pure convection (no time term) */
      cs_equation_param_t *eqp_yp = cs_field_get_equation_param(f_yp);
      eqp_yp->iconv = 1;
      eqp_yp->istat = 0;
      eqp_yp->idiff = 0;
      eqp_yp->idifft = 0;
      eqp_yp->relaxv = 1.; // No relaxation, even for steady algorithm.
      eqp_yp->blencv = 0.; // Pure upwind
      eqp_yp->epsilo = 1.e-5; // By default, not a high precision

      int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX
        + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;
      cs_field_set_key_int(f_yp, keydri, drift);
    }
  }

  if (   cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0
      && cs_glob_atmo_option->compute_z_ground > 0) {
    cs_field_t *f_ground = _add_variable_field("z_ground",
                                               "Z ground",
                                               1);
    cs_equation_param_t *eqp_ground = cs_field_get_equation_param(f_ground);
    eqp_ground->iconv = 1;
    eqp_ground->blencv = 0.; // Pure upwind
    eqp_ground->istat = 0;
    eqp_ground->nswrsm = 100;
    eqp_ground->epsrsm = 1.e-3;
    eqp_ground->idiff = 0;
    eqp_ground->idifft = 0;
    eqp_ground->relaxv = 1.; // No relaxation, even for steady algorithm.

    int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX
      + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;
    cs_field_set_key_int(f_ground, keydri, drift);
  }

  if (cs_glob_atmo_option->meteo_profile >= 2) {
    cs_field_t *f = _add_variable_field("meteo_pressure",
                                        "Meteo pressure",
                                        1);

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    eqp->iconv = 0;
    eqp->blencv = 0.; // Pure upwind
    eqp->istat = 0;
    eqp->idircl = 1;
    eqp->nswrsm = 100;
    eqp->nswrgr = 100;
    eqp->imrgra = 0;
    eqp->imligr = -1;
    eqp->epsilo = 0.000001;
    eqp->epsrsm = 1.e-3;
    eqp->epsrgr = 0.0001;
    eqp->climgr = 1.5;
    eqp->idiff = 1;
    eqp->idifft = 0;
    eqp->idften = 1;
    eqp->relaxv = 1.; // No relaxation, even for steady algorithm
    eqp->thetav = 1;

    _add_property_field("meteo_density",
                        "Meteo density",
                        1,
                        false);

    _add_property_field("meteo_temperature",
                        "Meteo temperature",
                        1,
                        false);

    _add_property_field("meteo_pot_temperature",
                        "Meteo pot temperature",
                        1,
                        false);

    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2) {
      _add_property_field("meteo_humidity",
                          "Meteo humidity",
                          1,
                          false);

      _add_property_field("meteo_drop_nb",
                          "Meteo drop nb",
                          1,
                          false);
    }

    _add_property_field("meteo_velocity",
                        "Meteo velocity",
                        3,
                        false);

    _add_property_field("meteo_tke",
                        "Meteo TKE",
                        1,
                        false);

    _add_property_field("meteo_eps",
                        "Meteo epsilon",
                        1,
                        false);

    /* DRSM models, store Rxz/k */
    if (cs_glob_turb_model->itytur == 3) {
      _add_property_field("meteo_shear_anisotropy",
                          "Meteo shear anisotropy",
                          1,
                          false);
    }
  }

  if (cs_glob_porosity_from_scan_opt->compute_porosity_from_scan) {
    cs_field_t *f;

    f = _add_property_field("nb_scan_points",
                            "nb scan points",
                            1,
                            false);
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);

    f = _add_property_field("solid_roughness",
                            "solid roughness",
                            1,
                            false);
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);

    f = _add_property_field("cell_scan_points_cog",
                            "Point centers",
                            3,
                            false);
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    cs_field_set_key_int(f, key_restart_id, CS_RESTART_IBM);

    f = _add_property_field("cell_scan_points_color",
                            "Cell color",
                            3,
                            false);
    cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
  }

  /* Additional postprocessing fields
     -------------------------------- */

  /* Fields used to save postprocessing data */

  /* Boundary efforts postprocessing for immersed boundaries, create field */
  if (cs_glob_porous_model == 3) {
    cs_field_t *f = cs_field_create("immersed_pressure_force",
                                    CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    false);

    if (cs_glob_turb_model->iturb != 0) {
      f = cs_field_create("immersed_boundary_uk",
                          CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
      f = cs_field_create("immersed_boundary_yplus",
                          CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
      cs_field_set_key_int(f, keylog, 1);
      f = cs_field_create("immersed_boundary_dplus",
                          CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS,
                          CS_MESH_LOCATION_CELLS,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    }
  }

  /* In case of ALE or postprocessing, ensure boundary forces are tracked */
  if (cs_glob_ale >= 1) {
    cs_field_find_or_create("boundary_forces",
                            CS_FIELD_EXTENSIVE | CS_FIELD_POSTPROCESS,
                            CS_MESH_LOCATION_BOUNDARY_FACES,
                            3,
                            false);
  }

  if (   cs_glob_wall_condensation->icondb >= 0
      || cs_glob_wall_condensation->icondv >= 0) {
    int f_id = cs_field_by_name_try("yplus")->id;
    cs_field_t *f = cs_field_find_or_create("yplus",
                                            CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                            CS_MESH_LOCATION_BOUNDARY_FACES,
                                            1,
                                            false);
    if (f_id < 0) { // Set some properties if the field is new
      cs_field_set_key_str(f, keylbl, "Yplus");
      cs_field_set_key_int(f, keylog, 1);
    }
  }

  /* Some mappings */
  cs_field_pointer_map_boundary();

  /* Cooling towers mapping */
  if (cs_glob_physical_model_flag[CS_COOLING_TOWERS] >= 0) {
    cs_ctwr_field_pointer_map();
  }

  cs_field_t *f_temp_b = cs_field_by_name_try("boundary_temperature");
  if (f_temp_b != NULL) {
    cs_field_t *f_temp = cs_field_by_name_try("temperature");
    if (f_temp != NULL)
      cs_field_set_key_str(f_temp_b, keylbl, cs_field_get_label(f_temp));
  }

  /* Set some field keys and number of previous values if needed
     ----------------------------------------------------------- */

  /* Density at the second previous time step for VOF algorithm
   * or dilatable algorithm */
  if (   cs_glob_vof_parameters->vof_model > 0
      || (   cs_glob_velocity_pressure_model->idilat > 1
          && cs_glob_velocity_pressure_param->ipredfl == 0)
      || cs_glob_fluid_properties->irovar == 1) {
    cs_field_set_n_time_vals(CS_F_(rho), 3);
    cs_field_set_n_time_vals(CS_F_(rho_b), 3);
  }
  else if (   cs_glob_velocity_pressure_param->icalhy > 0
           || cs_glob_fluid_properties->ipthrm == 1
           || cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
           || cs_glob_velocity_pressure_model->idilat > 1) {
    cs_field_set_n_time_vals(CS_F_(rho), 2);
    cs_field_set_n_time_vals(CS_F_(rho_b), 2);
  }

  /* Time extrapolation */

  /* Density */
  int t_ext = cs_field_get_key_int(CS_F_(rho), key_t_ext_id);
  if (t_ext == -1) {
    if (cs_glob_time_scheme->time_order == 1)
      t_ext = 0;
    else if (cs_glob_time_scheme->time_order == 2)
      t_ext = 0;
    cs_field_set_key_int(CS_F_(rho), key_t_ext_id, t_ext);
  }

  /* Molecular viscosity */
  t_ext = cs_field_get_key_int(CS_F_(mu), key_t_ext_id);
  if (t_ext == -1) {
    if (cs_glob_time_scheme->time_order == 1)
      t_ext = 0;
    else if (cs_glob_time_scheme->time_order == 2)
      t_ext = 0;
    cs_field_set_key_int(CS_F_(mu), key_t_ext_id, t_ext);
  }

  /* Turbulent viscosity */
  t_ext = cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id);
  if (t_ext == -1) {
    if (cs_glob_time_scheme->time_order == 1)
      t_ext = 0;
    else if (cs_glob_time_scheme->time_order == 2)
      t_ext = 0;
    cs_field_set_key_int(CS_F_(mu_t), key_t_ext_id, t_ext);
  }

  /* Specific heat */
  if (CS_F_(cp) != NULL) {
    t_ext = cs_field_get_key_int(CS_F_(cp), key_t_ext_id);
    if (t_ext == -1) {
      if (cs_glob_time_scheme->time_order == 1)
        t_ext = 0;
      else if (cs_glob_time_scheme->time_order == 2)
        t_ext = 0;
    }
    cs_field_set_key_int(CS_F_(cp), key_t_ext_id, t_ext);
  }

  /* Scalar diffusivity time extrapolation */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = (ks > -1) ? cs_field_get_key_int(f, ks) -1 : -1;
    if (scalar_id < 0)
      continue;
    int f_s_id = cs_field_get_key_int(f, kivisl);
    if (f_s_id >= 0) {
      cs_field_t *f_s = cs_field_by_id(f_s_id);
      t_ext = cs_field_get_key_int(f_s, key_t_ext_id);
      if (t_ext == -1) {
        if (cs_glob_time_scheme->time_order == 1)
          t_ext = 0;
        else if (cs_glob_time_scheme->time_order == 2)
          t_ext = 0;
      }
      cs_field_set_key_int(f_s, key_t_ext_id, t_ext);
    }
  }

  /* If time extrapolation, set previous values */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    t_ext = cs_field_get_key_int(f, key_t_ext_id);
    if (t_ext > 0) {
      int nprev = f->n_time_vals - 1;
      if (nprev < 1)
        cs_field_set_n_time_vals(f, nprev + 1);
    }
  }

  /* Stop the calculation if needed once all checks have been done */
  cs_parameters_error_barrier();

  /* Map some field pointers (could be deone "on the fly") */
  cs_field_pointer_map_base();
  cs_field_pointer_map_boundary();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read specific physical model data file
 */
/*----------------------------------------------------------------------------*/

static void
_read_specific_physics_data(void)
{
  /* Diffusion flame - 3-point chemistry
   * premix flame    - EBU model
   * premix flame    - LWC model */
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] != -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] != -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_LW] != -1)
    cs_f_colecd();

  /* Diffusion flame - steady laminar flamelet approach */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] != -1)
    cs_f_steady_laminar_flamelet_read_base();

  /* Pulverized coal combustion */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] != -1) {
    cs_gui_coal_model();
    cs_coal_read_data();
  }

  /* Joule effect, electric arc or ionic conduction */
  if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 1
      || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_electrical_model_param();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function calling the user-defined functions for the definition
 *        of computation parameters: we run through here in any case.
 */
/*----------------------------------------------------------------------------*/

static void
_init_user
(
  int   *nmodpp
)
{
  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();

  int icondb = cs_glob_wall_condensation->icondb;
  int condv = cs_glob_wall_condensation->icondv;

  /* Check for restart and read matching time steps and notebook values */
  cs_parameters_read_restart_info();

  /* Flow model selection through GUI */
  cs_gui_physical_model_select();

  /* Flow model selection through user Fortran subroutine */
  /* Warning : This is a user function maybe no need to translate it...
   * It is only called in the case 69_PARISFOG */
  cs_f_usppmo();
  cs_wall_condensation_set_onoff_state(icondb, condv);

  /* Other model selection through GUI */
  /* ALE parameters */
  cs_gui_ale_params();

  /* Thermal model */
  cs_gui_thermal_model();

  /* Turbulence model */
  cs_gui_turb_model();

  cs_gui_cp_params();
  cs_gui_dt();
  cs_gui_hydrostatic_pressure();

  /* Gravity and Coriolis
   * Presence or not of gravity may be needed to determine whether some fields
   * are created, so this is called before cs_user_model (to provide a
   * user call site allowing to modify GUI-defined value programatically
   * before property fields are created). */
  cs_gui_physical_constants();

  /* Activate radiative transfer model */
  cs_gui_radiative_transfer_parameters();
  cs_user_radiative_transfer_parameters();

  /* Flow and other model selection through user C routines */
  cs_user_model();

  /* Set additional members of turbulence model and RANS
     turbulence model strucutre */
  cs_turbulence_init_models();

  /* If CDO is active, initialize the context structures for models which
   * have been activated */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_OFF) {
    /* Groundwater flow module */
    if (cs_gwf_is_activated())
      cs_gwf_init_model_context();
  }

  /* Activate CDO for ALE */
  if (cs_glob_ale == CS_ALE_CDO)
    cs_ale_activate();

  if (cs_glob_ale == CS_ALE_LEGACY)
    cs_gui_mobile_mesh_structures_add();

  /* Read thermomechanical data for specific physics */
  _read_specific_physics_data();

  /* Other model parameters, including user-defined scalars */
  cs_gui_user_variables();
  cs_gui_user_arrays();
  cs_gui_calculator_functions();

  /* Solid zones */
  cs_velocity_pressure_set_solid();

  /* Initialize parameters for specific physics */
  cs_rad_transfer_options();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    _init_variable_fields();
    _create_variable_fields();
    cs_f_fldvar(nmodpp);

    /* Map pointers */
    cs_field_pointer_map_base();
    cs_field_pointer_map_boundary();

    /* Activate pressure correction model if CDO mode is not stand-alone */
    cs_pressure_correction_model_activate();
  }

  if (cs_glob_ale != CS_ALE_NONE)
    cs_gui_ale_diffusion_type();

  cs_gui_laminar_viscosity();

  /* Specific physics modules
   * Note: part of what is inside ppini1 could be moved here
   * so that usipsu / cs_user_parameters can be used by the user
   * to modify default settings */

  /* Atmospheric flows */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1)
    cs_f_atini1();

  /* Compressible flows */
  cs_gui_hydrostatic_equ_param();
  const cs_field_t *f_id = cs_field_by_name_try("velocity");
  if (f_id != NULL) {
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] != -1)
      cs_runaway_check_define_field_max(f_id->id, 1.e5);
    else
      cs_runaway_check_define_field_max(f_id->id, 1.e4);
  }

  /* Atmospheric module */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1) {
    /* Advanced init/allocation for the soil model */
    if (cs_glob_atmo_option->soil_cat >= 0)
      cs_f_solcat(1);
  }

  /* Initialization of global parameters */
  cs_gui_output_boundary();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    _create_property_fields();

  /* Initialization of additional user parameters */
  cs_gui_checkpoint_parameters();

  cs_gui_dt_param();

  /* Local numerical options */
  cs_gui_equation_parameters();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_gui_numerical_options();

  /* Physical properties */
  cs_gui_physical_properties();

  /* Turbulence reference values (uref, almax) */
  cs_gui_turb_ref_values();

  /* Set turbulence constants according to model choices.
   * This can be overwritten by the user in cs_user_parameters() */
  cs_turb_compute_constants(-1);

  /* Scamin, scamax, turbulent flux model, diffusivities
   * (may change physical properties for some scalars, so called
   * after cs_gui_physical_properties). */
  cs_gui_scalar_model_settings();

  /* Porosity model */
  cs_gui_porous_model();

  /* Init fan */
  cs_gui_define_fans();

  /* Init error estimator */
  cs_gui_error_estimator();

  /* Initialize base evaluation functions */
  cs_function_default_define();

  /* User functions */
  cs_f_usipsu(nmodpp);
  cs_user_parameters(cs_glob_domain);

  /* If time step is local or variable, pass information to C layer, as it
   * may be needed for some field (or moment) definitions */
  if (cs_glob_time_step_options->idtvar != 0)
    cs_time_step_define_variable(1);

  if (cs_glob_time_step_options->idtvar == 2
    || cs_glob_time_step_options->idtvar == -1)
    cs_time_step_define_local(1);

  /* Initialize Fortran restarted computation flag (isuite) */
  cs_f_indsui();

  /* Default value of physical properties for the compressible model */
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] != -1) {
    /* EOS has been set above with the GUI or in cs_user_model.
     * The variability of the thermal conductivity
     * (diffusivity_id for itempk) and the volume viscosity (iviscv) has
     * been set in fldprp.
     *
     * Compute cv0 according to chosen EOS */

    cs_real_t l_cp[1] = {fluid_props->cp0};
    cs_real_t l_xmasmr[1] = {fluid_props->xmasmr};
    cs_real_t l_cv[1] = {-1};
    cs_cf_thermo_cv(l_cp, l_xmasmr, l_cv, 1);
    fluid_props->cv0 = l_cv[0];
  }

  if (cs_glob_porosity_ibm_opt->porosity_mode > 0)
    cs_porous_model_set_model(3);

  /* Varpos
   * If CDO mode only, skip this stage */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    _additional_fields_stage_1();

  /* Internal coupling */
  cs_gui_internal_coupling();
  cs_user_internal_coupling();

  cs_internal_coupling_setup();

  /* Mobile structures
   * After call to cs_gui_mobile_mesh_structures_add possible
   * call by user to cs_mobile_structures_add_n_structures */
  if (cs_glob_ale != CS_ALE_NONE)
    cs_mobile_structures_setup();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create some additional fields which depend on main field options.
 */
/*----------------------------------------------------------------------------*/

static void
_additional_fields_stage_3(void)
{
  /* Get ids */
  const int k_log = cs_field_key_id("log");
  const int k_vis = cs_field_key_id("post_vis");
  const int k_sca = cs_field_key_id("scalar_id");

  // Keys not stored globally
  const int k_turt = cs_field_key_id("turbulent_flux_model");

  // Key id for mass flux
  const int k_imasf = cs_field_key_id("inner_mass_flux_id");
  const int k_bmasf = cs_field_key_id("boundary_mass_flux_id");

  // Key id for gradient weighting
  const int k_wgrec = cs_field_key_id("gradient_weighting_id");

  // Key id for limiter
  const int k_cvlim = cs_field_key_id("convection_limiter_id");
  const int k_dflim = cs_field_key_id("diffusion_limiter_id");

  // Key id for slope test
  const int k_slts = cs_field_key_id("slope_test_upwind_id");

  // Key id of the coal scalar class
  const int k_ccl = cs_field_key_id("scalar_class");

  // Key id for drift scalar
  const int k_dri = cs_field_key_id("drift_scalar_model");

  // Key id for restart file
  const int k_restart_id = cs_field_key_id("restart_file");

  // Get number of fields
  const int n_fld = cs_field_n_fields();

  /* Global param */

  cs_velocity_pressure_param_t *vp_param
    = cs_get_glob_velocity_pressure_param();

  /* Equation param */

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(CS_F_(vel));

  /* Set keywords and add some additional fields
     ------------------------------------------- */

  /* User variables */

  int idfm = 0, iggafm = 0;

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = cs_field_get_key_int(f, k_sca) -1;
    if (scalar_id < 0)
      continue;

    const cs_equation_param_t *eqp_f = cs_field_get_equation_param_const(f);

    const int turb_flux_model      = cs_field_get_key_int(f, k_turt);
    const int turb_flux_model_type = turb_flux_model / 10;

    if (turb_flux_model_type != CS_TURB_TYPE_NONE) {
      if (turb_flux_model_type == 3)
        idfm = 1;

      /* GGDH or AFM on current scalar and if DFM, GGDH on the scalar variance
       */
      iggafm = 1;
    }
    else if (eqp_f->idften & CS_ANISOTROPIC_DIFFUSION) {
      /* If the user has chosen a tensorial diffusivity */
      idfm = 1;
    }

    /* Additional fields for Drift scalars
       is done in _additional_fields_stage_2 */
  }

  /* Reserved fields whose ids are not saved (may be queried by name) */

  const int iphydr = cs_glob_velocity_pressure_param->iphydr;
  if (iphydr == 1) {
    cs_field_t *f_vf = cs_field_find_or_create("volume_forces",
                                               CS_FIELD_INTENSIVE,
                                               CS_MESH_LOCATION_CELLS,
                                               3,
                                               false);

    cs_field_set_key_int(f_vf, k_log, 1);
    cs_field_set_key_int(f_vf, k_vis, 0);
    cs_field_set_key_int(f_vf, k_restart_id, CS_RESTART_AUXILIARY);
  }
  else if (iphydr == 2) {
    cs_field_t *f_hp = cs_field_find_or_create("hydrostatic_pressure_prd",
                                               CS_FIELD_INTENSIVE,
                                               CS_MESH_LOCATION_CELLS,
                                               1,
                                               false);

    cs_field_set_key_int(f_hp, k_restart_id, CS_RESTART_AUXILIARY);
  }

  /* Hybrid blending field */

  if (eqp_u->ischcv == 3) {
    cs_field_find_or_create("hybrid_blend",
                            CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                            CS_MESH_LOCATION_CELLS,
                            1,
                            false);
  }

  /* Friction velocity at the wall, in the case of a LES calculation
   * with van Driest-wall damping (delayed here rather than placed in
   * addfld, as idries may be set in cs_parameters_*_complete). */

  if (cs_glob_turb_model->itytur == 4 && cs_glob_turb_les_model->idries == 1) {
    cs_field_find_or_create("boundary_ustar",
                            CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                            CS_MESH_LOCATION_BOUNDARY_FACES,
                            1,
                            false);
  }

  if (vp_param->staggered == 1) {

    /* Head loss on interior faces */

    cs_field_create("inner_face_head_loss",
                    CS_FIELD_PROPERTY,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    1,
                    false);

    /* Head loss on boundary faces */

    cs_field_create("boundary_face_head_loss",
                    CS_FIELD_PROPERTY,
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    1,
                    false);

    /* Source term on interior faces */

    cs_field_create("inner_face_source_term",
                    CS_FIELD_PROPERTY,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    1,
                    false);
  }

  /* Interior mass flux field */

  bool previous_val = false;
  if (cs_glob_time_scheme->istmpf != 1 || vp_param->staggered == 1) {
    previous_val = true;
  }

  cs_field_t *f_imf = cs_field_create("inner_mass_flux",
                                      CS_FIELD_EXTENSIVE,
                                      CS_MESH_LOCATION_INTERIOR_FACES,
                                      1,
                                      previous_val);

  cs_field_set_key_int(f_imf, k_log, 0);
  cs_field_set_key_int(f_imf, k_vis, 0);

  /* Same mass flux for every variable, an other mass flux
   * might be defined hereafterwards */

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if ((f->type & CS_FIELD_VARIABLE)) {
      cs_field_set_key_int(f, k_imasf, f_imf->id);
    }
  }

  /* Rusanov flux */

  if (cs_glob_turb_rans_model->irijnu == 2) {
    cs_field_create("i_rusanov_diff",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    1,
                    false);

    cs_field_create("b_rusanov_diff",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    1,
                    false);
  }

  /* Godunov scheme flux */

  if (cs_glob_turb_rans_model->irijnu == 3) {

    cs_field_create("i_velocity",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    3,
                    false);

    cs_field_create("i_reynolds_stress",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    6,
                    false);

    cs_field_create("b_velocity",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    3,
                    false);

    cs_field_create("b_reynolds_stress",
                    CS_FIELD_EXTENSIVE,
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    6,
                    false);

  }

  /* Boundary mass flux field */

  previous_val = false;
  if (cs_glob_time_scheme->istmpf != 1 || vp_param->staggered == 1) {
    previous_val = true;
  }

  cs_field_t *f_bmf = cs_field_create("boundary_mass_flux",
                                      CS_FIELD_EXTENSIVE,
                                      CS_MESH_LOCATION_BOUNDARY_FACES,
                                      1,
                                      previous_val);

  cs_field_set_key_int(f_bmf, k_log, 0);
  cs_field_set_key_int(f_bmf, k_vis, 0);

  /* The same mass flux for every variable, an other mass flux
   * might be defined here afterwards */

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if ((f->type & CS_FIELD_VARIABLE)) {
      cs_field_set_key_int(f, k_bmasf, f_bmf->id);
    }
  }

  /* Add mass flux for scalar with a drift (one mass flux per class) */

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    const int iscdri = cs_field_get_key_int(f, k_dri);

    if (iscdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX) {

      /* Mass flux for the class on interior faces */

      char f_imf_name[60];
      snprintf(f_imf_name, 60, "inner_mass_flux_%s", f->name);

      cs_field_t *f_imf_d = cs_field_create(f_imf_name,
                                            CS_FIELD_PROPERTY,
                                            CS_MESH_LOCATION_INTERIOR_FACES,
                                            1,
                                            false);

      cs_field_set_key_int(f_imf_d, k_log, 0);

      /* Set the inner mass flux index */
      cs_field_set_key_int(f, k_imasf, f_imf_d->id);

      /* Mass flux for the class on boundary faces */

      char f_bmf_name[60];
      snprintf(f_bmf_name, 60, "boundary_mass_flux_%s", f->name);

      cs_field_t *f_bmf_d = cs_field_create(f_bmf_name,
                                            CS_FIELD_PROPERTY,
                                            CS_MESH_LOCATION_BOUNDARY_FACES,
                                            1,
                                            false);

      cs_field_set_key_int(f_bmf_d, k_log, 0);

      /* Set the inner mass flux index */
      cs_field_set_key_int(f, k_bmasf, f_bmf_d->id);

      /* Index of the class, all member of the class share the same mass flux */
      const int icla = cs_field_get_key_int(f, k_ccl);

      /* If the scalar is the representant of a class, then
       * set the mass flux index to all members of the class */
      if (icla != 0) {
        for (int fj_id = 0; fj_id < n_fld; fj_id++) {
          cs_field_t *fj    = cs_field_by_id(fj_id);
          const int   iclap = cs_field_get_key_int(fj, k_ccl);

          if (icla == iclap
              && ((fj->type & CS_FIELD_VARIABLE) == CS_FIELD_VARIABLE)) {
            cs_field_set_key_int(fj, k_imasf, f_imf_d->id);
            cs_field_set_key_int(fj, k_bmasf, f_bmf_d->id);
          }
        }
      }

      /* Get the scalar's output options
       * (except non-reconstructed boundary output) */

      int       iopchr = cs_field_get_key_int(f, k_vis);
      const int ilog   = cs_field_get_key_int(f, k_log);

      if (iopchr & CS_POST_BOUNDARY_NR) {
        iopchr -= CS_POST_BOUNDARY_NR;
      }

      /* If the mass flux is imposed, no need of drift_tau nor drift_vel */
      if (!(iscdri & CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX)) {
        /* Relaxation time */

        cs_field_t *f_dt = cs_field_create_by_composite_name
                             (f->name,
                              "drift_tau",
                              CS_FIELD_PROPERTY,
                              CS_MESH_LOCATION_CELLS,
                              1,
                              false);

        /* Set the same visualization options as the scalar */
        cs_field_set_key_int(f_dt, k_log, ilog);
        cs_field_set_key_int(f_dt, k_vis, iopchr);

        /* Drift velocity */

        cs_field_t *f_dv = cs_field_create_by_composite_name
                             (f->name,
                              "drift_vel",
                              CS_FIELD_PROPERTY,
                              CS_MESH_LOCATION_CELLS,
                              3,
                              false);

        /* Set the same visualization options as the scalar */
        cs_field_set_key_int(f_dv, k_log, ilog);
        cs_field_set_key_int(f_dv, k_vis, iopchr);
      }

      /* Interaction time particle--eddies */
      if (iscdri & CS_DRIFT_SCALAR_TURBOPHORESIS) {
        cs_field_t *f_ddt = cs_field_create_by_composite_name
                              (f->name,
                               "drift_turb_tau",
                               CS_FIELD_PROPERTY,
                               CS_MESH_LOCATION_CELLS,
                               1,
                               false);

        /* Set the same visualization options as the scalar */
        cs_field_set_key_int(f_ddt, k_log, ilog);
        cs_field_set_key_int(f_ddt, k_vis, iopchr);
      }
    }
  }

  /* Add various associated fields for variables */

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;

    const cs_equation_param_t *eqp_f = cs_field_get_equation_param_const(f);

    /* Add weight fields for variables to compute gradient */

    if (eqp_f->iwgrec == 1 && eqp_f->idiff > 0) {
      int idimf = -1;
      if (eqp_f->idften & CS_ISOTROPIC_DIFFUSION)
        idimf = 1;
      else if (eqp_f->idften & CS_ANISOTROPIC_DIFFUSION)
        idimf = 6;

      const int fl_id = cs_field_get_key_int(f, k_wgrec);

      if (fl_id > -1) {
        const cs_field_t *fl = cs_field_by_id(fl_id);
        if (fl->dim != idimf) {
          cs_parameters_error(CS_ABORT_IMMEDIATE,
                              _("initial data setup"),
                              _("Variable %s should be assigned a\n"
                                "gradient_weighting_field of dimension %d\n"
                                "but has already been assigned field %s \n"
                                "of dimension %d.\n"),
                              f->name,
                              idimf,
                              fl->name,
                              fl->dim);
        }
      }
      else {
        char f_name[128];
        snprintf(f_name, 127, "gradient_weighting_%s", f->name);
        f_name[127] = '\0';

        cs_field_t *f_gw = cs_field_create(f_name,
                                           0,
                                           CS_MESH_LOCATION_CELLS,
                                           idimf,
                                           false);
        cs_field_set_key_int(f, k_wgrec, f_gw->id);
      }
    }

    /* Postprocessing of slope tests */

    const int ifctsl = cs_field_get_key_int(f, k_slts);
    if (ifctsl != -1
        && eqp_f->iconv > 0 && eqp_f->blencv > 0 && eqp_f->isstpc == 0) {
      char f_name[128];
      snprintf(f_name, 127, "%s_slope_upwind", f->name);
      f_name[127] = '\0';

      cs_field_t *f_su = cs_field_create(f_name,
                                         CS_FIELD_POSTPROCESS,
                                         CS_MESH_LOCATION_CELLS,
                                         1,
                                         false);

      cs_field_set_key_int(f_su, k_vis, CS_POST_ON_LOCATION);
      cs_field_set_key_int(f, k_slts, f_su->id);
    }

    /* Diffusion limiter */

    const int ikdf = cs_field_get_key_int(f, k_dflim);
    if (ikdf != -1) {
      char f_name[128];
      snprintf(f_name, 60, "%s_diff_lim", f->name);
      f_name[127] = '\0';

      cs_field_t *f_dflim = cs_field_create(f_name,
                                            CS_FIELD_PROPERTY,
                                            CS_MESH_LOCATION_CELLS,
                                            f->dim,
                                            false);

      cs_field_set_key_int(f_dflim, k_log, 1);
      cs_field_set_key_int(f_dflim, k_vis, CS_POST_ON_LOCATION);
      cs_field_set_key_int(f, k_cvlim, f_dflim->id);
    }

    /* Convection limiter */

    const int icv = cs_field_get_key_int(f, k_cvlim);
    if (eqp_f->isstpc == 2 || icv != -1) {
      char f_name[128];
      snprintf(f_name, 127, "%s_conv_lim", f->name);
      f_name[127] = '\0';

      cs_field_t *f_cvlim = cs_field_create(f_name,
                                            CS_FIELD_PROPERTY,
                                            CS_MESH_LOCATION_CELLS,
                                            f->dim,
                                            false);

      cs_field_set_key_int(f_cvlim, k_log, 1);
      cs_field_set_key_int(f, k_cvlim, f_cvlim->id);
    }

  } /* End of loop on fields */

  /* Fan id visualization */

  cs_fan_field_create();

  /* VOF */

  cs_vof_field_create();

  /* Turbulent anisotropic viscosity or user defined tensor diffusivity
   * for a scalar (exclusive or).*/

  if (idfm == 1 || iggafm == 1
      || (cs_glob_turb_model->order == CS_TURB_SECOND_ORDER
          && cs_glob_turb_rans_model->idirsm == 1)) {

    cs_field_t *f_atv = cs_field_create("anisotropic_turbulent_viscosity",
                                        CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_CELLS,
                                        6,
                                        false);

    cs_field_set_key_int(f_atv, k_log, 0);

    if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM && iggafm == 1) {
      cs_field_t *f_atvs
        = cs_field_create("anisotropic_turbulent_viscosity_scalar",
                          CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                          CS_MESH_LOCATION_CELLS,
                          6,
                          false);

      cs_field_set_key_int(f_atvs, k_log, 0);
    }
  }

  /* Change some field settings
     -------------------------- */

  if (cs_glob_ale > 0) {
    cs_field_t *f_imasf
      = cs_field_by_id(cs_field_get_key_int(CS_F_(p), k_imasf));
    cs_field_set_n_time_vals(f_imasf, 2);
    cs_field_t *f_bmasf
      = cs_field_by_id(cs_field_get_key_int(CS_F_(p), k_bmasf));
    cs_field_set_n_time_vals(f_bmasf, 2);
  }

  /* Set some field keys
     ------------------- */

  /* Copy imrgra into the field structure if still at default */

  for (int f_id = 0; f_id < n_fld; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;

    cs_equation_param_t *eqp_f = cs_field_get_equation_param(f);

    if (eqp_f->imrgra < 0) {
      eqp_f->imrgra = cs_glob_space_disc->imrgra;
    }
  }

  /* Check if scalars are buoyant and set n_buoyant_scal accordingly.
   * It is then used in tridim to update buoyant scalars and density
   * in U-P loop */

  cs_velocity_pressure_set_n_buoyant_scalars();

  /* For Low Mach and compressible (increment) algorithms, particular care
   * must be taken when dealing with density in the unsteady term in the
   * velocity pressure loop */

  if (cs_glob_fluid_properties->irovar == 1
      && (cs_glob_velocity_pressure_model->idilat > 1
          || cs_glob_vof_parameters->vof_model > 0
          || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)) {

    /* EOS density, imposed after the correction step, so we need
     * to keep the previous one, which is in balance with the mass */

    cs_field_t *f_dm = cs_field_create("density_mass",
                                       CS_FIELD_PROPERTY,
                                       CS_MESH_LOCATION_CELLS,
                                       1,
                                       false);
    cs_field_set_key_int(f_dm, k_log, 0);
    cs_field_set_key_int(f_dm, k_vis, 0);

    cs_field_t *f_bdm = cs_field_create("boundary_density_mass",
                                        CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_BOUNDARY_FACES,
                                        1,
                                        false);
    cs_field_set_key_int(f_bdm, k_log, 0);
    cs_field_set_key_int(f_bdm, k_vis, 0);
  }

  /* Update field pointer mappings
     ----------------------------- */
  cs_field_pointer_map_base();
  cs_field_pointer_map_boundary();
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup computation based on provided user data and functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_setup(void)
{
  /* Initialize modules before user has access */
  _init_setup();

  int nmodpp = 0;
  for (int i = 1; i < CS_N_PHYSICAL_MODEL_TYPES; i++) {
    if (cs_glob_physical_model_flag[i] > -1)
      nmodpp++;
  }

  /* User input, variable definitions */
  _init_user(&nmodpp);

  cs_f_ppini1();

  /* Map Fortran pointers to C global data */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1)
    cs_at_data_assim_initialize();

  /* Initialize lagr structures */
  cs_lagr_map_specific_physics();

  int have_thermal_model = 0;
  if (cs_thermal_model_field() != NULL)
    have_thermal_model = 1;

  int is_restart = cs_restart_present();
  cs_real_t dtref = cs_glob_time_step->dt_ref;
  cs_time_scheme_t *time_scheme = cs_get_glob_time_scheme();

  cs_lagr_options_definition(is_restart,
                             have_thermal_model,
                             dtref,
                             &time_scheme->iccvfg);
  cs_lagr_add_fields();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    /* Additional fields if not in CDO mode only */
    _additional_fields_stage_2();

    /* Changes after user initialization and additional fields dependent on
     * main fields options. */
    cs_parameters_global_complete();

    _additional_fields_stage_3();
  }

  cs_parameters_eqp_complete();

  /* Time moments called after additional creation */
  cs_gui_time_moments();
  cs_user_time_moments();

  /* GUI based boundary condition definitions */
  cs_gui_boundary_conditions_define(NULL);

  /* Some final settings */
  cs_gui_output();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    /* Warning: Used in 0 validation cases ? */
    cs_f_usipes(&nmodpp);

    /* Avoid a second spurious call to this function
     * called in the C part if CDO is activated */
    if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_OFF) {
      cs_user_boundary_conditions_setup(cs_glob_domain);
      cs_user_finalize_setup(cs_glob_domain);
    }
  }

  cs_parameters_output_complete();

  /* Coherency checks */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_parameters_check();

  cs_log_printf(CS_LOG_DEFAULT,
                "\n"
                "No error detected during the data verification.\n");

  /* Print output */
  cs_f_impini();
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
