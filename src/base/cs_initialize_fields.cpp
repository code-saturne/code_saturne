/*============================================================================
 * Initialize fields.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_cf_thermo.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_ibm.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_porosity_from_scan.h"
#include "cs_scalar_clipping.h"
#include "cs_time_step.h"
#include "cs_turbulence_init.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_initialize_fields.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_initialize_fields.c

  \brief Various field initializations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computed variable and property initialization.
 *
 * The time step, the indicator of wall distance computation are also
 * initialized just before reading a restart file or use the user
 * initializations.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_0(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const int n_fields = cs_field_n_fields();

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kdflim = cs_field_key_id("diffusion_limiter_id");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");

  /* Initialize temperature to reference temperature if present */
  if (CS_F_(t) != nullptr) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->t0,
                             CS_F_(t)->val);
  }

  /* Initialize boundary temperature to "marker" if present */
  if (CS_F_(t_b) != nullptr) {
    cs_array_real_set_scalar(m->n_b_faces,
                             -cs_math_big_r,
                             CS_F_(t_b)->val);
  }

  /* Time step
     --------- */

  /* dt might be used on the halo cells during the ALE initialization
     otherwise dt is synchronized in the pressure correction step. */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           cs_glob_time_step->dt_ref,
                           CS_F_(dt)->val);

  /* Initialize physical properties
     ------------------------------ */

  /* Density */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           fluid_props->ro0,
                           CS_F_(rho)->val);
  cs_array_real_set_scalar(m->n_b_faces,
                           fluid_props->ro0,
                           CS_F_(rho_b)->val);

  /* Boussinesq (set to impossible value at this stage) */
  if (cs_glob_velocity_pressure_model->idilat == 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             -1.,
                             cs_field_by_name("thermal_expansion")->val);
  }

  /* Molecular viscosity */
  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           fluid_props->viscl0,
                           CS_F_(mu)->val);

  /* Specific heat */
  if (fluid_props->icp >= 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->cp0,
                             CS_F_(cp)->val);
  }

  /* The total pressure will be initialized to P0 + rho.g.r in inivar
     if the user has not done a prior initialization;
     so mark to "uninit" value here */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             - cs_math_infinite_r,
                             cs_field_by_name("total_pressure")->val);
  }

  /* Initialization of mix_mol_mas with default values (air)
     (used in cs_cf_thermo_default_init) */
  if (cs_glob_physical_model_flag[CS_GAS_MIX] >= 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->xmasmr,
                             cs_field_by_name("mix_mol_mas")->val);
  }

  /* Default initialisations for the compressible model */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
    /* In compressible, for now, the temperature is not solved but is a field of
       type variable anyway. The reference value has to be taken into account. */
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->t0,
                             CS_F_(t)->val);

    /* Default isochoric specific heat (cv0), total energy and density */
    cs_cf_thermo_default_init();

    /* Default diffusivity for total energy */
    double visls_0 = cs_field_get_key_double(CS_F_(t), kvisl0);
    visls_0 /= fluid_props->cv0;
    cs_field_set_key_double(CS_F_(e_tot), kvisl0, visls_0);
  }

  /* Scalars diffusivity */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = cs_field_get_key_int(f, keysca) - 1;
    if (scalar_id < 0)
      continue;

    int ifcvsl = cs_field_get_key_int(f, kivisl);
    double visls_0 = cs_field_get_key_double(f, kvisl0);
    /* Diffusivity at current time step (and previous one if second order) */
    if (ifcvsl >= 0) {
      cs_field_t *f_vsl = cs_field_by_id(ifcvsl);
      for (int i = 0; i < f_vsl->n_time_vals; i++) {
        cs_array_real_set_scalar(m->n_cells_with_ghosts,
                                 visls_0,
                                 f_vsl->vals[i]);
      }
    }
  }

  /* Mesh viscosity for ALE */

  if (cs_glob_ale != CS_ALE_NONE) {
    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(CS_F_(mesh_u));
    const int idftnp = eqp->idften;

    if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
      cs_real_t ref_val[6] = {1., 1., 1., 0., 0., 0.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              6,
                              ref_val,
                              cs_field_by_name("mesh_viscosity")->val);
    }
    else if (idftnp & CS_ISOTROPIC_DIFFUSION) {
      cs_array_real_set_scalar(m->n_cells_with_ghosts,
                               1,
                               cs_field_by_name("mesh_viscosity")->val);
    }

    cs_gui_mesh_viscosity();
  }

  /* Porosity */

  if (cs_glob_porous_model > 0) {
    cs_real_t *porosity = cs_field_by_name("porosity")->val;

    if (cs_glob_porosity_from_scan_opt->compute_porosity_from_scan)
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 0., porosity);

    else if (cs_glob_porosity_ibm_opt->porosity_mode > 0) {
      cs_field_t *f_ifp = cs_field_by_name_try("i_face_porosity");
      cs_field_t *f_bfp = cs_field_by_name_try("b_face_porosity");
      if (f_ifp != nullptr)
        cs_array_real_set_scalar(m->n_i_faces, 1., f_ifp->val);
      if (f_bfp != nullptr)
        cs_array_real_set_scalar(m->n_b_faces, 1., f_bfp->val);
    }

    else
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 1., porosity);

    /* Tensorial porosity */
    if (cs_glob_porous_model == 2) {
      cs_real_t ref_val[6] = {1., 1., 1., 0., 0., 0.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              6,
                              ref_val,
                              cs_field_by_name("tensorial_porosity")->val);
    }
  }

  /* Reference value for P*
     ---------------------- */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           cs_glob_fluid_properties->pred0,
                           CS_F_(p)->val);

  /* VoF
     --- */

  /* Note that the initialization here should not be strictly
     necessary, as cs_scalar_clipping should handle it, but
     doing so here will avoid activating the clipping
     counters; maybe a better manner of doing things would
     be for the counters to be deactivated at this step,
     which comes before user settings */

  if (cs_glob_vof_parameters->vof_model > 0) {
    const int kscmin = cs_field_key_id_try("min_scalar_clipping");

    cs_field_t *f_vf = cs_field_by_name("void_fraction");
    cs_real_t clvfmn = cs_field_get_key_double(f_vf, kscmin);

    cs_array_real_set_scalar(m->n_cells_with_ghosts, clvfmn, f_vf->val);
  }

  /* Turbulence
     ---------- */

  cs_turbulence_init_by_ref_quantities();

  /* Clip scalar variables (except k-epsilon, see above)
     --------------------------------------------------- */

  /* Clipping of non-variance scalars */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;
    if (cs_field_get_key_int(f, kscavr) >= 0)  /* is a variance */
      continue;

    cs_scalar_clipping(f);
  }

  /* Clipping of variances  */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;

    if (cs_field_get_key_int(f, kscavr) >= 0)
      cs_scalar_clipping(f);
  }

  /* ALE
     --- */

  if (cs_glob_ale > 0) {
    cs_field_t *f_coord0 = cs_field_by_name("vtx_coord0");
    cs_array_real_copy(m->n_vertices*3, m->vtx_coord, f_coord0->val);
  }

  /* Specific field initializations
     ------------------------------ */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    /* Diffusion limiters;
       Note that allowing a default value initialization for 1D-fields
       would allow handling this upon field definition */

    const int dl_id = cs_field_get_key_int(f, kdflim);
    if (dl_id > -1) {
      cs_field_t *f_dl = cs_field_by_id(dl_id);
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 1.0, f_dl->val);
    }

  }

  /* Previous value initializations
     ------------------------------

     This is not restricted to variables
     (for example Cp, Cv, and mass fluxes may be handled here also,
     depending on the time stepping options). */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    for (int n_p = f->n_time_vals; n_p > 1; n_p--)
      cs_field_current_to_previous(f);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
