/*============================================================================
 * Thermodynamic pressure and density.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_ale.h"
#include "atmo/cs_atmo_variables.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "base/cs_reducers.h"
#include "cfbl/cs_cf_model.h"
#include "cdo/cs_domain.h"
#include "elec/cs_elec_model.h"
#include "cdo/cs_equation.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_gas_mix.h"
#include "gui/cs_gui.h"
#include "base/cs_ht_convert.h"
#include "turb/cs_les_mu_t.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_thermal_model.h"
#include "turb/cs_turbulence_ml.h"
#include "turb/cs_turbulence_model.h"
#include "turb/cs_turbulence_ke.h"
#include "turb/cs_turbulence_kw.h"
#include "turb/cs_turbulence_rij.h"
#include "turb/cs_turbulence_sa.h"
#include "turb/cs_turbulence_v2f.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_vof.h"
#include "base/cs_wall_condensation.h"
#include "base/cs_wall_condensation_1d_thermal.h"
#include "base/cs_1d_wall_thermal.h"

#include "pprt/cs_physical_model.h"
#include "cogz/cs_combustion_ebu.h"
#include "cogz/cs_combustion_slfm.h"
#include "cogz/cs_combustion_physical_properties.h"
#include "comb/cs_coal_physical_properties.h"
#include "ctwr/cs_ctwr_physical_properties.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_physical_properties_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_physical_properties_default.cpp

  Compute physical properties which are variable in time.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes.
 *============================================================================*/

void
cs_f_physical_properties2(void);

void
cs_combustion_lw_physical_prop(void);

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Modify eddy viscosity using the Reboud correction:
 *         \f[
 *         \mu_t'= \dfrac{\rho_v + (1-\alpha)^{mcav}(\rho_l-\rho_v)}{\rho}\mu_t.
 *         \f]
 *
 * \param[in]  n_cells    Number of cells
 * \param[in]  rho1       reference density of fluid 1
 * \param[in]  rho2       reference density of fluid 2
 * \param[in]  crom       Density array
 *
 */
/*----------------------------------------------------------------------------*/

static void
_cavitation_correct_visc_turb(const cs_lnum_t  n_cells,
                              const cs_real_t  rho1,
                              const cs_real_t  rho2,
                              const cs_real_t  crom[])
{
  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *cvar_voidf = CS_F_(void_f)->val;

  const cs_real_t mcav = cs_get_glob_cavitation_parameters()->mcav;

  cs_dispatch_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t drho = (rho1 - rho2);
    const cs_real_t p_void = pow(1.0 - cvar_voidf[c_id], mcav);
    const cs_real_t rho_max = cs::max(crom[c_id], cs_math_epzero);

    const cs_real_t frho = (rho2 + p_void*drho) / rho_max;
    visct[c_id] = frho*visct[c_id];
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return error if density or molecular viscosity are  not constant
 *
 * \param[in]  name      name of density or molecular viscosity
 * \param[in]  n_elts    number of elements (face or cells)
 * \param[in]  val       values of density or molecular viscosity
 * \param[in]  val_ref   value reference
 */
/*----------------------------------------------------------------------------*/

static void
_field_is_constant(const char       *name,
                   const cs_lnum_t  n_elts,
                   const cs_real_t  val[],
                   const cs_real_t  val_ref)
{
  bool is_constant = true;

  cs_lnum_t ii = 0;
  for (ii = 0; ii < n_elts; ii++) {
    if (fabs(val[ii] - val_ref) <= cs_math_epzero)
      continue;
    is_constant = false;
    break;
  }

  if (!is_constant)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: abort in the physical quantities computation\n"
                "=====\n"
                "Incoherency between parameters for %s.\n"
                "%s has been declared constant\n"
                "but the field value has been modified\n"
                "and is not equal to the ref value anymore.\n"
                "ref=%e, value=%e.\n"
                "The calculation will not be run.\n"
                "Check the interface, cs_user_parameters\n"
                "and cs_user_physical_properties."),
              name, name, val_ref, val[ii]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the eddy viscosity
 *
 * * \param[in]  n_cells    number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_compute_turbulence_mu(const cs_lnum_t  n_cells)
{
  // Laminar
  if (cs_glob_turb_model->model == CS_TURB_NONE)
    cs_array_real_fill_zero(n_cells, CS_F_(mu_t)->val);

  // Mixing length model
  else if (cs_glob_turb_model->model == CS_TURB_MIXING_LENGTH)
    cs_turbulence_ml_mu_t();

  // k-epsilon
  else if (cs_glob_turb_model->itytur == 2)
    cs_turbulence_ke_mu_t(-1);

  else if (cs_glob_turb_model->order == CS_TURB_SECOND_ORDER
      && cs_glob_turb_model->type == CS_TURB_RANS)
    cs_turbulence_rij_mu_t(-1);

  // LES (Smagorinsky, dynamic Smagorinsky, or Wale)
  else if (cs_glob_turb_model->type == CS_TURB_LES) {
    if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_CONST)
      cs_les_mu_t_smago_const();
    else if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN)
      cs_les_mu_t_smago_dyn();
    else if (cs_glob_turb_model->model == CS_TURB_LES_WALE)
      cs_les_mu_t_wale();
    else if (cs_glob_turb_model->model == CS_TURB_LES_KSGS)
      cs_les_mu_t_ksgs();
    else if (cs_glob_turb_model->model == CS_TURB_LES_TAUSGS)
      cs_les_mu_t_tausgs();
  }

  // v2f (phi-model and BL-v2/k)
  else if (cs_glob_turb_model->itytur == 5) {
    if (cs_glob_turb_model->model == CS_TURB_V2F_PHI)
      cs_turbulence_v2f_phi_mu_t();
     else if (cs_glob_turb_model->model == CS_TURB_V2F_BL_V2K)
       cs_turbulence_v2f_bl_v2k_mu_t();
  }

  // k-omega SST
  else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA)
    cs_turbulence_kw_mu_t(-1);

  // Spalart-Allmaras
  else if (cs_glob_turb_model->model == CS_TURB_SPALART_ALLMARAS)
    cs_turbulence_sa_mu_t();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute anisotropic turbulent viscosity
 *
 *  \param[in]  n_scal      number of scalar field
 *  \param[in]  scalar_idx  id of scalar field
 *  \param[in]  mq          mesh quantities
 *  \param[in]  n_cells     number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_compute_anisotropic_turbulent_viscosity
(cs_lnum_t                   n_cells,
 int                         n_scal,
 const int                   scalar_idx[],
 const cs_mesh_quantities_t  *mq)
{
  bool idfm = false;
  bool iggafm = false;
  bool iebdfm = false;

  const int kturt = cs_field_key_id("turbulent_flux_model");

  for (int f_id = 0; f_id < n_scal; f_id++) {

    const cs_field_t *f = cs_field_by_id(scalar_idx[f_id]);

    const int turb_flux_model = cs_field_get_key_int(f, kturt);
    const int turb_flux_model_type = turb_flux_model / 10;
    if (turb_flux_model == 31)
      iebdfm = true;
    if (turb_flux_model_type == 3)
      idfm = true;
    if (turb_flux_model_type > 0)
      iggafm = true;

  }

  if (!(   (cs_glob_turb_rans_model->idirsm == 1)
        && ((idfm) || (cs_glob_turb_model->itytur == 3))))
    return;

  cs_real_6_t *visten
    = (cs_real_6_t *)cs_field_by_name
                       ("anisotropic_turbulent_viscosity")->val;

  if (cs_glob_turb_model->itytur == 3) {
    cs_turbulence_rij_anisotropic_mu_t(mq,
                                       n_cells,
                                       -1,
                                       idfm,
                                       iggafm,
                                       iebdfm);
  }
  else {
    cs_array_real_fill_zero(6*n_cells, (cs_real_t *)visten);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  clipping of density, viscosity and specific heat
 *
 * \param[in]  first_pass  first time passing
 * \param[in]  n_cells     number of cells
 * \param[in]  n_b_faces   number of boundary faces
 * \param[in]  n_scal      number of scalar field
 * \param[in]  scalar_idx  id of scalar field
 * \param[in]  mq          mesh quantities
 * \param[in]  brom        boundary density
 */
/*----------------------------------------------------------------------------*/

static void
_clip_rho_mu_cp(bool                         first_pass,
                cs_lnum_t                    n_cells,
                cs_lnum_t                    n_b_faces,
                int                          n_scal,
                int                          scalar_idx[],
                const cs_mesh_quantities_t  *mq,
                const cs_real_t              brom[])
{
  int iscacp = 0;
  constexpr int n_fields = 5; // number of fields for log
  const cs_field_t *fields[n_fields]
    = {CS_F_(rho), CS_F_(mu), CS_F_(mu_t), CS_F_(cp), CS_F_(rho_b)};

  char tmp_s[64] = "";
  if (CS_F_(cp) != nullptr) {
    const int kscacp  = cs_field_key_id("is_temperature");
    for (int f_id = 0; f_id < n_scal; f_id++) {
      const cs_field_t *f_scal = cs_field_by_id(scalar_idx[f_id]);
      iscacp = cs_field_get_key_int(f_scal, kscacp);
    }
  }

  const cs_lnum_t *c_disable_flag = nullptr;
  if (mq->has_disable_flag)
    c_disable_flag = mq->c_disable_flag;

  cs_dispatch_context ctx;

  // Set turbulent viscosity to 0 in disabled cells
  cs_real_t *visct = CS_F_(mu_t)->val;
  if (mq->has_disable_flag) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      if (c_disable_flag[c_id])
        visct[c_id] = 0.;
    });
  }

  // No use for logging of turbulence viscosity without turbulence
  if (cs_glob_turb_model->model == CS_TURB_NONE)
    fields[2] = nullptr;

  /* Min max */

  struct cs_double_n<2> rd[n_fields];
  struct cs_reduce_min_max_nr<1> reducer[n_fields];

  cs_real_t varminmax[n_fields*2];

  for (int ii = 0; ii < n_fields; ii++) {
    if (fields[ii] == nullptr)
      continue;
    cs_lnum_t n_elts = cs_mesh_location_get_n_elts(fields[ii]->location_id)[0];
    const cs_real_t *cpro_var = fields[ii]->val;
    ctx.parallel_for_reduce
      (n_elts, rd[ii], reducer[ii],
       [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_double_n<2> &res) {
         res.r[0] = cpro_var[c_id];
         res.r[1] = cpro_var[c_id];
       });
  }
  ctx.wait();

  for (int ii = 0; ii < n_fields; ii++) {
    if (fields[ii] == nullptr)
      continue;
    if (fields[ii] == nullptr) {
      varminmax[ii] = 0;
      varminmax[ii+n_fields] = 0;
    }
    else {
      varminmax[ii]          =  rd[ii].r[0];
      varminmax[ii+n_fields] = -rd[ii].r[1];
    }
  }

  cs_parall_min(n_fields*2, CS_REAL_TYPE, varminmax);

  // Loggging
  if (first_pass) {
    bft_printf(" -----------------------------------------\n"
               " Property           Min. value  Max. value\n"
               " -----------------------------------------\n");
    for (int ii = 0; ii < n_fields; ii++) {
      if (fields[ii] == nullptr)
        continue;
      int width = 16;
      cs_log_strpad(tmp_s, fields[ii]->name, width, 64);
      bft_printf(" %s   %10.04e  %10.04e\n", tmp_s,
                 varminmax[ii], -varminmax[n_fields+ii]);
    }
    bft_printf(" -----------------------------------------\n");
  }

  // verification of physical values
  for (int ii = 0; ii < n_fields; ii++) {
    if (varminmax[ii] >= 0.0)
      continue;

    /* we do not clip turbulent viscosity
       in dynamic LES model, because we have
       done clipping on the total viscosity */
    if (ii == 2 &&  cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN)
      continue;

    // for specific heat
    if (ii == 3 && iscacp == 0)
      continue;

    bft_error(__FILE__, __LINE__, 0,
              _("Abort in the physical quantities computation.\n\n"
                "The physical property %s is variable and the"
                "minimum reached is %10.12e\n"
                "Verify that this property has been defined and"
                "that the chosen law leads to correct values.\n\n"
                "For turbulent viscosity you can modify"
                " cs_user_physical_properties_turb_viscosity\n"),
              fields[ii]->name, varminmax[ii]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log and check scalar diffusivity min/max
 *
 *  \param[in]  first_pass  first time passing
 *  \param[in]  n_scal      number of scalar field
 *  \param[in]  scalar_idx  id of scalar field
 *  \param[in]  n_cells     number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_check_log_scalar_diff(const bool        first_pass,
                       const int         n_scal,
                       const int         scalar_idx[],
                       const cs_lnum_t   n_cells)
{
  if (n_scal < 1)
    return;

  bool wrote_header = false;
  int error_count = 0;
  constexpr int block_size = 32;
  cs_real_t visminmax[block_size*2];

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");

  char tmp_s[64] = "";

  const int n_blocks = (n_scal % block_size) ?
    n_scal/block_size+1 : n_scal/block_size;

  struct cs_double_n<2> rd[block_size];
  struct cs_reduce_min_max_nr<1> reducer[block_size];

  cs_dispatch_context ctx;

  /* Loop on blocks of fields to regroup possible MPI operations */
  for (int b_id = 0; b_id < n_blocks; b_id++) {

    // 1st Loop on fields
    int s_idx = 0, reduce_idx = 0;
    for (int s_id = b_id*block_size;
         s_id < n_scal && s_idx < block_size;
         s_id++, s_idx++) {
      const cs_field_t *f = cs_field_by_id(scalar_idx[s_id]);
      const int ifcvsl = cs_field_get_key_int(f, kivisl);

      if (ifcvsl >= 0) {
        const cs_real_t *cpro_vis = cs_field_by_id(ifcvsl)->val;

        ctx.parallel_for_reduce
          (n_cells, rd[reduce_idx], reducer[reduce_idx],
           [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_double_n<2> &res) {
             res.r[0] = cpro_vis[c_id];
             res.r[1] = cpro_vis[c_id];
           });
        reduce_idx++;
      }
    } // End 1st loop on fields

    // Parallel reduction if needed
    int reduce_count = reduce_idx;
    if (reduce_count > 0) {
      ctx.wait();

      for (int i = 0; i < reduce_count; i++) {
        visminmax[i] = rd[i].r[0];
        visminmax[i + reduce_count] = -rd[i].r[1];
      }

      cs_parall_min(reduce_count*2, CS_REAL_TYPE, visminmax);
    }

    // 2nd Loop on fields
    s_idx = 0, reduce_idx = 0;
    for (int s_id = b_id*block_size;
         s_id < n_scal && s_idx < block_size;
         s_id++, s_idx++) {
      const cs_field_t *f = cs_field_by_id(scalar_idx[s_id]);
      const int ifcvsl = cs_field_get_key_int(f, kivisl);

      const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

      cs_real_t vismax = -cs_math_big_r, vismin = cs_math_big_r;

      if (ifcvsl >= 0) {
        vismin = visminmax[reduce_idx];
        vismax = visminmax[reduce_idx + reduce_count];
        reduce_idx++;
      }
      else {
        const cs_real_t visls_0 = cs_field_get_key_double(f, kvisl0);
        vismax = visls_0;
        vismin = visls_0;
      }

      if (eqp->verbosity > 0 || first_pass || vismin <= 0.0) {
        if (!wrote_header) {
          wrote_header = true;
          bft_printf("\n --- Diffusivity:\n"
                     " -----------------------------------------------\n"
                     " Scalar           Index   Min. value  Max. value\n"
                     " -----------------------------------------------\n");
        }
        int width = 16;
        cs_log_strpad(tmp_s, cs_field_get_label(f), width, 64);
        bft_printf(" %s     %d   %10.04e  %10.04e\n",
                   tmp_s, s_id, vismin, vismax);
      }

      if (vismin < 0)
        error_count += 1;

    } // 2nd loop on scalars

    assert(reduce_idx == reduce_count);

  } // Loop on blocks

  if (wrote_header)
    bft_printf(" -----------------------------------------------\n");

  // Physical value checks

  if (error_count > 0) {
      bft_error
        (__FILE__, __LINE__, 0,
         _("Abort in the physical quantities computation.\n\n"
           "The diffusivity of %d variables contain negative values "
           " (see log above).\n"
           "Verify that these properties have been defined and"
           "that the chosen law leads to correct values."),
         error_count);
  }

  const cs_field_t *th_f = cs_field_by_name_try("thermal_expansion");
  if (th_f != nullptr) {

    const cs_real_t *cpro_beta = th_f->val;

    ctx.parallel_for_reduce
      (n_cells, rd[0], reducer[0],
       [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_double_n<2> &res) {
         res.r[0] = cpro_beta[c_id];
         res.r[1] = cpro_beta[c_id];
       });
    ctx.wait();

    cs_real_t varmn = rd[0].r[0];
    cs_parall_min(1, CS_REAL_TYPE, &varmn);

    if (varmn < 0.0) {
      bft_error(__FILE__, __LINE__, 0,
                _("Abort in the physical quantities computation\n\n"
                  "The physical property %s is variable and the"
                  "minimum reached is %g, which is incorrect.\n"
                  "Verify that this property has been defined and"
                  "that the chosen law leads to correct values."),
                th_f->name, varmn);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  check and log mesh velocity diffusivity
 *
 *  \param[in]  first_pass  first time passing
 *  \param[in]  n_cells     number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_check_log_mesh_diff(bool       first_pass,
                     cs_lnum_t  n_cells)
{
  if (   cs_glob_ale == CS_ALE_NONE
      || cs_glob_time_step->nt_cur != cs_glob_time_step->nt_prev + 1)
    return;

  bool ok = false;
  char tmp_s[64] = "";
  cs_equation_param_t *eqp = cs_equation_param_by_name("mesh_velocity");

  if (eqp->idften & CS_ANISOTROPIC_LEFT_DIFFUSION) {
    const cs_real_6_t *cpro_visma_v = (const cs_real_6_t *)CS_F_(vism)->val;
    for (int ii = 0; ii < 6; ii++) {
      cs_real_t varmx = cpro_visma_v[0][ii];
      cs_real_t varmn = cpro_visma_v[0][ii];
      for (cs_lnum_t c_id = 1; c_id < n_cells; c_id++) {
        varmx = cs::max(varmx, cpro_visma_v[c_id][ii]);
        varmn = cs::min(varmn, cpro_visma_v[c_id][ii]);
      }

      cs_parall_max(1, CS_REAL_TYPE, &varmx);
      cs_parall_min(1, CS_REAL_TYPE, &varmn);

      if ( (eqp->verbosity > 0) || (first_pass) || (varmn < 0.0)) {
        if (!ok) {
          ok = true;
          bft_printf(" --- Mesh viscosity (ALE method)\n"
                     " -----------------------------------------\n"
                     " Property           Min. value  Max. value\n"
                     " -----------------------------------------\n");
        }
        int width = 16;
        cs_log_strpad(tmp_s, CS_F_(vism)->name, width, 64);
        bft_printf("  %s           %10.12e       %10.12e",
                   tmp_s, varmn, varmx);
      }

      // Physical value checks
      if (varmn < 0.0) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Warning: abort in the physical quantities computation\n"
                  "========\n"
                  "The mesh viscosity has not been correctly defined.\n"
                  "The calculation will not be run.\n"
                  "The minimum reached is %10.12e\n"
                  "Verify he definition of this property"), varmn);
      }
    }
  }
  else if (eqp->idften & CS_ISOTROPIC_DIFFUSION) {
    const cs_real_t *cpro_visma_s = CS_F_(vism)->val;

    cs_real_t varmx = cpro_visma_s[0];
    cs_real_t varmn = cpro_visma_s[0];

    for (cs_lnum_t c_id = 1; c_id < n_cells; c_id++) {
      varmx = cs::max(varmx, cpro_visma_s[c_id]);
      varmn = cs::min(varmn, cpro_visma_s[c_id]);
    }

    cs_parall_max(1, CS_REAL_TYPE, &varmx);
    cs_parall_min(1, CS_REAL_TYPE, &varmn);

    if ( (eqp->verbosity > 0) || (first_pass) || (varmn < 0.0)) {
      if (!ok) {
        ok = true;
        bft_printf(" --- Mesh viscosity (ALE method)\n"
                   " -----------------------------------------\n"
                   " Property           Min. value  Max. value\n"
                   " -----------------------------------------\n");
      }
      bft_printf("  %s           %10.12e       %10.12e",
                 cs_field_get_label(CS_F_(vism)), varmn, varmx);
    }

    // Physical value checks
    if (varmn < 0.0)
      bft_error(__FILE__, __LINE__, 0,
                _("Abort in the physical quantities computation\n\n"
                  "The mesh viscosity has not been correctly defined.\n"
                  "The calculation will not be run.\n"
                  "The minimum reached is %10.12e\n"
                  "Verify he definition of this property"), varmn);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize boundary temperature
 */
/*----------------------------------------------------------------------------*/

static void
_init_boundary_temperature(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  if (cs_field_by_name_try("boundary_temperature") == nullptr)
    return;

  cs_field_t *fld = cs_field_by_name_try("temperature");
  cs_real_t *field_s_b = cs_field_by_name_try("boundary_temperature")->val;

  cs_dispatch_context ctx;

  if (fld != nullptr) {

    const cs_real_t *field_s_v = fld->val;
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      if (field_s_b[face_id] <= -cs_math_big_r)
        field_s_b[face_id] = field_s_v[b_face_cells[face_id]];
    });

  }
  else if (   cs_glob_thermal_model->thermal_variable
           == CS_THERMAL_MODEL_ENTHALPY) {

    fld = cs_field_by_name_try("enthalpy");
    if (fld != nullptr) {

      const cs_real_t *field_s_v = fld->val;

      cs_real_t *ttmp = nullptr;   /* n_cells should be sufficient ? */
      CS_MALLOC(ttmp, n_cells_ext, cs_real_t);
      cs_ht_convert_h_to_t_cells(field_s_v, ttmp);

      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        if (field_s_b[face_id] <= -cs_math_big_r)
          field_s_b[face_id] = ttmp[b_face_cells[face_id]];
      });
      ctx.wait();

      CS_FREE(ttmp);

    }

  } /* Enthalpy */

  // Last resort
  if (fld != nullptr) {
    const cs_real_t t0 = cs_glob_fluid_properties->t0;
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      if (field_s_b[face_id] <= -cs_math_big_r)
        field_s_b[face_id] = t0;
    });
  }

  ctx.wait();

  // For wall condensation, initialize to user-prescribed value
  if (cs_glob_wall_condensation->icondb < 0)
    return;

  // TODO move this to cs_wall_condensation.cpp

  const cs_real_t *ztpar = cs_glob_wall_condensation->ztpar;

  const cs_lnum_t nfbpcd = cs_glob_wall_condensation->nfbpcd;
  const cs_lnum_t *ifbpcd = cs_glob_wall_condensation->ifbpcd;
  const cs_lnum_t *iztag1d = cs_glob_wall_condensation->iztag1d;
  const cs_lnum_t *izzftcd = cs_glob_wall_condensation->izzftcd;
  const cs_real_t *ztpar0 = cs_glob_wall_cond_1d_thermal->ztpar0;

  const cs_lnum_t nfpt1d = cs_glob_1d_wall_thermal->nfpt1d;
  const cs_lnum_t *ifpt1d = cs_glob_1d_wall_thermal->ifpt1d;
  const cs_real_t *tppt1d = cs_glob_1d_wall_thermal->tppt1d;

  cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

  h_ctx.parallel_for(nfbpcd, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    const cs_lnum_t iz = izzftcd[ii];
    const cs_lnum_t face_id = ifbpcd[ii];
    if (iztag1d[iz] == 0)
      field_s_b[face_id] = ztpar[iz];
    else if (iztag1d[iz] == 1)
     field_s_b[face_id] = ztpar0[iz];
    else {
      for (cs_lnum_t jj = 0; jj < nfpt1d; jj++)
        if (ifpt1d[jj] == face_id)
          field_s_b[face_id] = tppt1d[jj];
    }
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief First computation of physical properties for specific models.
 *
 * This is called before the user sets/modifies physical properties which
 * are variable in time
 *
 * \param[in]     iterns     Navier-Stokes sub-iterations indicator:
 *                           - if strictly negative, indicate that this
 *                                function is called outside Navier-Stokes loop
 *                           - if positive, Navier-Stokes iteration number.
 * \param[in, out]   mbrom   filling indicator of romb
 */
/*----------------------------------------------------------------------------*/

static void
_physical_properties_update_models_stage_1(int   iterns,
                                           int  *mbrom)
{
  const int *pm_flag = cs_glob_physical_model_flag;

  if (pm_flag[CS_PHYSICAL_MODEL_FLAG] <= 0)
    return;

  /* After this point, models considered are mutually exclusive */

  if (pm_flag[CS_COMBUSTION_3PT] >= 0) {
    cs_combustion_physical_properties_update_d3p();
    *mbrom = 1;
  }
  else if (pm_flag[CS_COMBUSTION_SLFM] >= 0) {
    cs_combustion_slfm_physical_properties(iterns);
    *mbrom = 1;
  }
  else if (pm_flag[CS_COMBUSTION_EBU] >= 0) {
    cs_combustion_ebu_physical_prop(mbrom);
  }
  else if (pm_flag[CS_COMBUSTION_LW] >= 0) {
    cs_combustion_lw_physical_prop();
  }

  if (pm_flag[CS_COMBUSTION_COAL] >= 0)
    cs_coal_physprop(mbrom);

  else if (   pm_flag[CS_JOULE_EFFECT] >= 1
           || pm_flag[CS_ELECTRIC_ARCS] >= 1) {
    /* - For Joule effect, the user must define property laws
     *   (density, ...).
     * - For electric arcs, properties are interpolated from tabulated data. */
    cs_elec_physical_properties(cs_glob_domain);
  }

  else if (pm_flag[CS_COOLING_TOWERS] != -1) {
    const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
    cs_ctwr_phyvar_update(fp->ro0, fp->t0, fp->p0);
  }

  /* Atmospheric Flows (except constant density) */
  else if (pm_flag[CS_ATMOSPHERIC] >= 1) {
    cs_atmo_physical_properties_update();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief First computation of physical properties for specific models.
 *
 * This is called before the user sets/modifies physical properties which
 * are variable in time
 */
/*----------------------------------------------------------------------------*/

static void
_physical_properties_update_models_stage_2(void)
{
  if (cs_glob_physical_model_flag[CS_GAS_MIX] >= 0)
    cs_gas_mix_physical_properties();

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0)
    cs_cf_physical_properties();
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fills physical properties which are variable in time
 *
 * \param[in]     iterns     Navier-Stokes sub-iterations indicator:
 *                           - if strictly negative, indicate that this
 *                                function is called outside Navier-Stokes loop
 *                           - if positive, Navier-Stokes iteration number.
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_properties_update(int   iterns)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_field_t *rho_b_f = CS_F_(rho_b);

  int mbrom = 0;
  static bool first_pass = true;

  // Densities at boundaries are computed in cs_vof_compute_linear_rho_mu for VoF
  if (cs_glob_vof_parameters->vof_model > 0)
    mbrom = 1;

  // First computation of physical properties for specific physics
  // BEFORE the user
  _physical_properties_update_models_stage_1(iterns,
                                             &mbrom);

  /* Interface code_saturne
     ---------------------- */
  cs_gui_physical_variable();

  if (mbrom == 0 && n_b_faces > 0 && rho_b_f != nullptr)
    rho_b_f->val[0] = -cs_math_big_r;

  if (cs_glob_thermal_model->thermal_variable == CS_THERMAL_MODEL_ENTHALPY)
    cs_ht_convert_h_to_t_cells_solid();

  cs_user_physical_properties(cs_glob_domain);

  if (mbrom == 0 && n_b_faces > 0 && rho_b_f != nullptr) {
    if (rho_b_f->val[0] > -cs_math_big_r)
      mbrom = 1;
  }

  // Finalization of physical properties for specific physics
  // AFTER the user
  _physical_properties_update_models_stage_2();

  // Boundary density based on adjacent cell value if not explicitly set.
  if (mbrom == 0 && rho_b_f != nullptr) {
    assert(CS_F_(rho) != nullptr);
    const cs_real_t  *crom = CS_F_(rho)->val;
    cs_real_t  *brom = rho_b_f->val;

    cs_dispatch_context ctx;
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      brom[face_id] = crom[b_face_cells[face_id]];
    });
    ctx.wait();
  }

  // Only density may be updated in Navier-Stokes loop
  if (iterns > 0)
    return;

  /* At the first time step of the calculation
   *     If we indicated that rho (visc) was constant
   *       and we modified it in cs_user_physical_properties, it doesn't work
   *     We use irovar (ivivar) to write and read
   *       rho (visc) in the file continued */
  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1) {
    if (cs_glob_fluid_properties->irovar == 0) {
      const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
      _field_is_constant("the density",
                          n_cells,
                          CS_F_(rho)->val,
                          ro0);

      _field_is_constant("the density at the boundary",
                          n_b_faces,
                          CS_F_(rho_b)->val,
                          ro0);
    }
    if (cs_glob_fluid_properties->ivivar == 0) {
      const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
      _field_is_constant("The molecular viscosity",
                         n_cells,
                         CS_F_(mu)->val,
                         viscl0);
    }
  }

  /* Compute the eddy viscosity
     -------------------------- */
  _compute_turbulence_mu(n_cells);

  /* storage id of scalar fields */
  const int keysca = cs_field_key_id("scalar_id");
  const int n_fields = cs_field_n_fields();

  int n_scal = 0;
  int *scalar_idx = nullptr;
  CS_MALLOC(scalar_idx, n_fields, int);

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    const int sc_id = cs_field_get_key_int(f, keysca) - 1;
    if (sc_id < 0)
      continue;
    scalar_idx[n_scal] = f_id;
    n_scal ++;
  }

  /* Anisotropic turbulent viscosity (symmetric)
     ------------------------------------------- */
  _compute_anisotropic_turbulent_viscosity(n_cells, n_scal, scalar_idx, mq);

  /* Eddy viscosity correction for cavitating flows
     ---------------------------------------------- */
  const int i_vof_mass_transfer
    = (cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER);
  const cs_cavitation_parameters_t *cavit_param
    = (const cs_cavitation_parameters_t *)cs_get_glob_cavitation_parameters();
  if (i_vof_mass_transfer != 0 && cavit_param->icvevm == 1) {
    if (   (cs_glob_turb_model->itytur == 2)
        || (cs_glob_turb_model->itytur == 5)
        || (cs_glob_turb_model->model  ==  CS_TURB_K_OMEGA)
        || (cs_glob_turb_model->model  ==  CS_TURB_SPALART_ALLMARAS) ) {
    _cavitation_correct_visc_turb(n_cells,
                                  cs_glob_vof_parameters->rho1,
                                  cs_glob_vof_parameters->rho2,
                                  CS_F_(rho)->val);
    }
  }

  /* User modification of the turbulent
     viscosity and symmetric tensor diffusivity
     ------------------------------------------ */
  cs_user_physical_properties_turb_viscosity(cs_glob_domain);

  /* Rusanov flux
     ------------ */
  cs_turbulence_rij_compute_rusanov();

  // Physical value checks
  if (rho_b_f != nullptr)
    _clip_rho_mu_cp(first_pass, n_cells, n_b_faces, n_scal, scalar_idx, mq,
                    rho_b_f->val);

  // Calculation of scalar limits and printing
  _check_log_scalar_diff(first_pass, n_scal, scalar_idx, n_cells);

  CS_FREE(scalar_idx);

  // Calculation of mesh viscosity bounds in ALE
  _check_log_mesh_diff(first_pass, n_cells);

  /* Initialize boundary temperature if present and not initialized yet
     ------------------------------------------------------------------ */

  _init_boundary_temperature();

  first_pass = false;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
