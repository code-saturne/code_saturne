/*============================================================================
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_restart_default.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "base/cs_theta_scheme.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 *
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Called at the very start of the time step
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var_stage1(void)
{
  /* Initialization */

  cs_dispatch_context ctx, ctx_i, ctx_b;

  const int key_t_ext_id = cs_field_key_id("time_extrapolated");

  cs_real_t *cpro_rho_mass = nullptr;
  cs_real_t *bpro_rho_mass = nullptr;

  if (cs_field_by_name_try("density_mass") != nullptr)
    cpro_rho_mass = cs_field_by_name_try("density_mass")->val;
  if (cs_field_by_name_try("boundary_density_mass") != nullptr)
    bpro_rho_mass = cs_field_by_name_try("boundary_density_mass")->val;

  const cs_velocity_pressure_param_t  *vp_param
    = cs_glob_velocity_pressure_param;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;
  cs_real_t *imasfl_pre = cs_field_by_id(iflmas)->val_pre;
  cs_real_t *bmasfl_pre = cs_field_by_id(iflmab)->val_pre;

  const cs_real_t *crom   = CS_F_(rho)->val;
  const cs_real_t *brom   = CS_F_(rho_b)->val;

  cs_field_t *f_cp = CS_F_(cp);

  const int time_order = cs_glob_time_scheme->time_order;

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;

  const int k_sca = cs_field_key_id("scalar_id");

  /*----------------------------------------------------------------------
   * Store the previous mass flux (n-1->n) in *_mass_flux_prev
   * Note: if there is no previous values, nothing is done for explicit
   * schemes (istmpf=0) a specific treatment is done because previous value
   * is used as a work array.
   *----------------------------------------------------------------------*/

  if (cs_glob_time_scheme->istmpf > 0 || vp_param->staggered == 1) {
    if (vp_param->itpcol == 0) {
      cs_field_current_to_previous(cs_field_by_id(iflmas));
      cs_field_current_to_previous(cs_field_by_id(iflmab));
    }
    else {
      ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        cs_real_t flux      = imasfl[f_id];
        imasfl[f_id]     = 2.*imasfl[f_id] - imasfl_pre[f_id];
        imasfl_pre[f_id] = flux;
      });
      ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        cs_real_t flux      = bmasfl[f_id];
        bmasfl[f_id]     = 2.*bmasfl[f_id] - bmasfl_pre[f_id];
        bmasfl_pre[f_id] = flux;
      });

      ctx_i.wait();
      ctx_b.wait();

    }
  }

  /*----------------------------------------------------------------------
   * If required, the density at time step n-1 is updated
   * Note that for VOF and dilatable algorithms, density at time step n-2
   * is also updated
   * Note that at the begining of the calculation, previous values have
   * been initialized by inivar
   *----------------------------------------------------------------------*/

  if (cs_glob_fluid_properties->irovar > 0) {
    if (time_order == 2 && vp_param->itpcol == 1) {
      const cs_real_t *croma  = CS_F_(rho)->val_pre;
      const cs_real_t *cromaa = CS_F_(rho)->vals[2];
      const cs_real_t *broma  = CS_F_(rho_b)->val_pre;
      const cs_real_t *bromaa = CS_F_(rho_b)->vals[2];

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_rho_mass[c_id] = 3.*crom[c_id] - 3.*croma[c_id]
                              + cromaa[c_id];
      });

      ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        bpro_rho_mass[f_id] = 3.*brom[f_id] - 3.*broma[f_id]
                              + bromaa[f_id];
      });
    }

    ctx.wait();
    ctx_b.wait();

    cs_field_current_to_previous(CS_F_(rho));
    cs_field_current_to_previous(CS_F_(rho_b));

    if (((cs_glob_velocity_pressure_model->idilat > 1
          || cs_glob_vof_parameters->vof_model > 0) && vp_param->itpcol == 0)
        || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3) {
      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_rho_mass[c_id] = crom[c_id];
      });
      ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        bpro_rho_mass[f_id] = brom[f_id];
      });

      ctx.wait();
      ctx_b.wait();

    }
  }

  /* Extrapolate viscosities and heat capacity if required */
  if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0)
    cs_field_current_to_previous(CS_F_(mu));
  if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0)
    cs_field_current_to_previous(CS_F_(mu_t));
  if (f_cp != nullptr) {
    if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0)
      cs_field_current_to_previous(f_cp);
  }

  /* Extrapolate scalars if required */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
      continue;
    int i_sca = cs_field_get_key_int(f, k_sca);
    if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0))
      cs_field_current_to_previous(f);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 *
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Called at just after cs_physical_properties_update (and thus before
 * cs_solve_navier_stokes)
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var_stage2(void)
{
  /* Initialization */

  cs_dispatch_context ctx;

  const int key_t_ext_id = cs_field_key_id("time_extrapolated");
  const int kthetvs = cs_field_key_id("diffusivity_extrapolated");

  cs_real_t *cpro_viscl  = CS_F_(mu)->val;
  cs_real_t *cproa_viscl = CS_F_(mu)->val_pre;
  cs_real_t *cpro_visct  = CS_F_(mu_t)->val;
  cs_real_t *cproa_visct = CS_F_(mu_t)->val_pre;

  cs_field_t *f_cp = CS_F_(cp);
  cs_real_t *cpro_cp = nullptr;
  cs_real_t *cproa_cp = nullptr;
  if (f_cp != nullptr) {
    cpro_cp  = f_cp->val;
    cproa_cp = f_cp->val_pre;
  }

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const int k_sca = cs_field_key_id("scalar_id");

  /*----------------------------------------------------------------------
   * Update previous values : only done if values are unavailable at the
   * first time step or when reading the restart file
   *----------------------------------------------------------------------*/

  if (cs_restart_present() == 1) {
    /* Update densities, viscosities and heat capacity if required */
    if (cs_glob_fluid_properties->irovar > 0) {
      if (cs_restart_get_field_read_status(CS_F_(rho)->id) == 0)
        cs_field_current_to_previous(CS_F_(rho));
      if (cs_restart_get_field_read_status(CS_F_(rho_b)->id) == 0)
        cs_field_current_to_previous(CS_F_(rho_b));
    }
    if (cs_restart_get_field_read_status(CS_F_(mu)->id) == 0)
      cs_field_current_to_previous(CS_F_(mu));
    if (cs_restart_get_field_read_status(CS_F_(mu_t)->id) == 0)
      cs_field_current_to_previous(CS_F_(mu_t));
    if (f_cp != nullptr) {
      if (cs_restart_get_field_read_status(f_cp->id) == 0)
        cs_field_current_to_previous(f_cp);
    }

    /* Update scalars diffusivity if required */
    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
        continue;
      int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
      if ((k_f_id > 0) && (cs_restart_get_field_read_status(k_f_id) == 0))
        cs_field_current_to_previous(cs_field_by_id(k_f_id));
    }
  }

  /*----------------------------------------------------------------------
   * Extrapolation of viscosities, heat capacity and scalars diffusivity
   * for the theta scheme. Values in the previous time step are kept before
   * the end the time iteration
   *----------------------------------------------------------------------*/

  /* Extrapolate viscosities and heat capacity */
  if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0) {
    const cs_real_t theta = cs_glob_time_scheme->thetvi;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t _viscol = cpro_viscl[c_id];
      cpro_viscl[c_id] = (1. + theta) * cpro_viscl[c_id]
                         - theta * cproa_viscl[c_id];
      cproa_viscl[c_id] = _viscol;
    });
  }
  if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0) {
    const cs_real_t theta = cs_glob_time_scheme->thetvi;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t _viscot = cpro_visct[c_id];
      cpro_visct[c_id] = (1. + theta) * cpro_visct[c_id]
                         - theta * cproa_visct[c_id];
      cproa_visct[c_id] = _viscot;
    });
  }
  if (f_cp != nullptr) {
    if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0) {
      const cs_real_t theta = cs_glob_time_scheme->thetcp;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t _cp = cpro_cp[c_id];
        cpro_cp[c_id] = (1. + theta) * cpro_cp[c_id]
                            - theta * cproa_cp[c_id];
        cproa_cp[c_id] = _cp;
      });
    }
  }

  ctx.wait();

  /* Extrapolate scalars diffusivity */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
      continue;
    int i_sca = cs_field_get_key_int(f, k_sca);

    /* Is a scalar and should be computed with theta scheme */
    if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0)) {

      const cs_real_t theta = cs_field_get_key_double(f, kthetvs);
      int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
      cs_real_t *cpro_visca = cs_field_by_id(k_f_id)->val;
      cs_real_t *cproa_visca = cs_field_by_id(k_f_id)->val_pre;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t _visca = cpro_visca[c_id];
        cpro_visca[c_id] = (1. + theta) * cpro_visca[c_id]
                           - theta * cproa_visca[c_id];
        cproa_visca[c_id] = _visca;
      });

      ctx.wait();

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 *
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Called just after cs_solve_navier_stokes, in loops for U/P and ALE
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var_stage3(void)
{
  /* Initialization */

  cs_dispatch_context ctx_i, ctx_b;

  const cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;

  const cs_velocity_pressure_param_t  *vp_param
    = cs_glob_velocity_pressure_param;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;
  cs_real_t *imasfl_pre = cs_field_by_id(iflmas)->val_pre;
  cs_real_t *bmasfl_pre = cs_field_by_id(iflmab)->val_pre;

  /*----------------------------------------------------------------------
   * Only update an unique mass flux.
   * For the second order scheme (istmpf=2), theta is set to 0.5.
   * For the explicit scheme (istmpf=0), the mass flux is set to the
   * previous one F(n) in the last iteration. Another modification will
   * be done in iappel = 4.
   * Nothing is done for standard scheme (istmpf=1).
   * In case of iterations in cs_solve_navier_stokes, the following process
   * is applied in all iterations except the last one if istmpf=0
   *----------------------------------------------------------------------*/

  if (cs_glob_time_scheme->istmpf == 2 && vp_param->itpcol == 1) {
    cs_real_t theta = 0.5;
    cs_real_t aa = 1./(2.-theta);
    cs_real_t bb = (1.-theta)/(2.-theta);

    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      imasfl[f_id] = aa * imasfl[f_id] + bb * imasfl_pre[f_id];
    });
    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      bmasfl[f_id] = aa * bmasfl[f_id] + bb * bmasfl_pre[f_id];
    });
  }
  else if (cs_glob_time_scheme->istmpf == 0) {
    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      imasfl[f_id] = imasfl_pre[f_id];
    });
    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      bmasfl[f_id] = bmasfl_pre[f_id];
    });
  }

  ctx_i.wait();
  ctx_b.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 *
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * After cs_solve_navier_stokes and loops for U/P and ALE
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var_stage4(void)
{
  /* Initialization */

  cs_dispatch_context ctx_i, ctx_b;

  const cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;
  cs_real_t *imasfl_pre = cs_field_by_id(iflmas)->val_pre;
  cs_real_t *bmasfl_pre = cs_field_by_id(iflmab)->val_pre;

  /* At the very start of the time step */

  /*----------------------------------------------------------------------
   * Only update an unique mass flux.
   * Nothing is done for standard and second order schemes (istmpf=1 and 2)
   * For the explicit scheme (istmpf=0), F(n+1) is saved in imasfl_pre but
   * computations are done with F(n) in imasfl
   * In case of iterations in cs_solve_navier_stokes, the following process
   * is applied in all iterations if istmpf differs from 0 and only in the
   * last iteration if istmpf=0. By doing so, from the second sub-iteration
   * the computation will be done with F(n+1) and no longer with F(n).
   * It is assumed here that the user chose to do sub-iterations also to
   * have an implicit mass flux
   *----------------------------------------------------------------------*/

  if (cs_glob_time_scheme->istmpf == 0) {
    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_real_t flux = imasfl[f_id];
      imasfl[f_id] = imasfl_pre[f_id];
      imasfl_pre[f_id] = flux;
    });
    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_real_t flux = bmasfl[f_id];
      bmasfl[f_id] = bmasfl_pre[f_id];
      bmasfl_pre[f_id] = flux;
    });

    ctx_i.wait();
    ctx_b.wait();

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 *
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Called just after cs_solve_transported_variables
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var_stage5(void)
{
  /* Initialization */

  cs_dispatch_context ctx, ctx_i, ctx_b;

  const int key_t_ext_id = cs_field_key_id("time_extrapolated");

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;
  cs_real_t *imasfl_pre = cs_field_by_id(iflmas)->val_pre;
  cs_real_t *bmasfl_pre = cs_field_by_id(iflmab)->val_pre;

  cs_real_t *cpro_viscl  = CS_F_(mu)->val;
  cs_real_t *cproa_viscl = CS_F_(mu)->val_pre;
  cs_real_t *cpro_visct  = CS_F_(mu_t)->val;
  cs_real_t *cproa_visct = CS_F_(mu_t)->val_pre;

  cs_field_t *f_cp = CS_F_(cp);
  cs_real_t *cpro_cp = nullptr;
  cs_real_t *cproa_cp = nullptr;
  if (f_cp != nullptr) {
    cpro_cp  = f_cp->val;
    cproa_cp = f_cp->val_pre;
  }

  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;

  const int k_sca = cs_field_key_id("scalar_id");


  /*----------------------------------------------------------------------
   * Previous modifications are reverted in order to be ready for the next
   * iteration in time.
   * Nothing is done for standard and second order schemes (istmpf=1 and 2)
   * For the explicit scheme (istmpf=0), imasfl is set to F(n+1)
   *----------------------------------------------------------------------*/

  /* Restore the mass flux */
  if (cs_glob_time_scheme->istmpf == 0) {
    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      imasfl[f_id] = imasfl_pre[f_id];
    });
    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      bmasfl[f_id] = bmasfl_pre[f_id];
    });
  }

  /* Restore physical properties */
  if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0)
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cpro_viscl[c_id] = cproa_viscl[c_id];
    });
  if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0)
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cpro_visct[c_id] = cproa_visct[c_id];
    });
  if (f_cp != nullptr) {
    if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_cp[c_id] = cproa_cp[c_id];
      });
    }
  }

  ctx.wait();
  ctx_i.wait();
  ctx_b.wait();

  /* Restore scalars diffusivity */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
      continue;
    int i_sca = cs_field_get_key_int(f, k_sca);

    /* Is a scalar and should be computed with theta scheme */
    if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0)) {

      int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
      cs_real_t *cpro_visca = cs_field_by_id(k_f_id)->val;
      cs_real_t *cproa_visca = cs_field_by_id(k_f_id)->val_pre;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_visca[c_id] = cproa_visca[c_id];
      });

      ctx.wait();

    }
  }
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
