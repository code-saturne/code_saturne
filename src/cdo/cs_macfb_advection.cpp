/*============================================================================
 * Build discrete convection operators for MAC schemes
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
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_bc.h"
#include "cs_flag.h"
#include "cs_macfb_builder.h"
#include "cs_math.h"
#include "cs_property.h"
#include "cs_scheme_geometry.h"

#if defined(DEBUG) && !defined(NDEBUG) /* For debugging purpose */
#include "cs_dbg.h"
#include "cs_log.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_advection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_macfb_advection.c

  \brief Build discrete advection operators for MAC face-based schemes

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_MACFB_ADVECTION_DBG 0

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Select face index for upwind scheme for a dual face
 *
 * \param[in]      f_idx   index of the face (inside the control volume)
 * \param[in]      fe_idx  index of the other face
 * \param[in]      flux    signed flux at the dual face
 */
/*----------------------------------------------------------------------------*/

static inline short int
_get_face_idx_upw(const short int f_idx,
                  const short int fe_idx,
                  const cs_real_t flux)
{
  return (flux >= 0.) ? f_idx : fe_idx;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform preprocessing such as the computation of the advection flux
 *         at the expected location in order to be able to build the advection
 *         matrix. Follow the prototype given by cs_macfb_adv_open_hook_t
 *         Default case.
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      macb     pointer to a cs_macfb_builder_t structure
 * \param[in]      csys     pointer to a cs_cell_sys_t structure
 * \param[in, out] input    nullptr or pointer to a structure cast on-the-fly
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_open_default(const cs_equation_param_t *eqp,
                                const cs_cell_mesh_t      *cm,
                                const cs_macfb_builder_t  *macb,
                                const cs_cell_sys_t       *csys,
                                void                      *input,
                                cs_cell_builder_t         *cb)
{
  CS_NO_WARN_IF_UNUSED(csys);
  CS_NO_WARN_IF_UNUSED(input);

  assert(eqp->adv_extrapol == CS_PARAM_ADVECTION_EXTRAPOL_NONE);

  /* Compute the flux across the primal faces. Store in cb->adv_fluxes */

  cs_advection_field_macb_dface_flux(
    macb, eqp->adv_field, cb->t_bc_eval, cb->adv_fluxes);

  if (eqp->adv_scaling_property != nullptr) {

    /* Loop on inner faces */
    for (short int fi = 0; fi < 6; fi++) {

      cs_real_t scaling = eqp->adv_scaling_property->ref_value;
      if (cs_property_is_uniform(eqp->adv_scaling_property) == false) {

        scaling = cs_property_get_face_value(
          cm->f_ids[fi], cb->t_pty_eval, eqp->adv_scaling_property);
      }

      cb->adv_fluxes[fi] *= scaling;

      /* Loop on outer faces */
      for (short int fj = 0; fj < 4; fj++) {
        const short int shift_j = 4 * fi + fj;

        cb->adv_fluxes[6 + shift_j] *= scaling;
      }
    }

  } /* Apply a scaling factor to the advective flux */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_macfb_adv_close_hook_t
 *         Default vector-valued case.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_close_default_vect(const cs_cell_mesh_t     *cm,
                                      const cs_macfb_builder_t *macb,
                                      const cs_cell_builder_t  *cb,
                                      cs_cell_sys_t            *csys)
{
  /* Add the local convection operator to the local system */

  cs_sdm_add_block_topleft(csys->mat, cm->n_fc, macb->n_dofs, cb->loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operation done after the matrix related to the advection term has
 *         been defined.
 *         Follow the prototype given by cs_macfb_adv_close_hook_t
 *         Explicit treatment without extrapolation for vector-valued DoFs.
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_close_exp_none_vect(const cs_cell_mesh_t     *cm,
                                       const cs_macfb_builder_t *macb,
                                       const cs_cell_builder_t  *cb,
                                       cs_cell_sys_t            *csys)
{
  CS_NO_WARN_IF_UNUSED(macb);
  CS_NO_WARN_IF_UNUSED(cm);
  CS_NO_WARN_IF_UNUSED(csys);
  CS_NO_WARN_IF_UNUSED(cb);

  bft_error(__FILE__, __LINE__, 0, "%s: not implemented.", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the cellwise advection operator for MAC-Fb schemes
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          Case of an advection term without a diffusion operator. In this
 *          situation, a numerical issue may arise if an internal or a border
 *          face is such that there is no advective flux. A special treatment
 *          is performed to tackle this issue.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_no_diffusion(const cs_equation_param_t *eqp,
                                const cs_cell_mesh_t      *cm,
                                const cs_macfb_builder_t  *macb,
                                cs_macfb_adv_scheme_t     *scheme_func,
                                cs_cell_sys_t             *csys,
                                cs_cell_builder_t         *cb)
{
  assert(eqp->space_scheme == CS_SPACE_SCHEME_MACFB);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ));

  /* Initialize the local matrix structure */

  cs_sdm_t *adv = cb->loc;

  if (cb->cell_flag & CS_FLAG_SOLID_CELL) {
    /* Nothing to do. No advection in the current cell volume */

    cs_sdm_init(cm->n_fc, macb->n_fc, adv);

    return;
  }

  /* Define the local operator for advection. Boundary conditions are also
     treated here since there are always weakly enforced */

  scheme_func(cm, macb, cb, adv, csys->rhs);

  /* Handle the specific case when there is no diffusion and no advection
   * flux. In this case, a zero row may appear leading to the divergence of
   * linear solver. To circumvent this issue, one set the boundary face value
   * to the cell value. This is equivalent to enforce a homogeneous Neumann
   * behavior. */

  assert(cs_equation_param_has_diffusion(eqp) == false);

  bft_error(__FILE__, __LINE__, 0, "%s: not implemented.", __func__);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, nullptr)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection fluxes");
    cs_log_printf(CS_LOG_DEFAULT, "\n beta_fluxes>>");
    for (int f = 0; f < cm->n_fc; f++)
      cs_log_printf(
        CS_LOG_DEFAULT, "f%d;% -5.3e|", cm->f_ids[f], cb->adv_fluxes[f]);
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, nullptr, nullptr, adv);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main function to build the cellwise advection operator for MAC
 *          face-based schemes.
 *          The local matrix related to this operator is stored in cb->loc
 *
 *          One assumes that a diffusion term is present so that there is no
 *          need to perform additional checkings on the well-posedness of the
 *          operator.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      scheme_func  pointer to the function building the system
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection(const cs_equation_param_t *eqp,
                   const cs_cell_mesh_t      *cm,
                   const cs_macfb_builder_t  *macb,
                   cs_macfb_adv_scheme_t     *scheme_func,
                   cs_cell_sys_t             *csys,
                   cs_cell_builder_t         *cb)
{
  assert(eqp->space_scheme == CS_SPACE_SCHEME_MACFB);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ));

  /* Initialize the local matrix structure */

  cs_sdm_t *adv = cb->loc;

  if (cb->cell_flag & CS_FLAG_SOLID_CELL) {
    /* Nothing to do. No advection in the current cell volume */

    cs_sdm_init(cm->n_fc, macb->n_fc, adv);

    return;
  }

  /* Remark: The flux across the dual faces is stored in cb->adv_fluxes and
     should have been computed previously in a function compliant with the
     cs_macfb_adv_open_hook_t prototype */

  /* Define the local operator for advection. Boundary conditions are also
     treated here since there are always weakly enforced */

  scheme_func(cm, macb, cb, adv, csys->rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, nullptr)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection fluxes");
    cs_log_printf(CS_LOG_DEFAULT, "\n beta_fluxes>>");
    for (int f = 0; f < macb->n_fc; f++)
      cs_log_printf(
        CS_LOG_DEFAULT, "f%d;% -5.3e|", macb->f_ids[f], cb->adv_fluxes[f]);
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, nullptr, nullptr, adv);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - non-conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - upwind scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_upwnoc(const cs_cell_mesh_t     *cm,
                          const cs_macfb_builder_t *macb,
                          cs_cell_builder_t        *cb,
                          cs_sdm_t                 *adv,
                          cs_real_t                *rhs)
{
  /* Compute the convection matrix */

  /* Sanity checks */
  assert(cm != nullptr && cb != nullptr && macb != nullptr);
  assert(adv != nullptr && rhs != nullptr);
  assert(cm->n_fc == 6);

  const cs_real_t *fluxes = cb->adv_fluxes;
  assert(fluxes != nullptr);

  /* Initialize objects */
  const cs_lnum_t n_cols = macb->n_dofs;
  cs_sdm_init(cm->n_fc, n_cols, adv);

  /* Loop on inner faces */
  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const short int fi_shift = fi * n_cols;
    const short int fo       = macb->f_opp_idx[fi];
    assert(fo >= 0 && fo < 6);

    /* Face info */
    const cs_real_t un_s_vol_cv = 1.0 / macb->f_vol_cv[fi];

    const cs_real_t flux_fi = fluxes[fi] * macb->f_sgn_axis[fo];
    const cs_real_t val_fi  = flux_fi * un_s_vol_cv;

    const short int fi_upw_idx = _get_face_idx_upw(fi, fo, flux_fi);
    /* entries */
    adv->val[fi_shift + fi]         = -val_fi;
    adv->val[fi_shift + fi_upw_idx] = val_fi;

    /* Loop on outer faces */
    for (short int fj = 0; fj < 4; fj++) {
      const short int shift_j = 4 * fi + fj;
      const short int fj_idx  = macb->f2f_idx[shift_j];

      const short int f0 = macb->f2fo_idx[2 * shift_j + 0];
      assert(f0 >= 0 && f0 < 6);

      const cs_real_t flux_fj = fluxes[6 + shift_j] * macb->f_sgn_axis[f0];
      const cs_real_t val_fj  = flux_fj * un_s_vol_cv;

      /* entry */
      adv->val[fi_shift + fi] -= val_fj;

      if (fj_idx >= 0) {
        const short int fj_upw_idx = _get_face_idx_upw(fi, fj_idx, flux_fj);

        /* entry */
        adv->val[fi_shift + fj_upw_idx] += val_fj;
      }
      else {
        /* To close convection reconstruction */
        rhs[fi] -= val_fj * macb->dir_values[shift_j];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - upwind scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_upwcsv(const cs_cell_mesh_t     *cm,
                          const cs_macfb_builder_t *macb,
                          cs_cell_builder_t        *cb,
                          cs_sdm_t                 *adv,
                          cs_real_t                *rhs)
{
  /* Compute the convection matrix */

  /* Sanity checks */
  assert(cm != nullptr && cb != nullptr && macb != nullptr);
  assert(adv != nullptr && rhs != nullptr);
  assert(cm->n_fc == 6);

  const cs_real_t *fluxes = cb->adv_fluxes;
  assert(fluxes != nullptr);

  /* Initialize objects */
  const cs_lnum_t n_cols = macb->n_dofs;
  cs_sdm_init(cm->n_fc, n_cols, adv);

  /* Loop on inner faces */
  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const short int fi_shift = fi * n_cols;
    const short int fo       = macb->f_opp_idx[fi];
    assert(fo >= 0 && fo < 6);

    /* Face info */
    const cs_real_t un_s_vol_cv = 1.0 / macb->f_vol_cv[fi];

    const cs_real_t flux_fi = fluxes[fi] * macb->f_sgn_axis[fo];
    const cs_real_t val_fi  = flux_fi * un_s_vol_cv;

    const short int fi_upw_idx = _get_face_idx_upw(fi, fo, flux_fi);
    /* entry */
    adv->val[fi_shift + fi_upw_idx] = val_fi;

    /* Loop on outer faces */
    for (short int fj = 0; fj < 4; fj++) {
      const short int shift_j = 4 * fi + fj;
      const short int fj_idx  = macb->f2f_idx[shift_j];

      const short int f0 = macb->f2fo_idx[2 * shift_j + 0];
      assert(f0 >= 0 && f0 < 6);

      const cs_real_t flux_fj = fluxes[6 + shift_j] * macb->f_sgn_axis[f0];
      const cs_real_t val_fj  = flux_fj * un_s_vol_cv;

      if (fj_idx >= 0) {
        const short int fj_upw_idx = _get_face_idx_upw(fi, fj_idx, flux_fj);

        /* extra-diagonal entry */
        adv->val[fi_shift + fj_upw_idx] += val_fj;
      }
      else {
        /* To close convection reconstruction */
        rhs[fi] -= val_fj * macb->dir_values[shift_j];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - non-conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - centered scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_cennoc(const cs_cell_mesh_t     *cm,
                          const cs_macfb_builder_t *macb,
                          cs_cell_builder_t        *cb,
                          cs_sdm_t                 *adv,
                          cs_real_t                *rhs)
{
  /* Compute the convection matrix */

  /* Sanity checks */
  assert(cm != nullptr && cb != nullptr && macb != nullptr);
  assert(adv != nullptr && rhs != nullptr);
  assert(cm->n_fc == 6);

  const cs_real_t *fluxes = cb->adv_fluxes;
  assert(fluxes != nullptr);

  /* Initialize objects */
  const cs_lnum_t n_cols = macb->n_dofs;
  cs_sdm_init(cm->n_fc, n_cols, adv);

  /* Loop on inner faces */
  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const short int fi_shift = fi * n_cols;
    const short int fo       = macb->f_opp_idx[fi];
    assert(fo >= 0 && fo < 6);

    /* Face info */
    const cs_real_t un_s2_vol_cv = 1.0 / (2.0 * macb->f_vol_cv[fi]);

    const cs_real_t flux_fi = fluxes[fi] * macb->f_sgn_axis[fo];
    const cs_real_t val_fi  = flux_fi * un_s2_vol_cv;

    /* diagonal entry */
    adv->val[fi_shift + fi] = -val_fi;
    /* extra-diagonal entry */
    adv->val[fi_shift + fo] = val_fi;

    /* Loop on outer faces */
    for (short int fj = 0; fj < 4; fj++) {
      const short int shift_j = 4 * fi + fj;
      const short int fj_idx  = macb->f2f_idx[shift_j];

      const short int f0 = macb->f2fo_idx[2 * shift_j + 0];
      assert(f0 >= 0 && f0 < 6);

      const cs_real_t flux_fj = fluxes[6 + shift_j] * macb->f_sgn_axis[f0];
      const cs_real_t val_fj  = flux_fj * un_s2_vol_cv;

      /* diagonal entry */
      adv->val[fi_shift + fi] -= val_fj;

      if (fj_idx >= 0) {
        /* extra-diagonal entry */
        adv->val[fi_shift + fj_idx] = val_fj;
      }
      else {
        /* To close convection reconstruction */
        rhs[fi] -= val_fj * macb->dir_values[shift_j];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a MAC
 *         face-based scheme
 *         - conservative formulation \f$ \vec{beta} \cdot \nabla u \f$
 *         - centered scheme
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      macb    pointer to a cs_macfb_builder_t structure
 * \param[in]      cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] adv     pointer to a local matrix to build
 * \param[in, out] rhs     pointer to a cs_real_t array. It is filled inside the
 *                        function. Have to preallocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_advection_cencsv(const cs_cell_mesh_t     *cm,
                          const cs_macfb_builder_t *macb,
                          cs_cell_builder_t        *cb,
                          cs_sdm_t                 *adv,
                          cs_real_t                *rhs)
{
  /* Compute the convection matrix */

  /* Sanity checks */
  assert(cm != nullptr && cb != nullptr && macb != nullptr);
  assert(adv != nullptr && rhs != nullptr);
  assert(cm->n_fc == 6);

  const cs_real_t *fluxes = cb->adv_fluxes;
  assert(fluxes != nullptr);

  /* Initialize objects */
  const cs_lnum_t n_cols = macb->n_dofs;
  cs_sdm_init(cm->n_fc, n_cols, adv);

  /* Loop on inner faces */
  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const short int fi_shift = fi * n_cols;
    const short int fo       = macb->f_opp_idx[fi];
    assert(fo >= 0 && fo < 6);

    /* Face info */
    const cs_real_t un_s2_vol_cv = 1.0 / (2.0 * macb->f_vol_cv[fi]);

    const cs_real_t flux_fi = fluxes[fi] * macb->f_sgn_axis[fo];
    const cs_real_t val_fi  = flux_fi * un_s2_vol_cv;

    /* diagonal entry */
    adv->val[fi_shift + fi] = val_fi;
    /* extra-diagonal entry */
    adv->val[fi_shift + fo] = val_fi;

    /* Loop on outer faces */
    for (short int fj = 0; fj < 4; fj++) {
      const short int shift_j = 4 * fi + fj;
      const short int fj_idx  = macb->f2f_idx[shift_j];

      const short int f0 = macb->f2fo_idx[2 * shift_j + 0];
      assert(f0 >= 0 && f0 < 6);

      const cs_real_t flux_fj = fluxes[6 + shift_j] * macb->f_sgn_axis[f0];
      const cs_real_t val_fj  = flux_fj * un_s2_vol_cv;

      /* diagonal entry */
      adv->val[fi_shift + fi] += val_fj;

      if (fj_idx >= 0) {
        /* extra-diagonal entry */
        adv->val[fi_shift + fj_idx] = val_fj;
      }
      else {
        /* To close convection reconstruction */
        rhs[fi] -= val_fj * macb->dir_values[shift_j];
      }
    }
  }
}

END_C_DECLS
