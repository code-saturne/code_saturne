/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
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
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_solve.h"
#include "cs_cdo_toolbox.h"
#include "cs_equation_bc.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_MACFB_VECTEQ_DBG 0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */

static cs_cell_sys_t     **cs_macfb_cell_sys = nullptr;
static cs_cell_builder_t **cs_macfb_cell_bld = nullptr;

/* Pointer to shared structures */

static const cs_cdo_quantities_t *cs_shared_quant;
static const cs_cdo_connect_t    *cs_shared_connect;
static const cs_time_step_t      *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the local builder structure used for building the system
 *         cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_cell_builder_t *
_cell_builder_create(const cs_cdo_connect_t *connect)
{
  const int n_fc = connect->n_max_fbyc;
  assert(n_fc == 6);
  const int n_dofs = 30;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  /* Since it relies on the scalar case, n_fc should be enough */

  BFT_MALLOC(cb->adv_fluxes, n_dofs, double);
  memset(cb->adv_fluxes, 0, n_dofs * sizeof(double));

  BFT_MALLOC(cb->ids, n_dofs, int);
  memset(cb->ids, 0, n_dofs * sizeof(int));

  int size = n_fc * n_dofs;
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size * sizeof(cs_real_t));

  size = 2 * n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size * sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */

  cb->aux = cs_sdm_square_create(n_dofs);
  cb->loc = cs_sdm_square_create(n_dofs);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system.
 *          Case of vector-valued MAC-Fb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] diff_hodge  pointer to a Hodge op. for diffusion and its pty
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vfb_apply_bc(const cs_equation_param_t   *eqp,
              const cs_macfb_vecteq_t     *eqc,
              const cs_equation_builder_t *eqb,
              const cs_cell_mesh_t        *cm,
              cs_face_mesh_t              *fm,
              cs_hodge_t                  *diff_hodge,
              cs_cell_sys_t               *csys,
              cs_cell_builder_t           *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann)
      for (short int i = 0; i < cm->n_fc; i++)
        csys->rhs[i] -= csys->neu_values[i];

    /* Weakly enforced Dirichlet BCs for cells attached to the boundary
       csys is updated inside (matrix and rhs) */

    if (cs_equation_param_has_diffusion(eqp)) {

      if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
          || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
        eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);
    }

    if (csys->has_sliding)
      eqc->enforce_sliding(eqp, cm, fm, diff_hodge, cb, csys);

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED
        || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

      /* Enforced Dirichlet BCs for cells attached to the boundary
       * csys is updated inside (matrix and rhs). This is close to a strong
       * way to enforce Dirichlet BCs */

      eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);
    }

  } /* Boundary cell */

  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    cs_equation_builder_enforce_block_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(
        "\n>> MAC-fb vecteq: Cell system after the internal enforcement", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the linear system arising from a vector steady-state
 *         diffusion equation with a MAC-Fb scheme
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_solve_system(const bool                 cur2prev,
              const cs_equation_param_t *eqp,
              cs_equation_builder_t     *eqb,
              cs_macfb_vecteq_t         *eqc,
              cs_field_t                *fld)
{
  cs_cdo_system_helper_t *sh = eqb->system_helper;

  const cs_lnum_t n_faces = cs_shared_quant->n_faces;

  cs_timer_t ts = cs_timer_time();

  if (cur2prev && eqc->face_values_pre != nullptr)
    cs_array_real_copy(n_faces, eqc->face_values, eqc->face_values_pre);

  /* Solve the linear system (treated as a scalar-valued system) */

  cs_real_t       normalization = 1.0; /* TODO */
  cs_sles_t *sles = cs_sles_find_or_add(eqp->sles_param->field_id, nullptr);
  cs_matrix_t    *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t *range_set = cs_cdo_system_get_range_set(sh, 0);

  cs_cdo_solve_scalar_system(n_faces,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             normalization,
                             true, /* rhs_redux */
                             sles,
                             eqc->face_values,
                             sh->rhs);

  cs_timer_t te = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &ts, &te);

  /* Update fields */

  cs_macfb_vecteq_update_fields(&(eqb->tce), fld, cur2prev);

  /* Free remaining buffers */

  cs_sles_free(sles);
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed. This is stored inside eqb
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_setup(cs_real_t                  t_eval,
                      const cs_mesh_t           *mesh,
                      const cs_equation_param_t *eqp,
                      cs_equation_builder_t     *eqb)
{
  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_cdo_connect_t    *connect = cs_shared_connect;

  /* Initialize and compute the values of the Dirichlet BC */

  BFT_MALLOC(eqb->dir_values, 3 * quant->n_b_faces, cs_real_t);
  cs_array_real_fill_zero(3 * quant->n_b_faces, eqb->dir_values);

  cs_equation_bc_dirichlet_at_faces(
    mesh, quant, connect, eqp, eqb->face_bc, t_eval, eqb->dir_values);

  /* Internal enforcement of DoFs  */

  if (cs_equation_param_has_internal_enforcement(eqp))
    eqb->enforced_values = cs_enforcement_define_at_faces(
      connect, eqp->n_enforcements, eqp->enforcement_params);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *          The algebraic system for time t^{n+1} is going to be built knowing
 *          previous field at time t^{n} and potentially the field at time
 *          t^{n-1}. Make sure to be consistent between the call to
 *          current_to_previous and the parameters vel_{f}_n/nm1 given
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      val_f_n     face DoFs at time step n
 * \param[in]      val_f_nm1   face DoFs at time step n-1 or nullptr
 * \param[in, out] macb        pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_init_cell_system(const cs_cell_mesh_t        *cm,
                                 const cs_equation_param_t   *eqp,
                                 const cs_equation_builder_t *eqb,
                                 const cs_real_t              val_f_n[],
                                 const cs_real_t              val_f_nm1[],
                                 cs_macfb_builder_t          *macb,
                                 cs_cell_sys_t               *csys,
                                 cs_cell_builder_t           *cb)
{
  /* Cell-wise view of the linear system to build */

  /* Sanity check*/
  if (cm->n_fc != 6) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _(" %s: The cell %d has not 6 faces for MAC-fb scheme\n"),
              __func__,
              cm->c_id);
  }

  /* Only face dofs */
  const int n_dofs = macb->n_dofs;

  csys->c_id   = cm->c_id;
  csys->n_dofs = n_dofs;

  assert(macb->n_fc == n_dofs);

  /* Initialize the local system */

  cs_cell_sys_reset(n_dofs, csys);

  cs_sdm_init(cm->n_fc, n_dofs, csys->mat);

  /* One has to keep the same numbering for faces between cell mesh and cell
     system */

  for (int f = 0; f < n_dofs; f++) {

    const cs_lnum_t f_id = macb->f_ids[f];
    csys->dof_ids[f]     = macb->dof_ids[f];
    if (val_f_n != nullptr) { /* Case of steady algo. */
      csys->val_n[f] = val_f_n[f_id];
    }
  }

  if (val_f_nm1
      != nullptr) { /* State at n-1 is given (2nd order time scheme) */

    for (int f = 0; f < n_dofs; f++) {

      const cs_lnum_t f_id = macb->f_ids[f];
      csys->val_nm1[f]     = val_f_nm1[f_id];
    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    cs_equation_bc_set_cw_macfb(
      cm, eqp, eqb->face_bc, eqb->dir_values, macb, csys);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif

  } /* Border cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys))
    cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize stuctures for a gven cell
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      c_id         cell id
 * \param[in]      vel_f_n      velocity face DoFs of the previous time step
 * \param[in, out] cm           pointer to a cellwise view of the mesh
 * \param[in, out] macb         pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_init_build(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_equation_param_t   *eqp,
                           const cs_equation_builder_t *eqb,
                           const cs_lnum_t              c_id,
                           const cs_real_t              vel_f_n[],
                           cs_cell_mesh_t              *cm,
                           cs_macfb_builder_t          *macb,
                           cs_cell_sys_t               *csys,
                           cs_cell_builder_t           *cb)
{

  /* Set the current cell flag */

  cb->cell_flag = connect->cell_flag[c_id];

  /* Set the local mesh structure for the current cell */

  cs_cell_mesh_build(c_id,
                     cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                     connect,
                     quant,
                     cm);

  /* Set the local builder structure for the current cell */

  cs_macfb_builder_cellwise_setup(cm, connect, quant, macb);

  /* For the problem, the global system writes:
   *
   *     |        |         |
   *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
   *     |        |         |  A is csys->mat in what follows
   *     |--------|---------|  The viscous part arising from the MAC-Fb
   *     |        |         |  schemes for vector-valued variables and
   *     |   B    |    0    |  additional terms as the linearized
   *     |        |         |  convective term
   *
   * Set the local (i.e. cellwise) structures for the current cell
   */

  cs_macfb_vecteq_init_cell_system(cm,
                                   eqp,
                                   eqb,
                                   vel_f_n,
                                   nullptr, /* no n-1 state is given */
                                   macb,
                                   csys,
                                   cb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term for a vector-valued MAC scheme
 *         and add it to the local rhs
 *
 * \param[in]      cm          pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      t_eval      time at which the source term is evaluated
 * \param[in]      coef        scaling of the time source (for theta schemes)
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] cb          pointer to a \ref cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_sourceterm(const cs_cell_mesh_t      *cm,
                           const cs_equation_param_t *eqp,
                           cs_macfb_builder_t        *macb,
                           const cs_real_t            t_eval,
                           const cs_real_t            coef,
                           cs_equation_builder_t     *eqb,
                           cs_cell_builder_t         *cb,
                           cs_cell_sys_t             *csys)
{

  /* Reset the local contribution */

  memset(csys->source, 0, csys->n_dofs * sizeof(cs_real_t));

  cs_source_term_compute_cellwise(eqp->n_source_terms,
                                  (cs_xdef_t *const *)eqp->source_terms,
                                  cm,
                                  eqb->source_mask,
                                  eqb->compute_source,
                                  t_eval,
                                  (void *)macb,
                                  cb,
                                  csys->source);

  for (short int f = 0; f < cm->n_fc; f++) {
    csys->rhs[f] += coef * csys->source[f];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion term in the
 *          vector-valued CDO-Fb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      diff_pty    pointer to a cs_property_t structure
 *                             for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_diffusion(const cs_equation_param_t *eqp,
                          const cs_cell_mesh_t      *cm,
                          const cs_macfb_builder_t  *macb,
                          const cs_property_t       *diff_pty,
                          cs_cell_sys_t             *csys,
                          cs_cell_builder_t         *cb)
{
  if (cs_equation_param_has_diffusion(eqp)) {

    /* Compute the diffusion matrix */

    cs_macfb_diffusion(cm, macb, diff_pty, cb->loc, csys->rhs);

    /* Add the local diffusion operator to the local system */

    cs_sdm_add_block_topleft(csys->mat, cm->n_fc, macb->n_dofs, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> MAC-fb vecteq: Cell system after diffusion", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the convection term in the
 *          vector-valued MAC-Fb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_advection(const cs_equation_param_t *eqp,
                          const cs_macfb_vecteq_t   *eqc,
                          const cs_cell_mesh_t      *cm,
                          const cs_macfb_builder_t  *macb,
                          cs_cell_sys_t             *csys,
                          cs_cell_builder_t         *cb)
{
  if (cs_equation_param_has_convection(eqp)
      && ((cb->cell_flag & CS_FLAG_SOLID_CELL) == 0)) {

    /* Open hook: Compute the advection flux for the numerical scheme and store
   the advection fluxes across primal faces */

    eqc->advection_open(eqp, cm, macb, csys, eqc->advection_input, cb);

    /* Compute the local advection matrix */

    eqc->advection_main(eqp, cm, macb, eqc->advection_scheme, csys, cb);

    /* Close hook: Modify if needed the computed advection matrix and update
       the local system */

    eqc->advection_close(cm, macb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> MAC-fb vecteq: Cell system after convection",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the convection, diffusion,
 *          reaction terms in vector-valued MAC-Fb schemes.
 *          mass_hodge could be set to nullptr if a Voronoi algo. is used.
 *          Otherwise, the mass matrix should be pre-computed.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      diff_pty    pointer to a cs_property_data_t structure
 *                              for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_conv_diff_reac(const cs_equation_param_t   *eqp,
                               const cs_equation_builder_t *eqb,
                               const cs_macfb_vecteq_t     *eqc,
                               const cs_cell_mesh_t        *cm,
                               const cs_macfb_builder_t    *macb,
                               const cs_property_data_t    *diff_pty,
                               cs_cell_sys_t               *csys,
                               cs_cell_builder_t           *cb)
{

  /* Diffusion term */

  cs_macfb_vecteq_diffusion(eqp, cm, macb, diff_pty->property, csys, cb);

  /* Convection term */

  cs_macfb_vecteq_advection(eqp, eqc, cm, macb, csys, cb);

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    /* Update the value of the reaction property(ies) if needed */

    cs_equation_builder_set_reaction_pty_cw(eqp, eqb, cm, cb);

    bft_error(__FILE__, __LINE__, 0, "Reaction is not implemented.");

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> MAC-fb vecteq: Cell system after reaction", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the matrix and rhs for a vector-valued MAC scheme
 *          and Euler implicit. Values are added in place
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      cb          pointer to a \ref cs_cell_builder_t structure
 * \param[in]      dt          value of the time step
 * \param[in, out] csys        pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_euler_implicit_term(const cs_equation_param_t *eqp,
                                    const cs_cell_mesh_t      *cm,
                                    const cs_macfb_builder_t  *macb,
                                    const cs_cell_builder_t   *cb,
                                    const cs_real_t            dt,
                                    cs_cell_sys_t             *csys)
{
  assert(csys != nullptr);

  cs_sdm_t *mat = csys->mat;

  assert(mat->n_rows >= cm->n_fc);

  const cs_lnum_t n_cols = mat->n_cols;

  const cs_real_t inv_dt = 1. / dt;

  /* Loop on inner faces */
  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const cs_real_t rho_f = cs_property_get_face_value(
      cm->f_ids[fi], cb->t_pty_eval, eqp->time_property);

    cs_real_t val_fi = rho_f * macb->f_vol_cv[fi] * inv_dt;

    /* if not a boundary face: divided by 2 since
     * two cells share this face */

    if (csys->bf_ids[fi] <= -1) {
      val_fi *= 0.5;
    }

    /* diagonal entry */
    mat->val[fi * n_cols + fi] += val_fi;

    /* rhs */
    csys->rhs[fi] += val_fi * csys->val_n[fi];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with MAC-fb scheme
 *
 * \param[in]      csys         pointer to a cs_cell_sys_t structure
 * \param[in, out] block        pointer to a block structure
 * \param[in, out] rhs          array of values for the rhs
 * \param[in, out] asb          pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_assembly(const cs_cell_sys_t   *csys,
                         cs_cdo_system_block_t *block,
                         cs_real_t             *rhs,
                         cs_cdo_assembly_t     *asb)
{
  assert(asb != nullptr && block != nullptr); /* Sanity check */
  assert(block->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  cs_cdo_system_dblock_t *db
    = static_cast<cs_cdo_system_dblock_t *>(block->block_pointer);

  /* Matrix assembly */

  db->assembly_func(csys->mat, csys->dof_ids, db->range_set, asb, db->mav);

  /* RHS assembly only on faces */

#pragma omp critical
  {
    for (short int f = 0; f < csys->n_dofs; f++)
      rhs[csys->dof_ids[f]] += csys->rhs[f];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables associated to cells in case of a MAC-Fb
 *         scheme. This has to be done after a resolution.
 *
 * \param[in, out] tce       pointer to a timer counter
 * \param[in, out] fld       pointer to a cs_field_t structure
 * \param[in]      cur2prev  true if one performs "current to previous" op.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_update_fields(cs_timer_counter_t *tce,
                              cs_field_t         *fld,
                              bool                cur2prev)
{
  cs_timer_t t0 = cs_timer_time();

  /* Copy current field values to previous values */

  if (cur2prev)
    cs_field_current_to_previous(fld);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(tce, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a MAC-Fb scheme:
 *           - steady scheme
 *           - implicit Euler scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_solve_steady_implicit(bool                       cur2prev,
                                      const cs_mesh_t           *mesh,
                                      const int                  field_id,
                                      const cs_equation_param_t *eqp,
                                      cs_equation_builder_t     *eqb,
                                      void                      *context)
{
  cs_timer_t t0 = cs_timer_time();

  const cs_cdo_connect_t    *connect = cs_shared_connect;
  const cs_cdo_quantities_t *quant   = cs_shared_quant;

  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)context;

  const cs_time_step_t *ts        = cs_shared_time_step;
  const cs_real_t       time_eval = ts->t_cur + ts->dt[0];

  cs_field_t             *fld = cs_field_by_id(field_id);
  cs_cdo_system_helper_t *sh  = eqb->system_helper;

  /* Build an array storing the Dirichlet values at faces.
   * First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */

  cs_macfb_vecteq_setup(time_eval, mesh, eqp, eqb);

  /* Initialize the local system: rhs, matrix and assembler values */

  cs_real_t *rhs = nullptr;

  cs_cdo_system_helper_init_system(sh, &rhs);

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    const int t_id = cs_get_thread_id();

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t     *fm   = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t     *cm   = cs_cdo_local_get_cell_mesh(t_id);
    cs_macfb_builder_t *macb = cs_macfb_get_builder(t_id);
    cs_cell_sys_t      *csys = cs_macfb_cell_sys[t_id];
    cs_cell_builder_t  *cb   = cs_macfb_cell_bld[t_id];
    cs_cdo_assembly_t  *asb  = cs_cdo_assembly_get(t_id);
    cs_hodge_t         *diff_hodge = (eqc->diffusion_hodge == nullptr)
                                       ? nullptr
                                       : eqc->diffusion_hodge[t_id];

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = time_eval;
    cb->t_bc_eval  = time_eval;
    cb->t_st_eval  = time_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* 1- Initialize structures */

      cs_macfb_vecteq_init_build(
        connect, quant, eqp, eqb, c_id, eqc->face_values, cm, macb, csys, cb);

      /**/
      cs_macfb_vecteq_conv_diff_reac(
        eqp, eqb, eqc, cm, macb, diff_hodge->pty_data, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) /* SOURCE TERM */
        cs_macfb_vecteq_sourceterm(
          cm, eqp, macb, time_eval, 1.0, eqb, cb, csys);

      if (cs_equation_param_has_time(eqp)) {
        assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);
        cs_macfb_vecteq_euler_implicit_term(eqp, cm, macb, cb, ts->dt[0], csys);
      }

      /* Apply BOUNDARY CONDITIONS */

      _vfb_apply_bc(eqp, eqc, eqb, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> MAC-fb vecteq: (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_macfb_vecteq_assembly(csys, sh->blocks[0], sh->rhs, asb);

    } /* Main loop on cells */

  } /* OpenMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system (treated as a scalar-valued system) */

  _solve_system(cur2prev, eqp, eqb, eqc, fld);

  cs_cdo_system_helper_reset(sh); /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a MAC-Fb scheme and an implicit/explicit theta scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_solve_theta(bool                       cur2prev,
                            const cs_mesh_t           *mesh,
                            const int                  field_id,
                            const cs_equation_param_t *eqp,
                            cs_equation_builder_t     *eqb,
                            void                      *context)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);
  CS_UNUSED(cur2prev);
  CS_UNUSED(mesh);
  CS_UNUSED(field_id);
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(context);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a MAC-Fb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_macfb_vecteq_is_initialized(void)
{
  if (cs_macfb_cell_sys == nullptr || cs_macfb_cell_bld == nullptr)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to MAC
 *         vector-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_init_sharing(const cs_cdo_quantities_t *quant,
                             const cs_cdo_connect_t    *connect,
                             const cs_time_step_t      *time_step)
{
  /* Assign static const pointers */

  cs_shared_quant     = quant;
  cs_shared_connect   = connect;
  cs_shared_time_step = time_step;

  /* Specific treatment for handling openMP */

  BFT_MALLOC(cs_macfb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_macfb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_macfb_cell_sys[i] = nullptr;
    cs_macfb_cell_bld[i] = nullptr;
  }

  const int n_max_dofs = 30;
  const int n_max_fbyc = 30;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t *cb   = _cell_builder_create(connect);
    cs_macfb_cell_bld[t_id] = cb;

    int *block_size = nullptr;
    cs_macfb_cell_sys[t_id]
      = cs_cell_sys_create(n_max_dofs, n_max_fbyc, 1, block_size);
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t *cb = _cell_builder_create(connect);
  cs_macfb_cell_bld[0]  = cb;

  int *block_size = nullptr;
  cs_macfb_cell_sys[0]
    = cs_cell_sys_create(n_max_dofs, n_max_fbyc, 1, block_size);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   double pointer to a \ref cs_cell_sys_t structure
 * \param[out]  cb     double pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_get(cs_cell_sys_t **csys, cs_cell_builder_t **cb)
{
  const int t_id = cs_get_thread_id();

  *csys = cs_macfb_cell_sys[t_id];
  *cb   = cs_macfb_cell_bld[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO face-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_finalize_sharing(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_macfb_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_macfb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_macfb_cell_sys[0]));
  cs_cell_builder_free(&(cs_macfb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_macfb_cell_sys);
  BFT_FREE(cs_macfb_cell_bld);
  cs_macfb_cell_bld = nullptr;
  cs_macfb_cell_sys = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_macfb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_macfb_vecteq_init_context(cs_equation_param_t   *eqp,
                             int                    var_id,
                             int                    bflux_id,
                             cs_equation_builder_t *eqb)
{
  assert(eqp != nullptr && eqb != nullptr);
  if (eqp->space_scheme != CS_SPACE_SCHEME_MACFB || eqp->dim != 3)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid type of equation.\n"
              " Expected: vector-valued MAC face-based equation.",
              __func__);

  const cs_cdo_connect_t *connect = cs_shared_connect;
  const cs_lnum_t         n_faces = connect->n_faces[CS_ALL_FACES];

  cs_macfb_vecteq_t *eqc = nullptr;

  BFT_MALLOC(eqc, 1, cs_macfb_vecteq_t);

  eqc->var_field_id   = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system - face dofs only */

  eqc->n_faces = n_faces;
  eqc->n_dofs  = n_faces;

  eqb->sys_flag = CS_FLAG_SYS_VECTOR;
  eqb->msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */

  eqb->bdy_flag
    = CS_FLAG_COMP_PV | CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ;

  BFT_MALLOC(eqc->face_values, n_faces, cs_real_t);
  cs_array_real_fill_zero(n_faces, eqc->face_values);
  BFT_MALLOC(eqc->face_values_pre, n_faces, cs_real_t);
  cs_array_real_fill_zero(n_faces, eqc->face_values_pre);

  bool need_eigen
    = (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
       || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
        ? true
        : false;

  /* Diffusion term */

  eqc->diffusion_hodge  = nullptr;
  eqc->enforce_robin_bc = nullptr;

  if (cs_equation_param_has_diffusion(eqp)) {

    eqc->diffusion_hodge = cs_hodge_init_context(connect,
                                                 eqp->diffusion_property,
                                                 &(eqp->diffusion_hodgep),
                                                 true,        /* tensor ? */
                                                 need_eigen); /* eigen ? */

    /* If necessary, enrich the mesh flag to account for a property defined
     * by an analytical expression. In this case, one evaluates the definition
     * as the mean value over the cell */

    const cs_xdef_t *diff_def = eqp->diffusion_property->defs[0];
    if (diff_def->type == CS_XDEF_BY_ANALYTIC_FUNCTION)
      eqb->msh_flag
        |= cs_quadrature_get_flag(diff_def->qtype, cs_flag_primal_cell);

  } /* Diffusion part */

  eqc->enforce_dirichlet = nullptr;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
  case CS_PARAM_BC_ENFORCE_PENALIZED:
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    /* Nothing for the moment */
    break;

  default:
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);
  }

  eqc->enforce_sliding = nullptr;
  if (eqb->face_bc->n_sliding_faces > 0) {

    /* There is at least one face with a sliding condition to handle */

    eqb->bdy_flag |= CS_FLAG_COMP_HFQ;
    eqc->enforce_sliding = cs_cdo_diffusion_vfb_wsym_sliding;
  }

  /* Advection part */

  cs_macfb_set_advection_function(eqp, eqb, (cs_macfb_priv_t *)eqc);

  /* Reaction term */

  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodgep.algo != CS_HODGE_ALGO_VORONOI)
      bft_error(__FILE__,
                __LINE__,
                0,
                "%s: Eq. %s: Invalid type of discretization for the reaction"
                " term\n",
                __func__,
                eqp->name);

    /* If necessary, enrich the mesh flag to account for a property defined
     * by an analytical expression. In this case, one evaluates the definition
     * as the mean value over the cell */

    for (short int ir = 0; ir < eqp->n_reaction_terms; ir++) {
      const cs_xdef_t *rea_def = eqp->reaction_properties[ir]->defs[0];
      if (rea_def->type == CS_XDEF_BY_ANALYTIC_FUNCTION)
        eqb->msh_flag
          |= cs_quadrature_get_flag(rea_def->qtype, cs_flag_primal_cell);
    }
  }

  /* Unsteady term */

  if (cs_equation_param_has_time(eqp)) {

    if (eqp->time_hodgep.algo == CS_HODGE_ALGO_VORONOI) {
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    }
    else {
      bft_error(__FILE__,
                __LINE__,
                0,
                "%s: Eq. %s: Invalid type of discretization for time"
                " term\n",
                __func__,
                eqp->name);
    }
  }

  /* Source term part */

  eqc->source_terms = nullptr;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_faces, cs_real_t);
    cs_array_real_fill_zero(n_faces, eqc->source_terms);

  } /* There is at least one source term */

  /* Mass matrix */

  eqc->mass_hodgep.inv_pty = false;
  eqc->mass_hodgep.type    = CS_HODGE_TYPE_FB;
  eqc->mass_hodgep.algo    = CS_HODGE_ALGO_VORONOI;
  eqc->mass_hodgep.coef    = 0.0;

  eqc->mass_hodge = nullptr;

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) {

    eqc->mass_hodge = cs_hodge_init_context(connect,
                                            nullptr,
                                            &(eqc->mass_hodgep),
                                            false,  /* tensor ? */
                                            false); /* eigen ? */

    if (eqp->verbosity > 1) {
      cs_log_printf(CS_LOG_SETUP,
                    "#### Parameters of the mass matrix of the equation %s\n",
                    eqp->name);
      cs_hodge_param_log("Mass matrix", nullptr, eqc->mass_hodgep);
    }
  }

  /* Helper structures (range set, interface set, matrix structure and all the
     assembly process) */

  cs_cdo_system_helper_t *sh             = nullptr;
  cs_lnum_t               col_block_size = n_faces;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_DEFAULT,
                                   1,               /* n_col_blocks */
                                   &col_block_size, /* col_block_size */
                                   1);              /* n_blocks */

  /* Choose the right class of matrix to avoid copy.
   * The way to perform the assembly may change if an external librairy is used
   * for solving the linear system */

  cs_cdo_system_matrix_class_t matclass;

  switch (eqp->sles_param->solver_class) {

  case CS_PARAM_SOLVER_CLASS_CS:
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  case CS_PARAM_SOLVER_CLASS_HYPRE:
#if defined(HAVE_HYPRE)
    matclass = CS_CDO_SYSTEM_MATRIX_HYPRE;
#else
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
#endif
    break;

  default:
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;
  }

  cs_cdo_system_add_dblock(sh,
                           0, /* block_id */
                           matclass,
                           cs_flag_mac_primal_face,
                           n_faces,
                           1,     /* stride */
                           true,  /* interlaced */
                           true); /* unrolled */

  cs_cdo_system_build_block(sh, 0);

  eqb->system_helper = sh;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_macfb_vecteq_t structure
 *
 * \param[in, out]  data   pointer to a cs_macfb_vecteq_t structure
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_macfb_vecteq_free_context(void *data)
{
  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)data;

  if (eqc == nullptr)
    return eqc;

  /* Free temporary buffers */

  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->face_values);
  BFT_FREE(eqc->face_values_pre);

  cs_hodge_free_context(&(eqc->diffusion_hodge));
  cs_hodge_free_context(&(eqc->mass_hodge));

  BFT_FREE(eqc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of vector-valued Mac-Fb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_init_values(cs_real_t                  t_eval,
                            const int                  field_id,
                            const cs_mesh_t           *mesh,
                            const cs_equation_param_t *eqp,
                            cs_equation_builder_t     *eqb,
                            void                      *context)
{
  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_cdo_connect_t    *connect = cs_shared_connect;

  cs_macfb_vecteq_t *eqc    = (cs_macfb_vecteq_t *)context;
  cs_field_t        *fld    = cs_field_by_id(field_id);
  cs_real_t         *c_vals = fld->val;
  cs_real_t         *f_vals = eqc->face_values;

  /* Check that a face interface has been defined */

  if (eqp->n_ic_defs > 0 && cs_glob_n_ranks > 1 && connect->face_ifs == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Interface set structure at faces not allocated.\n",
              __func__);

  /* By default, 0 is set as initial condition for the computational domain */

  /* 1 value per face and 3 values per cell (an interpolation of face
   * velocity)*/
  cs_array_real_fill_zero(quant->n_faces, f_vals);
  cs_array_real_fill_zero(3 * quant->n_cells, c_vals);

  if (eqp->n_ic_defs > 0) {

    cs_lnum_t *def2f_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
    cs_lnum_t *def2f_idx = nullptr;

    BFT_MALLOC(def2f_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_cdo_sync_vol_def_at_faces(
      eqp->n_ic_defs, eqp->ic_defs, def2f_idx, def2f_ids);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */

      const cs_xdef_t *def          = eqp->ic_defs[def_id];
      const cs_lnum_t  n_f_selected = def2f_idx[def_id + 1] - def2f_idx[def_id];
      const cs_lnum_t *selected_lst = def2f_ids + def2f_idx[def_id];

      switch (def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_at_faces_by_value(
          def, n_f_selected, selected_lst, f_vals);
        cs_evaluate_potential_at_cells_by_value(def, c_vals);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION: {
        const cs_param_dof_reduction_t red = eqp->dof_reduction;
        switch (red) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_evaluate_potential_at_faces_by_analytic(
            def, t_eval, n_f_selected, selected_lst, f_vals);
          cs_evaluate_potential_at_cells_by_analytic(def, t_eval, c_vals);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_evaluate_average_on_faces_by_analytic(
            def, t_eval, n_f_selected, selected_lst, f_vals);
          cs_evaluate_average_on_cells_by_analytic(def, t_eval, c_vals);
          break;

        default:
          bft_error(__FILE__,
                    __LINE__,
                    0,
                    " %s: Incompatible reduction for equation %s.\n",
                    __func__,
                    eqp->name);
          break;

        } /* Switch on possible reduction types */

      } break;

      case CS_XDEF_BY_QOV: /* TODO */
      default:
        bft_error(__FILE__,
                  __LINE__,
                  0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__,
                  eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

    /* Free */

    BFT_FREE(def2f_idx);

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */

  /* Allocate field for Dirichlet BC*/
  cs_real_t *f_diri_vals = nullptr;
  BFT_MALLOC(f_diri_vals, 3 * quant->n_b_faces, cs_real_t);

  cs_equation_bc_dirichlet_at_faces(
    mesh, quant, connect, eqp, eqb->face_bc, t_eval, f_diri_vals);

  for (int bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {

    const cs_lnum_t f_id = quant->n_i_faces + bf_id;
    f_vals[f_id]         = f_diri_vals[3 * bf_id + quant->face_axis[f_id]];
  }

  BFT_FREE(f_diri_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_current_to_previous(const cs_equation_param_t *eqp,
                                    cs_equation_builder_t     *eqb,
                                    void                      *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)context;
  cs_field_t        *fld = cs_field_by_id(eqc->var_field_id);

  /* Face values */

  if (eqc->face_values_pre != nullptr)
    cs_array_real_copy(eqc->n_faces, eqc->face_values, eqc->face_values_pre);

  /* Cell values */

  cs_field_current_to_previous(fld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_extra_post(const cs_equation_param_t *eqp,
                           cs_equation_builder_t     *eqb,
                           void                      *context)
{
  CS_UNUSED(eqp);

  cs_timer_t t0 = cs_timer_time();

  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)context;

  const cs_field_t *field        = cs_field_by_id(eqc->var_field_id);
  const cs_lnum_t   n_i_faces    = cs_shared_connect->n_faces[CS_INT_FACES];
  const cs_real_t  *bface_values = eqc->face_values + n_i_faces;

  /* In case of postprocessing of the border faces, one has to check if there
     is a mesh modification. In particular, a removal of 2D extruded border
     faces*/

  bool use_parent = (cs_shared_quant->remove_boundary_faces) ? false : true;

  /* Field post-processing */

  char *postlabel = nullptr;
  int   len       = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true, /* interlaced arrays */
                    use_parent,
                    CS_POST_TYPE_cs_real_t,
                    nullptr,              /* values on cells */
                    nullptr,              /* values at internal faces */
                    bface_values,         /* values at border faces */
                    cs_shared_time_step); /* time step management structure */

  BFT_FREE(postlabel);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at mesh cells from the inverse operation
 *         w.r.t. the static condensation (DoF used in the linear system are
 *         located at primal faces)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size 3*n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_macfb_vecteq_get_cell_values(void *context, bool previous)
{
  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)context;

  if (eqc == nullptr)
    return nullptr;

  cs_field_t *pot = cs_field_by_id(eqc->var_field_id);

  if (previous)
    return pot->val_pre;
  else
    return pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size 3*n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_macfb_vecteq_get_face_values(void *context, bool previous)
{
  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)context;

  if (eqc == nullptr)
    return nullptr;

  if (previous) {
    assert(eqc->face_values_pre != nullptr);
    return eqc->face_values_pre;
  }
  else
    return eqc->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_read_restart(cs_restart_t *restart,
                             const char   *eqname,
                             void         *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */

  if (restart == nullptr)
    return;
  if (eqname == nullptr)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is nullptr", __func__);
  if (scheme_context == nullptr)
    bft_error(
      __FILE__, __LINE__, 0, " %s: Scheme context is nullptr", __func__);

  int                retcode = CS_RESTART_SUCCESS;
  cs_macfb_vecteq_t *eqc     = (cs_macfb_vecteq_t *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */

  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Check section */

  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     i_ml_id,
                                     3, /* vector-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      i_ml_id,
                                      3, /* vector-valued */
                                      CS_TYPE_cs_real_t,
                                      eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id  = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  cs_real_t *b_values = eqc->face_values + cs_shared_quant->n_i_faces;

  /* Define the section name */

  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Check section */

  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     b_ml_id,
                                     3, /* vector-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      b_ml_id,
                                      3, /* vector-valued */
                                      CS_TYPE_cs_real_t,
                                      b_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_vecteq_write_restart(cs_restart_t *restart,
                              const char   *eqname,
                              void         *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */

  if (restart == nullptr)
    return;
  if (eqname == nullptr)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is nullptr", __func__);

  const cs_macfb_vecteq_t *eqc = (const cs_macfb_vecteq_t *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */

  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Write interior face section */

  cs_restart_write_section(restart,
                           sec_name,
                           i_ml_id,
                           3, /* vector-valued */
                           CS_TYPE_cs_real_t,
                           eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  const cs_real_t *b_values = eqc->face_values + cs_shared_quant->n_i_faces;

  /* Define the section name */

  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Write boundary face section */

  cs_restart_write_section(restart,
                           sec_name,
                           b_ml_id,
                           3, /* vector-valued */
                           CS_TYPE_cs_real_t,
                           b_values);
}

END_C_DECLS

/*----------------------------------------------------------------------------*/
