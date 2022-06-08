/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of scalar-valued equations with source terms
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_toolbox.h"
#include "cs_cdo_solve.h"
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdofb_scaleq.c

  \brief Build an algebraic CDO face-based system for unsteady
         convection-diffusion-reaction of scalar-valued equations with
         source terms
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_SCALEQ_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */

static cs_cell_sys_t      **cs_cdofb_cell_sys = NULL;
static cs_cell_builder_t  **cs_cdofb_cell_bld = NULL;

/* Pointer to shared structures */

static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_cell_builder_create(const cs_cdo_connect_t   *connect)
{
  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->adv_fluxes, n_fc, double);
  memset(cb->adv_fluxes, 0, n_fc*sizeof(double));

  BFT_MALLOC(cb->ids, n_fc, int);
  memset(cb->ids, 0, n_fc*sizeof(int));

  int  size = n_fc*(n_fc+1);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(double));

  size = 2*n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of operators */

  cb->loc = cs_sdm_square_create(n_fc + 1);
  cb->aux = cs_sdm_square_create(n_fc + 1);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings and build the
 *         list of DoFs associated to an (internal) enforcement
 *
 * \param[in]      t_eval    time at which one evaluates BCs
 * \param[in]      mesh      pointer to a cs_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_setup(cs_real_t                     t_eval,
       const cs_mesh_t              *mesh,
       const cs_equation_param_t    *eqp,
       cs_equation_builder_t        *eqb)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Compute the values of the Dirichlet BC */

  BFT_MALLOC(eqb->dir_values, quant->n_b_faces, cs_real_t);
  memset(eqb->dir_values, 0, quant->n_b_faces*sizeof(cs_real_t));

  cs_equation_compute_dirichlet_fb(mesh, quant, connect, eqp, eqb->face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0], /* static variable */
                                   eqb->dir_values);

  /* Internal enforcement of DoFs  */

  if (cs_equation_param_has_internal_enforcement(eqp))
    eqb->enforced_values =
      cs_enforcement_define_at_faces(connect,
                                     eqp->n_enforcements,
                                     eqp->enforcement_params);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      val_f_pre   face values used as the previous one
 * \param[in]      val_c_pre   cell values used as the previous one
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_sfb_init_cell_system(const cs_cell_mesh_t         *cm,
                      const cs_equation_param_t    *eqp,
                      const cs_equation_builder_t  *eqb,
                      const cs_real_t               val_f_pre[],
                      const cs_real_t               val_c_pre[],
                      cs_cell_sys_t                *csys,
                      cs_cell_builder_t            *cb)
{
  /* Cell-wise view of the linear system to build */

  const int  n_dofs = cm->n_fc + 1;

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;

  /* Initialize the local system */

  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_square_init(n_dofs, csys->mat);

  for (short int f = 0; f < cm->n_fc; f++) {
    csys->dof_ids[f] = cm->f_ids[f];
    csys->val_n[f] = val_f_pre[cm->f_ids[f]];
  }
  csys->dof_ids[cm->n_fc] = cm->c_id;
  csys->val_n[cm->n_fc] = val_c_pre[cm->c_id];

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    cs_equation_bc_set_cw_fb(cm,
                             eqp,
                             eqb->face_bc,
                             eqb->dir_values,
                             csys,
                             cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in scalar-valued CDO-Fb schemes.
 *          mass_hodge could be set to NULL if a Voronoi algo. is used.
 *          Otherwise, the mass matrix is computed.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_sfb_conv_diff_reac(const cs_equation_param_t     *eqp,
                    const cs_equation_builder_t   *eqb,
                    const cs_cdofb_scaleq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_hodge_t                    *mass_hodge,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */
    assert(mass_hodge != NULL);

    /* Build the mass matrix and store it in mass_hodge->matrix */

    eqc->get_mass_matrix(cm, mass_hodge, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Local mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids,
                  mass_hodge->matrix);
    }
#endif
  }

  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */
    assert(diff_hodge != NULL);

    /* Set the diffusion property */

    if (!(eqb->diff_pty_uniform))
      cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                     diff_hodge);

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */

    eqc->get_stiffness_matrix(cm, diff_hodge, cb);

    /* Add the local diffusion operator to the local system */

    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp) &&
      ((cb->cell_flag & CS_FLAG_SOLID_CELL) == 0)) {  /* ADVECTION TERM
                                                       * ============== */

    /* Open hook: Compute the advection flux for the numerical scheme and store
       the advection fluxes across primal faces */

    eqc->advection_open(eqp, cm, csys, eqc->advection_input, cb);

    /* Define the local advection matrix. Boundary conditions are treated at
       this stage since there are always weakly enforced. */

    eqc->advection_main(eqp, cm, csys, eqc->advection_scheme, cb);

    /* Close hook: Modify if needed the computed advection matrix and update
       the local system */

    eqc->advection_close(eqp, cm, csys, cb, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after advection", csys);
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) {  /* REACTION TERM
                                               * ============= */

    /* Update the value of the reaction property(ies) if needed */

    cs_equation_builder_set_reaction_pty_cw(eqp, eqb, cm, cb);

    if (eqp->reaction_hodgep.algo == CS_HODGE_ALGO_VORONOI) {

      /* Use a \mathbb{P}_0 reconstruction in the cell
       *
       * Update the local system with reaction term. Only the row attached to
       * the current cell is involved */

      assert(csys->mat->n_cols == csys->n_dofs);
      double  *c_row = csys->mat->val + cm->n_fc*csys->n_dofs;
      c_row[cm->n_fc] += cb->rpty_val * cm->vol_c;

    }
    else {

      assert(eqp->reaction_hodgep.algo == CS_HODGE_ALGO_COST);
      assert(eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX);

      /* Update local system matrix with the reaction term */

      cs_sdm_add_mult(csys->mat, cb->rpty_val, mass_hodge->matrix);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system after reaction", csys);
#endif
  }

  /* Advanced user function: enable to add any user terms */

  if (cs_equation_param_has_user_hook(eqp)) {

    eqb->hook_function(eqp, eqb, eqc, cm,
                       eqb->hook_context,
                       mass_hodge, diff_hodge,
                       csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after user hook",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the part of boundary conditions that should be done before
 *          the static condensation and the time scheme.
 *          Case of scalar-valued CDO-Fb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_sfb_apply_bc_partly(const cs_equation_param_t     *eqp,
                     const cs_cdofb_scaleq_t       *eqc,
                     const cs_cell_mesh_t          *cm,
                     cs_face_mesh_t                *fm,
                     cs_hodge_t                    *diff_hodge,
                     cs_cell_sys_t                 *csys,
                     cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed BEFORE the static condensation */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann)
      for (short int f  = 0; f < cm->n_fc; f++)
        csys->rhs[f] -= csys->neu_values[f];

    /* Weakly enforced Dirichlet BCs for cells attached to the boundary
       csys is updated inside (matrix and rhs) */

    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_robin) {
        assert(eqc->enforce_robin_bc != NULL);
        eqc->enforce_robin_bc(eqp, cm, fm, diff_hodge, cb, csys);
      }

      if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */

        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
          eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

    }

  } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
  if (cs_dbg_cw_test(eqp, cm, csys))
    cs_cell_sys_dump(">> Cell system matrix after BC & before condensation",
                     csys);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the remaining treatment for the boundary conditions to the
 *          local system.
 *          Case of scalar-valued CDO-Fb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_sfb_apply_remaining_bc(const cs_equation_param_t     *eqp,
                        const cs_equation_builder_t   *eqb,
                        const cs_cdofb_scaleq_t       *eqc,
                        const cs_cell_mesh_t          *cm,
                        cs_face_mesh_t                *fm,
                        cs_hodge_t                    *diff_hodge,
                        cs_cell_sys_t                 *csys,
                        cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed AFTER the static condensation */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

      /* Enforced Dirichlet BCs for cells attached to the boundary
       * csys is updated inside (matrix and rhs). This is close to a strong
       * way to enforce Dirichlet BCs */

      eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

    }

  } /* Boundary cell */

  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    cs_equation_builder_enforce_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the residual normalization at the cellwise level according
 *         to the requested type of renormalization
 *         This function should be called after the static condensation.
 *         Case of CDO face-based scalar-valued system.
 *
 * \param[in]  type        type of renormalization
 * \param[in]  cm          pointer to a cs_cell_mesh_t structure
 * \param[in]  csys        pointer to a cs_cell_sys_t structure
 *
 * \return the value of the cellwise contribution to the normalization of
 *         the residual
 */
/*----------------------------------------------------------------------------*/

static double
_sfb_cw_rhs_normalization(cs_param_resnorm_type_t     type,
                          const cs_cell_mesh_t       *cm,
                          const cs_cell_sys_t        *csys)
{
  double  _rhs_norm = 0;

  switch (type) {

  case CS_PARAM_RESNORM_WEIGHTED_RHS:
    for (short int i = 0; i < cm->n_fc; i++)
      _rhs_norm += cm->pvol_f[i] * csys->rhs[i]*csys->rhs[i];
    break;

  case CS_PARAM_RESNORM_FILTERED_RHS:
    for (short int i = 0; i < cm->n_fc; i++) {
      if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
        continue;
      else if (csys->dof_is_forced[i])
        continue;
      else
        _rhs_norm += csys->rhs[i]*csys->rhs[i];
    }
    break;

  case CS_PARAM_RESNORM_NORM2_RHS:
    for (short int i = 0; i < cm->n_fc; i++)
      _rhs_norm += csys->rhs[i]*csys->rhs[i];
    break;

  default:
    break; /* Nothing to do */

  } /* Type of residual normalization */

  return _rhs_norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step
 *
 * \param[in]      csys   pointer to a cellwise view of the system
 * \param[in, out] block  pointer to a block structure
 * \param[in, out] rhs    right-hand side array
 * \param[in, out] eqc    context for this kind of discretization
 * \param[in, out] asb    pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_sfb_assemble(const cs_cell_sys_t         *csys,
              cs_cdo_system_block_t       *block,
              cs_real_t                   *rhs,
              cs_cdofb_scaleq_t           *eqc,
              cs_cdo_assembly_t           *asb)
{
  assert(asb != NULL && block != NULL && rhs != NULL); /* Sanity check */
  assert(block->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  cs_cdo_system_dblock_t  *db = block->block_pointer;

  /* Matrix assembly */

  db->assembly_func(csys->mat, csys->dof_ids, db->range_set, asb, db->mav);

  /* RHS assembly (only on faces since a static condensation has been performed
     to reduce the size) so that n_dofs = n_fc */

# pragma omp critical
  {
    for (short int f = 0; f < csys->n_dofs; f++)
      rhs[csys->dof_ids[f]] += csys->rhs[f];
  }

  if (eqc->source_terms != NULL) { /* Source term */

    /* Reset the value of the source term for the cell DoF
       Source term is only hold by the cell DoF in face-based schemes */

    eqc->source_terms[csys->c_id] = csys->source[csys->n_dofs];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables related to CDO-Fb system after a resolution
 *
 * \param[in, out] tce       pointer to a timer counter
 * \param[in, out] fld       pointer to a cs_field_t structure
 * \param[in, out] eqc       pointer to a context structure
 * \param[in]      cur2prev  true if one performs "current to previous" op.
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_cell_fields(cs_timer_counter_t      *tce,
                    cs_field_t              *fld,
                    cs_cdofb_scaleq_t       *eqc,
                    bool                     cur2prev)
{
  cs_timer_t  t0 = cs_timer_time();

  /* Copy current field values to previous values */

  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */

  cs_static_condensation_recover_scalar(cs_shared_connect->c2f,
                                        eqc->rc_tilda,
                                        eqc->acf_tilda,
                                        eqc->face_values,
                                        fld->val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(tce, &t0, &t1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Fb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdofb_scaleq_is_initialized(void)
{
  if (cs_cdofb_cell_sys == NULL || cs_cdofb_cell_bld == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         scalar-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_init_sharing(const cs_cdo_quantities_t     *quant,
                             const cs_cdo_connect_t        *connect,
                             const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */

  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Specific treatment for handling openMP */

  BFT_MALLOC(cs_cdofb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdofb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdofb_cell_sys[i] = NULL;
    cs_cdofb_cell_bld[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdofb_cell_sys[t_id] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                                 connect->n_max_fbyc,
                                                 1, NULL);
    cs_cdofb_cell_bld[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdofb_cell_sys[0] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                            connect->n_max_fbyc,
                                            1, NULL);
  cs_cdofb_cell_bld[0] = _cell_builder_create(connect);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_cdofb_cell_sys[t_id];
  *cb = cs_cdofb_cell_bld[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO face-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_finalize_sharing(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdofb_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_cdofb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdofb_cell_sys[0]));
  cs_cell_builder_free(&(cs_cdofb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_cdofb_cell_sys);
  BFT_FREE(cs_cdofb_cell_bld);
  cs_cdofb_cell_bld = NULL;
  cs_cdofb_cell_sys = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO face-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  cs_lnum_t  n_cells = connect->n_cells;
  cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  cs_cdofb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdofb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */

  eqc->n_faces = n_faces;
  eqc->n_dofs = n_faces + n_cells;

  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PF | CS_FLAG_COMP_DEQ |
    CS_FLAG_COMP_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */

  eqb->bd_msh_flag = CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ;

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */

  /* Values at each face (interior and border) i.e. take into account BCs */

  BFT_MALLOC(eqc->face_values, n_faces, cs_real_t);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) eqc->face_values[i] = 0;

  eqc->face_values_pre = NULL;
  if (cs_equation_param_has_time(eqp)) {
    BFT_MALLOC(eqc->face_values_pre, n_faces, cs_real_t);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_faces; i++) eqc->face_values_pre[i] = 0;
  }

  /* Store the last computed values of the field at cell centers and the data
     needed to compute the cell values from the face values.
     No need to synchronize all these quantities since they are only cellwise
     quantities. */

  BFT_MALLOC(eqc->rc_tilda, n_cells, cs_real_t);
  BFT_MALLOC(eqc->acf_tilda, connect->c2f->idx[n_cells], cs_real_t);

  memset(eqc->rc_tilda, 0, sizeof(cs_real_t)*n_cells);
  memset(eqc->acf_tilda, 0, sizeof(cs_real_t)*connect->c2f->idx[n_cells]);

  bool  need_eigen =
    (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
     eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;

  /* Diffusion term */

  eqc->get_stiffness_matrix = NULL;
  eqc->diffusion_hodge = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    eqc->diffusion_hodge = cs_hodge_init_context(connect,
                                                 eqp->diffusion_property,
                                                 &(eqp->diffusion_hodgep),
                                                 true,        /* tensor ? */
                                                 need_eigen); /* eigen ? */

    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_COST:
      eqc->get_stiffness_matrix = cs_hodge_fb_cost_get_stiffness;
      break;

    case CS_HODGE_ALGO_BUBBLE:
      eqc->get_stiffness_matrix = cs_hodge_fb_bubble_get_stiffness;
      break;

    case CS_HODGE_ALGO_VORONOI:
      eqc->get_stiffness_matrix = cs_hodge_fb_voro_get_stiffness;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

    /* If necessary, enrich the mesh flag to account for a property defined
     * by an analytical expression. In this case, one evaluates the definition
     * as the mean value over the cell */

    const cs_xdef_t *diff_def = eqp->diffusion_property->defs[0];
    if (diff_def->type == CS_XDEF_BY_ANALYTIC_FUNCTION)
      eqb->msh_flag |= cs_quadrature_get_flag(diff_def->qtype,
                                              cs_flag_primal_cell);

  } /* Diffusion */

  eqc->enforce_robin_bc = cs_cdo_diffusion_sfb_cost_robin;
  eqc->enforce_sliding = NULL;  /* Not useful in the scalar-valued case */

  /* Dirichlet boundary condition enforcement */

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_FLAG_COMP_PFC | CS_FLAG_COMP_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_sfb_weak_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_FLAG_COMP_PFC | CS_FLAG_COMP_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_sfb_wsym_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Advection */

  cs_cdofb_set_advection_function(eqp, eqb, eqc);

  /* Reaction term */

  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodgep.algo == CS_HODGE_ALGO_COST) {
      eqb->msh_flag |= CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
      eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
    }

    /* If necessary, enrich the mesh flag to account for a property defined
     * by an analytical expression. In this case, one evaluates the definition
     * as the mean value over the cell */

    for (short int ir = 0; ir < eqp->n_reaction_terms; ir++) {
      const cs_xdef_t *rea_def = eqp->reaction_properties[ir]->defs[0];
      if (rea_def->type == CS_XDEF_BY_ANALYTIC_FUNCTION)
        eqb->msh_flag |= cs_quadrature_get_flag(rea_def->qtype,
                                                cs_flag_primal_cell);
    }

  }

  /* Unsteady term */

  if (cs_equation_param_has_time(eqp)) {

    if (eqp->time_hodgep.algo == CS_HODGE_ALGO_VORONOI) {
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    }
    else if (eqp->time_hodgep.algo == CS_HODGE_ALGO_COST) {
      if (eqp->do_lumping)
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
      else {
        eqb->msh_flag |= CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
      }
    }

  }

  /* Mass matrix */

  eqc->mass_hodgep.inv_pty  = false;
  eqc->mass_hodgep.type = CS_HODGE_TYPE_FB;
  eqc->mass_hodgep.algo = CS_HODGE_ALGO_COST;
  eqc->mass_hodgep.coef = cs_math_1ov3;

  eqc->get_mass_matrix = NULL;
  eqc->mass_hodge = NULL;

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) {

    eqc->get_mass_matrix = cs_hodge_fb_get;
    eqc->mass_hodge = cs_hodge_init_context(connect,
                                            NULL,
                                            &(eqc->mass_hodgep),
                                            false,  /* tensor ? */
                                            false); /* eigen ? */

    if (eqp->verbosity > 1) {
      cs_log_printf(CS_LOG_SETUP,
                    "#### Parameters of the mass matrix of the equation %s\n",
                    eqp->name);
      cs_hodge_param_log("Mass matrix", NULL, eqc->mass_hodgep);
    }

  }

  /* Source term */

  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cells, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_cells);

  } /* There is at least one source term */

  /* Helper structures (range set, interface set, matrix structure and all the
     assembly process) */

  cs_cdo_system_helper_t  *sh = NULL;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_DEFAULT,
                                   1,        /* n_col_blocks */
                                   &n_faces, /* col_block_sizes */
                                   1);       /* n_blocks */

  /* Choose the right class of matrix to avoid copy.
   * The way to perform the assembly may change if an external librairy is used
   * for solving the linear system */

  cs_cdo_system_matrix_class_t  matclass;

  switch (eqp->sles_param->solver_class) {

  case CS_PARAM_SLES_CLASS_CS:
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  case CS_PARAM_SLES_CLASS_HYPRE:
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

  cs_cdo_system_add_dblock(sh, 0,
                           matclass,
                           cs_flag_primal_face,
                           n_faces,
                           1,       /* stride */
                           false,   /* interlaced (useless for scalar) */
                           true);   /* unrolled (useless for scalar )*/

  cs_cdo_system_build_block(sh, 0); /* build/set structures */

  eqb->system_helper = sh;

  /* Renormalization of the residual */

  if (eqp->sles_param->resnorm_type == CS_PARAM_RESNORM_WEIGHTED_RHS)
    eqb->msh_flag |= CS_FLAG_COMP_PFC;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_scaleq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdofb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_free_context(void   *data)
{
  cs_cdofb_scaleq_t   *eqc  = (cs_cdofb_scaleq_t *)data;

  if (eqc == NULL)
    return eqc;

  cs_hodge_free_context(&(eqc->diffusion_hodge));
  cs_hodge_free_context(&(eqc->mass_hodge));

  /* Free temporary buffers */

  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->face_values);
  if (eqc->face_values_pre != NULL)
    BFT_FREE(eqc->face_values_pre);

  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->acf_tilda);

  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Fb schemes.
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
cs_cdofb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *c_vals = fld->val;
  cs_real_t  *f_vals = eqc->face_values;

  /* Check that a face interface has been defined */

  if (eqp->n_ic_defs > 0 && cs_glob_n_ranks > 1 && connect->face_ifs == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Interface set structure at faces not allocated.\n",
              __func__);

  /* By default, 0 is set as initial condition for the computational domain */

  memset(f_vals, 0, quant->n_faces*sizeof(cs_real_t));
  memset(c_vals, 0, quant->n_cells*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    cs_lnum_t  *def2f_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
    cs_lnum_t  *def2f_idx = NULL;

    BFT_MALLOC(def2f_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_cdo_sync_vol_def_at_faces(eqp->n_ic_defs, eqp->ic_defs,
                                 def2f_idx,
                                 def2f_ids);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */

      const cs_xdef_t  *def = eqp->ic_defs[def_id];
      const cs_lnum_t  n_f_selected = def2f_idx[def_id+1] - def2f_idx[def_id];
      const cs_lnum_t  *selected_lst = def2f_ids + def2f_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_at_faces_by_value(def,
                                                n_f_selected,
                                                selected_lst,
                                                f_vals);
        cs_evaluate_potential_at_cells_by_value(def, c_vals);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        {
          const cs_param_dof_reduction_t  red = eqp->dof_reduction;
          switch (red) {

          case CS_PARAM_REDUCTION_DERHAM:
            cs_evaluate_potential_at_faces_by_analytic(def,
                                                       t_eval,
                                                       n_f_selected,
                                                       selected_lst,
                                                       f_vals);
            cs_evaluate_potential_at_cells_by_analytic(def, t_eval, c_vals);
            break;

          case CS_PARAM_REDUCTION_AVERAGE:
            cs_evaluate_average_on_faces_by_analytic(def,
                                                     t_eval,
                                                     n_f_selected,
                                                     selected_lst,
                                                     f_vals);
            cs_evaluate_average_on_cells_by_analytic(def, t_eval, c_vals);
            break;

          default:
            bft_error(__FILE__, __LINE__, 0,
                      " %s: Incompatible reduction for equation %s.\n",
                      __func__, eqp->name);
            break;

          } /* Switch on possible reduction types */

        }
        break;

      case CS_XDEF_BY_QOV:      /* TODO */
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

    BFT_FREE(def2f_idx);

    if (fld->val_pre != NULL)
      memcpy(fld->val_pre, c_vals, quant->n_cells*sizeof(cs_real_t));

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */

  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0],
                                   f_vals + quant->n_i_faces);

  if (eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values,
           quant->n_faces*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Fb scheme
 *         Use for interpolation purpose from cell values to face values.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      cell_values  array of cell values
 * \param[in]      field_id     id of the variable field
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in, out] eqb          pointer to a cs_equation_builder_t structure
 * \param[in, out] context      pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_interpolate(const cs_mesh_t            *mesh,
                            const cs_real_t            *cell_values,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_cdo_system_helper_t  *sh = eqb->system_helper;

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces.
   * First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */

  _setup(time_eval, mesh, eqp, eqb);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities if needed */

    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set values used as the previous one (useless here) */

    cs_real_t  *val_f_pre = eqc->face_values;
    cs_real_t  *val_c_pre = fld->val;

    /* ---------------------------------------------
     * Main loop on cells to build the linear system
     * --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _sfb_init_cell_system(cm, eqp, eqb, val_f_pre, val_c_pre,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. */

      _sfb_conv_diff_reac(eqp, eqb, eqc, cm, mass_hodge, diff_hodge, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

        /* Reset the local contribution */

        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      } /* End of term source */

      /* BOUNDARY CONDITIONS + CONDENSATION
       * ================================== */

      /* Apply a part of BC before the condensation */

      _sfb_apply_bc_partly(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      { /* Reduce the system size since one has the knowledge of the cell
           value */

        /* Reshape the local system */

        for (short int i = 0; i < cm->n_fc; i++) {

          double  *old_i = csys->mat->val + csys->n_dofs*i;   /* Old "i" row  */
          double  *new_i = csys->mat->val + cm->n_fc*i;       /* New "i" row */

          for (short int j = 0; j < cm->n_fc; j++)
            new_i[j] = old_i[j];

          /* Update RHS: RHS = RHS - Afc*pc */

          csys->rhs[i] -= cell_values[csys->c_id] * old_i[cm->n_fc];

        }

        csys->n_dofs = cm->n_fc;
        csys->mat->n_rows = csys->mat->n_cols = cm->n_fc;

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after condensation",
                         csys);
#endif

      /* Remaining part of boundary conditions */

      _sfb_apply_remaining_bc(eqp, eqb, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Compute a cellwise norm of the RHS for the normalization of the
         residual during the resolution of the linear system */

      rhs_norm += _sfb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* ASSEMBLY PROCESS
       * ================ */

      _sfb_assemble(csys, sh->blocks[0], rhs, eqc, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             n_faces,
                             rhs,
                             &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  cs_cdo_solve_scalar_system(n_faces,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             rhs_norm,
                             true, /* rhs_redux */
                             sles,
                             eqc->face_values,
                             rhs);

  /* Update field (cell values ar known) */

  memcpy(fld->val, cell_values, quant->n_cells*sizeof(cs_real_t));

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Fb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_cdo_system_helper_t  *sh = eqb->system_helper;

  /* Build an array storing the Dirichlet values at faces
   * First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */

  _setup(time_eval, mesh, eqp, eqb);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set values used as the previous one (useless here) */

    cs_real_t  *val_f_pre = eqc->face_values;
    cs_real_t  *val_c_pre = fld->val;

    /* ---------------------------------------------
     * Main loop on cells to build the linear system
     * --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _sfb_init_cell_system(cm, eqp, eqb, val_f_pre, val_c_pre,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction terms to the local
         system. Mass matrix is computed inside if needed during the building */

      _sfb_conv_diff_reac(eqp, eqb, eqc, cm, mass_hodge, diff_hodge, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

        /* Reset the local contribution */

        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      } /* End of term source */

      /* BOUNDARY CONDITIONS + STATIC CONDENSATION
       * ========================================= */

      /* Apply a part of BC before the static condensation */

      _sfb_apply_bc_partly(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* STATIC CONDENSATION
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_scalar_eq(connect->c2f,
                                       eqc->rc_tilda, eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of boundary conditions */

      _sfb_apply_remaining_bc(eqp, eqb, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Compute a cellwise norm of the RHS for the normalization of the
         residual during the resolution of the linear system */

      rhs_norm += _sfb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* ASSEMBLY PROCESS
       * ================ */

      _sfb_assemble(csys, sh->blocks[0], rhs, eqc, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  if (cur2prev && eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values, sizeof(cs_real_t)*n_faces);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             n_faces,
                             rhs,
                             &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  cs_cdo_solve_scalar_system(n_faces,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             rhs_norm,
                             true, /* rhs_redux */
                             sles,
                             eqc->face_values,
                             rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Update field associated to cells */

  _update_cell_fields(&(eqb->tce), fld, eqc, cur2prev);

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar
 *         convection/diffusion/reaction equation with a CDO-Fb scheme and an
 *         implicit Euler scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_implicit(bool                        cur2prev,
                               const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_cdo_system_helper_t  *sh = eqb->system_helper;
  cs_field_t  *fld = cs_field_by_id(field_id);

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  /* Build an array storing the Dirichlet values at faces
   * Always evaluated at t_cur + dt */

  _setup(ts->t_cur + ts->dt[0], mesh, eqp, eqb);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  time_eval = t_cur + dt_cur;
    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set values used as the previous one */

    cs_real_t  *val_f_pre = NULL, *val_c_pre = NULL;
    if (cur2prev)
      val_f_pre = eqc->face_values, val_c_pre = fld->val;
    else  /* Non-linear process --> values at t^n in *_pre  */
      val_f_pre = eqc->face_values_pre, val_c_pre = fld->val_pre;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _sfb_init_cell_system(cm, eqp, eqb, val_f_pre, val_c_pre,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction terms to the local
         system. Mass matrix is computed inside if needed during the building */

      _sfb_conv_diff_reac(eqp, eqb, eqc, cm, mass_hodge, diff_hodge, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

        /* Reset the local contribution */

        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      } /* End of term source */

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _sfb_apply_bc_partly(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (!(eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm,
                                                 eqp->time_property,
                                                 time_eval);

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                      or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Simply add an entry in mat[cell, cell] */

        csys->rhs[cm->n_fc] += ptyc * csys->val_n[cm->n_fc];
        csys->mat->val[cm->n_fc*csys->n_dofs + cm->n_fc] += ptyc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = mass_hodge->matrix;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */

        /* Update rhs with csys->mat*p^n */

        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */

        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after time treatment",
                         csys);
#endif

      /* STATIC CONDENSATION
       * ===================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_scalar_eq(connect->c2f,
                                       eqc->rc_tilda,
                                       eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */

      _sfb_apply_remaining_bc(eqp, eqb, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Compute a cellwise norm of the RHS for the normalization of the
         residual during the resolution of the linear system */

      rhs_norm += _sfb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* ASSEMBLY PROCESS
       * ================ */

      _sfb_assemble(csys, sh->blocks[0], rhs, eqc, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  if (cur2prev && eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values, sizeof(cs_real_t)*n_faces);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             n_faces,
                             rhs,
                             &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  cs_cdo_solve_scalar_system(n_faces,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             rhs_norm,
                             true, /* rhs_redux */
                             sles,
                             eqc->face_values,
                             rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Update field associated to cells */

  _update_cell_fields(&(eqb->tce), fld, eqc, cur2prev);

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar
 *         convection/diffusion/reaction equation with a CDO-Fb scheme and an
 *         implicit/explicit theta scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_theta(bool                        cur2prev,
                            const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_cdo_system_helper_t  *sh = eqb->system_helper;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO ||
         eqp->time_scheme == CS_TIME_SCHEME_THETA);

  /* Store the current face values as previous */

  memcpy(eqc->face_values_pre, eqc->face_values,
         quant->n_faces*sizeof(cs_real_t));

  /* Detect the first call (in this case, we compute the initial source term)*/

  bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces
   * Should not be t_eval since one sets the Dirichlet values
   * Dirichlet boundary conditions are always evaluated at t_cur + dt
   */

  _setup(ts->t_cur + ts->dt[0], mesh, eqp, eqb);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  inv_dtcur = 1./dt_cur;
    const double  tcoef = 1 - eqp->theta;

    /* Set times at which one evaluates quantities when needed
     * t_pty_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
     * since t^(n+1) = t^n + dt
     */

    cb->t_pty_eval = t_cur + eqp->theta*dt_cur;
    cb->t_bc_eval = t_cur + dt_cur;
    cb->t_st_eval = t_cur + dt_cur;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set values used as the previous one */

    cs_real_t  *val_f_pre = NULL, *val_c_pre = NULL;
    if (cur2prev)
      val_f_pre = eqc->face_values, val_c_pre = fld->val;
    else  /* Non-linear process --> values at t^n in *_pre  */
      val_f_pre = eqc->face_values_pre, val_c_pre = fld->val_pre;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _sfb_init_cell_system(cm, eqp, eqb, val_f_pre, val_c_pre,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction terms to the local
         system. Mass matrix is computed inside if needed during the building */

      _sfb_conv_diff_reac(eqp, eqb, eqc, cm, mass_hodge, diff_hodge, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        if (compute_initial_source) { /* First time step */

          /* Reset the local contribution */

          memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

          cs_source_term_compute_cellwise(eqp->n_source_terms,
                      (cs_xdef_t *const *)eqp->source_terms,
                                          cm,
                                          eqb->source_mask,
                                          eqb->compute_source,
                                          t_cur,
                                          mass_hodge,
                                          cb,
                                          csys->source);

          csys->rhs[cm->n_fc] += tcoef * csys->source[cm->n_fc];

        }
        else { /* Add the contribution of the previous time step */

          csys->rhs[cm->n_fc] += tcoef * eqc->source_terms[cm->c_id];

        }

        /* Reset the local contribution */

        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        csys->rhs[cm->n_fc] += eqp->theta * csys->source[cm->n_fc];

      } /* End of term source */

       /* First part of BOUNDARY CONDITIONS
        *               ===================
        * Apply a part of BC before time (csys->mat is going to be multiplied
        * by theta when applying the time scheme) */

      _sfb_apply_bc_partly(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */

      double  *adr_pn = cb->values;
      cs_sdm_square_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++) /* n_dofs = n_vc */
        csys->rhs[i] -= tcoef * adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */

      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs mass_mat * p_n */

      if (!(eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm,
                                                 eqp->time_property,
                                                 cb->t_pty_eval);

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Only the cell row is involved in the time evolution */

        csys->rhs[cm->n_fc] += ptyc*csys->val_n[cm->n_fc];

        /* Simply add an entry in mat[cell, cell] */

        csys->mat->val[cm->n_fc*(csys->n_dofs + 1)] += ptyc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = mass_hodge->matrix;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */

        /* Update rhs with mass_mat*p^n */

        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */

        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after adding time", csys);
#endif

      /* STATIC CONDENSATION
       * ===================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_scalar_eq(connect->c2f,
                                       eqc->rc_tilda, eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * ===================================== */

      _sfb_apply_remaining_bc(eqp, eqb, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Compute a cellwise norm of the RHS for the normalization of the
         residual during the resolution of the linear system */

      rhs_norm += _sfb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* ASSEMBLY PROCESS
       * ================ */

      _sfb_assemble(csys, sh->blocks[0], rhs, eqc, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  if (cur2prev && eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values, sizeof(cs_real_t)*n_faces);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             n_faces,
                             rhs,
                             &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  cs_cdo_solve_scalar_system(n_faces,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             rhs_norm,
                             true, /* rhs_redux */
                             sles,
                             eqc->face_values,
                             rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Update field associated to cells */

  _update_cell_fields(&(eqb->tce), fld, eqc, cur2prev);

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *         Case of scalar-valued CDO face-based scheme
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdofb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const char *func_name = __func__;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  /* Allocate and initialize the structure storing the balance evaluation */

  cs_cdo_balance_t  *eb = cs_cdo_balance_create(cs_flag_primal_cell,
                                                quant->n_cells);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, ts, eqp, eqb, eqc, pot, eb, cs_cdofb_cell_bld, \
         func_name)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Set times at which one evaluates quantities if needed */

    switch (eqp->time_scheme) {

    case CS_TIME_SCHEME_EULER_EXPLICIT:
      cb->t_pty_eval = t_cur;
      cb->t_bc_eval = t_cur;
      cb->t_st_eval = t_cur;
      break;
    case CS_TIME_SCHEME_CRANKNICO:
      cb->t_pty_eval = t_cur + 0.5*dt_cur;
      cb->t_bc_eval = t_cur + dt_cur;
      cb->t_st_eval = t_cur + dt_cur;
      break;
    case CS_TIME_SCHEME_THETA:
      cb->t_pty_eval = t_cur + eqp->theta*dt_cur;
      cb->t_bc_eval = t_cur + dt_cur;
      cb->t_st_eval = t_cur + dt_cur;
      break;

    default: /* Implicit (Forward Euler or BDF2) */
      cb->t_pty_eval = t_cur + dt_cur;
      cb->t_bc_eval = t_cur + dt_cur;
      cb->t_st_eval = t_cur + dt_cur;
      break;

      } /* Switch on time scheme */

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set inside the OMP section so that each thread has its own value */

    cs_real_t  _p_cur[10], _p_prev[10], _p_theta[10];
    cs_real_t  *p_cur = NULL, *p_prev = NULL, *p_theta = NULL;

    if (connect->n_max_fbyc > 10) {
      BFT_MALLOC(p_cur, connect->n_max_fbyc, cs_real_t);
      BFT_MALLOC(p_prev, connect->n_max_fbyc, cs_real_t);
      BFT_MALLOC(p_theta, connect->n_max_fbyc, cs_real_t);
    }
    else {
      p_cur = _p_cur;
      p_prev = _p_prev;
      p_theta = _p_theta;
    }

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the value of the current potential */

      for (short int f = 0; f < cm->n_fc; f++)
        p_cur[f] = eqc->face_values[cm->f_ids[f]];
      p_cur[cm->n_fc] = pot->val[cm->c_id];

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                      * =========== */
        assert(mass_hodge != NULL);

        /* Build the mass matrix and store it in mass_hodge->matrix */

        eqc->get_mass_matrix(cm, mass_hodge, cb);

      }

      /* Unsteady term + time scheme */

      if (cs_equation_param_has_time(eqp)) {

        if (!(eqb->time_pty_uniform))
          cb->tpty_val = cs_property_value_in_cell(cm,
                                                   eqp->time_property,
                                                   cb->t_pty_eval);

        /* Set the value of the current potential */

        for (short int f = 0; f < cm->n_fc; f++)
          p_prev[f] = eqc->face_values_pre[cm->f_ids[f]];
        p_prev[cm->n_fc] = pot->val_pre[cm->c_id];

        /* Get the value of the time property */

        const double  tptyc = inv_dtcur * cb->tpty_val;

        /* Assign local matrix to a mass matrix to define */

        CS_CDO_OMP_ASSERT(cb->loc->n_rows == cb->loc->n_cols);
        CS_CDO_OMP_ASSERT(cb->loc->n_rows == cm->n_fc + 1);

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG)
          eb->unsteady_term[c_id] += tptyc * cm->vol_c *
            (p_cur[cm->n_fc] - p_prev[cm->n_fc]);
        else
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Not implemented yet.", func_name);

      } /* End of time contribution */

      /* Set p_theta */

      switch (eqp->time_scheme) {

      case CS_TIME_SCHEME_EULER_EXPLICIT:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = p_prev[i];
        break;
      case CS_TIME_SCHEME_CRANKNICO:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = 0.5*(p_cur[i] + p_prev[i]);
        break;
      case CS_TIME_SCHEME_THETA:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = eqp->theta*p_cur[i] + (1-eqp->theta)*p_prev[i];
        break;

      default:                  /* Implicit (Euler or BDF2) */
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = p_cur[i];
        break;

      } /* Switch on time scheme */

      /* Reaction term */

      if (cs_equation_param_has_reaction(eqp)) {

        /* Update the value of the reaction property(ies) if needed */

        cs_equation_builder_set_reaction_pty_cw(eqp, eqb, cm, cb);

        /* Define the local reaction property */

        const double  rpty_val = cb->rpty_val * cm->vol_c;
        eb->reaction_term[c_id] += rpty_val * p_theta[cm->n_fc];

      } /* Reaction */

      /* Diffusion term */

      if (cs_equation_param_has_diffusion(eqp)) {

        assert(diff_hodge != NULL);

        /* Set the diffusion property */

        if (!(eqb->diff_pty_uniform))
          cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                         diff_hodge);

        /* Define the local stiffness matrix: local matrix owned by the cellwise
           builder (store in cb->loc) */

        eqc->get_stiffness_matrix(cm, diff_hodge, cb);

        cs_real_t  *res = cb->values;
        memset(res, 0, (cm->n_fc + 1)*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->loc, p_theta, res);

        eb->diffusion_term[cm->c_id] += res[cm->n_fc];

      } /* End of diffusion */

      /* Advection term */

      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix and store the advection fluxes
           across primal faces (Boundary conditions are treated at this stage
           since there are always weakly enforced) */

        /* TODO: Boundary condition and csys --> set to NULL up to now */

        eqc->advection_main(eqp, cm, NULL, eqc->advection_scheme, cb);

        cs_real_t  *res = cb->values;
        memset(res, 0, (cm->n_fc + 1)*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->loc, p_theta, res);

        eb->advection_term[cm->c_id] += res[cm->n_fc];

      } /* End of advection */

      /* Source term */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */

        cs_real_t  *src = cb->values;
        memset(src, 0, (cm->n_fc + 1)*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        src);

        eb->source_term[cm->c_id] += src[cm->n_fc];

      } /* End of term source contribution */

    } /* Main loop on cells */

    if (p_cur != _p_cur) {
      BFT_FREE(p_cur);
      BFT_FREE(p_prev);
      BFT_FREE(p_theta);
    }

  } /* OPENMP Block */

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    eb->balance[c_id] =
      eb->unsteady_term[c_id]  + eb->reaction_term[c_id]  +
      eb->diffusion_term[c_id] + eb->advection_term[c_id] +
      eb->source_term[c_id];

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  return eb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux accross primal faces.
 *         Interior faces first and then boundary faces.
 *         Values at faces are recovered thanks to the equation builder
 *         Case of scalar-valued CDO-Fb schemes
 *
 * \param[in]       c_values    values for the potential at cells
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   value of the diffusive flux at primal faces
  */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_diff_flux_faces(const cs_real_t             *c_values,
                                const cs_equation_param_t   *eqp,
                                cs_real_t                    t_eval,
                                cs_equation_builder_t       *eqb,
                                void                        *context,
                                cs_real_t                   *diff_flux)
{
  if (diff_flux == NULL)
    return;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* If no diffusion, return after resetting */

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, quant->n_faces*sizeof(cs_real_t));
    return;
  }

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  assert(eqp != NULL && eqb != NULL && eqc!= NULL);
  assert(eqc->diffusion_hodge != NULL);

  cs_timer_t  t0 = cs_timer_time();

  cs_hodge_compute_t  *get_diffusion_hodge =
    cs_hodge_get_func(__func__, eqp->diffusion_hodgep);

  const cs_real_t  *f_values = eqc->face_values;

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(t_eval, quant, connect, eqp, eqb, diff_flux, c_values,         \
         f_values, get_diffusion_hodge, cs_cdofb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */

    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_hodge_t  *diff_hodge = eqc->diffusion_hodge[t_id];
    cs_eflag_t  msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ |
      CS_FLAG_COMP_DEQ;
    cs_cdo_diffusion_cw_flux_t  *compute_flux =
      cs_cdo_diffusion_sfb_get_face_flux;

    /* Set inside the OMP section so that each thread has its own value */

    double  *pot = cs_cdo_local_get_d_buffer(t_id);
    double  *flx = pot + cm->n_fc + 1;

    /* Set times at which one evaluates quantities if needed */

    cb->t_pty_eval = cb->t_bc_eval = cb->t_st_eval = t_eval;

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_hodge_set_property_value(0, cb->t_pty_eval, 0, diff_hodge);

    /* Define the flux by cellwise contributions */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform)
        cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, 0, diff_hodge);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */

      for (short int f = 0; f < cm->n_fc; f++)
        pot[f] = f_values[cm->f_ids[f]];
      pot[cm->n_fc] = c_values[cm->c_id];

      /* Compute the local Hodge operator */

      get_diffusion_hodge(cm, diff_hodge, cb);

      /* Compute and store the flux */

      compute_flux(cm, pot, diff_hodge, cb, flx);

      /* Store the fluxes (the flux for the face f seen from the cell c_id is
       * equal to minus the flux from the same face but from the adjacent cell.
       * So, one only keeps the flux whose the face normal is aligned with the
       * arbitrary direction chosen for this face */

      for (short int f = 0; f < cm->n_fc; f++)
        diff_flux[cm->f_ids[f]] = flx[f];

    } /* Loop on cells */

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an approximation of the the diffusive flux across each
 *         boundary face.
 *         Case of scalar-valued CDO-Fb schemes
 *
 * \param[in]       t_eval    time at which one performs the evaluation
 * \param[in]       eqp       pointer to a cs_equation_param_t structure
 * \param[in]       pot_f     array of values at faces
 * \param[in]       pot_c     array of values at cells
 * \param[in, out]  eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  context   pointer to a scheme builder structure
 * \param[in, out]  bflux     pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                   const cs_equation_param_t   *eqp,
                                   const cs_real_t             *pot_f,
                                   const cs_real_t             *pot_c,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context,
                                   cs_real_t                   *bflux)
{
  if (bflux == NULL)
    return;

  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(bflux, 0, quant->n_b_faces*sizeof(cs_real_t));

    cs_timer_t  t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
    return;
  }

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  assert(eqc->diffusion_hodge != NULL);

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)                   \
  shared(quant, connect, eqp, eqb, eqc, bflux, pot_c, pot_f,            \
         cs_cdofb_cell_bld)                                             \
  firstprivate(t_eval)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;
    const cs_adjacency_t  *f2c = connect->f2c;
    const cs_lnum_t  fidx_shift = f2c->idx[quant->n_i_faces];

    cs_real_t  *pot = NULL;
    BFT_MALLOC(pot, connect->n_max_fbyc + 1, cs_real_t); /* +1 for cell */

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */

    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];

    assert(diff_hodge != NULL);

    /* Set times at which one evaluates quantities if needed */

    cb->t_pty_eval = cb->t_bc_eval = cb->t_st_eval = t_eval;

    cs_eflag_t  msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ;
    cs_eflag_t  add_flag = CS_FLAG_COMP_DEQ;

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_hodge_set_property_value(0, cb->t_pty_eval, 0, diff_hodge);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {

      const cs_lnum_t  f_id = bf_id + quant->n_i_faces;
      const cs_lnum_t  c_id = f2c->ids[bf_id + fidx_shift];

      switch (face_bc->flag[bf_id]) {

      case CS_CDO_BC_HMG_NEUMANN:
        bflux[bf_id] = 0.;
        break;

      case CS_CDO_BC_NEUMANN:
        {
          /* Set the local mesh structure for the current cell */

          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

          cs_equation_compute_neumann_sfb(cb->t_bc_eval,
                                          face_bc->def_ids[bf_id],
                                          f,
                                          eqp,
                                          cm,
                                          bflux + bf_id);
        }
        break;

      default:
        { /* Reconstruct a normal flux at the boundary face */

          /* Set the local mesh structure for the current cell */

          cs_cell_mesh_build(c_id, msh_flag | add_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
          if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif
          /* Set the diffusion property */

          if (!(eqb->diff_pty_uniform))
            cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, 0, diff_hodge);

          /* Define a local buffer keeping the value of the discrete potential
             for the current cell */

          for (short int ff = 0; ff < cm->n_fc; ff++)
            pot[ff] = pot_f[cm->f_ids[ff]];
          pot[cm->n_fc] = pot_c[c_id];

          /* Compute the boundary flux and store it */

          cs_cdo_diffusion_sfb_cost_flux(f, cm, pot, diff_hodge,
                                         cb, bflux + bf_id);
        }
        break;

      } /* End of switch */

    } /* End of loop on boundary faces */

    BFT_FREE(pot);

  } /* End of OpenMP block */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(eqc->var_field_id);

  /* Face values */

  if (eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values,
           sizeof(cs_real_t)*eqc->n_faces);

  /* Cell values */

  cs_field_current_to_previous(fld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
  CS_UNUSED(eqp);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  const cs_field_t  *field = cs_field_by_id(eqc->var_field_id);
  const cs_lnum_t  n_i_faces = cs_shared_connect->n_faces[CS_INT_FACES];
  const cs_real_t  *bface_values = eqc->face_values + n_i_faces;

  /* In case of postprocessing of the border faces, one has to check if there
     is a mesh modification. In particular, a removal of 2D extruded border
     faces*/

  bool  use_parent = (cs_shared_quant->remove_boundary_faces) ? false : true;

  /* Field post-processing */

  char *postlabel = NULL;
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,                  /* interlaced arrays */
                    use_parent,
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  /* values on cells */
                    NULL,                  /* values at internal faces */
                    bface_values,          /* values at border faces */
                    cs_shared_time_step);  /* time step management structure */


  BFT_FREE(postlabel);

  cs_timer_t  t1 = cs_timer_time();
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
 * \return  a pointer to an array of cs_real_t (size n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_scaleq_get_cell_values(void      *context,
                                bool       previous)
{
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;

  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

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
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_scaleq_get_face_values(void    *context,
                                bool     previous)
{
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;

  if (previous) {
    assert(eqc->face_values_pre != NULL);
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
cs_cdofb_scaleq_read_restart(cs_restart_t    *restart,
                             const char      *eqname,
                             void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */

  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);
  if (scheme_context == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme context is NULL", __func__);

  int retcode = CS_RESTART_SUCCESS;
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */

  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Check section */

  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     i_ml_id,
                                     1, /* scalar-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      i_ml_id,
                                      1, /* scalar-valued */
                                      CS_TYPE_cs_real_t,
                                      eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  cs_real_t  *b_values = eqc->face_values + cs_shared_quant->n_i_faces;

  /* Define the section name */

  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Check section */

  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     b_ml_id,
                                     1, /* scalar-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      b_ml_id,
                                      1, /* scalar-valued */
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
cs_cdofb_scaleq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */

  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);

  const cs_cdofb_scaleq_t  *eqc = (const cs_cdofb_scaleq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */

  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Write interior face section */

  cs_restart_write_section(restart,
                           sec_name,
                           i_ml_id,
                           1,   /* scalar-valued */
                           CS_TYPE_cs_real_t,
                           eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  const cs_real_t  *b_values = eqc->face_values + cs_shared_quant->n_i_faces;

  /* Define the section name */

  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Write boundary face section */

  cs_restart_write_section(restart,
                           sec_name,
                           b_ml_id,
                           1,   /* scalar-valued */
                           CS_TYPE_cs_real_t,
                           b_values);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
