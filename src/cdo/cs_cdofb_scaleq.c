/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of scalar-valued equations with source terms
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
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
 *  Header for the current file
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
static const cs_matrix_structure_t  *cs_shared_ms;

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

  BFT_MALLOC(cb->ids, n_fc, short int);
  memset(cb->ids, 0, n_fc*sizeof(short int));

  int  size = n_fc*(n_fc+1);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(double));

  size = 2*n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_fc + 1); /* +1 -> if used as a mass matrix */
  cb->loc = cs_sdm_square_create(n_fc + 1);
  cb->aux = cs_sdm_square_create(n_fc + 1);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *
 * \param[in]      t_eval        time at which one evaluates BCs
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      eqp           pointer to a cs_equation_param_t structure
 * \param[in, out] eqb           pointer to a cs_equation_builder_t structure
 * \param[in, out] p_dir_values  pointer to the Dirichlet values to set
 */
/*----------------------------------------------------------------------------*/

static void
_setup_bc(cs_real_t                     t_eval,
          const cs_mesh_t              *mesh,
          const cs_equation_param_t    *eqp,
          cs_equation_builder_t        *eqb,
          cs_real_t                    *p_dir_values[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, quant->n_b_faces, cs_real_t);
  memset(dir_values, 0, quant->n_b_faces*sizeof(cs_real_t));

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   cs_shared_connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0], /* static variable */
                                   dir_values);
  *p_dir_values = dir_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdofb_scaleq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_fb_cell_system(const cs_flag_t               cell_flag,
                     const cs_cell_mesh_t         *cm,
                     const cs_equation_param_t    *eqp,
                     const cs_equation_builder_t  *eqb,
                     const cs_cdofb_scaleq_t      *eqc,
                     const cs_real_t               dir_values[],
                     const cs_real_t               field_tn[],
                     cs_real_t                     t_eval,
                     cs_cell_sys_t                *csys,
                     cs_cell_builder_t            *cb)
{
  /* Cell-wise view of the linear system to build */
  const int  n_dofs = cm->n_fc + 1;

  csys->c_id = cm->c_id;
  csys->cell_flag = cell_flag;
  csys->n_dofs = n_dofs;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_square_init(n_dofs, csys->mat);

  for (short int f = 0; f < cm->n_fc; f++) {
    csys->dof_ids[f] = cm->f_ids[f];
    csys->val_n[f] = eqc->face_values[cm->f_ids[f]];
  }
  csys->dof_ids[cm->n_fc] = cm->c_id;
  csys->val_n[cm->n_fc] = field_tn[cm->c_id];

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    cs_equation_fb_set_cell_bc(cm,
                               cs_shared_connect,
                               cs_shared_quant,
                               eqp,
                               eqb->face_bc,
                               dir_values,
                               t_eval,
                               csys,
                               cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

  /* Set the properties for this cell if not uniform */
  cs_equation_init_properties_cw(eqp, eqb, t_eval, cell_flag, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in CDO-Fb schemes.
 *
 * \param[in]      time_eval   time at which analytic function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_fb_advection_diffusion_reaction(double                         time_eval,
                                 const cs_equation_param_t     *eqp,
                                 const cs_equation_builder_t   *eqb,
                                 const cs_cdofb_scaleq_t       *eqc,
                                 const cs_cell_mesh_t          *cm,
                                 cs_face_mesh_t                *fm,
                                 cs_cell_sys_t                 *csys,
                                 cs_cell_builder_t             *cb)
{
  CS_UNUSED(eqb);
  CS_UNUSED(fm);
  CS_UNUSED(time_eval);

  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */
    eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

    /* Add the local diffusion operator to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp)) {  /* ADVECTION TERM
                                                 * ============== */

    /* Define the local advection matrix and store the advection
       fluxes across primal faces */
    cs_cdofb_advection_build(eqp, cm, time_eval, eqc->adv_func, cb);

    /* Add it to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Local system after advection", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */

    /* Build the mass matrix adn store it in cb->hdg */
    eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Local mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, cb->hdg);
    }
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) {  /* REACTION TERM
                                               * ============= */

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) {

      /* Use a \mathbb{P}_0 reconstruction in the cell
       *
       * Update the local system with reaction term. Only the row attached to
       * the current cell is involved */
      assert(csys->mat->n_cols == csys->n_dofs);
      double  *c_row = csys->mat->val + cm->n_fc*csys->n_dofs;
      c_row[cm->n_fc] += cb->rpty_val * cm->vol_c;

    }
    else {

      assert(eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_COST);
      assert(eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX);

      /* Update local system matrix with the reaction term
         cb->hdg corresponds to the current mass matrix */
      cs_sdm_add_mult(csys->mat, cb->rpty_val, cb->hdg);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Local system after reaction", csys);
#endif
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the part of boundary conditions that should be done before
 *          the static condensation and the time scheme (case of CDO-Fb schemes)
 *
 * \param[in]      time_eval   time at which analytical function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_fb_apply_bc_partly(cs_real_t                      time_eval,
                    const cs_equation_param_t     *eqp,
                    const cs_cdofb_scaleq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  CS_UNUSED(time_eval);

  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed BEFORE the static condensation */
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann)
      for (short int f  = 0; f < cm->n_fc; f++)
        csys->rhs[f] += csys->neu_values[f];

    /* Weakly enforced Dirichlet BCs for cells attached to the boundary
       csys is updated inside (matrix and rhs) */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
          eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
        eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

    if (cs_equation_param_has_convection(eqp)) { /* Always weakly enforced */
      eqc->adv_func_bc(eqp, cm, cb, csys);
    }

  } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
  if (cs_dbg_cw_test(eqp, cm, csys))
    cs_cell_sys_dump(">> Local system matrix after BC & before condensation",
                     csys);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system in CDO-Fb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_fb_apply_remaining_bc(const cs_equation_param_t     *eqp,
                       const cs_cdofb_scaleq_t       *eqc,
                       const cs_cell_mesh_t          *cm,
                       cs_face_mesh_t                *fm,
                       cs_cell_sys_t                 *csys,
                       cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed AFTER the static condensation */
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
        eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

      /* Enforced Dirichlet BCs for cells attached to the boundary
       * csys is updated inside (matrix and rhs). This is close to a strong
       * way to enforce Dirichlet BCs */
      eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Fb scheme
 *
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] x        solution of the linear system (in: initial guess)
 * \param[in, out] b        right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

static int
_solve_fbs_system(cs_sles_t                    *sles,
                  const cs_matrix_t            *matrix,
                  const cs_equation_param_t    *eqp,
                  cs_real_t                    *x,
                  cs_real_t                    *b)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  cs_range_set_t  *rset = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_real_t  *xsol = NULL;

  const cs_lnum_t  n_scatter_elts = n_faces;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

  if (n_cols > n_scatter_elts) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_elts*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,        /* stride */
                                              n_scatter_elts,
                                              matrix,
                                              rset,
                                              xsol, b);

  /* Solve the linear solver */
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  sles_param = eqp->sles_param;

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    sles_param.eps,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, n_iters, residual, nnz);

  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, x);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_SCALEQ_DBG,
                        x, b, n_faces);
#endif

  /* Free what can be freed at this stage */
  cs_sles_free(sles);

  if (n_cols > n_scatter_elts)
    BFT_FREE(xsol);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables related to CDO-Fb system after a resolution
 *
 * \param[in, out] tce     pointer to a timer counter
 * \param[in, out] fld     pointer to a cs_field_t structure
 * \param[in, out] eqc     pointer to a context structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_fields(cs_timer_counter_t      *tce,
               cs_field_t              *fld,
               cs_cdofb_scaleq_t       *eqc)
{
  cs_timer_t  t0 = cs_timer_time();

  /* Copy current field values to previous values */
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
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_init_common(const cs_cdo_quantities_t     *quant,
                            const cs_cdo_connect_t        *connect,
                            const cs_time_step_t          *time_step,
                            const cs_matrix_structure_t   *ms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms = ms;

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
 * \brief  Get the pointer to the related cs_matrix_structure_t
 *
 * \return a  pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_structure_t *
cs_cdofb_scaleq_matrix_structure(void)
{
  return cs_shared_ms;
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
cs_cdofb_scaleq_finalize_common(void)
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
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO face-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];

  cs_cdofb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdofb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */
  eqc->n_dofs = n_faces + n_cells;

  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PF | CS_CDO_LOCAL_DEQ |
    CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ;

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

  /* Diffusion part */
  /* -------------- */

  eqc->get_stiffness_matrix = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqc->get_stiffness_matrix = cs_hodge_fb_cost_get_stiffness;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->get_stiffness_matrix = cs_hodge_fb_voro_get_stiffness;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Diffusion part */

  /* Dirichlet boundary condition enforcement */
  eqc->enforce_dirichlet = NULL;
  switch (eqp->enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_sfb_weak_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_sfb_wsym_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Advection part */
  eqc->adv_func = NULL;
  eqc->adv_func_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    cs_xdef_type_t  adv_deftype =
      cs_advection_field_get_deftype(eqp->adv_field);

    if (adv_deftype == CS_XDEF_BY_ANALYTIC_FUNCTION)
      eqb->msh_flag |= CS_CDO_LOCAL_FEQ;

    /* Boundary conditions for advection */
    eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FEQ;

    switch (eqp->adv_formulation) {

    case CS_PARAM_ADVECTION_FORM_CONSERV:
      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        if (cs_equation_param_has_diffusion(eqp)) {
          eqc->adv_func = cs_cdo_advection_fb_upwcsv_di;
          eqc->adv_func_bc = cs_cdo_advection_fb_bc_wdi;
        }
        else {
          eqc->adv_func = cs_cdo_advection_fb_upwcsv;
          eqc->adv_func_bc = cs_cdo_advection_fb_bc;
        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for face-based discretization",
                  __func__);

      } /* Scheme */
      break; /* Conservative formulation */

    case CS_PARAM_ADVECTION_FORM_NONCONS:
      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        if (cs_equation_param_has_diffusion(eqp)) {
          eqc->adv_func = cs_cdo_advection_fb_upwnoc_di;
          eqc->adv_func_bc = cs_cdo_advection_fb_bc_wdi;
        }
        else {
          eqc->adv_func = cs_cdo_advection_fb_upwnoc;
          eqc->adv_func_bc = cs_cdo_advection_fb_bc;
        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for face-based discretization",
                  __func__);

      } /* Scheme */
      break; /* Non-conservative formulation */

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of formulation for the advection term",
                __func__);

    } /* Switch on the formulation */

  }

  /* Reaction part */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {
      eqb->msh_flag |= CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
      eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
    }

  } /* Reaction */

  /* Time part */
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) {
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    }
    else if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {
      if (eqp->do_lumping)
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
      else {
        eqb->msh_flag |= CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
      }
    }

  }

  /* Source term part */
  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cells, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_cells);

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t struct. */
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_FB;
  eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_COST;
  eqc->hdg_mass.coef = 1.0; /* not useful in this case */

  eqc->get_mass_matrix = cs_hodge_fb_get_mass;

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

  /* By default, 0 is set as initial condition for the computational domain */
  memset(f_vals, 0, quant->n_faces*sizeof(cs_real_t));
  memset(c_vals, 0, quant->n_cells*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */
    cs_flag_t  f_dof_flag = CS_FLAG_SCALAR | cs_flag_primal_face;
    cs_flag_t  c_dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(f_dof_flag, def, f_vals);
        cs_evaluate_potential_by_value(c_dof_flag, def, c_vals);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        {
          const cs_param_dof_reduction_t  red = eqp->dof_reduction;
          switch (red) {

          case CS_PARAM_REDUCTION_DERHAM:
            cs_evaluate_potential_by_analytic(f_dof_flag, def, t_eval, f_vals);
            cs_evaluate_potential_by_analytic(c_dof_flag, def, t_eval, c_vals);
            break;
          case CS_PARAM_REDUCTION_AVERAGE:
            cs_evaluate_average_on_faces_by_analytic(def, t_eval, f_vals);
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

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */
  const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;

  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0],
                                   f_vals + quant->n_i_faces);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Fb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_real_t  time_eval = 0; /* dummy variable */

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  _setup_bc(time_eval, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_faces, cs_real_t);
# pragma omp parallel for if  (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, rs,           \
         dir_values, fld, cs_cdofb_cell_sys, cs_cdofb_cell_bld)
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

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_fb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                           dir_values, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      _fb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

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
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      } /* End of term source */

      /* BOUNDARY CONDITIONS + STATIC CONDENSATION
       * ========================================= */

      /* Apply a part of BC before the static condensation */
      _fb_apply_bc_partly(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* STATIC CONDENSATION
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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of boundary conditions */
      _fb_apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int f = 0; f < cm->n_fc; f++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->f_ids[f]] += csys->rhs[f];

      if (eqc->source_terms != NULL) { /* Source term */

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        eqc->source_terms[c_id] = csys->source[cm->n_fc];

      }

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system */
  _solve_fbs_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    eqc->face_values, rhs);

  /* Update field */
  _update_fields(&(eqb->tce), fld, eqc);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar
 *         convection/diffusion/reaction equation with a CDO-Fb scheme and an
 *         implicit Euler scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_implicit(const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Sanity checks */
  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t0 = cs_timer_time();

  /* Store the current face values as previous */
  memcpy(eqc->face_values_pre, eqc->face_values,
         quant->n_faces*sizeof(cs_real_t));

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  _setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_faces, cs_real_t);
# pragma omp parallel for if  (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, rs,           \
         dir_values, fld, cs_cdofb_cell_sys, cs_cdofb_cell_bld)
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

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_fb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                           dir_values, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      _fb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

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
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      } /* End of term source */

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      _fb_apply_bc_partly(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                      or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Simply add an entry in mat[cell, cell] */
        csys->rhs[cm->n_fc] += ptyc * csys->val_n[cm->n_fc];
        csys->mat->val[cm->n_fc*csys->n_dofs + cm->n_fc] += ptyc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = cb->hdg;

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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */

      _fb_apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int f = 0; f < cm->n_fc; f++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->f_ids[f]] += csys->rhs[f];

      if (eqc->source_terms != NULL) { /* Source term */

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        eqc->source_terms[c_id] = csys->source[cm->n_fc];

      }

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system */
  _solve_fbs_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    eqc->face_values, rhs);

  /* Update field */
  _update_fields(&(eqb->tce), fld, eqc);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar
 *         convection/diffusion/reaction equation with a CDO-Fb scheme and an
 *         implicit/explicit theta scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_solve_theta(const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;
  const double  tcoef = 1 - eqp->theta;

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Sanity checks */
  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO ||
         eqp->time_scheme == CS_TIME_SCHEME_THETA);

  cs_timer_t  t0 = cs_timer_time();

  /* Store the current face values as previous */
  memcpy(eqc->face_values_pre, eqc->face_values,
         quant->n_faces*sizeof(cs_real_t));

  /* Detect the first call (in this case, we compute the initial source term)*/
  _Bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  /* Should not be t_eval since one sets the Dirichlet values */
  _setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_faces, cs_real_t);
# pragma omp parallel for if  (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, rs,           \
         dir_values, fld, cs_cdofb_cell_sys, cs_cdofb_cell_bld,         \
         compute_initial_source)
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

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_fb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                           dir_values, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      _fb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

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
                                          NULL,  /* No input structure */
                                          cb,    /* mass matrix is cb->hdg */
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
                                        t_cur + dt_cur,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        csys->rhs[cm->n_fc] += eqp->theta * csys->source[cm->n_fc];

      } /* End of term source */

       /* First part of BOUNDARY CONDITIONS
        *               ===================
        * Apply a part of BC before time (csys->mat is going to be multiplied
        * by theta when applying the time scheme) */
      _fb_apply_bc_partly(time_eval, eqp, eqc, cm, fm, csys, cb);

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
      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Only the cell row is involved in the time evolution */
        csys->rhs[cm->n_fc] += ptyc*csys->val_n[cm->n_fc];

        /* Simply add an entry in mat[cell, cell] */
        csys->mat->val[cm->n_fc*(csys->n_dofs + 1)] += ptyc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = cb->hdg;

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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */

      _fb_apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int f = 0; f < cm->n_fc; f++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->f_ids[f]] += csys->rhs[f];

      /* Reset the value of the source term for the cell DoF
         Source term is only hold by the cell DoF in face-based schemes */
      if (cs_equation_param_has_sourceterm(eqp))
        eqc->source_terms[c_id] = csys->source[cm->n_fc];  /* Source term */

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system. Overwrite face_values with the computed solution */
  _solve_fbs_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    eqc->face_values, rhs);

  /* Update field */
  _update_fields(&(eqb->tce), fld, eqc);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
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
 * \return a pointer to a \ref cs_equation_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_cdofb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];
  const cs_real_t  inv_dtcur = 1./dt_cur;
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  /* Allocate and initialize the structure storing the balance evaluation */
  cs_equation_balance_t  *eb = cs_equation_balance_create(cs_flag_primal_cell,
                                                          quant->n_cells);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, pot, eb, cs_cdofb_cell_bld)
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

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Initialize the properties at a cellwise level */
      cs_equation_init_properties_cw(eqp, eqb, time_eval, cell_flag, cm, cb);

      /* Set the value of the current potential */
      for (short int f = 0; f < cm->n_fc; f++)
        p_cur[f] = eqc->face_values[cm->f_ids[f]];
      p_cur[cm->n_fc] = pot->val[cm->c_id];

      /* Unsteady term + time scheme */
      if (cs_equation_param_has_time(eqp)) {

        /* Set the value of the current potential */
        for (short int f = 0; f < cm->n_fc; f++)
          p_prev[f] = eqc->face_values_pre[cm->f_ids[f]];
        p_prev[cm->n_fc] = pot->val_pre[cm->c_id];

        /* Get the value of the time property */
        const double  tptyc = inv_dtcur * cb->tpty_val;

        /* Assign local matrix to a mass matrix to define */
        cs_sdm_t  *mass_mat = cb->loc;
        assert(mass_mat->n_rows == mass_mat->n_cols);
        assert(mass_mat->n_rows == cm->n_fc + 1);

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG)
          eb->unsteady_term[c_id] += tptyc * cm->vol_c *
            (p_cur[cm->n_fc] - p_prev[cm->n_fc]);
        else
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Not implemented yet.", __func__);

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

      default:                  /* Implicit */
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = p_cur[i];
        break;

      } /* Switch on time scheme */

       /* Reaction term */
      if (cs_equation_param_has_reaction(eqp)) {

        /* Define the local reaction property */
        const double  rpty_val = cb->rpty_val * cm->vol_c;
        eb->reaction_term[c_id] += rpty_val * p_theta[cm->n_fc];

      } /* Reaction */

      /* Diffusion term */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        cs_real_t  *res = cb->values;
        memset(res, 0, (cm->n_fc + 1)*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->loc, p_theta, res);

        eb->diffusion_term[cm->c_id] += res[cm->n_fc];

      } /* End of diffusion */

      /* Advection term */
      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix */
        cs_cdofb_advection_build(eqp, cm, time_eval, eqc->adv_func, cb);

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
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
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
 * \brief  Compute the diffusive flux across each boundary face.
 *         Case of scalar-valued CDO-Fb schemes
 *
 * \param[in]       t_eval    time at which one performs the evaluation
 * \param[in]       eqp       pointer to a cs_equation_param_t structure
 * \param[in]       pot_f     array of values at faces
 * \param[in]       pot_c     array of values at cells
 * \param[in, out]  eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  bflux     pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                   const cs_equation_param_t   *eqp,
                                   const cs_real_t             *pot_f,
                                   const cs_real_t             *pot_c,
                                   cs_equation_builder_t       *eqb,
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

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, bflux, pot_c, pot_f,                 \
         cs_cdofb_cell_bld)
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

#if defined(DEBUG) && !defined(NDEBUG)
    cs_cell_mesh_reset(cm);
#endif

    cs_flag_t  msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ;
    cs_flag_t  add_flag = CS_CDO_LOCAL_DEQ;

    if (eqb->diff_pty_uniform) /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

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
          cs_real_t  *neu_values = cb->values;

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

          cs_equation_compute_neumann_fb(t_eval,
                                         face_bc->def_ids[bf_id],
                                         f,
                                         quant,
                                         eqp,
                                         cm,
                                         neu_values);

          bflux[bf_id] = neu_values[f];
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
          if (!eqb->diff_pty_uniform)
            cs_property_tensor_in_cell(cm,
                                       eqp->diffusion_property,
                                       t_eval,
                                       eqp->diffusion_hodge.inv_pty,
                                       cb->dpty_mat);

          /* Define a local buffer keeping the value of the discrete potential
             for the current cell */
          for (short int ff = 0; ff < cm->n_fc; ff++)
            pot[ff] = pot_f[cm->f_ids[ff]];
          pot[cm->n_fc] = pot_c[c_id];

          /* Compute the boundary flux and store it */
          cs_cdo_diffusion_sfb_cost_flux(f, eqp, cm, pot, cb, bflux + bf_id);
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
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname); /* avoid a compilation warning */
  CS_UNUSED(eqp);

  char *postlabel = NULL;
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_real_t  *face_pdi = cs_cdofb_scaleq_get_face_values(data);

  /* Field post-processing */
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,
                    true,                  /* true = original mesh */
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  /* values on cells */
                    NULL,                  /* values at internal faces */
                    face_pdi + n_i_faces,  /* values at border faces */
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
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_scaleq_get_cell_values(void      *context)
{
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_scaleq_get_face_values(void    *context)
{
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;
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
