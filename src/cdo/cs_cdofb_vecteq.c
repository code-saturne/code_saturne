/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
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
 * Local headers
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
#include "cs_cdofb_priv.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_VECTEQ_DBG      0

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
 * \brief  Initialize the local builder structure used for building the system
 *         cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_cell_builder_t *
_cell_builder_create(const cs_cdo_connect_t   *connect)
{
  const int  n_fc = connect->n_max_fbyc;
  const int  n_dofs = n_fc + 1;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_dofs, short int);
  memset(cb->ids, 0, n_dofs*sizeof(short int));

  int  size = CS_MAX(n_fc*n_dofs, 6*n_dofs);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_dofs; i++)
    block_sizes[i] = 3;

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_dofs);
  cb->aux = cs_sdm_square_create(n_dofs);

  cb->loc = cs_sdm_block_create(n_dofs, n_dofs, block_sizes, block_sizes);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the part of boundary conditions that should be done before
 *          the static condensation and the time scheme (case of CDO-Fb schemes)
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
_apply_bc_partly(const cs_equation_param_t     *eqp,
                 const cs_cdofb_vecteq_t       *eqc,
                 const cs_cell_mesh_t          *cm,
                 cs_face_mesh_t                *fm,
                 cs_cell_sys_t                 *csys,
                 cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed BEFORE the static condensation */
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann)
      for (short int f  = 0; f < 3*cm->n_fc; f++)
        csys->rhs[f] += csys->neu_values[f];

    /* Weakly enforced Dirichlet BCs for cells attached to the boundary
       csys is updated inside (matrix and rhs) */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
          eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
        eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

    if (csys->has_sliding)
      eqc->enforce_sliding(eqp, cm, fm, cb, csys);

  } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
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
_apply_remaining_bc(const cs_equation_param_t     *eqp,
                    const cs_cdofb_vecteq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed AFTER the static condensation */
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

      /* Enforced Dirichlet BCs for cells attached to the boundary
       * csys is updated inside (matrix and rhs). This is close to a strong
       * way to enforce Dirichlet BCs */
      eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

  } /* Boundary cell */
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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

void
cs_cdofb_vecteq_setup_bc(cs_real_t                     t_eval,
                         const cs_mesh_t              *mesh,
                         const cs_equation_param_t    *eqp,
                         cs_equation_builder_t        *eqb,
                         cs_real_t                    *p_dir_values[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, 3*quant->n_b_faces, cs_real_t);
  memset(dir_values, 0, 3*quant->n_b_faces*sizeof(cs_real_t));

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
 * \param[in]      eqc         pointer to a cs_cdofb_vecteq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_init_cell_system(const cs_flag_t               cell_flag,
                                 const cs_cell_mesh_t         *cm,
                                 const cs_equation_param_t    *eqp,
                                 const cs_equation_builder_t  *eqb,
                                 const cs_cdofb_vecteq_t      *eqc,
                                 const cs_real_t               dir_values[],
                                 const cs_real_t               field_tn[],
                                 cs_real_t                     t_eval,
                                 cs_cell_sys_t                *csys,
                                 cs_cell_builder_t            *cb)
{
  /* Cell-wise view of the linear system to build */
  const int  n_blocks = cm->n_fc + 1;
  const int  n_dofs = 3*n_blocks;

  csys->cell_flag = cell_flag;
  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_block_init(csys->mat, n_blocks, n_blocks, block_sizes, block_sizes);

  /* One has to keep the same numbering for faces between cell mesh and cell
     system */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_lnum_t  f_id = cm->f_ids[f];
    for (int k = 0; k < 3; k++) {
      csys->dof_ids[3*f + k] = 3*f_id + k;
      csys->val_n[3*f + k] = eqc->face_values[3*f_id + k];
    }

  }

  for (int k = 0; k < 3; k++) {

    const cs_lnum_t  dof_id = 3*cm->c_id+k;
    const cs_lnum_t  _shift = 3*cm->n_fc + k;

    csys->dof_ids[_shift] = dof_id;
    csys->val_n[_shift] = field_tn[dof_id];

  }

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 2
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

void
cs_cdofb_vecteq_diffusion(double                         time_eval,
                          const cs_equation_param_t     *eqp,
                          const cs_equation_builder_t   *eqb,
                          const cs_cdofb_vecteq_t       *eqc,
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

    if (eqp->diffusion_hodge.is_iso == false)
      bft_error(__FILE__, __LINE__, 0, " %s: Case not handle yet\n",
                __func__);

    /* Add the local diffusion operator to the local system */
    const cs_real_t  *sval = cb->loc->val;
    for (int bi = 0; bi < cm->n_fc + 1; bi++) {
      for (int bj = 0; bj < cm->n_fc + 1; bj++) {

        /* Retrieve the 3x3 matrix */
        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
        assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

        const cs_real_t  _val = sval[(cm->n_fc+1)*bi+bj];
        bij->val[0] += _val;
        bij->val[4] += _val;
        bij->val[8] += _val;

      }
    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Fb scheme
 *
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      field_id id related to the variable field of this equation
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] x        solution of the linear system (in: initial guess)
 * \param[in, out] b        right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_vecteq_solve_system(cs_sles_t                    *sles,
                             const cs_matrix_t            *matrix,
                             const int                     field_id,
                             const cs_equation_param_t    *eqp,
                             cs_real_t                    *x,
                             cs_real_t                    *b)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* solving info */
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_solving_info_t sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  cs_range_set_t  *rset = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  cs_real_t  *xsol = NULL;

  const cs_lnum_t  n_scatter_elts = 3*n_faces;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

  if (n_cols > n_scatter_elts) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_elts*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,            /* stride */
                                              n_scatter_elts,
                                              matrix,
                                              rset,
                                              xsol, b);

  /* Solve the linear solver */
  sinfo.rhs_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  sles_param = eqp->sles_param;

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    sles_param.eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, sinfo.n_it, sinfo.res_norm, nnz);

  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, x);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_VECTEQ_DBG,
                        x, b, 3*n_faces);
#endif

  /* Free what can be freed at this stage */
  if (n_cols > n_scatter_elts)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector steady-state
 *         diffusion equation with a CDO-Fb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  /* First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */
  cs_cdofb_vecteq_setup_bc(time_eval, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

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
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    mab->n_x_dofs = 3;  /* vector-valued equation */

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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                                       dir_values, fld->val, time_eval,
                                       csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      cs_cdofb_vecteq_diffusion(time_eval, eqp, eqb, eqc, cm, fm, csys, cb);

      const _Bool has_sourceterm = cs_equation_param_has_sourceterm(eqp);
      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */

        cs_cdofb_vecteq_sourceterm(cm, eqp,
                                   time_eval, 1., /* time, scaling */
                                   cb, eqb, csys);

      } /* End of term source */

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      _apply_bc_partly(eqp, eqc, cm, fm, csys, cb);


      /* STATIC CONDENSATION
       * ===================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       eqc->rc_tilda,
                                       eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */
      _apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm,
                               eqc, mab, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values
   * Steady, but let us suppose we have an initial condition */
  cs_field_current_to_previous(fld);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t1, &t2);

  /* Now solve the system. face_val is overwritten during the solve step */
  cs_real_t *face_val = eqc->face_values;
  cs_sles_t *sles = cs_sles_find_or_add(field_id, NULL);

  cs_cdofb_vecteq_solve_system(sles, matrix, field_id, eqp, face_val, rhs);

  /* Update field */
  cs_timer_t  t3 = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                        eqc->rc_tilda, eqc->acf_tilda,
                                        face_val, fld->val);

  cs_timer_t  t4 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t3, &t4);

  /* Free remaining buffers */
  cs_sles_free(sles);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a CDO-Fb scheme and an implicit Euler scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_solve_implicit(const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Sanity checks */
  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

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
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    mab->n_x_dofs = 3;  /* vector-valued equation */

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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                                       dir_values, fld->val, time_eval,
                                       csys, cb);

      const short int  n_f = cm->n_fc;

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      cs_cdofb_vecteq_diffusion(time_eval, eqp, eqb, eqc, cm, fm, csys, cb);

      const _Bool has_sourceterm = cs_equation_param_has_sourceterm(eqp);
      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */

        cs_cdofb_vecteq_sourceterm(cm, eqp,
                                   time_eval, 1., /* time, scaling */
                                   cb, eqb, csys);

      } /* End of term source */

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      _apply_bc_partly(eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                      or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_f, n_f);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_f + k] += ptyc * csys->val_n[3*n_f+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: Only diagonal time treatment "
            "available so far.\n", __func__);

      /* STATIC CONDENSATION
       * ===================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       eqc->rc_tilda,
                                       eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */
      _apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm,
                               eqc, mab, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t1, &t2);

  /* Now solve the system. face_val is overwritten during the solve step */
  cs_real_t *face_val = eqc->face_values;
  cs_sles_t *sles = cs_sles_find_or_add(field_id, NULL);

  cs_cdofb_vecteq_solve_system(sles, matrix, field_id, eqp, face_val, rhs);

  cs_timer_t  t3 = cs_timer_time();

  /* Update field.
   * Compute values at cells pc from values at faces pf
   * pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                        eqc->rc_tilda, eqc->acf_tilda,
                                        face_val, fld->val);

  cs_timer_t  t4 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t3, &t4);

  /* Free remaining buffers */
  cs_sles_free(sles);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a CDO-Fb scheme and an implicit/explicit theta scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_solve_theta(const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  const double  tcoef = 1 - eqp->theta;

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Sanity checks */
  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO ||
         eqp->time_scheme == CS_TIME_SCHEME_THETA);

  cs_timer_t  t0 = cs_timer_time();

  /* Detect the first call (in this case, we compute the initial source term)*/
  _Bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;

  /* Should be t_cur + dt_cur since one sets the Dirichlet values */
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, rs, dir_values, \
         fld, cs_cdofb_cell_sys, cs_cdofb_cell_bld, compute_initial_source)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    mab->n_x_dofs = 3;  /* vector-valued equation */

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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                                       dir_values, fld->val, time_eval,
                                       csys, cb);

      const short int  n_f = cm->n_fc;

      /* Build and add the diffusion/advection/reaction term to the local
         system. */
      cs_cdofb_vecteq_diffusion(time_eval, eqp, eqb, eqc, cm, fm, csys, cb);

      const _Bool has_sourceterm = cs_equation_param_has_sourceterm(eqp);
      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */
        if (compute_initial_source) { /* First time step */

          cs_cdofb_vecteq_sourceterm(cm, eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     cb, eqb, csys);

        }
        else { /* Add the contribution of the previous time step */

          for (short int k = 0; k < 3; k++)
            csys->rhs[3*n_f + k] += tcoef * eqc->source_terms[3*c_id + k];

        }

        cs_cdofb_vecteq_sourceterm(cm, eqp,
                                   t_cur+dt_cur, eqp->theta, /* time, scaling */
                                   cb, eqb, csys);

      } /* End of term source */

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      _apply_bc_partly(eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_block_matvec(csys->mat, csys->val_n, adr_pn);
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

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_f, n_f);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_f + k] += ptyc * csys->val_n[3*n_f+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Only diagonal time treatment available so far.\n",
                  __func__);

      /* STATIC CONDENSATION
       * ===================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       eqc->rc_tilda,
                                       eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* Remaining part of BOUNDARY CONDITIONS
       * =================================== */
      _apply_remaining_bc(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm,
                               eqc, mab, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t1, &t2);

  /* Now solve the system. face_val is overwritten during the solve step */
  cs_real_t  *face_val = eqc->face_values;
  cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

  cs_cdofb_vecteq_solve_system(sles, matrix, field_id, eqp, face_val, rhs);

  cs_timer_t  t3 = cs_timer_time();

  /* Update field
   * Compute values at cells pc from values at faces pf
   * pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                        eqc->rc_tilda, eqc->acf_tilda,
                                        face_val, fld->val);

  cs_timer_t  t4 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t3, &t4);

  /* Free remaining buffers */
  cs_sles_free(sles);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Fb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdofb_vecteq_is_initialized(void)
{
  if (cs_cdofb_cell_sys == NULL || cs_cdofb_cell_bld == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vector-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_init_common(const cs_cdo_quantities_t     *quant,
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

  const short int  n_blocks = connect->n_max_fbyc + 1;
  const short int  n_max_dofs = 3*n_blocks;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _cell_builder_create(connect);
    short int  *block_sizes = cb->ids;
    for (int i = 0; i < n_blocks; i++)
      block_sizes[i] = 3;

    cs_cdofb_cell_sys[t_id] = cs_cell_sys_create(n_max_dofs,
                                                 n_blocks - 1,
                                                 n_blocks,
                                                 block_sizes);
    cs_cdofb_cell_bld[t_id] = cb;
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _cell_builder_create(connect);
  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  cs_cdofb_cell_sys[0] =  cs_cell_sys_create(n_max_dofs,
                                             n_blocks - 1,
                                             n_blocks,
                                             block_sizes);
  cs_cdofb_cell_bld[0] = cb;
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
cs_cdofb_vecteq_matrix_structure(void)
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
cs_cdofb_vecteq_get(cs_cell_sys_t       **csys,
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
cs_cdofb_vecteq_finalize_common(void)
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
 * \brief  Initialize a cs_cdofb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB || eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid type of equation.\n"
              " Expected: vector-valued CDO face-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];

  cs_cdofb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdofb_vecteq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */
  eqc->n_dofs = 3*(n_faces + n_cells);

  eqb->sys_flag = CS_FLAG_SYS_VECTOR;
  eqb->msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE |
    CS_CDO_LOCAL_FEQ;

  BFT_MALLOC(eqc->face_values, 3*n_faces, cs_real_t);
  BFT_MALLOC(eqc->rc_tilda, 3*n_cells, cs_real_t);
# pragma omp parallel if (3*n_cells > CS_THR_MIN)
  {
    /* Values at each face (interior and border) i.e. take into account BCs */
#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < 3*n_faces; i++) eqc->face_values[i] = 0;

    /* Store the last computed values of the field at cell centers and the data
       needed to compute the cell values from the face values.
       No need to synchronize all these quantities since they are only cellwise
       quantities. */
#   pragma omp for
    for (cs_lnum_t i = 0; i < 3*n_cells; i++) eqc->rc_tilda[i] = 0;
  }

  /* Assume the 3x3 matrix is diagonal */
  BFT_MALLOC(eqc->acf_tilda, 3*connect->c2f->idx[n_cells], cs_real_t);
  memset(eqc->acf_tilda, 0, 3*connect->c2f->idx[n_cells]*sizeof(cs_real_t));

  /* Diffusion */
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

  } /* Diffusion */

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_block_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_block_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_vfb_weak_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_vfb_wsym_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  eqc->enforce_sliding = NULL;
  if (eqb->face_bc->n_sliding_faces > 0) {
    /* There is at least one face with a sliding condition to handle */
    eqb->bd_msh_flag |= CS_CDO_LOCAL_HFQ;
    eqc->enforce_sliding = cs_cdo_diffusion_vfb_wsym_sliding;
  }

  /* Advection */
  eqc->adv_func = NULL;
  eqc->adv_func_bc = NULL;

  /* Time */
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

  /* Source term */
  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, 3*n_cells, cs_real_t);
#   pragma omp parallel for if (3*n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < 3*n_cells; i++) eqc->source_terms[i] = 0;

  } /* There is at least one source term */

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_block_matrix;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_vecteq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdofb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_vecteq_free_context(void   *data)
{
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)data;

  if (eqc == NULL)
    return eqc;

  /* Free temporary buffers */
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->face_values);
  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->acf_tilda);

  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of vector-valued CDO-Fb schemes.
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
cs_cdofb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *c_vals = fld->val;
  cs_real_t  *f_vals = eqc->face_values;

  /* By default, 0 is set as initial condition for the computational domain */
  memset(f_vals, 0, 3*quant->n_faces*sizeof(cs_real_t));
  memset(c_vals, 0, 3*quant->n_cells*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */
    cs_flag_t  f_dof_flag = CS_FLAG_VECTOR | cs_flag_primal_face;
    cs_flag_t  c_dof_flag = CS_FLAG_VECTOR | cs_flag_primal_cell;

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
                                   f_vals + 3*quant->n_i_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname);    /* avoid a compilation warning */
  CS_UNUSED(eqp);

  char *postlabel = NULL;
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_real_t  *face_pdi = cs_cdofb_vecteq_get_face_values(data);

  /* Field post-processing */
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_pdi + n_i_faces,  // values at border faces
                    cs_shared_time_step);  // time step management structure


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
cs_cdofb_vecteq_get_cell_values(void      *context)
{
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return  pot->val;
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
cs_cdofb_vecteq_get_face_values(void    *context)
{
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;

  if (eqc == NULL)
    return NULL;
  else
    return  eqc->face_values;
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
cs_cdofb_vecteq_read_restart(cs_restart_t    *restart,
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
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t  *)scheme_context;

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

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  cs_real_t  *b_values = eqc->face_values + 3*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     b_ml_id,
                                     1, /* vector-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */
  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      b_ml_id,
                                      1, /* vector-valued */
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
cs_cdofb_vecteq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);

  const cs_cdofb_vecteq_t  *eqc = (const cs_cdofb_vecteq_t  *)scheme_context;

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
                           3,   /* vector-valued */
                           CS_TYPE_cs_real_t,
                           eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  const cs_real_t  *b_values = eqc->face_values + 3*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Write boundary face section */
  cs_restart_write_section(restart,
                           sec_name,
                           b_ml_id,
                           3,   /* vector-valued */
                           CS_TYPE_cs_real_t,
                           b_values);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
