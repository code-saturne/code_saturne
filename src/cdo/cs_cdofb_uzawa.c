/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an Augmented Lagrangian-Uzawa algorithm
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_blas.h"
#include "cs_cdo_bc.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_navsto.h"
#include "cs_dbg.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_uzawa.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_uzawa.c
 *
 * \brief Build an algebraic CDO face-based system for the Navier-Stokes
 * equations and solved it with an Augmented Lagrangian-Uzawa algorithm
 *
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_uzawa_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         vector-valued unknowns
 */

typedef struct {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_uzawa_t (owned by \ref cs_navsto_system_t)
   *  containing the settings related to a Uzawa algorithm
   */

  cs_navsto_uzawa_t   *coupling_context;

  /*!
   * @name Main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   */

  /*! \var velocity
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the velocity
   */

  cs_field_t  *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t  *pressure;

  /*! \var divergence
   *  Pointer to \ref cs_real_t containing the values of the divergence on the
   *  cells
   */

  cs_field_t  *divergence;

 /*!
   * @}
   * @name Parameters of the algorithm
   * Easy access to useful features and parameters of the algorithm
   * @{
   */

  /*! \var is_gdscale_uniform
   *  Bool telling if the scaling in front of the grad-div operator is uniform.
   */

  bool  is_gdscale_uniform;

  /*! \var residual
   *  Last known value of the residual of the algorithm. It is insightful if
   *  using a limit-steady or steady scheme
   */

  cs_real_t residual;

  /*! \var last_iter
   *  Last actually computed iteration. It is useful when a limit-steady scheme
   *  has reached convergence
   */

  cs_lnum_t last_iter;

  /*!
   * @}
   * @name Performance monitoring
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   * @{
   */

  /*! \var timer
   *  Cumulated elapsed time for building and solving the Navier--Stokes system
   */
  cs_timer_counter_t  timer;

  /*! @} */

} cs_cdofb_uzawa_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_UZAWA_DBG      3

#define _dp3  cs_math_3_dot_product
#define _n3   cs_math_3_norm

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t     *cs_shared_quant;
static const cs_cdo_connect_t        *cs_shared_connect;
static const cs_time_step_t          *cs_shared_time_step;
static const cs_matrix_structure_t   *cs_shared_scal_ms;
static const cs_matrix_structure_t   *cs_shared_vect_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence of a cell using the \ref cs_cdo_quantities_t
 *         structure
 *
 * \param[in]     c_id         cell ID
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_dof        values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_cell_divergence(const cs_lnum_t               c_id,
                     const cs_cdo_quantities_t    *quant,
                     const cs_adjacency_t         *c2f,
                     const cs_real_t              *f_dof)
{
  cs_real_t  div = 0.0;

  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

    const cs_lnum_t  f_id = c2f->ids[f];
    const cs_real_t  *_val = f_dof + 3*f_id;

    if (f_id < quant->n_i_faces) {
      const cs_real_t *_nuf = quant->i_face_normal + 3*f_id;

      div += c2f->sgn[f]*quant->i_face_surf[f_id]*_dp3(_val, _nuf)/_n3(_nuf);

    }
    else {

      const cs_lnum_t  bf_id = f_id - quant->n_i_faces;
      const cs_real_t  *_nuf = quant->b_face_normal + 3*bf_id;

      div += c2f->sgn[f]*quant->b_face_surf[bf_id]*_dp3(_val, _nuf)/_n3(_nuf);

    } /* Boundary face */
  } /* Loop on cell faces */

  div /= quant->cell_vol[c_id];

  return div;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence vector associated to the current cell.
 *         WARNING: mind that, differently form the original definition, the
 *         result here is not divided by the cell volume
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out] div        array related to the divergence operator
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_divergence_vect(const cs_cell_mesh_t  *cm,
                     cs_real_t              div[])
{
  /* D(\hat{u}) = \frac{1}{|c|} \sum_{f_c} \iota_{fc} u_f.f
   * But, when integrating [[ p, q ]]_{P_c} = |c| p_c q_c
   * Thus, the volume in the divergence drops
   */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_real_t  i_f = cm->f_sgn[f] * pfq.meas;

    cs_real_t  *_div_f = div + 3*f;
    _div_f[0] = i_f * pfq.unitv[0];
    _div_f[1] = i_f * pfq.unitv[1];
    _div_f[2] = i_f * pfq.unitv[2];

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the grad-div part to the local matrix (i.e. for the current
 *         cell)
 *
 * \param[in]      n_fc       local number of faces for the current cell
 * \param[in]      zeta       scalar coefficient for the grad-div operator
 * \param[in]      div        divergence
 * \param[in, out] mat        local system matrix to update
 */
/*----------------------------------------------------------------------------*/

static void
_add_grad_div(short int          n_fc,
              const cs_real_t    zeta,
              const cs_real_t    div[],
              cs_sdm_t          *mat)
{
  cs_sdm_t  *b = NULL;

  /* Avoid dealing with cell DoFs which are not impacted */
  for (short int bi = 0; bi < n_fc; bi++) {

    const cs_real_t  *divi = div + 3*bi;
    const cs_real_t  zt_di[3] = {zeta*divi[0], zeta*divi[1], zeta*divi[2]};

    /* Begin with the diagonal block */
    b = cs_sdm_get_block(mat, bi, bi);
    assert(b->n_rows == b->n_cols && b->n_rows == 3);
    for (short int l = 0; l < 3; l++) {
      cs_real_t *m_l = b->val + 3*l;
      for (short int m = 0; m < 3; m++)
        m_l[m] += zt_di[l] * divi[m];
    }

    /* Continue with the extra-diag. blocks */
    for (short int bj = bi+1; bj < n_fc; bj++) {

      b = cs_sdm_get_block(mat, bi, bj);
      assert(b->n_rows == b->n_cols && b->n_rows == 3);
      cs_real_t *mij  = b->val;
      b = cs_sdm_get_block(mat, bj, bi);
      assert(b->n_rows == b->n_cols && b->n_rows == 3);
      cs_real_t *mji  = b->val;

      const cs_real_t *divj = div + 3*bj;

      for (short int l = 0; l < 3; l++) {

        /* Diagonal: 3*l+l = 4*l */
        const cs_real_t  gd_coef_ll = zt_di[l]*divj[l];
        mij[4*l] += gd_coef_ll;
        mji[4*l] += gd_coef_ll;

        /* Extra-diagonal: Use the symmetry of the grad-div */
        for (short int m = l+1; m < 3; m++) {
          const short int  lm = 3*l+m, ml = 3*m+l;
          const cs_real_t  gd_coef_lm = zt_di[l]*divj[m];
          mij[lm] += gd_coef_lm;
          mji[ml] += gd_coef_lm;
          const cs_real_t  gd_coef_ml = zt_di[m]*divj[l];
          mij[ml] += gd_coef_ml;
          mji[lm] += gd_coef_ml;
        }
      }

    } /* Loop on column blocks: bj */
  } /* Loop on row blocks: bi */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a Stokes equation with a CDO
 *         face-based scheme solved with Uzawa-Augmented Lagrangian algo.
 *         One works cellwise and then proceed to the assembly
 *         Remark: The only difference between this algorithm and the standard
 *         artificial compressibility algorithm is that, in the former, the
 *         scaling of the grad.div is not divided by the time step
 *
 * \param[in]      mesh      pointer to a \ref cs_mesh_t structure
 * \param[in]      vel_vals  pointer to the current value of the velocity field
 *                           on cells
 * \param[in]      pr_vals   pointer to the current value of the pressure field
 * \param[in]      dt_cur    current value of the time step
 * \param[in]      zeta      scaling of the grad-div
 * \param[in, out] sc        scheme context for Uzawa algorithm
 * \param[in, out] eq        pointer to the momentum \ref cs_equation_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_system_uzawa(const cs_mesh_t       *mesh,
                    const cs_real_t       *vel_vals,
                    const cs_real_t       *pr_vals,
                    double                 dt_cur,
                    const cs_property_t   *zeta,
                    cs_cdofb_uzawa_t      *sc,
                    cs_equation_t         *eq)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_equation_param_t *eqp = eq->param;
  const cs_real_t  odt = cs_equation_param_has_time(eqp) ? 1./dt_cur : 1.0;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_equation_builder_t  *eqb = eq->builder;
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t*)eq->scheme_context;
  cs_real_t  *rhs = eq->rhs;
  cs_matrix_t  *matrix = eq->matrix;

  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* Dirichlet values at boundary faces are first computed */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, 3*quant->n_b_faces, cs_real_t);
  memset(dir_values, 0, 3*quant->n_b_faces*sizeof(cs_real_t));

  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_cur + dt_cur,
                                   NULL, /* Should be a cell builder but it is
                                            not used*/
                                   dir_values);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,      \
         dir_values, zeta, vel_vals, pr_vals, sc)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  time_eval = t_cur + dt_cur;

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    cs_real_t  zeta_c = cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* For the stationary Stokes problem:
       *
       *     |        |         |
       *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
       *     |        |         |  A is csys->mat in what follows
       *     |--------|---------|  The viscous part arising from the CDO-Fb
       *     |        |         |  schemes for vector-valued variables
       *     |   B    |    0    |
       *     |        |         |
       */

      /* Set the local (i.e. cellwise) structures for the current cell */
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                                       dir_values, vel_vals, time_eval, /* in */
                                       csys, cb);                      /* out */

      const short int  n_fc = cm->n_fc, f_dofs = 3*n_fc;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform)) {
          cs_property_tensor_in_cell(cm,
                                     eqp->diffusion_property,
                                     time_eval,
                                     eqp->diffusion_hodge.inv_pty,
                                     cb->dpty_mat);

          if (eqp->diffusion_hodge.is_iso)
            cb->dpty_val = cb->dpty_mat[0][0];
        }

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        if (eqp->diffusion_hodge.is_iso == false)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Case not handle yet\n", __func__);

        /* Add the local diffusion operator to the local system
           sval stores the stiffness matrix corresponding to the scalar case
         */
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

          } /* Loop on columns blocks: bj */
        }  /* Loop on row blocks: bi */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
        if (cs_dbg_cw_test(eqp, cm, csys))
          cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
      } /* END OF DIFFUSION */

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */
      _get_divergence_vect(cm, cb->aux->val);

      const cs_real_t *div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, csys)) {
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT,
                          "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, div[3*f]*ovol, div[3*f+1]*ovol, div[3*f+2]*ovol);
        } /* critical section */
      }
#endif

      /* Update the property */
      if ( !(sc->is_gdscale_uniform) )
        zeta_c = cs_property_value_in_cell(cm, zeta, time_eval);

      _add_grad_div(n_fc, zeta_c*ovol, div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL, /* No input structure */
                                        cb,   /* mass matrix is cb->hdg */
                                        csys->source);

        for (int k = 0; k < 3; k++) /* DoFs related to the cell velocity */
          csys->rhs[3*cm->n_fc + k] += csys->source[3*cm->n_fc + k];

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        for (int k = 0; k < 3; k++)
          eqc->source_terms[3*c_id + k] = csys->source[3*n_fc + k];

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
        if (cs_dbg_cw_test(eqp, cm, csys))
          cs_cell_sys_dump(">> Local system matrix after"
                           " diffusion and source term", csys);
#endif

      } /* If source term */

      /* 3b- OTHER RHS CONTRIBUTIONS */
      /* =========================== */

      /* Apply the operator gradient to the pressure field and add it to the
         rhs */
      cs_sdm_add_scalvect(f_dofs, pr_vals[c_id], div, csys->rhs);

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
       * Operations that have to be performed BEFORE the static condensation
       */
      if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

        /* Neumann boundary conditions */
        if (csys->has_nhmg_neumann) {
          for (short int f = 0; f < 3*cm->n_fc; f++)
            csys->rhs[f] += csys->neu_values[f];
        }

      } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after "
                         "weak BC and pressure contribution", csys);
#endif

      /* 4- TIME CONTRIBUTION */
      /* ==================== */

      if (cs_equation_param_has_time(eqp)) {
        /* !!! ATTENTION !!!
         * Forcing diagonal implicit
         */
        assert(eqp->time_scheme == CS_TIME_SCHEME_IMPLICIT);
        assert(eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG);
        const short int nf = cm->n_fc;
        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, nf, nf);

        const double  ptyc = cb->tpty_val * cm->vol_c * odt;

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*nf + k] += ptyc * csys->val_n[3*nf+k];
        /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      } /* Time block */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after time contribution and"
                         " before condensation", csys);
#endif

      /* 5- (STATIC CONDENSATION AND) UPDATE OF THE LOCAL SYSTEM */
      /* ========================================================= */
      /* Static condensation of the local system stored inside a block matrix
       * of size n_fc + 1 into a block matrix of size n_fc.
       * Store information in the context structure in order to be able to
       * compute the values at cell centers.
       */

      cs_static_condensation_vector_eq(connect->c2f,
                                       eqc->rc_tilda, eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after condensation", csys);
#endif

      /* 6- BOUNDARY CONDITIONS */
      /* ====================== */
      if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

          /* Weakly enforced Dirichlet BCs for cells attached to the boundary
             csys is updated inside (matrix and rhs) */
          eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

        }

      } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* 7- ASSEMBLY */
      /* =========== */

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];

      /* Matrix assembly */
      cs_equation_assemble_block_matrix(csys, rs, 3, mav);

      /* Assemble RHS */
      for (short int i = 0; i < 3*cm->n_fc; i++) {
#       pragma omp atomic
        rhs[csys->dof_ids[i]] += csys->rhs[i];
      }

    } /* Main loop on cells */

  } /* End of th OpenMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Performs the updates after the first iteration of the Uzawa algo and
 *         stores the divergence by cell
 *
 * \param[in]      relax     scaling of the div operator
 * \param[in]      eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in]      dt_cur    current value of the time step
 * \param[in]      vel_f     velocity DoFs on faces
 * \param[in, out] pr        pressure DoFs (on cells)
 * \param[in, out] div       divergence operator
 * \param[in, out] rhs       rhs for the next iteration
 */
/*----------------------------------------------------------------------------*/

static void
_update_pr_div_rhs(const cs_property_t          *relax,
                   const cs_equation_builder_t  *eqb,
                   cs_real_t                     dt_cur,
                   cs_real_t                    *vel_f,
                   cs_real_t                    *pr,
                   cs_real_t                    *div,
                   cs_real_t                    *rhs)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  time_eval = cs_shared_time_step->t_cur + dt_cur;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;
  const bool  rlx_n_unif = !(cs_property_is_uniform(relax));
  const cs_flag_t *bc_flag = eqb->face_bc->flag;

  /* Resetting new rhs */
  memset(rhs, 0, 3*quant->n_faces*sizeof(cs_real_t));

  /* Get the value of the relaxation parameter for the first cell */
  cs_real_t  rlx = cs_property_get_cell_value(0, time_eval, relax);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Update pressure value: p^{n+1} = p^n - relax div.u^n*/
    if (rlx_n_unif) rlx = cs_property_get_cell_value(c_id, time_eval, relax);

    /* Compute divergence and store it */
    const cs_real_t  div_c = _get_cell_divergence(c_id, quant, c2f, vel_f);

    /* Compute the increment for the pressure */
    const cs_real_t  delta_pc = rlx * div_c;

    div[c_id] = div_c;
    pr[c_id] -= delta_pc;

    for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

      const cs_lnum_t f_id  = c2f->ids[f];
      const cs_lnum_t bf_id = f_id - quant->n_i_faces;

      /* If (hmg or nhmg) Dirichlet boundary face, we should leave the rhs
       * null --> velocity is already known and thus increment is zero */
      if ((bf_id > -1) && /* Boundary face and... */
          cs_cdo_bc_is_dirichlet(bc_flag[bf_id]))
          continue;

      const cs_nvec3_t pfq = cs_quant_set_face_nvec(f_id, quant);

      /* Manually computing the divergence */
      /* No divide by volume since it is then integrated */
      const cs_real_t ifdv = c2f->sgn[f] * pfq.meas * delta_pc;
      cs_real_t  *_rhs = rhs + 3*f_id;

      /* Update the RHS */
#     pragma omp atomic
      _rhs[0] -= ifdv * pfq.unitv[0];
#     pragma omp atomic
      _rhs[1] -= ifdv * pfq.unitv[1];
#     pragma omp atomic
      _rhs[2] -= ifdv * pfq.unitv[2];

    } /* Loop on cell faces */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the residual related to the divergence of the velocity field
 *
 * \param[in]  div    array of divergence inside each cell (already computed)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_compute_residual(int              iter,
                  const cs_real_t *div)
{
  cs_real_t  res = 0.;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Compute the residual related to the divergence of the velocity field */
  cs_real_t  div_res = cs_dot_wxx(quant->n_cells, quant->cell_vol, div);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &div_res);
#endif

  res = sqrt(div_res);
  cs_log_printf(CS_LOG_DEFAULT, "  Uzawa iteration #%4d >> Residual: %8.6e",
                iter, res);
  cs_log_printf(CS_LOG_DEFAULT, "\n");

  return res;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare the incremental resolution of the ULA algo
 *
 * \param[in, out] eq     pointer to a \ref cs_equation_t structure
 * \param[in, out] x_f    unknown for the linear system (face Dofs)
 * \param[in, out] b      rhs for the linear
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_incremental_solve(cs_equation_t      *eq,
                           cs_real_t          *x_f,
                           cs_real_t          *b)
{
   const int  stride = 1;  /* Since the global numbering is adapted in each
                             case (scalar-, vector-valued equations) */

   /* Force the initial guess to zero which should be the expected solution.
    * WARNING: Necessary if Dirichlet algebraic enforcement!
    * Sometimes, it avoids solving the system another time if the solution
    * of the previous iteration is fine enough
    */

  memset(x_f, 0, eq->n_sles_scatter_elts * sizeof(cs_real_t));

  if (cs_glob_n_ranks > 1) { /* Parallel mode */
#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      b[i] = eq->rhs[i];
#else
    memcpy(b, eq->rhs, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    cs_interface_set_sum(eq->rset->ifs,
                         eq->n_sles_scatter_elts, stride, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE,/* type */
                        stride,      /* stride */
                        b,           /* in: size = n_sles_scatter_elts */
                        b);          /* out: size = n_sles_gather_elts */
  }
  else { /* Serial mode *** without periodicity *** */

    /* Nothing to do for the right-hand side */
    assert(eq->n_sles_gather_elts == eq->n_sles_scatter_elts);
    assert(b == eq->rhs);

  }
}

/*----------------------------------------------------------------------------*/
 /*!
 * \brief  Solve the Uzawa linear system (adaptation of cs_equation_solve())
 *
 * \param[in, out]  eq    pointer to the momentum cs_equation_t structure
 * \param[in, out]  x     pointer to temporary velocity on faces
 * \param[in, out]  b     pointer to auxiliary rhs, should be NULL
 */
/*----------------------------------------------------------------------------*/

static int
_uzawa_solve(cs_equation_t   *eq,
             cs_real_t       *x,
             cs_real_t       *b)
{
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);

  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_itsol_t  itsol_info = eq->param->itsol_info;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   eq->matrix,
                                                   CS_HALO_ROTATION_IGNORE,
                                                   itsol_info.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   b,
                                                   x,
                                                   0,      /* aux. size */
                                                   NULL);  /* aux. buffers */

  if (eq->param->sles_verbosity > 0) {

    const cs_lnum_t  size = eq->n_sles_gather_elts;
    const cs_lnum_t  *row_index, *col_id;
    const cs_real_t  *d_val, *x_val;

    cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, &d_val, &x_val);

    cs_gnum_t  nnz = row_index[size];
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&nnz, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  "  <%s/sles_cvg> code  %-d n_iters  %d residual  % -8.4e"
                  " nnz %lu\n",
                  eq->param->name, code, n_iters, residual, nnz);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG
    cs_dbg_dump_linear_system(eq->param->name, size, CS_CDOFB_UZAWA_DBG,
                              x, b,
                              row_index, col_id, x_val, d_val);
#endif
  }

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(eq->rset, CS_REAL_TYPE, 1, /* type and stride */
                         x, x);                     /* in/out */

    cs_range_set_scatter(eq->rset, CS_REAL_TYPE, 1, /* type and stride */
                         b, eq->rhs);               /* in/out */

  }

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  return  n_iters; /* Number of inner iterations for solving the linear system
                      at this step */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_init_common(const cs_cdo_quantities_t     *quant,
                           const cs_cdo_connect_t        *connect,
                           const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /*
    Matrix structure related to the algebraic system for scalar-valued equation
  */
  cs_shared_scal_ms = cs_cdofb_scaleq_matrix_structure();

  /*
    Matrix structure related to the algebraic system for vector-valued equation
  */
  cs_shared_vect_ms = cs_cdofb_vecteq_matrix_structure();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_cdofb_uzawa_t structure
 *
 * \param[in] nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a \ref cs_navsto_uzawa_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_uzawa_init_scheme_context(const cs_navsto_param_t     *nsp,
                                   void                        *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  /* Cast the coupling context (CC) */
  cs_navsto_uzawa_t  *cc = (cs_navsto_uzawa_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_uzawa_t  *sc = NULL;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n",
              __func__);

  BFT_MALLOC(sc, 1, cs_cdofb_uzawa_t);

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");
  sc->divergence = cs_field_by_name("velocity_divergence");

  sc->is_gdscale_uniform = true;
  sc->residual = DBL_MAX;
  sc->last_iter = -1;

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_uzawa_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_uzawa_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t  *)scheme_context;

  if (sc == NULL)
    return sc;

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the velocity values
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_init_velocity(const cs_navsto_param_t     *nsp,
                             void                        *scheme_context)
{
  CS_UNUSED(nsp);
  CS_UNUSED(scheme_context);

  /* Nothing to do. All is already done during the initialization of the
     momentum equation */
  return;
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_init_pressure(const cs_navsto_param_t     *nsp,
                             void                        *scheme_context)
{
  /* Sanity checks */
  assert(nsp != NULL && scheme_context != NULL);

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t *)scheme_context;

  /* Initial conditions for the pressure */
  if (nsp->n_pressure_ic_defs > 0) {

    assert(nsp->pressure_ic_defs != NULL);
    assert(sc != NULL);

    const cs_time_step_t *ts = cs_shared_time_step;
    cs_field_t *pr  = sc->pressure;
    cs_real_t  *values = pr->val;

    const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;
    const cs_flag_t   dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;
    const cs_real_t  t_cur = ts->t_cur;

    for (int def_id = 0; def_id < nsp->n_pressure_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      cs_xdef_t  *def = nsp->pressure_ic_defs[def_id];

      /* Initialize face-based array */
      switch (def->type) {

        /* Evaluating the integrals: the averages will be taken care of at the
         * end when ensuring zero-mean valuedness */
      case CS_XDEF_BY_VALUE:
        cs_evaluate_density_by_value(dof_flag, def, values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          /* Forcing BARY so that it is equivalent to DeRham (JB?)*/
          cs_xdef_set_quadrature(def, CS_QUADRATURE_BARY);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          /* Restoring the original */
          cs_xdef_set_quadrature(def, nsp->qtype);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Incompatible reduction for the field %s.\n"),
                    __func__, pr->name);

        }  /* Switch on possible reduction types */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Incompatible way to initialize the field %s.\n"),
                  __func__, pr->name);
        break;

      }  /* Switch on possible type of definition */

    }  /* Loop on definitions */

  /* We should ensure that the mean of the pressure is zero. Thus we compute
   * it and subtract it from every value.
   * NOTES:
   *  - It could be useful to stored this average somewhere
   *  - The procedure is not optimized (we can avoid setting the average if
   *    it's a value), but it is the only way to allow multiple definitions
   *    and definitions that do not cover all the domain. Moreover, we need
   *    information (e.g. cs_cdo_quant) which we do not know here
   */
  cs_cdofb_navsto_set_zero_mean_pressure(values);

  }  /* Not the default */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         a Uzawa-Lagrangian Augmented approach.
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute(const cs_mesh_t              *mesh,
                       const cs_navsto_param_t      *nsp,
                       void                         *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t0 = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t *)scheme_context;
  cs_navsto_uzawa_t  *cc = (cs_navsto_uzawa_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const size_t  rsize = sizeof(cs_real_t);
  const cs_property_t  *zeta = cc->zeta;
  const cs_lnum_t  n_cells = quant->n_cells, n_faces = quant->n_faces;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];

  cs_real_t  *pr = sc->pressure->val;
  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *vel_f = mom_eq->get_face_values(mom_eqc);
  cs_real_t  *div = sc->divergence->val;

  /* Update pressure and velocity fields */
  /* ----------------------------------- */

  /* Copy current field values to previous values */
  cs_field_current_to_previous(sc->velocity);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(sc->pressure);

  /* If residual normalization is needed */
  /* ----------------------------------- */
  cs_real_t norm_res = 1.;
  if (cs_shared_time_step->nt_cur > 0 ||
      nsp->n_pressure_ic_defs > 0) {
    cs_real_t l2_p = sqrt(cs_dot_wxx(quant->n_cells,quant->cell_vol,pr));

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      cs_parall_sum(1, CS_REAL_TYPE, &l2_p);
#endif
    if (l2_p > 10*mom_eq->param->itsol_info.eps)
      norm_res =  1. / l2_p;
  } /* If needs a normalization */

  /**********  INNER ITERATIONS - START  ***********/
  /* Convergence code:
   *  1 = OK,
   * -1 = max iter,
   * -2 = algo stagnated,
   * -3 = divergence,
   *  2 = default
   */
  short int  cvg_code = 2;
  cs_lnum_t  iter = 1, loc_solv_iter = 0, solv_iter = 0;
  double  res = DBL_MAX;

  /* ---------------------------------------- */
  /* 1) Build the linear system once for all  */
  /* ---------------------------------------- */

  mom_eq->initialize_system(mom_eq->param, mom_eq->builder, mom_eqc,
                            &(mom_eq->matrix),
                            &(mom_eq->rhs));

  _build_system_uzawa(mesh, vel_c, pr, dt_cur, zeta, sc, mom_eq);

  /* Sanity checks (up to now, only scalar field are handled) */
  assert(mom_eq->n_sles_gather_elts <= mom_eq->n_sles_scatter_elts);
  assert(mom_eq->n_sles_gather_elts == cs_matrix_get_n_rows(mom_eq->matrix));

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                " n_sles_gather_elts:  %d\n"
                " n_sles_scatter_elts: %d\n"
                " n_matrix_rows:       %d\n"
                " n_matrix_columns:    %d\n",
                mom_eq->n_sles_gather_elts, mom_eq->n_sles_scatter_elts,
                cs_matrix_get_n_rows(mom_eq->matrix),
                cs_matrix_get_n_columns(mom_eq->matrix));
#endif

  /* Prepare the call to the linear solver:
   *  - x_f is allocated inside (size 3*n_faces and set to the current value
   *    of the field related to mom_eq (i.e. the velocity). x_f = u_{f,k=0}
   *  - Handle parallelism (if // --> b is allocated since it gathers
   *    contribution from ranks sharing faces)
   */
  cs_real_t *x_f = NULL, *b = NULL;
  mom_eq->prepare_solving(mom_eq, &x_f, &b);

  /* First iteration: The Uzawa system is not formulated with an increment */
  solv_iter += _uzawa_solve(mom_eq, x_f, b);

  /* Update the velocity field: u_{c,k=1} and u_{f,k=1} */
  cs_cdofb_vecteq_update_field(x_f, mom_eq->rhs, mom_eq->param,
                               mom_eq->builder, mom_eqc, vel_c);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
  cs_dbg_darray_to_listing("FACE_VELOCITY_k=1", 3*n_faces, x_f, 9);
  cs_dbg_darray_to_listing("CELL_VELOCITY_k=1", 3*n_cells, vel_c, 9);
#endif

  /* Updates after the first resolution:
   *  the divergence: div = B.u_{f,k=1}
   *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
   *  rhs: -zeta Bt.B.u_{f,k=1}
   */
  _update_pr_div_rhs(zeta, mom_eq->builder, dt_cur, x_f, pr, div, mom_eq->rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE_k=1", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV_k=1", n_cells, div, 9);
#endif

  /* Compute residual */
  res = _compute_residual(iter, div) * norm_res;

  /*********** FIRST ITERATION - END ***************/

  if (res > nsp->residual_tolerance) { /* Iterate one more time */

    /* Resetting local cell-defined rhs:
     * In the cell part of the rhs, we have only the source term and/or the
     * term related to the time scheme: both are eliminated in the incremental
     * point of view */
    memset(mom_eqc->rc_tilda, 0, 3*n_cells*rsize);

    /* We now solve the incremental formulation
     * A := Laplacian + \zeta \grad\div
     * du^n_{k+1}] := u^n_{k+1} - u^n_k
     * A du^n_{k+1} = - \grad( p_{k-1}-p_{k-2} )
     *              = - \grad( \zeta \div(u^n_{k-1}) )
     */

    cs_real_t  *delta_vel_c = NULL, *delta_vel_f = x_f;

    BFT_MALLOC(delta_vel_c, 3*n_cells, cs_real_t);
    memset(delta_vel_c, 0, 3*n_cells*rsize);

    while (res > nsp->residual_tolerance && iter < nsp->max_algo_iter) {

      iter++;

      _prepare_incremental_solve(mom_eq, delta_vel_f, b);

      /* Solve */
      solv_iter += (loc_solv_iter = _uzawa_solve(mom_eq, delta_vel_f, b));
      if (loc_solv_iter == 0) {
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n  The inner iterations stagnated. Stopping.\n");
        cvg_code = -2;
        break;
      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      cs_dbg_darray_to_listing("FACE_DELTA_VELOCITY",
                               3*n_faces, delta_vel_f, 9);
#endif
      /* Compute delta_vel_c from the knowledge of delta_vel_f */
      cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                            mom_eqc->rc_tilda,
                                            mom_eqc->acf_tilda,
                                            delta_vel_f,
                                            delta_vel_c);

#     pragma omp parallel if (n_cells > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_faces; i++)
          vel_f[i] += delta_vel_f[i];

#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_cells; i++)
          vel_c[i] += delta_vel_c[i];

      } /* End of the OpenMP region */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      cs_dbg_darray_to_listing("CELL_DELTA_VELOCITY",
                               3*n_cells, delta_vel_c, 9);
#endif

      /* Updates after the first resolution:
       *  the divergence: div = B.u_{f,k=1}
       *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
       *  rhs: -zeta Bt.B.u_{f,k=1}
       */
      _update_pr_div_rhs(zeta, mom_eq->builder, dt_cur,
                         vel_f, pr, div, mom_eq->rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      cs_dbg_darray_to_listing("FACE_VELOCITY", 3*n_faces, vel_f, 9);
      cs_dbg_darray_to_listing("CELL_VELOCITY", 3*n_cells, vel_c, 9);
      cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
      cs_dbg_darray_to_listing("DIVERGENCE", n_cells, div, 9);
#endif

      /* Compute residual */
      res = _compute_residual(iter, div) * norm_res;

      if (res > 1e8) {
        cvg_code = -3;
        break;
      }

    } /* while */

    BFT_FREE(delta_vel_c);

  } /* If more than one iteration */

  /**************  INNER ITERATIONS - END  *************/

  if (res > nsp->residual_tolerance) {
    if (cvg_code == 2)
      cvg_code = -1;
  }
  else
    cvg_code = 1;

  cs_log_printf(CS_LOG_DEFAULT,
                "\n <Uzawa Summary>\n"
                "  Convergence.Code             %-d\n"
                "  Final.Residual               %7.6e\n"
                "  Uzawa.Iterations             %d\n"
                "  Cumulalted.Solver.Iterations %d, mean: %6.1f\n",
                cvg_code, res, iter, solv_iter, (float)solv_iter/iter);

  if (cvg_code < 0) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n ATTENTION: Uzawa algorithm did NOT converge.\n");
    if (cvg_code < -2)
      bft_error(__FILE__, __LINE__, 0, " Uzawa algorithm DIVERGED.\n");
  }

  sc->last_iter = iter;
  sc->residual = res;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
  cs_dbg_darray_to_listing("FINAL_OUTER_SOLUTION", 3*n_faces, vel_f, 9);
#endif

  /* -- Frees - Taken from cs_equation_solve */
  BFT_FREE(x_f);
  if (b != mom_eq->rhs)
    BFT_FREE(b);
  BFT_FREE(mom_eq->rhs);
  cs_sles_free(cs_sles_find_or_add(mom_eq->field_id, NULL));
  cs_matrix_destroy(&(mom_eq->matrix));

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values of the velocity on the faces
 *
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 *
 * \return a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_uzawa_get_face_velocity(void    *scheme_context)
{
  CS_UNUSED(scheme_context);

  return NULL;  /* Not in the scheme context */
}

#undef _dp3
#undef _n3

/*----------------------------------------------------------------------------*/

END_C_DECLS
