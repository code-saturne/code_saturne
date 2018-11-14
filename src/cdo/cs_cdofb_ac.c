/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an artificial compressibility algorithm
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

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

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
#include "cs_sles.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_navsto_utilities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_ac.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_ac.c
 *
 * \brief Build an algebraic CDO face-based system for the Navier-Stokes
 * equations and solved it with an artificial compressibility algorithm
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_ac_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         vector-valued unknowns
 */

typedef struct {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_ac_t (owned by \ref cs_navsto_system_t)
   *  containing the settings related to an artificial compressibility (AC)
   *  algorithm or vector penalty projection (VPP)
   */

  cs_navsto_ac_t   *coupling_context;

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

  /*! \var is_zeta_uniform
   *  Bool telling if the auxiliary parameter zeta is uniform. Not always
   *  necessary: zeta is typically used in Artificial Compressibility algos
   */

  bool is_zeta_uniform;

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

} cs_cdofb_ac_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_AC_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_scal_ms;
static const cs_matrix_structure_t  *cs_shared_vect_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a Stokes equation with a CDO
 *         face-based scheme solved with Artificial Compressibility algo.
 *         One works cellwise and then proceed to the assembly
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
_build_system_ac(const cs_mesh_t       *mesh,
                 const cs_real_t       *vel_vals,
                 const cs_real_t       *pr_vals,
                 double                 dt_cur,
                 const cs_property_t   *zeta,
                 cs_cdofb_ac_t         *sc,
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

    cs_real_t  o_zeta_c  = 1. / cs_property_get_cell_value(0, time_eval, zeta);

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
      cs_navsto_get_divergence_vect(cm, cb->aux->val);

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
      if ( !(sc->is_zeta_uniform) )
        o_zeta_c = 1. / cs_property_value_in_cell(cm, zeta, time_eval);

      cs_navsto_add_grad_div(n_fc,
                             cb->tpty_val * dt_cur * o_zeta_c * ovol,
                             div, csys->mat);

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
 * \brief  Performs the updates after the resolution of the momentum equation
 *         of the AC algo, update the pressure, and compute the divergence
 *
 * \param[in]      zeta      scaling of the div operator
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in]      dt_cur    current value of the time step
 * \param[in]      vel_f     velocity DoFs on faces
 * \param[in, out] pr        pressure DoFs (on cells)
 * \param[in, out] div       divergence operator
 */
/*----------------------------------------------------------------------------*/

static void
_update_pr_div(const cs_property_t          *zeta,
               const cs_equation_param_t    *eqp,
               const cs_equation_builder_t  *eqb,
               cs_real_t                     dt_cur,
               const cs_real_t               vel_f[],
               cs_real_t                    *pr,
               cs_real_t                    *div)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  time_eval = cs_shared_time_step->t_cur + dt_cur;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;
  const bool  zt_n_unif = !(cs_property_is_uniform(zeta));

  cs_real_t  o_zeta_c  = 1. / cs_property_get_cell_value(0, time_eval, zeta),
             t_pty     = cs_property_get_cell_value(0, time_eval,
                                                    eqp->time_property);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN) CS_CDO_OMP_SCHEDULE
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (zt_n_unif)
      o_zeta_c = 1. /cs_property_get_cell_value(c_id, time_eval, zeta);
    if (!eqb->time_pty_uniform)
      t_pty = cs_property_get_cell_value(c_id, time_eval, eqp->time_property);

    /* Compute divergence and store it */
    const cs_real_t  div_c =
      cs_navsto_get_cell_divergence(c_id, quant, c2f, vel_f);

    /* Compute the increment for the pressure */
    const cs_real_t  delta_pc = t_pty * dt_cur * o_zeta_c * div_c;

    div[c_id] = div_c;
    pr[c_id] -= delta_pc;

  } /* Loop on cells */
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
cs_cdofb_ac_init_common(const cs_cdo_quantities_t     *quant,
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
 * \brief  Initialize a \ref cs_cdofb_ac_t structure
 *
 * \param[in] nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a \ref cs_navsto_ac_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_ac_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_init_scheme_context(const cs_navsto_param_t     *nsp,
                                void                        *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  /* Cast the coupling context (CC) */
  cs_navsto_ac_t  *cc = (cs_navsto_ac_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_ac_t  *sc = NULL;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n",
              __func__);

  BFT_MALLOC(sc, 1, cs_cdofb_ac_t);

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");
  sc->divergence = cs_field_by_name("velocity_divergence");

  /* Only one vector equation */
  sc->is_zeta_uniform = true;//cs_property_is_uniform(cc->zeta);

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_ac_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t  *)scheme_context;

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
cs_cdofb_ac_init_velocity(const cs_navsto_param_t     *nsp,
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
cs_cdofb_ac_init_pressure(const cs_navsto_param_t     *nsp,
                          void                        *scheme_context)
{
  /* Sanity checks */
  assert(nsp != NULL && scheme_context != NULL);

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;

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
 *         a Ac-Lagrangian Augmented approach.
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute(const cs_mesh_t              *mesh,
                    const cs_navsto_param_t      *nsp,
                    void                         *scheme_context)
{
  CS_UNUSED(dt_cur);
  CS_UNUSED(nsp);

  cs_timer_t  t0 = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;
  cs_navsto_ac_t *cc = (cs_navsto_ac_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;

  const cs_property_t  *zeta = cc->zeta;

  cs_real_t  *pr = sc->pressure->val;
  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *div = sc->divergence->val;

  /* Update pressure and velocity fields */
  /* ----------------------------------- */

  /* Copy current field values to previous values */
  cs_field_current_to_previous(sc->velocity);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(sc->pressure);

  /* -------------------------- */
  /* 1) Build the linear system */
  /* -------------------------- */

  mom_eq->initialize_system(mom_eq->param, mom_eq->builder, mom_eqc,
                            &(mom_eq->matrix),
                            &(mom_eq->rhs));

  _build_system_ac(mesh, vel_c, pr, dt_cur, zeta, sc, mom_eq);

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

  cs_navsto_solve(mom_eq, x_f, b);

  /* Update the velocity field: u_{c} and u_{f} */
  cs_cdofb_vecteq_update_field(x_f, mom_eq->rhs, mom_eq->param,
                               mom_eq->builder, mom_eqc, vel_c);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
  cs_dbg_darray_to_listing("FACE_VELOCITY", 3*n_faces, x_f, 9);
  cs_dbg_darray_to_listing("CELL_VELOCITY", 3*n_cells, vel_c, 9);
#endif

  /* Updates after the first resolution:
   *  the divergence: div = B.u_f
   *  the pressure field: pr -= zeta * div(u_f)
   */
  _update_pr_div(zeta, mom_eq->param, mom_eq->builder, dt_cur, x_f, pr, div);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", n_cells, div, 9);
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

END_C_DECLS
