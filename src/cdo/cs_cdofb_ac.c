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
static const cs_matrix_structure_t  *cs_shared_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Performs the updates after the resolution of the momentum equation
 *         of the AC algo, update the pressure, and compute the divergence
 *
 * \param[in]      zeta      scaling of the div operator
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in]      eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in]      time_eval time at which the property/functions are evaluated
 * \param[in]      dt_cur    current value of the time step
 * \param[in]      vel_f     velocity DoFs on faces
 * \param[in, out] pr        pressure DoFs (on cells)
 * \param[in, out] div       divergence operator
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_pr_div(const cs_property_t          *zeta,
               const cs_equation_param_t    *eqp,
               const cs_equation_builder_t  *eqb,
               const cs_real_t               time_eval,
               const cs_real_t               dt_cur,
               const cs_real_t               vel_f[],
               cs_real_t                     pr[],
               cs_real_t                     div[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
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
    Matrix structure related to the algebraic system for vector-valued equation
  */
  cs_shared_ms = cs_cdofb_vecteq_matrix_structure();
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
  assert(nsp != NULL);

  /* Initial conditions for the pressure */
  if (nsp->n_pressure_ic_defs == 0)
    return; /* Nothing to do */

  assert(scheme_context != NULL && nsp->pressure_ic_defs != NULL);

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_flag_t   dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;

  cs_field_t *pr  = sc->pressure;
  cs_real_t  *values = pr->val;

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
      {
        const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;

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

      }
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
   *
   * NOTES:
   *  - It could be useful to stored this average somewhere
   *  - The procedure is not optimized (we can avoid setting the average if
   *    it's a value), but it is the only way to allow multiple definitions
   *    and definitions that do not cover all the domain. Moreover, we need
   *    information (e.g. cs_cdo_quantities_t) which we do not know here
   */
  cs_cdofb_navsto_set_zero_mean_pressure(values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and an Euler time scheme
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute_implicit(const cs_mesh_t              *mesh,
                             const cs_navsto_param_t      *nsp,
                             void                         *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t_cmpt = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;
  cs_navsto_ac_t *cc = (cs_navsto_ac_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_property_t  *zeta = cc->zeta;

  cs_real_t  *pr = sc->pressure->val;
  cs_field_t  *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;
  cs_real_t  *div = sc->divergence->val;

  /*----------------------------------------------------------------------------
   *
   *                      BUILD: START
   *
   *--------------------------------------------------------------------------*/

  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;

  /* Sanity checks */
  assert(cs_equation_param_has_time(mom_eqp) == true);
  assert(mom_eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces.
     Evaluation should be performed at t_cur + dt_cur */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix,        \
         mav, rs, dir_values, zeta, vel_c, pr, sc)
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
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

    cs_real_t  o_zeta_c  = 1. / cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, mom_eqb),
                         connect, quant, cm);

      /* Starts from the stationary Stokes problem where
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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, mom_eqp, mom_eqb, mom_eqc,
                                       dir_values, vel_c, time_eval,
                                       csys, cb);

      const short int  n_fc = cm->n_fc, f_dofs = 3*n_fc;

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */
      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqb, mom_eqc,
                                cm, fm, csys, cb);

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */
      cs_navsto_get_divergence_vect(cm, cb->aux->val);

      const cs_real_t *l_div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, l_div[3*f]*ovol, l_div[3*f+1]*ovol,
                          l_div[3*f+2]*ovol);
        } /* Critical section */
      }
#endif

      /* Update the property */
      if ( !(sc->is_zeta_uniform) )
        o_zeta_c = 1. / cs_property_value_in_cell(cm, zeta, time_eval);

      cs_navsto_add_grad_div(n_fc,
                             cb->tpty_val * dt_cur * o_zeta_c * ovol,
                             l_div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   time_eval, 1., /* time, scaling */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS */
      /* =========================== */

      /* Apply the operator gradient to the pressure field and add it to the
         rhs */
      cs_sdm_add_scalvect(f_dofs, pr[c_id], l_div, csys->rhs);

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      cs_cdofb_vecteq_apply_bc_partly(time_eval, mom_eqp, mom_eqc,
                                      cm, fm, csys, cb);

      /* 4- UNSTEADY TERM + TIME SCHEME
       * ============================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                      or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get the cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_fc, n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_fc + k] += ptyc * csys->val_n[3*n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: Only diagonal time treatment "
            "available so far.\n", __func__);

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda,
                                       mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      cs_cdofb_vecteq_apply_remaining_bc(mom_eqp, mom_eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm,
                               mav, rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*----------------------------------------------------------------------------
   *
   *                      BUILD: END
   *
   *--------------------------------------------------------------------------*/

  /* Copy current field values to previous values */
  cs_timer_t t_upd = cs_timer_time();

  cs_field_current_to_previous(vel_fld);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  /* Now solve the system */
  cs_real_t *vel_f = mom_eqc->face_values;
  cs_sles_t *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, vel_f, rhs);

  /* Update pressure, velocity and divergence fields */
  t_upd = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  /* Updates after the resolution:
   *  the divergence: div = B.u_f
   *  the pressure field: pr -= dt / zeta * div(u_f)
   */
  _update_pr_div(zeta, mom_eqp, mom_eqb, time_eval, dt_cur, vel_f, pr, div);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif

  /* Frees */
  cs_sles_free(sles);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and a theta time scheme
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute_theta(const cs_mesh_t              *mesh,
                          const cs_navsto_param_t      *nsp,
                          void                         *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t_cmpt = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;
  cs_navsto_ac_t *cc = (cs_navsto_ac_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  const cs_property_t  *zeta = cc->zeta;

  cs_real_t  *pr = sc->pressure->val;
  cs_field_t *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;
  cs_real_t  *div = sc->divergence->val;

  /*----------------------------------------------------------------------------
   *
   *                      BUILD: START
   *
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;
  const double  tcoef = 1 - mom_eqp->theta;

  /* Sanity checks */
  assert(cs_equation_param_has_time(mom_eqp) == true);
  assert(mom_eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         mom_eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);

  cs_timer_t  t_bld = cs_timer_time();

  /* Detect the first call (in this case, we compute the initial source term)*/
  _Bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix,        \
         mav, rs, dir_values, zeta, vel_c, pr, sc, compute_initial_source)
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
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

    cs_real_t  o_zeta_c  = 1. / cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, mom_eqb),
                         connect, quant, cm);

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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, mom_eqp, mom_eqb, mom_eqc,
                                       dir_values, vel_c, time_eval,
                                       csys, cb);

      const short int  n_fc = cm->n_fc, f_dofs = 3*n_fc;

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */
      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqb, mom_eqc,
                                cm, fm, csys, cb);

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */
      cs_navsto_get_divergence_vect(cm, cb->aux->val);

      const cs_real_t *l_div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT,
                          "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, l_div[3*f]*ovol, l_div[3*f+1]*ovol,
                          l_div[3*f+2]*ovol);
        } /* critical section */
      }
#endif

      /* Update the property */
      if ( !(sc->is_zeta_uniform) )
        o_zeta_c = 1. / cs_property_value_in_cell(cm, zeta, time_eval);

      cs_navsto_add_grad_div(n_fc,
                             cb->tpty_val * dt_cur * o_zeta_c * ovol,
                             l_div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      const _Bool has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */
        if (compute_initial_source) { /* First time step */

          cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     cb, mom_eqb, csys);

        }
        else { /* Add the contribution of the previous time step */

          for (short int k = 0; k < 3; k++)
            csys->rhs[3*n_fc + k] += tcoef * mom_eqc->source_terms[3*c_id + k];

        }

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   /* time,      scaling */
                                   t_cur+dt_cur, mom_eqp->theta,
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS */
      /* =========================== */

      /* Apply the operator gradient to the pressure field and add it to the
         rhs */
      cs_sdm_add_scalvect(f_dofs, pr[c_id], l_div, csys->rhs);

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      cs_cdofb_vecteq_apply_bc_partly(time_eval, mom_eqp, mom_eqc,
                                      cm, fm, csys, cb);

      /* 4- UNSTEADY TERM + TIME SCHEME
       * ============================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_block_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++) /* n_dofs = n_vc */
        csys->rhs[i] -= tcoef * adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */
      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= mom_eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs mass_mat * p_n */
      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        const double  ptyc = cb->tpty_val * cm->vol_c / dt_cur;

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_fc, n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_fc + k] += ptyc * csys->val_n[3*n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: Only diagonal time treatment "
            "available so far.\n", __func__);

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda,
                                       mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      cs_cdofb_vecteq_apply_remaining_bc(mom_eqp, mom_eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm,
                               mav, rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*----------------------------------------------------------------------------
   *
   *                      BUILD: END
   *
   *--------------------------------------------------------------------------*/

  /* Copy current field values to previous values */
  cs_timer_t t_upd = cs_timer_time();

  cs_field_current_to_previous(vel_fld);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  /* Now solve the system */
  cs_real_t *vel_f = mom_eqc->face_values;
  cs_sles_t *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp,
                               vel_f, rhs);

  /* Update pressure, velocity and divergence fields */
  t_upd = cs_timer_time();
  /* ----------------------------------------------- */

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  /* Updates after the resolution:
   *  the divergence: div = B.u_f
   *  the pressure field: pr -= dt / zeta * div(u_f)
   */
  _update_pr_div(zeta, mom_eqp, mom_eqb, time_eval, dt_cur, vel_f, pr, div);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif

  /* -- Frees -- */
  cs_sles_free(sles);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
