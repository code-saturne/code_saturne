/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an Augmented Lagrangian-Uzawa algorithm
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

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
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
#include "cs_log.h"
#include "cs_navsto_coupling.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sles.h"
#include "cs_sles_petsc.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

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

#define CS_CDOFB_UZAWA_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t     *cs_shared_quant;
static const cs_cdo_connect_t        *cs_shared_connect;
static const cs_time_step_t          *cs_shared_time_step;
static const cs_matrix_structure_t   *cs_shared_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a CG
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_amg_block_hook(void     *context,
                Mat       a,
                KSP       ksp)
{
  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  KSPSetType(ksp, KSPFCG);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* Apply modifications to the KSP structure */
  PetscInt  id, n_split;
  KSP  *uvw_subksp, _ksp;
  PC pc, _pc;

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);

  PCFieldSplitSetBlockSize(pc, 3);
  id = 0;
  PCFieldSplitSetFields(pc, "u", 1, &id, &id);
  id = 1;
  PCFieldSplitSetFields(pc, "v", 1, &id, &id);
  id = 2;
  PCFieldSplitSetFields(pc, "w", 1, &id, &id);

  PCSetFromOptions(pc);
  PCSetUp(pc);

  PCFieldSplitGetSubKSP(pc, &n_split, &uvw_subksp);
  assert(n_split == 3);

  for (id = 0; id < 3; id++) {

    _ksp = uvw_subksp[id];
    KSPSetType(_ksp, KSPPREONLY);
    KSPGetPC(_ksp, &_pc);
    PCSetType(_pc, PCHYPRE);
    PCHYPRESetType(_pc, "boomeramg");

    PCSetFromOptions(_pc);
    PCSetUp(_pc);

  }

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(uvw_subksp);
}
#endif

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
                   const cs_real_t               time_eval,
                   const cs_real_t               vel_f[],
                   cs_real_t                     pr[],
                   cs_real_t                     div[],
                   cs_real_t                     rhs[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
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
    const cs_real_t  div_c =
      cs_cdofb_navsto_cell_divergence(c_id, quant, c2f, vel_f);

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

static inline cs_real_t
_compute_residual(int              iter,
                  const cs_real_t *div)
{
  cs_real_t  res = 0.;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Compute the residual related to the divergence of the velocity field */
  cs_real_t  div_res = cs_dot_wxx(quant->n_cells, quant->cell_vol, div);

  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &div_res);

  res = sqrt(div_res);
  cs_log_printf(CS_LOG_DEFAULT, "  Uzawa iteration #%4d >> Residual: %8.6e",
                iter, res);
  cs_log_printf(CS_LOG_DEFAULT, "\n");

  return res;
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
    Matrix structure related to the algebraic system for vector-valued equation
  */
  cs_shared_ms = cs_cdofb_vecteq_matrix_structure();
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
 * \brief  Start setting-up the Navier-Stokes equations when an ALU algorithm
 *         is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_set_sles(const cs_navsto_param_t    *nsp,
                        void                       *context)
{
  cs_navsto_uzawa_t  *nsc = (cs_navsto_uzawa_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  switch (nsp->sles_strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp, field_id);
    break;

#if defined(HAVE_PETSC)
  case CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _amg_block_hook,
                         (void *)mom_eqp);
    break;
#else
  case CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please build a version of Code_Saturne with the PETSc support.",
              __func__, mom_eqp->name);
    break;
#endif /* HAVE_PETSC */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

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
cs_cdofb_uzawa_compute_steady(const cs_mesh_t              *mesh,
                              const cs_navsto_param_t      *nsp,
                              void                         *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t_cmp = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t *)scheme_context;
  cs_navsto_uzawa_t  *cc = (cs_navsto_uzawa_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const size_t  rsize = sizeof(cs_real_t);
  const cs_lnum_t  n_cells = quant->n_cells, n_faces = quant->n_faces;
  const cs_property_t  *zeta = cc->zeta;

  cs_real_t  *pr = sc->pressure->val;
  cs_field_t *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;
  cs_real_t  *vel_f = mom_eq->get_face_values(mom_eqc);
  cs_real_t  *div = sc->divergence->val;

  /* If residual normalization is needed */
  /* ----------------------------------- */
  cs_real_t o_norm_res = 1.;
  if (nsp->n_pressure_ic_defs > 0) {
    cs_real_t l2_p = sqrt(cs_dot_wxx(n_cells, quant->cell_vol, pr));

    if (cs_glob_n_ranks > 1)
      cs_parall_sum(1, CS_REAL_TYPE, &l2_p);

    if (l2_p > 10*mom_eq->param->sles_param.eps)
      o_norm_res =  1. / l2_p;
  } /* If needs a normalization */

  /*---------------------------------------------------------------------------
   *
   *                      BUILD: START
   *
   *--------------------------------------------------------------------------*/

  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  time_eval = t_cur; /* dummy variable */

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces (t_cur is a dummy
     argument) */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (n_cells > CS_THR_MIN) default(none)           \
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

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

    cs_real_t  zeta_c = cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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
      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      const cs_real_t *l_div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
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
      if ( !(sc->is_gdscale_uniform) )
        zeta_c = cs_property_value_in_cell(cm, zeta, time_eval);

      cs_cdofb_navsto_add_grad_div(n_fc, zeta_c*ovol, l_div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      const _Bool has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
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

      /* 4- TIME CONTRIBUTION */
      /* ==================== */
      /* Not applicable */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix before condensation", csys);
#endif

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
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

  } /* End of th OpenMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*---------------------------------------------------------------------------
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

  /* Prepare the call to the linear solver:
   *  - x_f is allocated inside (size 3*n_faces and set to the current value
   *    of the field related to mom_eq (i.e. the velocity). x_f = u_{f,k=0}
   *  - Handle parallelism (if // --> b is allocated since it gathers
   *    contribution from ranks sharing faces)
   */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  solv_iter += cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, vel_f, rhs);

  /* Update field */
  t_upd = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("FACE_VELOCITY_k=1", 3*n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("CELL_VELOCITY_k=1", 3*n_cells, vel_c, 9);
#endif

  /* Updates after the first resolution:
   *  the divergence: div = B.u_{f,k=1}
   *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
   *  rhs: -zeta Bt.B.u_{f,k=1}
   */
  _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("PRESSURE_k=1", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV_k=1", n_cells, div, 9);
#endif

  /* Compute residual */
  res = _compute_residual(iter, div) * o_norm_res;

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

    cs_real_t *delta_vel_c = NULL, *delta_vel_f = NULL; /* Auxiliary buffers */
    BFT_MALLOC(delta_vel_f, 3*n_faces, cs_real_t);
    BFT_MALLOC(delta_vel_c, 3*n_cells, cs_real_t);
    memset(delta_vel_c, 0, 3*n_cells*rsize);
    /* Initialization of delta_vel_f inside the loop */

    while (res > nsp->residual_tolerance && iter < nsp->max_algo_iter) {

      iter++;

      /* Null initial guess */
      memset(delta_vel_f, 0, 3*n_faces*rsize);

      solv_iter += (loc_solv_iter =
        cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, delta_vel_f, rhs));

      t_upd = cs_timer_time();
      /* Compute delta_vel_c from the knowledge of delta_vel_f */
      cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                            mom_eqc->rc_tilda,
                                            mom_eqc->acf_tilda,
                                            delta_vel_f,
                                            delta_vel_c);

      t_tmp = cs_timer_time();
      cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

      if (loc_solv_iter == 0) {
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n  The inner iterations stagnated. Stopping.\n");
        cvg_code = -2;
        break;
      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("FACE_DELTA_VELOCITY",
                               3*n_faces, delta_vel_f, 9);
#endif

#     pragma omp parallel if (3*n_cells > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_faces; i++)
          vel_f[i] += delta_vel_f[i];

#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_cells; i++)
          vel_c[i] += delta_vel_c[i];

      } /* End of the OpenMP region */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("CELL_DELTA_VELOCITY",
                               3*n_cells, delta_vel_c, 9);
#endif

      /* Updates after the first resolution:
       *  the divergence: div = B.u_{f,k=1}
       *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
       *  rhs: -zeta Bt.B.u_{f,k=1}
       */
      _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("FACE_VELOCITY", 3*n_faces, vel_f, 9);
      cs_dbg_darray_to_listing("CELL_VELOCITY", 3*n_cells, vel_c, 9);
      cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
      cs_dbg_darray_to_listing("DIVERGENCE", n_cells, div, 9);
#endif

      /* Compute residual */
      res = _compute_residual(iter, div) * o_norm_res;

      if (res > 1e8) {
        cvg_code = -3;
        break;
      }

    } /* while */

    BFT_FREE(delta_vel_c);
    BFT_FREE(delta_vel_f);

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
                "  Cumulated.Solver.Iterations %d, mean: %6.1f\n",
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
  cs_dbg_fprintf_system(mom_eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_UZAWA_DBG,
                        vel_f, NULL, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmp, &t_tmp);
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach and an Euler time scheme
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_implicit(const cs_mesh_t              *mesh,
                                const cs_navsto_param_t      *nsp,
                                void                         *scheme_context)
{
  cs_timer_t  t_cmp = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t *)scheme_context;
  cs_navsto_uzawa_t  *cc = (cs_navsto_uzawa_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const size_t  rsize = sizeof(cs_real_t);
  const cs_lnum_t  n_cells = quant->n_cells, n_faces = quant->n_faces;
  const cs_property_t  *zeta = cc->zeta;

  cs_real_t  *pr = sc->pressure->val;
  cs_field_t *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;
  cs_real_t  *vel_f = mom_eq->get_face_values(mom_eqc);
  cs_real_t  *div = sc->divergence->val;

  /* If residual normalization is needed */
  /* ----------------------------------- */
  cs_real_t o_norm_res = 1.;
  if (nsp->n_pressure_ic_defs > 0) {
    cs_real_t l2_p = sqrt(cs_dot_wxx(n_cells, quant->cell_vol, pr));

    if (cs_glob_n_ranks > 1)
      cs_parall_sum(1, CS_REAL_TYPE, &l2_p);

    if (l2_p > 10*mom_eq->param->sles_param.eps)
      o_norm_res =  1. / l2_p;

  } /* If needs a normalization */

  /*---------------------------------------------------------------------------
   *
   *                      BUILD: START
   *
   *--------------------------------------------------------------------------*/

  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
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

# pragma omp parallel if (n_cells > CS_THR_MIN) default(none)           \
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

    cs_real_t  zeta_c = cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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
      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      const cs_real_t *l_div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
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
      if ( !(sc->is_gdscale_uniform) )
        zeta_c = cs_property_value_in_cell(cm, zeta, time_eval);

      cs_cdofb_navsto_add_grad_div(n_fc, zeta_c*ovol, l_div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm)
        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   time_eval, 1., /* time, scaling */
                                   cb, mom_eqb, csys);

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

      /* 4- TIME CONTRIBUTION */
      /* ==================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                          or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix before condensation", csys);
#endif

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
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

  } /* End of th OpenMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*---------------------------------------------------------------------------
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

  /* Prepare the call to the linear solver:
   *  - x_f is allocated inside (size 3*n_faces and set to the current value
   *    of the field related to mom_eq (i.e. the velocity). x_f = u_{f,k=0}
   *  - Handle parallelism (if // --> b is allocated since it gathers
   *    contribution from ranks sharing faces)
   */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  solv_iter += cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, vel_f, rhs);

  /* Update field */
  t_upd = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("FACE_VELOCITY_k=1", 3*n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("CELL_VELOCITY_k=1", 3*n_cells, vel_c, 9);
#endif

  /* Updates after the first resolution:
   *  the divergence: div = B.u_{f,k=1}
   *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
   *  rhs: -zeta Bt.B.u_{f,k=1}
   */
  _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("PRESSURE_k=1", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV_k=1", n_cells, div, 9);
#endif

  /* Compute residual */
  res = _compute_residual(iter, div) * o_norm_res;

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

    cs_real_t *delta_vel_c = NULL, *delta_vel_f = NULL; /* Auxiliary buffers */
    BFT_MALLOC(delta_vel_f, 3*n_faces, cs_real_t);
    BFT_MALLOC(delta_vel_c, 3*n_cells, cs_real_t);
    memset(delta_vel_c, 0, 3*n_cells*rsize);
    /* Initialization of delta_vel_f inside the loop */

    while (res > nsp->residual_tolerance && iter < nsp->max_algo_iter) {

      iter++;

      /* Null initial guess */
      memset(delta_vel_f, 0, 3*n_faces*rsize);

      solv_iter += (loc_solv_iter =
        cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, delta_vel_f, rhs));

      t_upd = cs_timer_time();
      /* Compute delta_vel_c from the knowledge of delta_vel_f */
      cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                            mom_eqc->rc_tilda,
                                            mom_eqc->acf_tilda,
                                            delta_vel_f,
                                            delta_vel_c);

      t_tmp = cs_timer_time();
      cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

      if (loc_solv_iter == 0) {
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n  The inner iterations stagnated. Stopping.\n");
        cvg_code = -2;
        break;
      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("FACE_DELTA_VELOCITY",
                               3*n_faces, delta_vel_f, 9);
#endif

#     pragma omp parallel if (3*n_cells > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_faces; i++)
          vel_f[i] += delta_vel_f[i];

#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_cells; i++)
          vel_c[i] += delta_vel_c[i];

      } /* End of the OpenMP region */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("CELL_DELTA_VELOCITY",
                               3*n_cells, delta_vel_c, 9);
#endif

      /* Updates after the first resolution:
       *  the divergence: div = B.u_{f,k=1}
       *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
       *  rhs: -zeta Bt.B.u_{f,k=1}
       */
      _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("FACE_VELOCITY", 3*n_faces, vel_f, 9);
      cs_dbg_darray_to_listing("CELL_VELOCITY", 3*n_cells, vel_c, 9);
      cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
      cs_dbg_darray_to_listing("DIVERGENCE", n_cells, div, 9);
#endif

      /* Compute residual */
      res = _compute_residual(iter, div) * o_norm_res;

      if (res > 1e8) {
        cvg_code = -3;
        break;
      }

    } /* while */

    BFT_FREE(delta_vel_c);
    BFT_FREE(delta_vel_f);

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
                "  Cumulated.Solver.Iterations %d, mean: %6.1f\n",
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
  cs_dbg_fprintf_system(mom_eqp->name, ts->nt_cur, CS_CDOFB_UZAWA_DBG,
                        vel_f, NULL, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmp, &t_tmp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Uzawa-Lagrangian Augmented approach and a theta time scheme
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_uzawa_compute_theta(const cs_mesh_t              *mesh,
                             const cs_navsto_param_t      *nsp,
                             void                         *scheme_context)
{
  cs_timer_t  t_cmp = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_uzawa_t  *sc = (cs_cdofb_uzawa_t *)scheme_context;
  cs_navsto_uzawa_t  *cc = (cs_navsto_uzawa_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_cells = quant->n_cells, n_faces = quant->n_faces;
  const cs_property_t  *zeta = cc->zeta;
  const size_t  rsize = sizeof(cs_real_t);

  cs_real_t  *pr = sc->pressure->val;
  cs_field_t *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;
  cs_real_t  *vel_f = mom_eq->get_face_values(mom_eqc);
  cs_real_t  *div = sc->divergence->val;

  /* If residual normalization is needed */
  /* ----------------------------------- */
  cs_real_t o_norm_res = 1.;
  if (nsp->n_pressure_ic_defs > 0) {
    cs_real_t l2_p = sqrt(cs_dot_wxx(n_cells, quant->cell_vol, pr));

    if (cs_glob_n_ranks > 1)
      cs_parall_sum(1, CS_REAL_TYPE, &l2_p);

    if (l2_p > 10*mom_eq->param->sles_param.eps)
      o_norm_res =  1. / l2_p;

  } /* If needs a normalization */

  /*---------------------------------------------------------------------------
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

# pragma omp parallel if (n_cells > CS_THR_MIN) default(none)           \
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

    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

    cs_real_t  zeta_c = cs_property_get_cell_value(0, time_eval, zeta);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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
      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      const cs_real_t *l_div = cb->aux->val, ovol = 1. / cm->vol_c;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT,
                          "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, l_div[3*f]*ovol, l_div[3*f+1]*ovol,
                          l_div[3*f+2]*ovol);
        } /* Critical section */
      }
#endif

      /* Update the property */
      if ( !(sc->is_gdscale_uniform) )
        zeta_c = cs_property_value_in_cell(cm, zeta, time_eval);

      cs_cdofb_navsto_add_grad_div(n_fc, zeta_c*ovol, l_div, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion and grad-div (lhs)",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);

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

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix before condensation", csys);
#endif

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 1
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

  } /* End of th OpenMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*---------------------------------------------------------------------------
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

  /* Prepare the call to the linear solver:
   *  - x_f is allocated inside (size 3*n_faces and set to the current value
   *    of the field related to mom_eq (i.e. the velocity). x_f = u_{f,k=0}
   *  - Handle parallelism (if // --> b is allocated since it gathers
   *    contribution from ranks sharing faces)
   */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  solv_iter += cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, vel_f, rhs);

  /* Update field */
  t_upd = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("FACE_VELOCITY_k=1", 3*n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("CELL_VELOCITY_k=1", 3*n_cells, vel_c, 9);
#endif

  /* Updates after the first resolution:
   *  the divergence: div = B.u_{f,k=1}
   *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
   *  rhs: -zeta Bt.B.u_{f,k=1}
   */
  _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
  cs_dbg_darray_to_listing("PRESSURE_k=1", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV_k=1", n_cells, div, 9);
#endif

  /* Compute residual */
  res = _compute_residual(iter, div) * o_norm_res;

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

    cs_real_t *delta_vel_c = NULL, *delta_vel_f = NULL; /* Auxiliary buffers */
    BFT_MALLOC(delta_vel_f, 3*n_faces, cs_real_t);
    BFT_MALLOC(delta_vel_c, 3*n_cells, cs_real_t);
    memset(delta_vel_c, 0, 3*n_cells*rsize);
    /* Initialization of delta_vel_f inside the loop */

    while (res > nsp->residual_tolerance && iter < nsp->max_algo_iter) {

      iter++;

      /* Null initial guess */
      memset(delta_vel_f, 0, 3*n_faces*rsize);

      solv_iter += (loc_solv_iter =
        cs_cdofb_vecteq_solve_system(sles, matrix, mom_eqp, delta_vel_f, rhs));

      t_upd = cs_timer_time();
      /* Compute delta_vel_c from the knowledge of delta_vel_f */
      cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                            mom_eqc->rc_tilda,
                                            mom_eqc->acf_tilda,
                                            delta_vel_f,
                                            delta_vel_c);

      t_tmp = cs_timer_time();
      cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

      if (loc_solv_iter == 0) {
        cs_log_printf(CS_LOG_DEFAULT,
                      "\n  The inner iterations stagnated. Stopping.\n");
        cvg_code = -2;
        break;
      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("FACE_DELTA_VELOCITY",
                               3*n_faces, delta_vel_f, 9);
#endif

#     pragma omp parallel if (3*n_cells > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_faces; i++)
          vel_f[i] += delta_vel_f[i];

#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < 3*n_cells; i++)
          vel_c[i] += delta_vel_c[i];

      } /* End of the OpenMP region */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 3
      cs_dbg_darray_to_listing("CELL_DELTA_VELOCITY",
                               3*n_cells, delta_vel_c, 9);
#endif

      /* Updates after the first resolution:
       *  the divergence: div = B.u_{f,k=1}
       *  the pressure field: pr_{k=1} -= zeta * div(u_{f,k=1}
       *  rhs: -zeta Bt.B.u_{f,k=1}
       */
      _update_pr_div_rhs(zeta, mom_eqb, time_eval, vel_f, pr, div, rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_UZAWA_DBG > 2
      cs_dbg_darray_to_listing("FACE_VELOCITY", 3*n_faces, vel_f, 9);
      cs_dbg_darray_to_listing("CELL_VELOCITY", 3*n_cells, vel_c, 9);
      cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
      cs_dbg_darray_to_listing("DIVERGENCE", n_cells, div, 9);
#endif

      /* Compute residual */
      res = _compute_residual(iter, div) * o_norm_res;

      if (res > 1e8) {
        cvg_code = -3;
        break;
      }

    } /* while */

    BFT_FREE(delta_vel_c);
    BFT_FREE(delta_vel_f);

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
                "  Cumulated.Solver.Iterations %d, mean: %6.1f\n",
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
  cs_dbg_fprintf_system(mom_eqp->name, ts->nt_cur, CS_CDOFB_UZAWA_DBG,
                        vel_f, NULL, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmp, &t_tmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
