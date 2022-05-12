/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an artificial compressibility algorithm
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
#include "cs_cdo_blas.h"
#include "cs_cdo_solve.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_navsto.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_iter_algo.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_sles.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

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
 *         the Navier-Stokes equations with an Artificial compressibility
 *         algorithm
 */

typedef struct {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_ac_t (owned by \ref cs_navsto_system_t)
   *  containing the settings related to an artificial compressibility (AC)
   *  algorithm
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

  cs_field_t       *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t       *pressure;

  /*! \var divergence
   *  Pointer to \ref cs_real_t containing the values of the divergence on the
   *  cells
   */

  cs_field_t       *divergence;

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

  bool              is_zeta_uniform;

  /*!
   * @}
   * @name Advection quantities
   * Members related to the advection
   * @{
   *
   *  \var adv_field
   *  Pointer to the cs_adv_field_t related to the Navier-Stokes eqs (Shared)
   */

  cs_adv_field_t           *adv_field;

  /*! \var mass_flux_array
   *  Current values of the mass flux at primal faces (Shared)
   */

  cs_real_t                *mass_flux_array;

  /*! \var mass_flux_array_pre
   *  Previous values of the mass flux at primal faces (Shared)
   */

  cs_real_t                *mass_flux_array_pre;

  /*!
   * @}
   * @name Boundary conditions (BC) management
   * Functions and elements used for enforcing the BCs
   * @{
   *
   *  \var bf_type
   *  Array of boundary type for each boundary face. (Shared)
   */

  const cs_boundary_type_t       *bf_type;

  /*!
   * \var pressure_bc
   * Structure storing the metadata after processing the user-defined boundary
   * conditions related to the pressure field
   */

  cs_cdo_bc_face_t               *pressure_bc;

  /*! \var apply_fixed_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (no slip boundary)
   *
   *  \var apply_sliding_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (a tangential velocity is specified at the wall)
   *
   *  \var apply_velocity_inlet
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  boundary with a fixed velocity at the inlet
   *
   *  \var apply_symmetry
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  symmetry boundary
   */

  cs_cdo_apply_boundary_t        *apply_fixed_wall;
  cs_cdo_apply_boundary_t        *apply_sliding_wall;
  cs_cdo_apply_boundary_t        *apply_velocity_inlet;
  cs_cdo_apply_boundary_t        *apply_symmetry;

  /*!
   * \var add_gravity_term
   * Compute and add the source term related to the gravity vector
   * This can be the Boussinesq term or the hydrostatic term (rho*g)
   */

  cs_cdofb_navsto_source_t      *add_gravity_term;

  /*!
   * @}
   * @name Convergence monitoring
   * Structure used to drive the convergence of high-level iterative algorithms
   * @{
   */

  /*!
   * \var nl_algo
   * Structure driving the convergence of the non-linear algorithm
   */

  cs_iter_algo_t                *nl_algo;

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

  cs_timer_counter_t             timer;

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

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the Navier-Stokes related variables-fields
 *         to previous values
 *         Case of the artificial compressibility coupling algorithm.
 *
 * \param[in, out] sc    scheme context (\ref cs_cdofb_ac_t struct.)
 * \param[in, out] cc    coupling context (\ref cs_navsto_ac_t struct.)
 */
/*----------------------------------------------------------------------------*/

static inline void
_ac_fields_to_previous(cs_cdofb_ac_t        *sc,
                       cs_navsto_ac_t       *cc)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_field_current_to_previous(sc->velocity);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);

  /* Face velocity arrays */
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)cc->momentum->scheme_context;

  if (eqc->face_values_pre != NULL)
    memcpy(eqc->face_values_pre, eqc->face_values,
           3 * quant->n_faces * sizeof(cs_real_t));

  /* Mass flux arrays */
  memcpy(sc->mass_flux_array_pre, sc->mass_flux_array,
         quant->n_faces * sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence of the velocity in each cell
 *
 * \param[in]      vel_f     velocity DoFs on faces
 * \param[in, out] div       divergence of the given velocity in each cell
 */
/*----------------------------------------------------------------------------*/

static void
_ac_compute_div(const cs_real_t               vel_f[],
                cs_real_t                     div[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Update the divergence of the velocity field */

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    div[c_id] = cs_cdofb_navsto_cell_divergence(c_id,
                                                quant,
                                                cs_shared_connect->c2f,
                                                vel_f);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Performs the update of the pressure after the resolution of the
 *         momentum equation
 *
 * \param[in]      t_eval   time at which the property/functions are evaluated
 * \param[in]      dt_cur   current value of the time step
 * \param[in]      zeta     scaling of the div operator
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in]      eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in]      div_c    mean velocity divergence in each cell
 * \param[in, out] pr_c     pressure DoFs (on cells)
 */
/*----------------------------------------------------------------------------*/

static inline void
_ac_update_pr(const cs_real_t               t_eval,
              const cs_real_t               dt_cur,
              const cs_property_t          *zeta,
              const cs_equation_param_t    *eqp,
              const cs_equation_builder_t  *eqb,
              const cs_real_t               div_c[],
              cs_real_t                     pr_c[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  if (cs_property_is_uniform(zeta)) {

    const cs_real_t o_zeta_c = 1./cs_property_get_cell_value(0, t_eval, zeta);

      if (eqb->time_pty_uniform) {

        /* Best case scenario: everything is uniform */
        const cs_real_t t_pty = cs_property_get_cell_value(0, t_eval,
                                                           eqp->time_property);
        const cs_real_t  coef_mult = t_pty * dt_cur * o_zeta_c;

#       pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
          pr_c[c_id] -= coef_mult * div_c[c_id];

      }
      else { /* time_pty non uniform */

        const cs_real_t  coef_mult = dt_cur * o_zeta_c;

#       pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
          cs_real_t  t_pty = cs_property_get_cell_value(c_id, t_eval,
                                                        eqp->time_property);
          pr_c[c_id] -= coef_mult * t_pty * div_c[c_id];
        }

      }

  }
  else { /* zeta not uniform */

    if (eqb->time_pty_uniform) {

      const cs_real_t coef_mult =
        dt_cur * cs_property_get_cell_value(0, t_eval, eqp->time_property);

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_real_t  o_zeta_c = 1./cs_property_get_cell_value(c_id, t_eval, zeta);
        pr_c[c_id] -= coef_mult * o_zeta_c * div_c[c_id];
      }

    }
    else { /* Neither zeta nor time_pty is uniform */

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_real_t  o_zeta_c = 1./cs_property_get_cell_value(c_id, t_eval, zeta);
        cs_real_t  t_pty = cs_property_get_cell_value(c_id, t_eval,
                                                      eqp->time_property);
        pr_c[c_id] -= t_pty * dt_cur * o_zeta_c * div_c[c_id];
      }

    } /* time_pty non uniform */

  } /* zeta non uniform */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done before the static condensation
 *          Case of CDO Fb schemes with AC coupling
 *
 * \param[in]      sc          pointer to a cs_cdofb_ac_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary faces
 * \param[in]      diff_pty    pointer to \ref cs_property_data_t for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_ac_apply_bc_partly(const cs_cdofb_ac_t           *sc,
                    const cs_equation_param_t     *eqp,
                    const cs_cell_mesh_t          *cm,
                    const cs_boundary_type_t       bf_type[],
                    const cs_property_data_t      *diff_pty,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_diffusion(eqp) == true);

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the velocity-block and the right-hand side (part related to the
     * momentum equation). */

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann)
      for (short int i  = 0; i < 3*cm->n_fc; i++)
        csys->rhs[i] -= csys->neu_values[i];

    const cs_real_t  pc = sc->pressure->val[cm->c_id];

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];
      const cs_quant_t  pfq = cm->face[f];
      const cs_real_t  f_prs = pfq.meas * pc;
      cs_real_t  *f_rhs = csys->rhs + 3*f;

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_velocity_inlet(f, eqp, cm, diff_pty, cb, csys);
          f_rhs[0] -= f_prs * pfq.unitv[0];
          f_rhs[1] -= f_prs * pfq.unitv[1];
          f_rhs[2] -= f_prs * pfq.unitv[2];
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, diff_pty, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, diff_pty, cb, csys);
          f_rhs[0] -= f_prs * pfq.unitv[0];
          f_rhs[1] -= f_prs * pfq.unitv[1];
          f_rhs[2] -= f_prs * pfq.unitv[2];
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {
        /* Always weakly enforce the symmetric constraint on the
           velocity-block */
        sc->apply_symmetry(f, eqp, cm, diff_pty, cb, csys);
        f_rhs[0] -= f_prs * pfq.unitv[0];
        f_rhs[1] -= f_prs * pfq.unitv[1];
        f_rhs[2] -= f_prs * pfq.unitv[2];
      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the boundary conditions to the local system when this should be
 *         done after the static condensation. Apply the internal enforcement
 *         if needed.
 *         Case of CDO-Fb schemes with AC coupling
 *
 * \param[in]      sc          pointer to a cs_cdofb_ac_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary faces
 * \param[in]      diff_pty    pointer to \ref cs_property_data_t for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_ac_apply_remaining_bc(const cs_cdofb_ac_t           *sc,
                       const cs_equation_param_t     *eqp,
                       const cs_equation_builder_t   *eqb,
                       const cs_cell_mesh_t          *cm,
                       const cs_boundary_type_t      *bf_type,
                       const cs_property_data_t      *diff_pty,
                       cs_cell_sys_t                 *csys,
                       cs_cell_builder_t             *cb)
{
  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the divergence operator and the right-hand side related to the
     * mass equation.
     * Enforcement of Dirichlet BC in a stronger way if this is the choice
     */

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */

      const short int  f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */

        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          sc->apply_velocity_inlet(f, eqp, cm, diff_pty, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */

        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, diff_pty, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, diff_pty, cb, csys);
        }
      }

#if 0
      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {
        /* Weak-enforcement for the velocity-block (cf. _ac_apply_bc_partly) */
      }
#endif

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */

  } /* Boundary cell */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

    cs_equation_builder_enforce_block_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of an implicit Euler time scheme
 *
 * \param[in]      nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_pre   velocity face DoFs of the previous time step
 * \param[in]      vel_c_pre   velocity cell DoFs of the previous time step
 * \param[in]      prs_c_pre   pressure cell DoFs of the previous time step
 * \param[in, out] sc          void cast into to a \ref cs_cdofb_ac_t pointer
 */
/*----------------------------------------------------------------------------*/

static void
_implicit_euler_build(const cs_navsto_param_t  *nsp,
                      const cs_real_t           vel_f_pre[],
                      const cs_real_t           vel_c_pre[],
                      const cs_real_t           prs_c_pre[],
                      cs_cdofb_ac_t            *sc)
{
  /* Retrieve high-level structures */

  cs_navsto_ac_t  *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_cdo_system_helper_t  *mom_sh = mom_eqb->system_helper;

  /* Retrieve shared structures */

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  assert(cs_equation_param_has_time(mom_eqp) == true);
  assert(mom_eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(mom_eqb->dir_values != NULL);
#endif

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(nsp,
                                                                    connect);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (mom_eqc->diffusion_hodge == NULL) ? NULL:mom_eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (mom_eqc->mass_hodge == NULL) ? NULL:mom_eqc->mass_hodge[t_id];

    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  inv_dtcur = 1./dt_cur;

    cb->t_pty_eval = t_cur + dt_cur;
    cb->t_bc_eval = t_cur + dt_cur;
    cb->t_st_eval = t_cur + dt_cur;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(mom_eqp, mom_eqb, diff_hodge, cb);

    cs_real_t  o_zeta_c = 1./cs_property_get_cell_value(0, cb->t_pty_eval,
                                                        cc->zeta);

    /* ---------------------------------------------
     * Main loop on cells to build the linear system
     * --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag,
                                                            mom_eqb),
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

      cs_cdofb_vecteq_init_cell_system(cm, mom_eqp, mom_eqb,
                                       vel_f_pre, vel_c_pre,
                                       NULL, NULL, /* no n-1 state is given */
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define the divergence operator used in the linear system (div_op is
       *   equal to minus the divergence)
       */

      cs_cdofb_navsto_define_builder(cb->t_bc_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type, &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(mom_eqp, mom_eqb, mom_eqc, cm,
                                     mass_hodge, diff_hodge, csys, cb);

      /* Update the property related to grad-div */

      if ( !(sc->is_zeta_uniform) )
        o_zeta_c = 1./cs_property_value_in_cell(cm, cc->zeta, cb->t_pty_eval);

      /* Update the property related to the unsteady term */

      if (!(mom_eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm, mom_eqp->time_property,
                                                 cb->t_pty_eval);

      const short int  n_fc = cm->n_fc;

      cs_cdofb_navsto_add_grad_div(n_fc,
                                   cb->tpty_val * dt_cur * o_zeta_c / cm->vol_c,
                                   nsb.div_op, csys->mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after conv/diff/reac and grad-div"
                         " terms", csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm)
        cs_cdofb_vecteq_sourceterm(cm, mom_eqp, cb->t_st_eval,
                                   1., /* time scaling */
                                   mass_hodge,
                                   cb, mom_eqb, csys);

      /* 3b- OTHER RHS CONTRIBUTIONS
       * ===========================
       * Apply the operator gradient to the pressure field and add it to the
       * rhs */

      cs_sdm_add_scalvect(3*n_fc, -prs_c_pre[c_id], nsb.div_op, csys->rhs);

      /* Gravity effects and/or Boussinesq approximation rely on another
         strategy than classical source term. The treatment is more compatible
         with the pressure gradient by doing so. */

      if (sc->add_gravity_term != NULL)
        sc->add_gravity_term(nsp, cm, &nsb, csys);

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _ac_apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type,
                          diff_hodge->pty_data, csys, cb);

      /* 4- TIME CONTRIBUTION */
      /* ==================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

        /* Mass lumping or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get cell-cell block */

        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_fc, n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_fc + k] += ptyc * csys->val_n[3*n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        }

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Only diagonal time treatment available so far.");

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _ac_apply_remaining_bc(sc, mom_eqp, mom_eqb, cm, nsb.bf_type,
                             diff_hodge->pty_data, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys,
                               mom_sh->blocks[0], mom_sh->rhs, mom_eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_cdofb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  cs_cdo_system_helper_finalize_assembly(mom_sh);
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_ac_t structure
 *
 * \param[in] nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in] adv_field   pointer to \ref cs_adv_field_t structure
 * \param[in] mflux       current values of the mass flux across primal faces
 * \param[in] mflux_pre   current values of the mass flux across primal faces
 * \param[in] fb_type     type of boundary for each boundary face
 * \param[in] nsc_input   pointer to a \ref cs_navsto_ac_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_ac_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_init_scheme_context(const cs_navsto_param_t   *nsp,
                                cs_adv_field_t            *adv_field,
                                cs_real_t                 *mflux,
                                cs_real_t                 *mflux_pre,
                                cs_boundary_type_t        *fb_type,
                                void                      *nsc_input)
{
  assert(nsp != NULL && nsc_input != NULL); /* Sanity checks */

  /* Cast the coupling context (CC) */

  cs_navsto_ac_t  *cc = (cs_navsto_ac_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */

  cs_cdofb_ac_t  *sc = NULL;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n",
              __func__);

  BFT_MALLOC(sc, 1, cs_cdofb_ac_t);

  /* Quantities shared with the cs_navsto_system_t structure */

  sc->coupling_context = cc;
  sc->adv_field = adv_field;
  sc->mass_flux_array = mflux;
  sc->mass_flux_array_pre = mflux_pre;

  /* Quick access to the main fields */

  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE)
    sc->divergence = cs_field_by_name("velocity_divergence");
  else
    sc->divergence = NULL;

  /* Parameters related to the AC algorithm */

  sc->is_zeta_uniform = true; /* s_property_is_uniform(cc->zeta); */

  /* Boundary treatment */

  sc->bf_type = fb_type;

  /* Processing of the pressure boundary condition */

  sc->pressure_bc = cs_cdo_bc_face_define(CS_PARAM_BC_HMG_NEUMANN, /* Default */
                                          true, /* Steady BC up to now */
                                          1,    /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          cs_shared_quant->n_b_faces);

  cs_equation_param_t *mom_eqp = cc->momentum->param;
  cs_equation_builder_t  *mom_eqb = cc->momentum->builder;

  mom_eqb->bd_msh_flag |= CS_FLAG_COMP_PFC;

  /* Set the way to enforce the Dirichlet BC on the velocity
   * "fixed_wall" means a no-slip BC */

  sc->apply_symmetry = cs_cdofb_symmetry;

  switch (mom_eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_alge;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_alge;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_alge;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_pena;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_pena;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_pena;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_weak;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_weak;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_weak;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_wsym;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_wsym;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_wsym;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Take into account the gravity effect if needed */

  cs_cdofb_navsto_set_gravity_func(nsp, &(sc->add_gravity_term));

  /* Iterative algorithm to handle the non-linearity (Picard by default) */

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  sc->nl_algo = cs_iter_algo_create(nslesp->nl_algo_param);

  if (nslesp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    sc->nl_algo->context = cs_iter_algo_aa_create(nslesp->anderson_param,
                                                  cs_shared_quant->n_faces);

  /* Monitoring */

  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_ac_t structure
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

  /* Free BC structure */

  sc->pressure_bc = cs_cdo_bc_free(sc->pressure_bc);

  /* If the context is not NULL, this means that an Anderson algorithm has been
     activated otherwise nothing to do */

  cs_iter_algo_aa_free(sc->nl_algo);

  BFT_FREE(sc->nl_algo);

  /* Other pointers are only shared (i.e. not owner) */

  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an AC algorithm
 *         is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_set_sles(const cs_navsto_param_t    *nsp,
                     void                       *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  assert(nsp != NULL && nsc != NULL);

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  mom_eqp->sles_param->field_id = field_id;

  switch (nslesp->strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG:
#if defined(HAVE_PETSC)
    if (mom_eqp->sles_param->amg_type == CS_PARAM_AMG_NONE) {
#if defined(PETSC_HAVE_HYPRE)
      mom_eqp->sles_param->amg_type = CS_PARAM_AMG_HYPRE_BOOMER_V;
#else
      mom_eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_GAMG_V;
#endif
    }

    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         cs_navsto_sles_amg_block_hook,
                         (void *)nsp);
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please build a version of code_saturne with the PETSc support.",
              __func__, mom_eqp->name);
#endif /* HAVE_PETSC */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and an implicit Euler
 *         time scheme
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
  cs_timer_t  t_cmpt = cs_timer_time();

  /* Retrieve high-level structures */

  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;
  cs_navsto_ac_t *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_cdo_system_helper_t  *mom_sh = mom_eqb->system_helper;

  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve fields */

  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *vel_f = mom_eqc->face_values;
  cs_real_t  *div = sc->divergence->val;
  cs_real_t  *pr = sc->pressure->val;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces.
     Evaluation should be performed at t_cur + dt_cur */

  cs_cdofb_vecteq_setup(ts->t_cur + ts->dt[0], mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(mom_sh, &rhs);

  /* Main function for the building stage */

  _implicit_euler_build(nsp,
                        vel_f,  /* previous values for the velocity at faces */
                        vel_c,  /* previous values for the velocity at cells */
                        pr,     /* previous values for the pressure */
                        sc);

  /* Free temporary buffers and structures */

  cs_equation_builder_reset(mom_eqb);

  /* End of the system building */

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Copy current field values to previous values */

  cs_timer_t t_upd = cs_timer_time();

  /* Current to previous for main variable fields */

  _ac_fields_to_previous(sc, cc);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  /* Solve the linear system (treated as a scalar-valued system
   * with 3 times more DoFs) */

  cs_real_t  normalization = 1.0; /* TODO */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(mom_sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(mom_sh, 0);

  int  n_solver_iter = cs_cdo_solve_scalar_system(3*n_faces,
                                                  mom_eqp->sles_param,
                                                  matrix,
                                                  range_set,
                                                  normalization,
                                                  true, /* rhs_redux */
                                                  sles,
                                                  vel_f, /* updated here */
                                                  rhs);

  /* Update pressure, velocity at cells and divergence velocity fields */

  t_upd = cs_timer_time();

  /* Updates after the resolution:
   *  1. the divergence: div = B.u_f
   *  2. the cell velocity field
   *  2. the mass flux
   *  2. the pressure field: pr -= dt / zeta * div(u_f)
   */

  _ac_compute_div(vel_f, div);

  /* Compute values at cells vel_c from values at faces vel_f
     vel_c = acc^-1*(RHS - Acf*vel_f) */

  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  /* Compute the new mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, vel_f, sc->mass_flux_array);

  /* Update the pressure knowing the new divergence of the velocity */

  _ac_update_pr(ts->t_cur, ts->dt[0], cc->zeta, mom_eqp, mom_eqb, div, pr);

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  n_solver_iter);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", n_cells, div, 9);
#endif

  /* Frees */

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(mom_sh);      /* free rhs and matrix */

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a Artificial Compressibility approach and an implicit Euler
 *         time scheme
 *         For Picard - Navier--Stokes problems
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute_implicit_nl(const cs_mesh_t              *mesh,
                                const cs_navsto_param_t      *nsp,
                                void                         *scheme_context)
{
  cs_timer_t  t_cmpt = cs_timer_time();

  /* Retrieve high-level structures */

  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;
  cs_navsto_ac_t *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_cdo_system_helper_t  *mom_sh = mom_eqb->system_helper;
  cs_iter_algo_t  *nl_algo = sc->nl_algo;
  cs_param_nl_algo_t  nl_algo_type = nsp->sles_param->nl_algo_type;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve fields */

  cs_real_t  *vel_f = mom_eqc->face_values;
  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *div = sc->divergence->val;
  cs_real_t  *pr = sc->pressure->val;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces.
     Evaluation should be performed at t_cur + dt_cur */

  cs_cdofb_vecteq_setup(ts->t_cur + ts->dt[0], mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(mom_sh, &rhs);

  /* Main function for the building stage */

  _implicit_euler_build(nsp,
                        vel_f,  /* previous values for the velocity at faces */
                        vel_c,  /* previous values for the velocity at cells */
                        pr,     /* previous values for the pressure */
                        sc);

  /* End of the system building */

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  cs_timer_t t_upd = cs_timer_time();

  /* Copy current field values to previous values */

  _ac_fields_to_previous(sc, cc);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  cs_timer_t  t_solve_start = cs_timer_time();

  cs_iter_algo_reset_nl(nl_algo_type, nl_algo);

  /* Solve the linear system (treated as a scalar-valued system
   * with 3 times more DoFs) */

  cs_real_t  normalization = 1.0; /* TODO */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(mom_sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(mom_sh, 0);

  nl_algo->n_inner_iter = (nl_algo->last_inner_iter =
                           cs_cdo_solve_scalar_system(3*n_faces,
                                                      mom_eqp->sles_param,
                                                      matrix,
                                                      range_set,
                                                      normalization,
                                                      true, /* rhs_redux */
                                                      sles,
                                                      vel_f,
                                                      rhs));

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  t_upd = cs_timer_time();

  /* Updates after the resolution:
   *  1. the divergence: div = B.u_f
   *  2. the mass flux
   */

  _ac_compute_div(vel_f, div);

  /* Compute the new mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, vel_f, sc->mass_flux_array);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  const cs_real_t  *pr_c_pre = sc->pressure->val_pre;
  const cs_real_t  *vel_c_pre = sc->velocity->val_pre;
  const cs_real_t  *vel_f_pre = mom_eqc->face_values_pre;
  assert(vel_c_pre != NULL && vel_f_pre != NULL && pr_c_pre != NULL);

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  /* Set the normalization of the non-linear algo to the value of the first
     mass flux norm */

  nl_algo->normalization =
    sqrt(cs_cdo_blas_square_norm_pfsf(sc->mass_flux_array));

  /* Check the convergence status and update the nl_algo structure related
   * to the convergence monitoring */

  while (cs_cdofb_navsto_nl_algo_cvg(nl_algo_type,
                                     sc->mass_flux_array_pre,
                                     sc->mass_flux_array,
                                     nl_algo) == CS_SLES_ITERATING) {

    /* Main loop on cells to define the linear system to solve */

    cs_timer_t  t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */

    cs_cdo_system_helper_init_system(mom_sh, &rhs);
    cs_sles_free(sles), sles = NULL;

    /* Main loop on cells to define the linear system to solve */

    _implicit_euler_build(nsp,
                          vel_f_pre,  /* velocity at faces: previous values */
                          vel_c_pre,  /* velocity at cells: previous values */
                          pr_c_pre,   /* pressure at cells: previous values */
                          sc);

    /* End of the system building */

    cs_timer_t t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_build_start, &t_build_end);

    /* Solve the new system */

    t_solve_start = cs_timer_time();

    /* Redo the setup since the matrix is modified */

    sles = cs_sles_find_or_add(mom_eqp->sles_param->field_id, NULL);
    cs_sles_setup(sles, matrix);

    nl_algo->n_inner_iter += (nl_algo->last_inner_iter =
                         cs_cdo_solve_scalar_system(3*n_faces,
                                                    mom_eqp->sles_param,
                                                    matrix,
                                                    range_set,
                                                    normalization,
                                                    true, /* rhs_redux */
                                                    sles,
                                                    vel_f,
                                                    rhs));

    t_solve_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

    /* Compute the velocity divergence and retrieve its L2-norm */

    _ac_compute_div(vel_f, div);

    /* Compute the new mass flux used as the advection field */

    memcpy(sc->mass_flux_array_pre, sc->mass_flux_array,
           n_faces*sizeof(cs_real_t));

    cs_cdofb_navsto_mass_flux(nsp, quant, vel_f, sc->mass_flux_array);

  } /* Loop on non linear iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  nl_algo->n_inner_iter);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

  cs_iter_algo_post_check(__func__,
                          mom_eqp->name,
                          cs_param_get_nl_algo_label(nl_algo_type),
                          nl_algo);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_aa_free_arrays(nl_algo->context);

  /* Update pressure and the cell velocity */

  t_upd = cs_timer_time();

  /* Updates after the resolution the pressure field:
   * pr -= dt / zeta * div(u_f) */

  _ac_update_pr(ts->t_cur, ts->dt[0], cc->zeta, mom_eqp, mom_eqb, div, pr);

  /* Compute values at cells vel_c from values at faces vel_f
     vel_c = acc^-1*(RHS - Acf*vel_f) */

  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_AC_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", n_cells, pr, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", n_cells, div, 9);
#endif

  cs_equation_builder_reset(mom_eqb);
  cs_sles_free(sles);
  cs_cdo_system_helper_reset(mom_sh);      /* free rhs and matrix */

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
