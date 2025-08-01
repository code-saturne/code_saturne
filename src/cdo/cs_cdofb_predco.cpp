/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with a prediction/correction algorithm
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "base/cs_mem.h"

#include "base/cs_array.h"
#include "alge/cs_blas.h"
#include "cdo/cs_cdo_assembly.h"
#include "cdo/cs_cdo_bc.h"
#include "cdo/cs_cdo_solve.h"
#include "cdo/cs_cdofb_priv.h"
#include "cdo/cs_cdofb_scaleq.h"
#include "cdo/cs_cdofb_vecteq.h"
#include "cdo/cs_cdofb_navsto.h"
#include "cdo/cs_equation_bc.h"
#include "cdo/cs_equation_priv.h"
#include "cdo/cs_evaluate.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_post.h"
#include "cdo/cs_reco.h"
#include "cdo/cs_source_term.h"
#include "cdo/cs_static_condensation.h"
#include "base/cs_timer.h"

#if defined(HAVE_PETSC)
#include "alge/cs_sles_petsc.h"
#endif

#if defined(DEBUG) && !defined(NDEBUG)
#include "cdo/cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdofb_predco.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_predco.cpp
 *
 * \brief Build an algebraic CDO face-based system for the Navier-Stokes
 * equations and solved it with a prediction/correction algorithm. A first
 * equation related to the velocity prediction is solved and then a second
 * equation related to the pressure correction is solved to project the velocity
 * field into the space of divergence-free field.
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_predco_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         Navier-Stokes equations with a prediction/correction algorithm
 */

typedef struct {

  /*! \var coupling_context

   *  Pointer to a \ref cs_navsto_projection_t_t (owned by
   *  \ref cs_navsto_system_t) containing the settings related to a prjection
   *  or prediction/correction algorithm.
   */

  cs_navsto_projection_t   *coupling_context;

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

  /*! \var predicted_velocity_f
   * Values of the predicted velocity at faces.
   * This values may not be divergence-free
   */
  cs_real_t   *predicted_velocity_f;

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
  int                             pressure_rescaling;

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

} cs_cdofb_predco_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_PREDCO_DBG      0

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
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done before the static condensation
 *
 * \param[in]      sc          pointer to a cs_cdofb_predco_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in]      diff_pty    pointer to \ref cs_property_data_t for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_predco_apply_bc_partly(const cs_cdofb_predco_t       *sc,
                        const cs_equation_param_t     *eqp,
                        const cs_cell_mesh_t          *cm,
                        const cs_boundary_type_t      *bf_type,
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
      const cs_quant_t pfq = cm->face[f];
      const cs_real_t f_prs = pfq.meas * pc;
      cs_real_t *f_rhs = csys->rhs + 3*f;

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done after the static condensation
 *          Case of CDO-Fb schemes with Prediction-Correction coupling
 *
 * \param[in]      sc          pointer to a cs_cdofb_predco_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in]      diff_pty    pointer to \ref cs_property_data_t for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_predco_apply_remaining_bc(const cs_cdofb_predco_t       *sc,
                           const cs_equation_param_t     *eqp,
                           const cs_equation_builder_t   *eqb,
                           const cs_cell_mesh_t          *cm,
                           const cs_boundary_type_t      *bf_type,
                           const cs_property_data_t      *diff_pty,
                           cs_cell_sys_t                 *csys,
                           cs_cell_builder_t             *cb)
{
  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    cs_equation_builder_enforce_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }

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
        /* Weak-enforcement for the velocity-block
           (cf. _predco_apply_bc_partly) */
      }
#endif

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop boundary faces */

  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare and solve the pressure correction step
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] sc         pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_solve_pressure_correction(const cs_mesh_t              *mesh,
                           const cs_navsto_param_t      *nsp,
                           cs_cdofb_predco_t            *sc)
{
  cs_navsto_projection_t *cc = sc->coupling_context;
  cs_equation_t  *pre_eq = cc->correction;
  cs_equation_param_t  *pre_eqp = cs_equation_get_param(pre_eq);
  void  *pre_eqc = cs_equation_get_scheme_context(pre_eq);
  cs_equation_builder_t  *pre_eqb = pre_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t *velp_f = sc->predicted_velocity_f;

  /* Boundary conditions are always evaluated at t + dt */

  cs_timer_t  t_bld = cs_timer_time();

  /* Compute the source term */
# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Compute the divergence of the predicted velocity */

    const cs_real_t  div_c
      = cs_cdofb_navsto_cell_divergence(c_id, quant, connect->c2f, velp_f);

    cc->div_st[c_id] = - div_c * quant->cell_vol[c_id];
  }

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(pre_eqb->tcb), &t_bld, &t_tmp);

  /* Solve the equation related to the pressure increment */

  pre_eq->solve_steady_state(true,
                             mesh,
                             pre_eq->field_id,
                             pre_eqp,
                             pre_eqb,
                             pre_eqc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update velocity and pressure-related variables after the correction
 *         step
 *
 * \param[in, out] sc      pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_variables(const cs_navsto_param_t           *nsp,
                  const cs_cdofb_predco_t           *sc)
{
  cs_navsto_projection_t *cc = sc->coupling_context;
  cs_equation_t  *pre_eq = cc->correction;
  void  *pre_eqc = cs_equation_get_scheme_context(pre_eq);
  cs_equation_t  *mom_eq = cc->prediction;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_real_t *const  velp_f = sc->predicted_velocity_f;
  cs_real_t *velp_c = mom_eq->get_cell_values(mom_eqc, false);

  const cs_real_t *const dp_c = pre_eq->get_cell_values(pre_eqc, false);
  const cs_adjacency_t  *c2f = connect->c2f;

  /* Variables to update */

  cs_real_t  *pr_c = sc->pressure->val; /* cell DoFs for the pressure */
  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *vel_f = mom_eqc->face_values;

  cs_real_t  *div = nullptr;
  if (sc->divergence != nullptr)
    div = sc->divergence->val;

  cs_real_3_t *grad_dp = nullptr;
  if (cc->pressure_incr_gradient != nullptr)
    grad_dp = (cs_real_3_t*) cc->pressure_incr_gradient->val;

  cs_real_3_t  *grad_p = nullptr;
  if (nsp->post_flag & CS_NAVSTO_POST_PRESSURE_GRADIENT)
    grad_p = (cs_real_3_t*) cs_field_by_name_try("pressure_gradient")->val;

  cs_real_t  *mass_flux_array = sc->mass_flux_array;

  cs_timer_t  t_upd = cs_timer_time();

  cs_field_current_to_previous(sc->velocity);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);

  cs_real_t *diff_flux = nullptr;
  CS_MALLOC(diff_flux, n_faces, cs_real_t);

  cs_equation_compute_diffusive_flux(pre_eq,
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     cs_flag_primal_face,
                                     cs_glob_time_step->t_cur,
                                     diff_flux);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(cs_get_thread_id());
    cs_cell_sys_t     *csys = nullptr;
    cs_cell_builder_t *cb   = nullptr;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Reset the velocity at faces */

#   pragma omp for
    for (cs_lnum_t i = 0; i < 3*n_faces; i++) vel_f[i] = 0;

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         CS_FLAG_COMP_PF |  CS_FLAG_COMP_PFQ,
                         connect, quant, cm);

      /* Update the cell pressure */

      pr_c[c_id] += dp_c[c_id];

      if (grad_p != nullptr || grad_dp != nullptr) {
        /* Reconstruct the cell gradient
           ----------------------------- */

        const cs_lnum_t  s = c2f->idx[c_id],
                         e = c2f->idx[c_id+1];
        const cs_lnum_t  *c2f_ids = c2f->ids + s;
        const short int  *c2f_sgn = c2f->sgn + s;

        const cs_lnum_t n_fc = e - s;
        const cs_real_t vol_c = quant->cell_vol[c_id];
        const cs_real_t *cell_center = quant->cell_centers[c_id];
        cs_real_3_t grddp_reco = {0.0, 0.0, 0.0};

        for (short int f_id = 0; f_id < n_fc; f_id++) {

          cs_lnum_t ff = c2f_ids[f_id];

          const cs_real_t *face_center  = (ff < quant->n_i_faces) ?
            quant->i_face_center + 3*ff:
            quant->b_face_center + 3*(ff - quant->n_i_faces);

          for (short int k = 0; k < 3; k++)
            grddp_reco[k] += - c2f_sgn[f_id]/vol_c
                             *diff_flux[ff]*(face_center[k] - cell_center[k]);
        }

        if (grad_p != nullptr)
          for (short int k = 0; k < 3; k++)
            grad_p[c_id][k] += grddp_reco[k];

        if (grad_dp != nullptr)
          for (short int k = 0; k < 3; k++)
            grad_dp[c_id][k] = grddp_reco[k];
      }

      /* Partial update of the face velocity:
       * v_f^(n+1) = vp_f^(n+1) + dt*grd(incr_p)
       * Now: 0.5*dt_cur*grd_cell(incr_p) for an interior face
       * or       dt_cur*grd_cell(incr_p)     for a border face  */

      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_lnum_t  f_id = cm->f_ids[f];
        cs_real_3_t dtgradf_dp = {0.0, 0.0, 0.0};

        const cs_quant_t pfq = cm->face[f];
        for (short int k = 0; k < 3; k++)
          dtgradf_dp[k] = diff_flux[f_id]/pfq.meas*pfq.unitv[k];

        cs_real_t  f_coef = 1.0;
        if (f_id < quant->n_i_faces)
          f_coef = 0.5;
        cs_real_t  *_vel = vel_f + 3*f_id;
        for (int k = 0; k < 3; k++)
          _vel[k] += f_coef*dtgradf_dp[k];
      }

    } /* Loop on cells */

  } /* OpenMP block */

  CS_FREE(diff_flux);

  /* Parallel or periodic sum */

  if (connect->face_ifs != nullptr)
    cs_interface_set_sum(connect->face_ifs,
                         n_faces,
                         3,
                         true,
                         CS_REAL_TYPE,
                         vel_f);

  /* Update face-related unknowns */
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    /* v_f^(n+1) = vp_f^(n+1) + dt*grd(incr_p) */

    for (int k = 0; k < 3; k++)
      vel_f[3*f+k] += velp_f[3*f+k];

  } /* Loop on faces */

  /* Update the mass flux
   * Mass flux updating must be consistent
   * with face velocity updating
   * ----------------------------------
   * */

  cs_cdofb_navsto_mass_flux(nsp, quant, vel_f, mass_flux_array);

  /* Update the cell velocity */

  /* Compute values at cells pc from values at faces pf
     vc = acc^-1*(RHS - Acf*vf) */

  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f, vel_c);

  cs_array_real_copy(3*quant->n_cells, vel_c, velp_c);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    div[c_id] = 0.0;
    for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

      cs_lnum_t f_id = c2f->ids[f];
      div[c_id] += c2f->sgn[f]*mass_flux_array[f_id];

    } /* Loop on cell faces */

    div[c_id] /= quant->cell_vol[c_id];

  } /* Loop on cells */

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 2
  cs_dbg_darray_to_listing("PRED_VELOCITY_FACE", n_faces, velp_f, 9);
  cs_dbg_darray_to_listing("VELOCITY_FACE", n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_c, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif
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
cs_cdofb_predco_init_common(const cs_cdo_quantities_t     *quant,
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
 * \brief  Initialize a cs_cdofb_predco_t structure
 *
 * \param[in] nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in] adv_field   pointer to \ref cs_adv_field_t structure
 * \param[in] mflux       current values of the mass flux across primal faces
 * \param[in] mflux_pre   previous values of the mass flux across primal faces
 * \param[in] fb_type     type of boundary for each boundary face
 * \param[in] nsc_input   pointer to a \ref cs_navsto_projection_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_predco_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_predco_init_scheme_context(const cs_navsto_param_t   *nsp,
                                    cs_adv_field_t            *adv_field,
                                    cs_real_t                 *mflux,
                                    cs_real_t                 *mflux_pre,
                                    cs_boundary_type_t        *fb_type,
                                    void                      *nsc_input)
{
  assert(nsp != nullptr && nsc_input != nullptr);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Cast the coupling context (CC) */

  cs_navsto_projection_t  *cc = (cs_navsto_projection_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */

  cs_cdofb_predco_t *sc = nullptr;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  CS_MALLOC(sc, 1, cs_cdofb_predco_t);

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
    sc->divergence = nullptr;

  /* Values of the predicted velocity at faces */

  CS_MALLOC(sc->predicted_velocity_f, 3 * quant->n_faces, cs_real_t);
  cs_array_real_fill_zero(3*quant->n_faces, sc->predicted_velocity_f);

  /* Boundary treatment */

  sc->bf_type = fb_type;
  sc->pressure_bc =
    cs_cdo_bc_face_define(CS_BC_SYMMETRY,          /* Default */
                          true,                    /* Steady BC up to now */
                          1,                       /* Dimension */
                          nsp->n_pressure_bc_defs,
                          nsp->pressure_bc_defs,
                          quant->n_b_faces);

  sc->pressure_rescaling =
    cs_boundary_need_pressure_rescaling(quant->n_b_faces, fb_type);

  cs_equation_param_t  *mom_eqp = cc->prediction->param;
  cs_equation_builder_t  *mom_eqb = cc->prediction->builder;

  mom_eqb->bdy_flag |= CS_FLAG_COMP_PFC;

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

  /* Monitoring */

  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_predco_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_predco_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_predco_t  *sc = (cs_cdofb_predco_t  *)scheme_context;

  if (sc == nullptr)
    return sc;

  /* Free BC structure */

  sc->pressure_bc = cs_cdo_bc_free(sc->pressure_bc);

  CS_FREE(sc->predicted_velocity_f);

  /* Other pointers are only shared (i.e. not owner) */

  CS_FREE(sc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a prediction/correction approach and an implicit Euler time
 *         scheme
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_predco_compute_implicit(const cs_mesh_t              *mesh,
                                 const cs_navsto_param_t      *nsp,
                                 void                         *scheme_context)
{
  cs_timer_t  t_cmpt = cs_timer_time();

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve high-level structures */

  cs_cdofb_predco_t  *sc = (cs_cdofb_predco_t *)scheme_context;
  cs_navsto_projection_t *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->prediction;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_cdo_system_helper_t  *mom_sh = mom_eqb->system_helper;

  assert(cs_equation_param_has_time(mom_eqp) == true);
  assert(mom_eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  /* Retrieve fields */

  cs_real_t  *pr_c = sc->pressure->val; /* cell DoFs for the pressure */
  cs_real_t  *vel_c = mom_eq->get_cell_values(mom_eqc, false);

//  cs_array_real_fill_zero(quant->n_cells, pr_c);
  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces and if needed one
   * computes the enforced values. Evaluation should be performed at t_cur +
   * dt_cur
   */

  cs_cdofb_vecteq_setup(ts->t_cur + ts->dt[0], mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t *rhs = nullptr; /* Since it is nullptr, sh get sthe ownership */

  cs_cdo_system_helper_init_system(mom_sh, &rhs);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    const int  t_id = cs_get_thread_id();

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(nsp,
                                                                    connect);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t *diff_hodge = (mom_eqc->diffusion_hodge == nullptr)
                            ? nullptr
                            : mom_eqc->diffusion_hodge[t_id];
    cs_hodge_t *mass_hodge
      = (mom_eqc->mass_hodge == nullptr) ? nullptr : mom_eqc->mass_hodge[t_id];

    cs_cell_sys_t     *csys = nullptr;
    cs_cell_builder_t *cb   = nullptr;

    cs_cdofb_vecteq_get(&csys, &cb);

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  t_eval = t_cur + dt_cur;
    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = t_eval;
    cb->t_bc_eval = t_eval;
    cb->t_st_eval = t_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(mom_eqp, mom_eqb, diff_hodge, cb);

    /* ---------------------------------------------
     * Main loop on cells to build the linear system
     * --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_eflag_t  cell_mesh_flag =
        cs_equation_builder_cell_mesh_flag(cb->cell_flag, mom_eqb);

      cs_cell_mesh_build(c_id, cell_mesh_flag, connect, quant, cm);

      /* Starts from the steady Stokes problem where
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

      cs_cdofb_vecteq_init_cell_system(cm,
                                       mom_eqp,
                                       mom_eqb,
                                       mom_eqc->face_values,
                                       vel_c,
                                       nullptr,
                                       nullptr, /* no n-1 state is given */
                                       csys,
                                       cb);

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system after conv/diff/reac", csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm)
        cs_cdofb_vecteq_sourceterm(cm, mom_eqp, cb->t_st_eval,
                                   1., /* scaling */
                                   mass_hodge,
                                   cb, mom_eqb, csys);

      /* OTHER RHS CONTRIBUTIONS
       * =======================
       * Apply the operator gradient to the pressure field and add it to the
       * rhs */

      const short int  n_fc = cm->n_fc, nf_dofs = 3*n_fc;

      /* Pressure flux at the interior faces */
      cs_sdm_add_scalvect(nf_dofs, -pr_c[c_id], nsb.div_op, csys->rhs);

      /* Pressure flux at the boundary faces */
      for (short int i = 0; i < csys->n_bc_faces; i++) {
        if (nsb.bf_type[i] & CS_BOUNDARY_IMPOSED_P) {
          const short int  f = csys->_f_ids[i];
          assert(nsb.bf_type[i] & (CS_BOUNDARY_INLET | CS_BOUNDARY_OUTLET));
          for (short int k = 0; k < 3; k++)
            csys->rhs[3*f + k] += nsb.pressure_bc_val[i]*nsb.div_op[3*f + k];
        }
      }

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _predco_apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type,
                              diff_hodge->pty_data, csys, cb);

      /* 4- UNSTEADY TERM + TIME SCHEME
       * ============================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

        /* Mass lumping or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get the cell-cell block */

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

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _predco_apply_remaining_bc(sc, mom_eqp, mom_eqb, cm, nsb.bf_type,
                                 diff_hodge->pty_data, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys,
                               mom_sh->blocks[0], mom_sh->rhs, mom_eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffers */

    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(mom_sh);
  cs_equation_builder_reset(mom_eqb);

  /* End of the system building */

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Copy current field values to previous values */

  cs_field_t  *velp_fld = cc->predicted_velocity;

  cs_field_current_to_previous(velp_fld);

  cs_real_t  *velp_f = sc->predicted_velocity_f;

  /* Solve the linear system (treated as a scalar-valued system
   * with 3 times more DoFs) */

  cs_real_t  normalization = 1.0; /* TODO */
  cs_sles_t *sles = cs_sles_find_or_add(mom_eqp->sles_param->field_id, nullptr);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(mom_sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(mom_sh, 0);

  cs_cdo_solve_scalar_system(3*n_faces,
                             mom_eqp->sles_param,
                             matrix,
                             range_set,
                             normalization,
                             true, /* rhs_redux */
                             sles,
                             velp_f,
                             rhs);

  cs_timer_t  t_upd = cs_timer_time();

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(mom_sh);      /* free rhs and matrix */

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);

  /* SECOND MAJOR STEP: Solve the equation related to the correction step
   * Correction step = Solve a diffusion equation for the pressure increment */

  _solve_pressure_correction(mesh, nsp, sc);

  /* LAST MAJOR STEP: Update the pressure and the velocity */

  _update_variables(nsp, sc);

  if (sc->pressure_rescaling == CS_BOUNDARY_PRESSURE_RESCALING)
    cs_cdofb_navsto_rescale_pressure_to_ref(nsp, quant, pr_c);
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
