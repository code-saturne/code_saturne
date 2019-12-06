/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with a prediction/correction algorithm
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

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

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
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_sles.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_predco.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_predco.c
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

  /*! \var pressure_f
   * Values of the pressure at faces.
   */
  cs_real_t   *pressure_f;

  /*!
   * @}
   * @name Boundary conditions (BC) management
   * Routines and elements used for enforcing the BCs
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
static const cs_matrix_structure_t  *cs_shared_mom_ms;
static const cs_matrix_structure_t  *cs_shared_pre_ms;

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
 * \param[in]      prs_c       value of the pressure at the cell
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_apply_bc_partly(const cs_cdofb_predco_t       *sc,
                 const cs_equation_param_t     *eqp,
                 const cs_cell_mesh_t          *cm,
                 const cs_boundary_type_t      *bf_type,
                 const cs_real_t                prs_c,
                 cs_cell_sys_t                 *csys,
                 cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_diffusion(eqp) == true);

  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the velocity-block and the right-hand side (part related to the
     * momentum equation). */

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann)
      for (short int f  = 0; f < 3*cm->n_fc; f++)
        csys->rhs[f] += csys->neu_values[f];

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];
      const cs_quant_t pfq = cm->face[f];
      const cs_real_t f_prs = pfq.meas * prs_c;
      cs_real_t *f_rhs = csys->rhs + 3*f;

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);
          f_rhs[0] -= f_prs * pfq.unitv[0];
          f_rhs[1] -= f_prs * pfq.unitv[1];
          f_rhs[2] -= f_prs * pfq.unitv[2];
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, cb, csys);
          f_rhs[0] -= f_prs * pfq.unitv[0];
          f_rhs[1] -= f_prs * pfq.unitv[1];
          f_rhs[2] -= f_prs * pfq.unitv[2];
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {
        /* Always weakly enforce the symmetric constraint on the
           velocity-block */
        sc->apply_symmetry(f, eqp, cm, cb, csys);
        f_rhs[0] -= f_prs * pfq.unitv[0];
        f_rhs[1] -= f_prs * pfq.unitv[1];
        f_rhs[2] -= f_prs * pfq.unitv[2];
      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Local system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done after the static condensation
 *
 * \param[in]      sc          pointer to a cs_cdofb_predco_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_apply_remaining_bc(const cs_cdofb_predco_t       *sc,
                    const cs_equation_param_t     *eqp,
                    const cs_cell_mesh_t          *cm,
                    const cs_boundary_type_t      *bf_type,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

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
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {
        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, cb, csys);
        }
      }

#if 0
      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {
        /* Weak-enforcement for the velocity-block (cf. _apply_bc_partly) */
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
 * \param[in]      time_eval  time at which one evaluates quantities
 * \param[in, out] sc         pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_solve_pressure_correction(const cs_mesh_t              *mesh,
                           const cs_navsto_param_t      *nsp,
                           const cs_real_t               time_eval,
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
  const cs_real_t  *const pre_f = sc->pressure_f;

  cs_timer_t  t_bld = cs_timer_time();

  /* Compute the source term */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Compute the divergence of the predicted velocity */
    const cs_real_t  div_c
      = cs_cdofb_navsto_cell_divergence(c_id, quant, connect->c2f, velp_f);

    cc->div_st[c_id] = -div_c * quant->cell_vol[c_id];

  }

  /* Set the boundary conditions on the pressure increment if needed */
  for (int id = 0; id < nsp->n_pressure_bc_defs; id++) {

    const cs_xdef_t  *pdef = nsp->pressure_bc_defs[id];
    const cs_zone_t  *z = cs_boundary_zone_by_id(pdef->z_id);

    if (pdef->meta & CS_CDO_BC_DIRICHLET) {

      switch (pdef->type) {

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        /* Evaluate the boundary condition at each boundary face */
        switch(pre_eqp->dof_reduction) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_eval_at_b_faces_by_analytic(z->n_elts,
                                              z->elt_ids,
                                              false, /* compact output */
                                              mesh,
                                              connect,
                                              quant,
                                              time_eval,
                                              pdef->input,
                                              cc->bdy_pressure_incr);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_eval_avg_at_b_faces_by_analytic(z->n_elts,
                                                  z->elt_ids,
                                                  false, /* compact output */
                                                  mesh,
                                                  connect,
                                                  quant,
                                                  time_eval,
                                                  pdef->input,
                                                  pdef->qtype,
                                                  pdef->dim,
                                                  cc->bdy_pressure_incr);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Invalid type of reduction.\n"
                      " Stop computing the Dirichlet value.\n"), __func__);

        } /* switch on reduction */

        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          const cs_lnum_t  f_id = z->elt_ids[i];
          cc->bdy_pressure_incr[f_id] -= pre_f[f_id + quant->n_i_faces];
        }
        break;

      case CS_XDEF_BY_VALUE:    /* This is an homogeneous Dirichlet BC */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Not implemented yet.", __func__);
        break;
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0, "%s: Not implemented yet.", __func__);

  } /* Loop on pressure definitions */

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(pre_eqb->tcb), &t_bld, &t_tmp);

  /* Solve the equation related to the pressure increment */
  cs_cdofb_scaleq_solve_steady_state(mesh,
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
 * \param[in]      dt_cur  current value of the time step
 * \param[in, out] sc      pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_variables(const cs_real_t              dt_cur,
                  cs_cdofb_predco_t           *sc)
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
  const cs_field_t  *velp_fld = cc->predicted_velocity;
  const cs_real_t *const  velp_c = velp_fld->val;
  const cs_real_t *const  velp_f = sc->predicted_velocity_f;
  const cs_real_t *const  dp_f = cs_cdofb_scaleq_get_face_values(pre_eqc);
  const cs_real_t *const  dp_c = cs_cdofb_scaleq_get_cell_values(pre_eqc);

  /* Variables to update */
  cs_real_t  *pre_f = sc->pressure_f;
  cs_real_t  *pr_c = sc->pressure->val; /* cell DoFs for the pressure */
  cs_real_t  *vel_c = sc->velocity->val;
  cs_real_t  *vel_f = mom_eqc->face_values;
  cs_real_t  *div = sc->divergence->val;

  cs_timer_t  t_upd = cs_timer_time();

  cs_field_current_to_previous(sc->velocity);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Retrieve the cell mesh structure w.r.t. the mode */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(omp_get_thread_num());
#else
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);
#endif
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;
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
      pr_c[c_id] += dt_cur*dp_c[c_id];

      /* Evaluate the cell gradient for the pressure increment */
      cs_real_t  grd_dp[3] = {0, 0, 0};
      for (short int f = 0; f < cm->n_fc; f++) {
        const cs_real_t  f_coef = cm->f_sgn[f]*cm->face[f].meas;
        for (int k = 0; k < 3; k++)
          grd_dp[k] += f_coef*dp_f[cm->f_ids[f]]*cm->face[f].unitv[k];
      }

      /* Update the cell velocity */
      for (int k = 0; k < 3; k++)
        vel_c[3*c_id + k] = velp_c[3*c_id + k] + dt_cur*grd_dp[k];

      /* Partial update of the face velocity:
       * v_f^(n+1) = vp_f^(n+1) + dt*grd(incr_p)
       * Now: 0.5*dt_cur*grd_cell(incr_p) for an interior face
       * or       dt_cur*grd_cell(incr_p)     for a border face  */
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_lnum_t  f_id = cm->f_ids[f];

        cs_real_t  f_coef = dt_cur;
        if (f_id < quant->n_i_faces)
          f_coef *= 0.5;

        cs_real_t  *_vel = vel_f + 3*f_id;
        for (int k = 0; k < 3; k++)
          _vel[k] += f_coef*grd_dp[k];

      }

    } /* Loop on cells */

  } /* OpenMP block */

  /* Parallel sum */
  if (cs_glob_n_ranks > 1) {
    assert(connect->interfaces[CS_CDO_CONNECT_FACE_SP0] != NULL);
    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_FACE_SP0],
                         n_faces,
                         3,
                         true,
                         CS_REAL_TYPE,
                         vel_f);
  }

  /* Update face-related unknowns */
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    /* p^(n+1) = p^n + delta_p */
    pre_f[f] += dp_f[f]*dt_cur;

    /* v_f^(n+1) = vp_f^(n+1) + dt*grd(incr_p) */
    for (int k = 0; k < 3; k++)
      vel_f[3*f+k] += velp_f[3*f+k];

  } /* Loop on faces */

  const cs_adjacency_t  *c2f = connect->c2f;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Reset the divergence of the velocity before computing the updated
       value (vel_f has to be updated first) */
    div[c_id] = 0;
    for (cs_lnum_t jf = c2f->idx[c_id]; jf < c2f->idx[c_id+1]; jf++) {

      const cs_lnum_t  f_id = c2f->ids[jf];
      const cs_real_t  *_vel = vel_f + 3*f_id;
      const cs_real_t  *nf = cs_quant_get_face_vector_area(f_id, quant);

      div[c_id] += c2f->sgn[jf]*cs_math_3_dot_product(_vel, nf);

    } /* Loop on cell faces */

    const cs_real_t  ovc = 1./quant->cell_vol[c_id];
    div[c_id] *= ovc;

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
 * \brief  Retrieve the values of the pressure at faces
 *
 * \param[in]  context     pointer to a scheme context structure
 *
 * \return a pointer to the pressure values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_predco_get_face_pressure(void     *context)
{
  cs_cdofb_predco_t  *sc = (cs_cdofb_predco_t  *)context;

  if (sc == NULL)
    return NULL;

  return sc->pressure_f;
}

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

  /*
   * Matrix structure related to the algebraic system for
   * the momentum equation --> vector-valued equation
   * the pressure equation --> scalar-valued equation
   */
  cs_shared_mom_ms = cs_cdofb_vecteq_matrix_structure();
  cs_shared_pre_ms = cs_cdofb_scaleq_matrix_structure();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_cdofb_predco_t structure
 *
 * \param[in] nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in] fb_type     type of boundary for each boundary face
 * \param[in] nsc_input   pointer to a \ref cs_navsto_predco_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_predco_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_predco_init_scheme_context(const cs_navsto_param_t    *nsp,
                                    cs_boundary_type_t         *fb_type,
                                    void                       *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Cast the coupling context (CC) */
  cs_navsto_projection_t  *cc = (cs_navsto_projection_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_predco_t  *sc = NULL;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  BFT_MALLOC(sc, 1, cs_cdofb_predco_t);

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  /* Quick access to the main fields */
  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE)
    sc->divergence = cs_field_by_name("velocity_divergence");
  else
    sc->divergence = NULL;

  /* Values of the predicted velocity at faces */
  BFT_MALLOC(sc->predicted_velocity_f, 3*quant->n_faces, cs_real_t);
  memset(sc->predicted_velocity_f, 0, 3*quant->n_faces*sizeof(cs_real_t));

  /* Values of the pressure at faces */
  BFT_MALLOC(sc->pressure_f, quant->n_faces, cs_real_t);
  memset(sc->pressure_f, 0, quant->n_faces*sizeof(cs_real_t));

  /* Boundary treatment */
  sc->bf_type = fb_type;

  /* Processing of the pressure boundary condition */
  sc->pressure_bc = cs_cdo_bc_face_define(CS_CDO_BC_HMG_NEUMANN, /* Default */
                                          true, /* Steady BC up to now */
                                          1,    /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          cs_shared_quant->n_b_faces);

  cs_equation_param_t  *mom_eqp = cc->prediction->param;
  cs_equation_builder_t  *mom_eqb = cc->prediction->builder;

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

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_predco_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_predco_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_predco_t  *sc = (cs_cdofb_predco_t  *)scheme_context;

  if (sc == NULL)
    return sc;

  /* Free BC structure */
  sc->pressure_bc = cs_cdo_bc_free(sc->pressure_bc);

  BFT_FREE(sc->predicted_velocity_f);
  BFT_FREE(sc->pressure_f);

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a projection
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context   pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_predco_set_sles(const cs_navsto_param_t    *nsp,
                         void                       *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_navsto_param_sles_t  nslesp = nsp->sles_param;
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->prediction);
  int  field_id = cs_equation_get_field_id(nsc->prediction);

  mom_eqp->sles_param.field_id = field_id;

  switch (nslesp.strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG:
#if defined(HAVE_PETSC)
    if (mom_eqp->sles_param.amg_type == CS_PARAM_AMG_NONE) {
#if defined(PETSC_HAVE_HYPRE)
      mom_eqp->sles_param.amg_type = CS_PARAM_AMG_HYPRE_BOOMER;
#else
      mom_eqp->sles_param.amg_type = CS_PARAM_AMG_PETSC_GAMG;
#endif
    }

    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         cs_navsto_sles_amg_block_hook,
                         (void *)mom_eqp);
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please build a version of Code_Saturne with the PETSc support.",
              __func__, mom_eqp->name);
#endif /* HAVE_PETSC */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

  /* For the correction step, use the generic way to setup the SLES */
  cs_equation_param_t  *corr_eqp = cs_equation_get_param(nsc->correction);

  corr_eqp->sles_param.field_id = cs_equation_get_field_id(nsc->correction);
  cs_equation_param_set_sles(corr_eqp);

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
  CS_UNUSED(nsp);

  cs_timer_t  t_cmpt = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_predco_t  *sc = (cs_cdofb_predco_t *)scheme_context;
  cs_navsto_projection_t *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->prediction;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *mom_rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  cs_real_t  *pr_c = sc->pressure->val; /* cell DoFs for the pressure */
  cs_field_t  *vel_fld = sc->velocity;
  cs_real_t  *vel_c = vel_fld->val;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  /* Sanity checks */
  assert(cs_equation_param_has_time(mom_eqp) == true);
  assert(mom_eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces and ids of DoFs if
   * an enforcement of (internal) DoFs is requested
   * Evaluation should be performed at t_cur + dt_cur */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *enforced_ids = NULL;

  cs_cdofb_vecteq_setup(t_cur + dt_cur, mesh, mom_eqp, mom_eqb,
                        &dir_values, &enforced_ids);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_mom_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix, nsp,   \
         mav, mom_rs, dir_values, enforced_ids, vel_c, pr_c, sc)        \
  firstprivate(time_eval, inv_dtcur)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(connect);
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

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
                                       dir_values, enforced_ids,
                                       vel_c, time_eval,
                                       csys, cb);

      const short int  n_fc = cm->n_fc, nf_dofs = 3*n_fc;

       /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define the divergence operator used in the linear system (div_op is
       *  equal to minus the divergence)
       */
      cs_cdofb_navsto_define_builder(time_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */
      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqc, cm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */
      const bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   time_eval, 1., /* time, scaling */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* OTHER RHS CONTRIBUTIONS
       * =======================
       * Apply the operator gradient to the pressure field and add it to the
       * rhs
       */
      cs_sdm_add_scalvect(nf_dofs, -pr_c[c_id], nsb.div_op, csys->rhs);

      /* First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */
      _apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, pr_c[c_id], csys, cb);

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
        bft_error(__FILE__, __LINE__, 0,
                  "Only diagonal time treatment available so far.");

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */
      _apply_remaining_bc(sc, mom_eqp, cm, nsb.bf_type, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_PREDCO_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      cs_cdofb_vecteq_assembly(csys, mom_rs, cm, has_sourceterm,
                               mom_eqc, eqa, mav, rhs);

    } /* Main loop on cells */

    /* Free temporary buffers */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(enforced_ids);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Copy current field values to previous values */
  cs_field_t  *velp_fld = cc->predicted_velocity;

  cs_field_current_to_previous(velp_fld);

  cs_real_t  *velp_c = velp_fld->val;
  cs_real_t  *velp_f = sc->predicted_velocity_f;

  /* Solve the linear system (treated as a scalar-valued system
   * with 3 times more DoFs) */
  cs_real_t  normalization = 1.0; /* TODO */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eqp->sles_param.field_id, NULL);

  cs_equation_solve_scalar_system(3*n_faces,
                                  mom_eqp,
                                  matrix,
                                  mom_rs,
                                  normalization,
                                  true, /* rhs_redux */
                                  sles,
                                  velp_f,
                                  rhs);

  /* Update pressure, velocity and divergence fields */
  cs_timer_t  t_upd = cs_timer_time();

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        velp_f, velp_c);

  /* Frees */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);

  /* SECOND MAJOR STEP: Solve the equation related to the correction step
   * Correction step = Solve a diffusion equation for the pressure increment */

  _solve_pressure_correction(mesh, nsp, time_eval, sc);

  /* LAST MAJOR STEP: Update the pressure and the velocity */

  _update_variables(dt_cur, sc);

  if (nsp->n_pressure_bc_defs == 0)
    cs_cdofb_navsto_set_zero_mean_pressure(quant, pr_c);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
