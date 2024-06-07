/*============================================================================
 * Build an algebraic MAC face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
 * coupling)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_blas.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_blas.h"
#include "cs_cdofb_monolithic.h"
#include "cs_cdofb_monolithic_sles.h"
#include "cs_cdofb_navsto.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_macfb_builder.h"
#include "cs_macfb_monolithic_sles.h"
#include "cs_macfb_navsto.h"
#include "cs_macfb_priv.h"
#include "cs_macfb_vecteq.h"
#include "cs_navsto_coupling.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sdm.h"
#include "cs_solid_selection.h"
#include "cs_source_term.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_monolithic.h"
#include "cs_macfb_monolithic_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_macfb_monolithic.c
 *
 * \brief Build an algebraic MAC face-based system for the Navier-Stokes
 *        equations and solved it with a monolithic approach
 */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_MACFB_MONOLITHIC_DBG 0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_mesh_t           *cs_shared_mesh;
static const cs_cdo_quantities_t *cs_shared_quant;
static const cs_cdo_connect_t    *cs_shared_connect;
static const cs_time_step_t      *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the Navier-Stokes related variables-fields
 *         to previous values
 *         Case of the monolithic coupling algorithm.
 *
 * \param[in, out] sc    scheme context (\ref cs_macfb_monolithic_t struct.)
 * \param[in, out] cc    coupling context (\ref cs_navsto_monolithic_t struct.)
 */
/*----------------------------------------------------------------------------*/

static inline void
_mono_fields_to_previous(cs_macfb_monolithic_t *sc, cs_navsto_monolithic_t *cc)
{
  const cs_cdo_quantities_t *cdoq = cs_shared_quant;

  /* Cell unknows: pressure */

  cs_field_current_to_previous(sc->pressure);

  /* Face unknows: mass flux and face velocity */

  cs_array_real_copy(
    cdoq->n_faces, sc->mass_flux_array, sc->mass_flux_array_pre);

  cs_macfb_vecteq_t *eqc = (cs_macfb_vecteq_t *)cc->momentum->scheme_context;

  if (eqc->face_values_pre != nullptr)
    cs_array_real_copy(cdoq->n_faces, eqc->face_values, eqc->face_values_pre);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Make sure that the enforcement is done (remove tolerance that may
 *         arise from the resolution)
 *         Case of a monolithic coupling algorithm.
 *
 * \param[in, out]  vel_f    velocity at faces to enforce
 */
/*----------------------------------------------------------------------------*/

static void
_mono_enforce_solid_face_velocity(cs_real_t *vel_f)
{
  /* Enforcement of solid cells is always defined as follows for the momentum
   * equation:
   * CS_EQUATION_ENFORCE_BY_CELLS | CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE */

  cs_solid_selection_t *solid = cs_solid_selection_get();

  if (solid->n_g_cells > 0) {

    for (cs_lnum_t f = 0; f < cs_shared_connect->n_faces[0]; f++) {

      if (solid->face_is_solid[f]) {
        vel_f[f] = 0.;
      }

    } /* Loop on faces and search for faces with a tag */

  } /* There is at least one solid cell globally */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Rescale the pressure field
 *         Case of a monolithic coupling algorithm.
 *
 * \param[in]       nsp      set of parameters for the Navier-Stokes system
 * \param[in, out]  sc       scheme context
 * \param[in, out]  eqc  context of the momentum equation
 */
/*----------------------------------------------------------------------------*/

static void
_mono_update_related_cell_fields(const cs_navsto_param_t *nsp,
                                 cs_macfb_monolithic_t   *sc,
                                 cs_macfb_vecteq_t       *eqc)
{
  CS_UNUSED(eqc);
  const cs_cdo_quantities_t *quant = cs_shared_quant;

  /* Rescale pressure if needed */

  cs_field_t *pr_fld = sc->pressure;

  if (sc->pressure_rescaling == CS_BOUNDARY_PRESSURE_RESCALING)
    cs_cdofb_navsto_rescale_pressure_to_ref(nsp, quant, pr_fld->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_MONOLITHIC_DBG > 2
  const cs_real_t *vel_f = eqc->face_values;
  cs_dbg_darray_to_listing("FACE_VELOCITY", quant->n_faces, vel_f, 9);
#endif
#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_fld->val, 9);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize stuctures for a gven cell
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      sc           pointer to the scheme context
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      c_id         cell id
 * \param[in]      vel_f_n      velocity face DoFs of the previous time step
 * \param[in, out] cm           pointer to a cellwise view of the mesh
 * \param[in, out] macb         pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_build(const cs_navsto_param_t     *nsp,
            const cs_macfb_monolithic_t *sc,
            const cs_cdo_connect_t      *connect,
            const cs_cdo_quantities_t   *quant,
            const cs_equation_param_t   *eqp,
            const cs_equation_builder_t *eqb,
            const cs_lnum_t              c_id,
            const cs_real_t              vel_f_n[],
            cs_macfb_navsto_builder_t   *nsb,
            cs_cell_mesh_t              *cm,
            cs_macfb_builder_t          *macb,
            cs_cell_sys_t               *csys,
            cs_cell_builder_t           *cb)
{

  /* 1- Initialize common structures */

  cs_macfb_vecteq_init_build(
    connect, quant, eqp, eqb, c_id, vel_f_n, cm, macb, csys, cb);

  /* 1- SETUP THE NAVSTO LOCAL BUILDER
   * =================================
   * - Set the type of boundary
   * - Set the pressure boundary conditions (if required)
   * - Define the divergence operator used in the linear system (div_op is
   *   equal to minus the divergence)
   */

  cs_macfb_navsto_define_builder(
    cb->t_bc_eval, nsp, cm, macb, csys, sc->pressure_bc, sc->bf_type, nsb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes
 *         Common part between steady and unsteady case
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      sc           pointer to the scheme context
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      macb         pointer to a cs_macfb_builder_t structure
 * \param[in]      diff_pty     pointer to a cs_property_data_t structure
 *                              for diffusion
 * \param[in, out] nsb          pointer to a cs_macfb_navsto_builder_t structure
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_common_build(const cs_navsto_param_t     *nsp,
              const cs_macfb_monolithic_t *sc,
              const cs_equation_param_t   *eqp,
              const cs_macfb_vecteq_t     *eqc,
              cs_equation_builder_t       *eqb,
              const cs_cell_mesh_t        *cm,
              cs_macfb_builder_t          *macb,
              const cs_property_data_t    *diff_pty,
              const cs_real_t              time_scaling,
              cs_macfb_navsto_builder_t   *nsb,
              cs_cell_sys_t               *csys,
              cs_cell_builder_t           *cb)
{

  /* 2- VELOCITY (VECTORIAL) EQUATION */
  /* ================================ */

  cs_macfb_vecteq_conv_diff_reac(eqp, eqb, eqc, cm, macb, diff_pty, csys, cb);

  /* 3- SOUR_solve_monolithicCE TERM COMPUTATION (for the momentum equation) */
  /* ====================================================== */

  /* add source term (see cs_source_term.c) */
  if (cs_equation_param_has_sourceterm(eqp)) {
    cs_macfb_vecteq_sourceterm(
      cm, eqp, macb, cb->t_st_eval, time_scaling, eqb, cb, csys);
  }

  /* Gravity effects and/or Boussinesq approximation rely on another
     strategy than classical source term. The treatment is more compatible
     with the pressure gradient by doing so. */

  if (sc->add_gravity_term != nullptr) {
    sc->add_gravity_term(nsp, cm, nsb, csys);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system from a linear monolithic approach
 *         and update field
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] sc              pointer to a cs_macfb_monolithic_t strucutre
 */
/*----------------------------------------------------------------------------*/

static void
_solve_monolithic(const cs_navsto_param_t *nsp, cs_macfb_monolithic_t *sc)
{
  /* Retrieve high-level structures */

  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_builder_t  *eqb = eq->builder;

  /* Get fields */

  cs_real_t *u_f = eqc->face_values;
  cs_real_t *p_c = sc->pressure->val;

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);

  /* Solve the linear system */

  cs_timer_t t_solve_start = cs_timer_time();

  int iter = sc->solve(nsp, sc->saddle_solver, u_f, p_c);

  cs_timer_t t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(u_f);

  /* Now update the velocity and pressure fields associated to cells */

  _mono_update_related_cell_fields(nsp, sc, eqc);

  /* Compute the new mass flux used as the advection field */

  cs_macfb_navsto_mass_flux(nsp, cs_shared_quant, u_f, sc->mass_flux_array);

  if (nsp->verbosity > 1 && cs_log_default_is_active()) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- %s: iters: %d\n", __func__, iter);
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }

  /* Free */

  cs_saddle_solver_clean(sc->saddle_solver);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system from a nonlinear monolithic approach
 *         and update field
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] sc              pointer to a cs_macfb_monolithic_t strucutre
 */
/*----------------------------------------------------------------------------*/

static void
_solve_monolithic_nl(const cs_navsto_param_t *nsp, cs_macfb_monolithic_t *sc)
{
  /* Retrieve high-level structures */

  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_builder_t  *eqb = eq->builder;
  cs_iter_algo_t         *nl_algo = sc->nl_algo;

  /* Get fields */

  cs_real_t *u_f = eqc->face_values;
  cs_real_t *p_c = sc->pressure->val;

  /* Solve the linear system */

  cs_timer_t t_solve_start = cs_timer_time();

  cs_iter_algo_reset(nl_algo);

  /* Solve the new system: * Update the value of u_f and p_c */

  int last_inner_iter = sc->solve(nsp, sc->saddle_solver, u_f, p_c);

  cs_iter_algo_update_inner_iters(nl_algo, last_inner_iter);

  cs_timer_t t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(u_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the cell system.
 *          Case of MAC-Fb schemes with a monolithic velocity-pressure coupling
 *
 * \param[in]      sc        pointer to a cs_macfb_monolithic_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in]      cm        pointer to a cellwise view of the mesh
 * \param[in]      diff_pty  pointer to a \cs_property_data_t struct. for diff.
 * \param[in, out] csys      pointer to a cellwise view of the system
 * \param[in, out] cb        pointer to a cellwise builder
 * \param[in, out] nsb       builder structure for the NavSto system
 */
/*----------------------------------------------------------------------------*/

static void
_mono_apply_bc(const cs_macfb_monolithic_t *sc,
               const cs_equation_param_t   *eqp,
               const cs_equation_builder_t *eqb,
               const cs_cell_mesh_t        *cm,
               const cs_property_data_t    *diff_pty,
               cs_cell_sys_t               *csys,
               cs_cell_builder_t           *cb,
               cs_macfb_navsto_builder_t   *nsb)
{
  assert(cs_equation_param_has_diffusion(eqp) == true);

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    const cs_boundary_type_t *bf_type = nsb->bf_type;
    cs_real_t                *div_op  = nsb->div_op;
    cs_real_t *mass_rhs = sc->system_helper->rhs + cs_shared_quant->n_faces;

    /* Update the velocity-block and the right-hand side (part related to the
     * momentum equation). */

    /* Update the divergence operator and the right-hand side (part related to
     * the mass equation). */

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann) {
      for (short int i = 0; i < cm->n_fc; i++)
        csys->rhs[i] -= csys->neu_values[i];
    }

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */

      const short int f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {

        if (f < cm->n_fc) {

          /* Update mass RHS (constrain on the velocity divergence) from the
           * knowledge of the boundary face velocity */

          mass_rhs[cm->c_id] -= csys->dir_values[f] * div_op[f];

          /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */

          div_op[f] = 0;
        }

        /* Enforcement of the velocity DoFs for the velocity-block
         * Dirichlet BCs on the normal components of the velocity field. */

        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          sc->apply_velocity_inlet(f, eqp, cm, diff_pty, cb, csys);
        }
      }
      else if (bf_type[i] & CS_BOUNDARY_WALL) {

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */

        if (f < cm->n_fc) {
          div_op[f] = 0;
        }

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the normal components of the velocity field */

        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED
            || eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, diff_pty, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, diff_pty, cb, csys);
        }
      }
      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {

        /* Weak-enforcement for the velocity-block
         * Strong enforcement of u.n (--> dp/dn = 0) on the divergence */

        if (f < cm->n_fc) {
          div_op[f] = 0;
        }

        /* Always weakly enforce the symmetry constraint on the
           velocity-block */

        sc->apply_symmetry(f, eqp, cm, diff_pty, cb, csys);
      }
      else if (bf_type[i] & CS_BOUNDARY_IMPOSED_P) {

        /* Close the definition of the pressure gradient for this face */
        if (f < cm->n_fc) {
          csys->rhs[f] += div_op[f] * nsb->pressure_bc_val[i];
        }
        break;
      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */
  } /* Boundary cell */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

    cs_equation_builder_enforce_block_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_MONOLITHIC_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after BC enforcement", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with MAC-Fb schemes when the GKB or ALU algorithm is used as solver.
 *         Rely on cs_macfb_vecteq_assembly()
 *
 * \param[in]       csys              pointer to a cs_cell_sys_t structure
 * \param[in]       cm                pointer to a cs_cell_mesh_t structure
 * \param[in]       nsb               pointer to a navsto builder structure
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  asb               pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_velocity_full_assembly(const cs_cell_sys_t             *csys,
                        const cs_cell_mesh_t            *cm,
                        const cs_macfb_navsto_builder_t *nsb,
                        cs_macfb_monolithic_t           *sc,
                        cs_macfb_vecteq_t               *eqc,
                        cs_cdo_assembly_t               *asb)
{
  CS_UNUSED(eqc);

  const short int         n_f     = cm->n_fc;
  const cs_cdo_connect_t *connect = cs_shared_connect;

  cs_cdo_system_helper_t *sh = sc->system_helper;

  cs_real_t *_div = sc->block21_op + connect->c2f->idx[cm->c_id];

  /* 1. Store divergence operator in non assembly
   *    Take into account solid zone where DoF is set to zero */
  /* ======================================================== */

  if (csys->has_internal_enforcement) {

    for (int i = 0; i < n_f; i++) {

      if (csys->dof_is_forced[i]) {
        _div[i] = 0.; /* The velocity-block set the value of this DoF */
      }
      else {
        _div[i] = nsb->div_op[i];
      }
    }
  }
  else {
    memcpy(_div, nsb->div_op, n_f * sizeof(cs_real_t));
  }

  /* 1. Matrix assembly
   * ================== */

  const double gamma
    = cs_param_saddle_get_augmentation_coef(sc->saddle_solver->param);

  if (gamma > 0.) {
    cs_macfb_navsto_add_grad_div(cm->n_fc, gamma / cm->vol_c, _div, csys->mat);
  }

  cs_macfb_vecteq_assembly(csys, sh->blocks[0], sh->rhs, asb);

  /* 2. RHS assembly (mass eq. only the cell DoFs)
   * ============================================= */
  sh->rhs_array[1][cm->c_id] += nsb->mass_rhs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with MAC-Fb schemes
 *         Shares similarities with cs_cdo_assembly_block_matrix()
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb              pointer to a navsto builder structure
 * \param[in, out] sc               pointer to scheme context structure
 * \param[in, out] eqc              context structure for a vector-valued Fb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_full_system_assembly(const cs_cell_sys_t             *csys,
                      const cs_cell_mesh_t            *cm,
                      const cs_macfb_navsto_builder_t *nsb,
                      cs_macfb_monolithic_t           *sc,
                      cs_macfb_vecteq_t               *eqc,
                      cs_cdo_assembly_t               *asb)
{

  bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  CS_UNUSED(csys);
  CS_UNUSED(cm);
  CS_UNUSED(nsb);
  CS_UNUSED(sc);
  CS_UNUSED(eqc);
  CS_UNUSED(asb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         steady-state case
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_pre    velocity face DoFs of the previous time step
 * \param[in]      vel_f_nm1    nullptr (for unsteady computations)
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_navsto_param_t *nsp,
              const cs_real_t          vel_f_pre[],
              const cs_real_t          vel_f_nm1[],
              cs_macfb_monolithic_t   *sc)
{
  CS_NO_WARN_IF_UNUSED(vel_f_nm1);

  /* Retrieve shared structures */

  const cs_cdo_connect_t    *connect = cs_shared_connect;
  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_time_step_t      *ts      = cs_shared_time_step;

  /* Retrieve high-level structures */

  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(eqb->dir_values != nullptr);
#endif

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    const int t_id = cs_get_thread_id();

    cs_macfb_navsto_builder_t nsb
      = cs_macfb_navsto_create_builder(nsp, connect);
    cs_cell_mesh_t     *cm   = cs_cdo_local_get_cell_mesh(t_id);
    cs_macfb_builder_t *macb = cs_macfb_get_builder(t_id);
    cs_cdo_assembly_t  *asb  = cs_cdo_assembly_get(t_id);

    /* Hodge operator used only to compute diffusion property */
    cs_hodge_t *diff_hodge = (eqc->diffusion_hodge == nullptr)
                               ? nullptr
                               : eqc->diffusion_hodge[t_id];

    cs_cell_sys_t     *csys = nullptr;
    cs_cell_builder_t *cb   = nullptr;

    cs_macfb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = ts->t_cur; /* Dummy parameter if really steady */
    cb->t_bc_eval  = ts->t_cur; /* Dummy parameter if really steady */
    cb->t_st_eval  = ts->t_cur; /* Dummy parameter if really steady */

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    const cs_property_data_t *diff_pty = diff_hodge->pty_data;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* 1- Initialize structures */

      _init_build(nsp,
                  sc,
                  connect,
                  quant,
                  eqp,
                  eqb,
                  c_id,
                  vel_f_pre,
                  &nsb,
                  cm,
                  macb,
                  csys,
                  cb);

      /* 2- Compute steady part */

      _common_build(
        nsp, sc, eqp, eqc, eqb, cm, macb, diff_pty, 1.0, &nsb, csys, cb);

      /* 3- Apply boundaries condition */

      _mono_apply_bc(sc, eqp, eqb, cm, diff_pty, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> MAC-fb mono: (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, &nsb, sc, eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_macfb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  cs_cdo_system_helper_finalize_assembly(sc->system_helper);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of an implicit Euler time scheme
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_n      velocity face DoFs at time step n
 * \param[in]      vel_f_nm1    nullptr (not needed for this time scheme)
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_implicit_euler_build(const cs_navsto_param_t *nsp,
                      const cs_real_t          vel_f_n[],
                      const cs_real_t          vel_f_nm1[],
                      cs_macfb_monolithic_t   *sc)
{
  CS_NO_WARN_IF_UNUSED(vel_f_nm1);

  /* Retrieve shared structures */

  const cs_cdo_connect_t    *connect = cs_shared_connect;
  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_time_step_t      *ts      = cs_shared_time_step;

  /* Retrieve high-level structures */

  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;

  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(eqb->dir_values != nullptr);
#endif

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    const int t_id = cs_get_thread_id();

    cs_macfb_navsto_builder_t nsb
      = cs_macfb_navsto_create_builder(nsp, connect);
    cs_cell_mesh_t     *cm   = cs_cdo_local_get_cell_mesh(t_id);
    cs_macfb_builder_t *macb = cs_macfb_get_builder(t_id);
    cs_cdo_assembly_t  *asb  = cs_cdo_assembly_get(t_id);

    /* Hodge operator used only to compute diffusion property */
    cs_hodge_t *diff_hodge = (eqc->diffusion_hodge == nullptr)
                               ? nullptr
                               : eqc->diffusion_hodge[t_id];

    cs_cell_sys_t     *csys = nullptr;
    cs_cell_builder_t *cb   = nullptr;

    cs_macfb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    const cs_real_t t_eval = ts->t_cur + ts->dt[0];

    cb->t_pty_eval = t_eval;
    cb->t_bc_eval  = t_eval;
    cb->t_st_eval  = t_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    const cs_property_data_t *diff_pty = diff_hodge->pty_data;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* 1- Initialize structures */

      _init_build(nsp,
                  sc,
                  connect,
                  quant,
                  eqp,
                  eqb,
                  c_id,
                  vel_f_n,
                  &nsb,
                  cm,
                  macb,
                  csys,
                  cb);

      /* 2- Compute steady part */

      _common_build(
        nsp, sc, eqp, eqc, eqb, cm, macb, diff_pty, 1.0, &nsb, csys, cb);

      /* 3- compute unsteady terms */

      cs_macfb_vecteq_euler_implicit_term(eqp, cm, macb, cb, ts->dt[0], csys);

      /* 4- Apply boundaries condition */

      _mono_apply_bc(sc, eqp, eqb, cm, diff_pty, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> MAC-fb mono: (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, &nsb, sc, eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_macfb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  cs_cdo_system_helper_finalize_assembly(sc->system_helper);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of a theta time scheme
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_n      velocity face DoFs at time step n
 * \param[in]      vel_f_nm1    velocity face DoFs at time step n-1 or nullptr
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_theta_scheme_build(const cs_navsto_param_t *nsp,
                    const cs_real_t          vel_f_n[],
                    const cs_real_t          vel_f_nm1[],
                    cs_macfb_monolithic_t   *sc)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);
  CS_UNUSED(nsp);
  CS_UNUSED(vel_f_n);
  CS_UNUSED(vel_f_nm1);
  CS_UNUSED(sc);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in] eqp    equation parameter settings
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] quant      additional mesh quantities struct.
 * \param[in] connect    pointer to a \ref cs_cdo_connect_t struct.
 * \param[in] time_step  pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_monolithic_init_sharing(const cs_equation_param_t *eqp,
                                 const cs_mesh_t           *mesh,
                                 const cs_cdo_quantities_t *quant,
                                 const cs_cdo_connect_t    *connect,
                                 const cs_time_step_t      *time_step)
{
  assert(eqp->saddle_param != nullptr);
  assert(eqp->saddle_param->solver != CS_PARAM_SADDLE_SOLVER_NONE);

  /* Assign static const pointers */

  cs_shared_mesh      = mesh;
  cs_shared_quant     = quant;
  cs_shared_connect   = connect;
  cs_shared_time_step = time_step;

  /* Need to build special range set and interfaces ? */

  cs_param_sles_t *slesp = eqp->sles_param;

  if (slesp->precond_block_type != CS_PARAM_PRECOND_BLOCK_NONE) {
    if (slesp->solver_class != CS_PARAM_SOLVER_CLASS_PETSC) {
      bft_error(__FILE__, __LINE__, 0, " %s: Not implemented.\n", __func__);
    }
  }

  /* SLES needs these structures for advanced PETSc hooks */

  cs_macfb_monolithic_sles_init_sharing(mesh, connect, quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free shared pointers with lifecycle dedicated to this file
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_monolithic_finalize_common(void) {};

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_macfb_monolithic_t structure
 *
 * \param[in] nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in] adv_field    pointer to \ref cs_adv_field_t structure
 * \param[in] mflux        current values of the mass flux across primal faces
 * \param[in] mflux_pre    current values of the mass flux across primal faces
 * \param[in] bf_type      type of boundary for each boundary face
 * \param[in] cc_context   pointer to a \ref cs_navsto_monolithic_t structure
 *
 * \return a pointer to a new allocated \ref cs_macfb_monolithic_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_macfb_monolithic_init_scheme_context(const cs_navsto_param_t *nsp,
                                        cs_adv_field_t          *adv_field,
                                        cs_real_t               *mflux,
                                        cs_real_t               *mflux_pre,
                                        cs_boundary_type_t      *bf_type,
                                        void                    *cc_context)
{
  assert(nsp != nullptr && cc_context != nullptr);
  if (nsp->space_scheme != CS_SPACE_SCHEME_MACFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_cdo_connect_t    *connect = cs_shared_connect;

  /* Navier-Stokes scheme context (SC) */

  cs_macfb_monolithic_t *sc = nullptr;

  BFT_MALLOC(sc, 1, cs_macfb_monolithic_t);

  /* Cast the coupling context (CC) */

  cs_navsto_monolithic_t *cc      = (cs_navsto_monolithic_t *)cc_context;
  cs_equation_t          *eq      = cc->momentum;
  cs_equation_param_t    *eqp     = eq->param;
  cs_equation_builder_t  *eqb     = eq->builder;

  /* Quantities shared with the cs_navsto_system_t structure */

  sc->coupling_context    = cc;
  sc->adv_field           = adv_field;
  sc->mass_flux_array     = mflux;
  sc->mass_flux_array_pre = mflux_pre;

  /* Quick access to the main fields */

  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE)
    sc->divergence = cs_field_by_name("velocity_divergence");
  else
    sc->divergence = nullptr;

  /* Boundary treatment */

  sc->bf_type = bf_type;

  /* Processing of the pressure boundary condition */

  sc->pressure_bc = cs_cdo_bc_face_define(CS_BC_SYMMETRY, /* Default */
                                          true, /* Steady BC up to now */
                                          1,    /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          quant->n_b_faces);

  sc->pressure_rescaling
    = cs_boundary_need_pressure_rescaling(quant->n_b_faces, bf_type);

  /* Set the way to enforce the Dirichlet BC on the velocity
   * "fixed_wall" means a no-slip BC */

  eqb->bdy_flag |= CS_FLAG_COMP_PFC;

  sc->apply_symmetry     = cs_macfb_symmetry;
  sc->apply_sliding_wall = cs_macfb_block_dirichlet_alge;
  sc->apply_fixed_wall   = cs_macfb_block_dirichlet_alge;

  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    sc->apply_velocity_inlet = cs_macfb_block_dirichlet_alge;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    sc->apply_velocity_inlet = cs_macfb_block_dirichlet_pena;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    sc->apply_velocity_inlet = cs_macfb_block_dirichlet_weak;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    sc->apply_velocity_inlet = cs_macfb_block_dirichlet_wsym;
    break;

  default:
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);
  }

  /* Source terms induced by the gravity effect */

  cs_macfb_navsto_set_gravity_func(nsp, &(sc->add_gravity_term));

  /* Set the build function */

  sc->steady_build = _steady_build;

  switch (eqp->time_scheme) {

  case CS_TIME_SCHEME_STEADY:
    sc->build = _steady_build;
    break;

  case CS_TIME_SCHEME_EULER_IMPLICIT:
    sc->build = _implicit_euler_build;
    break;

  case CS_TIME_SCHEME_EULER_EXPLICIT:
  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
    sc->build = _theta_scheme_build;
    break;

  case CS_TIME_SCHEME_BDF2:
  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid time scheme.", __func__);

  } /* Switch on time schme */

  /* Linear algebra */
  /* -------------- */

  const cs_param_saddle_t *saddlep = eqp->saddle_param;

  /* Some saddle-point solver needs the (2,1)-block stored in an unassembled
     way. This corresponds to the -|c|.divergence operator */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_GKB:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    BFT_MALLOC(sc->block21_op, connect->c2f->idx[quant->n_cells], cs_real_t);
    break;

  default:
    /* Nothing to do */
    sc->block21_op = nullptr;
    break;
  }

  /* Define the layout of the system and how to assemble the system. It depends
    on the strategy to solve the saddle-point problem */

  cs_macfb_monolithic_sles_init_system_helper(nsp, saddlep, sc);

  /* Set the function pointer to assemble the linear system */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_GKB:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    sc->assemble = _velocity_full_assembly;
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM: /* Experimental */
    bft_error(__FILE__, __LINE__, 0, "%s: NOTAY is not implemented.", __func__);
    break;

  default:
    /* CS_PARAM_SADDLE_SOLVER_FGMRES
     * CS_PARAM_SADDLE_SOLVER_MUMPS */
    if (nsp->model_flag & CS_NAVSTO_MODEL_WITH_SOLIDIFICATION) {
      bft_error(__FILE__,
                __LINE__,
                0,
                "%s: Solidification is not implemented.",
                __func__);
    }
    else {
      sc->assemble = _full_system_assembly;
    }
    break;

  } /* Switch on saddle-point solver */

  /* Handle the resolution of a saddle-point system */

  cs_macfb_monolithic_sles_init_solver(nsp, saddlep, sc);

  /* Iterative algorithm to handle the non-linearity (Picard by default) */

  cs_iter_algo_type_t algo_type = CS_ITER_ALGO_TWO_LEVEL;

  if (nsp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    algo_type |= CS_ITER_ALGO_ANDERSON;
  else
    algo_type |= CS_ITER_ALGO_DEFAULT;

  sc->nl_algo = cs_iter_algo_create_with_settings(
    algo_type, nsp->verbosity, nsp->nl_cvg_param);

  if (nsp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_set_anderson_param(
      sc->nl_algo, nsp->anderson_param, quant->n_faces);

  /* Monitoring */

  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_macfb_monolithic_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_macfb_monolithic_free_scheme_context(void *scheme_context)
{
  cs_macfb_monolithic_t *sc = (cs_macfb_monolithic_t *)scheme_context;

  if (sc == nullptr)
    return sc;

  /* Free BC structure */

  sc->pressure_bc = cs_cdo_bc_free(sc->pressure_bc);

  /* Free the systemm helper */

  cs_cdo_system_helper_free(&(sc->system_helper));

  /* Block (2,1) may be allocated */

  BFT_FREE(sc->block21_op);

  /* Free the context structure for solving saddle-point system */

  cs_saddle_solver_free(&(sc->saddle_solver));

  /* Free the structure handling the non-linear algorithm */

  cs_iter_algo_free(&(sc->nl_algo));

  /* Other pointers are only shared (i.e. not owner) */

  BFT_FREE(sc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Stokes or Oseen system with a MAC face-based scheme
 *         using a monolithic approach
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_monolithic_steady(const cs_mesh_t         *mesh,
                           const cs_navsto_param_t *nsp,
                           void                    *scheme_context)
{
  cs_timer_t t_start = cs_timer_time();

  /* Retrieve high-level structures */

  cs_macfb_monolithic_t  *sc = (cs_macfb_monolithic_t *)scheme_context;
  cs_cdo_system_helper_t *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t      *ts    = cs_shared_time_step;
  const cs_real_t            t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces and DoF enforcement */

  cs_macfb_vecteq_setup(t_cur, mesh, eqp, eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t *rhs = nullptr; /* Since it is nullptr, sh get the ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  sc->steady_build(nsp,
                   eqc->face_values,
                   nullptr, /* no value at time step n-1 */
                   sc);

  /* Free temporary buffers and structures */

  cs_equation_builder_reset(eqb);

  /* End of the system building */

  cs_timer_t t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Solve linear system and update fields */

  _solve_monolithic(nsp, sc);

  /* Frees */

  cs_equation_builder_reset(eqb);
  cs_cdo_system_helper_reset(sh); /* free rhs and matrix */

  cs_timer_t t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the steady Navier-Stokes system with a MAC face-based scheme
 *        using a monolithic approach and a non-linear algorithm (Picard or
 *        Anderson) to solve the non-linearities arising from the advection
 *        term
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_monolithic_steady_nl(const cs_mesh_t         *mesh,
                              const cs_navsto_param_t *nsp,
                              void                    *scheme_context)
{
  cs_timer_t t_start = cs_timer_time();

  /* Retrieve high-level structures */

  cs_macfb_monolithic_t  *sc = (cs_macfb_monolithic_t *)scheme_context;
  cs_cdo_system_helper_t *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;
  cs_iter_algo_t         *nl_algo = sc->nl_algo;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_lnum_t            n_faces = quant->n_faces;
  const cs_time_step_t      *ts      = cs_shared_time_step;
  const cs_real_t            t_cur   = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces and DoF enforcement */

  cs_macfb_vecteq_setup(t_cur, mesh, eqp, eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t *rhs = nullptr; /* Since it is nullptr, sh get the ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  cs_real_t *u_f = eqc->face_values;

  sc->steady_build(nsp,
                   u_f,
                   nullptr, /* no value at time step n-1 */
                   sc);

  /* End of the system building */

  cs_timer_t t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);
  cs_iter_algo_reset(nl_algo);

  /* Solve the linear system */

  _solve_monolithic_nl(nsp, sc);

  /* Compute the new current mass flux used as the advection field */

  cs_macfb_navsto_mass_flux(nsp, quant, u_f, sc->mass_flux_array);

  /* Set the normalization of the non-linear algo to the value of the first
     mass flux norm */

  double normalization
    = sqrt(cs_cdo_blas_square_norm_pfsf(sc->mass_flux_array));

  cs_iter_algo_set_normalization(nl_algo, normalization);

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  /* Check the convergence status and update the nl_algo structure related
   * to the convergence monitoring */

  while (
    cs_macfb_navsto_nl_algo_cvg(
      nsp->nl_algo_type, sc->mass_flux_array_pre, sc->mass_flux_array, nl_algo)
    == CS_SLES_ITERATING) {

    /* Main loop on cells to define the linear system to solve */

    cs_timer_t t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */

    cs_cdo_system_helper_init_system(sh, &rhs);
    cs_saddle_solver_clean(sc->saddle_solver);

    sc->steady_build(nsp,
                     /* A current to previous op. has been done */
                     eqc->face_values_pre,
                     nullptr, /* no value at time step n-1 */
                     sc);

    /* End of the system building */

    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(eqb->tcb), &t_build_start, &t_build_end);

    /* Solve the linear system */

    _solve_monolithic_nl(nsp, sc);

    /* Compute the new mass flux used as the advection field */

    cs_array_real_copy(n_faces, sc->mass_flux_array, sc->mass_flux_array_pre);

    cs_macfb_navsto_mass_flux(nsp, quant, u_f, sc->mass_flux_array);

  } /* Loop on non-linear iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1 && cs_log_default_is_active()) {

    cs_log_printf(CS_LOG_DEFAULT,
                  " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  cs_iter_algo_get_n_inner_iter(nl_algo));
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }

  cs_iter_algo_check_warning(__func__,
                             eqp->name,
                             cs_param_get_nl_algo_label(nsp->nl_algo_type),
                             nl_algo);

  if (nsp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_release_anderson_arrays(
      (cs_iter_algo_aac_t *)nl_algo->context);

  /* Now compute/update the velocity and pressure fields */

  _mono_update_related_cell_fields(nsp, sc, eqc);

  /* Frees */

  cs_saddle_solver_clean(sc->saddle_solver);
  cs_cdo_system_helper_reset(sh); /* free rhs and matrix */

  cs_timer_t t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a MAC face-based scheme
 *         using a monolithic approach.
 *         According to the settings, this function can handle either an
 *         implicit Euler time scheme or more generally a theta time scheme.
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_monolithic(const cs_mesh_t         *mesh,
                    const cs_navsto_param_t *nsp,
                    void                    *scheme_context)
{
  const cs_timer_t t_start = cs_timer_time();

  /* Retrieve high-level structures */

  cs_macfb_monolithic_t  *sc = (cs_macfb_monolithic_t *)scheme_context;
  cs_cdo_system_helper_t *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t      *ts     = cs_shared_time_step;
  const cs_real_t            t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces and DoF enforcement */

  cs_macfb_vecteq_setup(t_eval, mesh, eqp, eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t *rhs = nullptr; /* Since it is nullptr, sh get sthe ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  cs_real_t *u_f     = eqc->face_values;
  cs_real_t *u_f_pre = eqc->face_values_pre;

  sc->build(nsp, u_f, u_f_pre, sc);

  /* End of the system building */

  cs_timer_t t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Solve linear system and update fields */

  _solve_monolithic(nsp, sc);

  /* Frees */

  cs_equation_builder_reset(eqb);
  cs_cdo_system_helper_reset(sh); /* free rhs and matrix */

  cs_timer_t t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a MAC face-based scheme
 *         using a monolithic approach.
 *         According to the settings, this function can handle either an
 *         implicit Euler time scheme or more generally a theta time scheme.
 *         Rely on Picard iterations to solve the non-linearities arising from
 *         the advection term
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_monolithic_nl(const cs_mesh_t         *mesh,
                       const cs_navsto_param_t *nsp,
                       void                    *scheme_context)
{
  cs_timer_t t_start = cs_timer_time();

  /* Retrieve high-level structures */

  cs_macfb_monolithic_t  *sc = (cs_macfb_monolithic_t *)scheme_context;
  cs_cdo_system_helper_t *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t          *eq  = cc->momentum;
  cs_macfb_vecteq_t      *eqc = (cs_macfb_vecteq_t *)eq->scheme_context;
  cs_equation_param_t    *eqp = eq->param;
  cs_equation_builder_t  *eqb = eq->builder;
  cs_iter_algo_t         *nl_algo = sc->nl_algo;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t *quant   = cs_shared_quant;
  const cs_lnum_t            n_faces = quant->n_faces;
  const cs_time_step_t      *ts      = cs_shared_time_step;
  const cs_real_t            t_eval  = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces and DoF enforcement */

  cs_macfb_vecteq_setup(t_eval, mesh, eqp, eqb);

  /* Initialize the rhs */

  cs_real_t *rhs = nullptr;

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  cs_real_t *u_f     = eqc->face_values;     /* cur.  velocity at faces */
  cs_real_t *u_f_pre = eqc->face_values_pre; /* prev. velocity at faces */

  sc->build(nsp, u_f, u_f_pre, sc);

  /* End of the system building */

  cs_timer_t t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);
  cs_iter_algo_reset(nl_algo);

  /* Solve the linear system */

  _solve_monolithic_nl(nsp, sc);

  /* Compute the new mass flux used as the advection field */

  cs_macfb_navsto_mass_flux(nsp, quant, u_f, sc->mass_flux_array);

  /* Set the normalization of the non-linear algo to the value of the first
     mass flux norm */

  double normalization
    = sqrt(cs_cdo_blas_square_norm_pfsf(sc->mass_flux_array));

  cs_iter_algo_set_normalization(nl_algo, normalization);

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  /* Since a current to previous op. has been done:
   *   sc->mass_flux_array_pre -> flux at t^n= t^n,0 (= t^(n-1)
   *   sc->mass_flux_array     -> flux at t^n,1 (call to .._navsto_mass_flux */

  cs_sles_convergence_state_t cvg_status = cs_macfb_navsto_nl_algo_cvg(
    nsp->nl_algo_type, sc->mass_flux_array_pre, sc->mass_flux_array, nl_algo);

  cs_real_t *mass_flux_array_k   = nullptr;
  cs_real_t *mass_flux_array_kp1 = sc->mass_flux_array;

  while (cvg_status == CS_SLES_ITERATING) {

    /* Start of the system building */

    cs_timer_t t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */

    cs_cdo_system_helper_init_system(sh, &rhs);
    cs_saddle_solver_clean(sc->saddle_solver);

    /* Main loop on cells to define the linear system to solve */

    sc->build(nsp,
              /* A current to previous op. has been done */
              u_f_pre,
              nullptr, /* no n-1 state is given */
              sc);

    /* End of the system building */

    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(eqb->tcb), &t_build_start, &t_build_end);

    /* Solve the linear system */

    _solve_monolithic_nl(nsp, sc);

    /* mass_flux_array_k <-- mass_flux_array_kp1; update mass_flux_array_kp1 */

    if (mass_flux_array_k == nullptr)
      BFT_MALLOC(mass_flux_array_k, n_faces, cs_real_t);
    cs_array_real_copy(n_faces, mass_flux_array_kp1, mass_flux_array_k);

    cs_macfb_navsto_mass_flux(nsp, quant, u_f, mass_flux_array_kp1);

    /* Check the convergence status and update the nl_algo structure related
     * to the convergence monitoring */

    cvg_status = cs_macfb_navsto_nl_algo_cvg(
      nsp->nl_algo_type, mass_flux_array_k, mass_flux_array_kp1, nl_algo);

  } /* Loop on non-linear iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1 && cs_log_default_is_active()) {

    cs_log_printf(CS_LOG_DEFAULT,
                  " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  cs_iter_algo_get_n_inner_iter(nl_algo));
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }

  cs_iter_algo_check_warning(__func__,
                             eqp->name,
                             cs_param_get_nl_algo_label(nsp->nl_algo_type),
                             nl_algo);

  if (nsp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_release_anderson_arrays(
      (cs_iter_algo_aac_t *)nl_algo->context);

  /* Now compute/update the velocity and pressure fields */

  _mono_update_related_cell_fields(nsp, sc, eqc);

  /* Frees */

  cs_saddle_solver_clean(sc->saddle_solver);
  cs_cdo_system_helper_reset(sh); /* free rhs and matrix */
  cs_equation_builder_reset(eqb);
  if (mass_flux_array_k != nullptr)
    BFT_FREE(mass_flux_array_k);

  cs_timer_t t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
