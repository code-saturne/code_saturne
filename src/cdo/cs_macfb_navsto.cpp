/*============================================================================
 * Functions shared for MAC face-based schemes for the discretization of the
 * Navier-Stokes system
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_blas.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_blas.h"
#include "cs_cdo_toolbox.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_cdofb_navsto.h"
#include "cs_equation_bc.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_macfb_builder.h"
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_reco.h"
#include "cs_sdm.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_navsto.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_macfb_navsto.c
 *
 * \brief Shared functions among all face-based schemes for building and
 *        solving Stokes and Navier-Stokes problem
 *
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_MACFB_NAVSTO_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */

#define _dp3 cs_math_3_dot_product

/*============================================================================
 * Private variables
 *============================================================================*/

static cs_macfb_navsto_boussinesq_type_t cs_macfb_navsto_boussinesq_type
  = CS_CDOFB_NAVSTO_BOUSSINESQ_FACE_DOF;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the way to compute the Boussinesq approximation
 *
 * \param[in] type     type of algorithm to use
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_set_boussinesq_algo(cs_macfb_navsto_boussinesq_type_t type)
{
  cs_macfb_navsto_boussinesq_type = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and allocate a local NavSto builder when MAC-Fb schemes are
 * used
 *
 * \param[in] nsp         set of parameters to define the NavSto system
 * \param[in] connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a cs_cdofb_navsto_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_navsto_builder_t
cs_macfb_navsto_create_builder(const cs_navsto_param_t *nsp,
                               const cs_cdo_connect_t  *connect)
{
  cs_macfb_navsto_builder_t nsb = { .rho_c           = 1.,
                                    .div_op          = nullptr,
                                    .mass_rhs        = 0.,
                                    .bf_type         = nullptr,
                                    .pressure_bc_val = nullptr };

  assert(nsp != nullptr);
  nsb.rho_c = nsp->mass_density->ref_value;

  if (connect == nullptr)
    return nsb;

  /* Number of faces (with full stencil) linked to a cell*/
  const int n_max_fbyc_x = 30;

  BFT_MALLOC(nsb.div_op, connect->n_max_fbyc, cs_real_t);
  BFT_MALLOC(nsb.bf_type, n_max_fbyc_x, cs_boundary_type_t);
  BFT_MALLOC(nsb.pressure_bc_val, n_max_fbyc_x, cs_real_t);

  return nsb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy the given cs_macfb_navsto_builder_t structure
 *
 * \param[in, out] nsb   pointer to the cs_macfb_navsto_builder_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_free_builder(cs_macfb_navsto_builder_t *nsb)
{
  if (nsb != nullptr) {
    BFT_FREE(nsb->div_op);
    BFT_FREE(nsb->bf_type);
    BFT_FREE(nsb->pressure_bc_val);
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_NAVSTO_DBG > 3
  cs_log_printf(CS_LOG_DEFAULT, ">> Free Navsto builder.\n");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the members of the cs_macfb_navsto_builder_t structure
 *
 * \param[in]      t_eval     time at which one evaluates the pressure BC
 * \param[in]      nsp        set of parameters to define the NavSto system
 * \param[in]      cm         cellwise view of the mesh
 * \param[in]      macb       macfb builder
 * \param[in]      csys       cellwise view of the algebraic system
 * \param[in]      pr_bc      set of definitions for the presuure BCs
 * \param[in]      bf_type    type of boundaries for all boundary faces
 * \param[in, out] nsb        builder to update
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_define_builder(cs_real_t                  t_eval,
                               const cs_navsto_param_t   *nsp,
                               const cs_cell_mesh_t      *cm,
                               const cs_macfb_builder_t  *macb,
                               const cs_cell_sys_t       *csys,
                               const cs_cdo_bc_face_t    *pr_bc,
                               const cs_boundary_type_t  *bf_type,
                               cs_macfb_navsto_builder_t *nsb)
{
  assert(cm != nullptr && macb != nullptr && csys != nullptr
         && nsp != nullptr); /* sanity checks */

  nsb->mass_rhs = 0; /* Reset the mass rhs */

  /* Update the value of the mass density for the current cell if needed */
  /* TODO: Case of a uniform but not constant in time */

  if (!cs_property_is_uniform(nsp->mass_density))
    nsb->rho_c = cs_property_value_in_cell(cm, nsp->mass_density, t_eval);

  /* Build the divergence operator:
   *        Div(u) = \frac{1}{|c|} \sum_{f_c} |fc| u_f nu_fc.e^i_f
   * We return the opposite -Div(u) because of the mass conservation
   * -div(u) = g
   */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t pfq = cm->face[f];

    nsb->div_op[f] = -(pfq.meas / cm->vol_c) * macb->f_sgn_axis[f];

  } /* Loop on cell faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_MACFB_NAVSTO_DBG > 2
  if (cs_dbg_cw_test(nullptr, cm, csys)) {
#pragma omp critical
    {
      cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
      for (short int f = 0; f < cm->n_fc; f++)
        cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e\n", f, nsb->div_op[f]);
    } /* Critical section */
  }
#endif

  /* Build local arrays related to the boundary conditions */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering and the boundary face id in
       the mesh numbering */

    const short int f     = csys->_f_ids[i];
    const cs_lnum_t bf_id = macb->f_ids[f] - cm->bface_shift;

    /* Set the type of boundary */

    nsb->bf_type[i] = bf_type[bf_id];

    /* Set the pressure BC if required */

    if (nsb->bf_type[i] & CS_BOUNDARY_IMPOSED_P) {

      assert(nsb->bf_type[i] & (CS_BOUNDARY_INLET | CS_BOUNDARY_OUTLET));

      /* Add a Dirichlet for the pressure field */

      const short int  def_id = pr_bc->def_ids[bf_id];
      const cs_xdef_t *def    = nsp->pressure_bc_defs[def_id];
      assert(pr_bc != nullptr);

      switch (def->type) {
      case CS_XDEF_BY_VALUE: {
        const cs_real_t *constant_val = (cs_real_t *)def->context;
        nsb->pressure_bc_val[i]       = constant_val[0];
      } break;

      case CS_XDEF_BY_ARRAY: {
        cs_xdef_array_context_t *c
          = static_cast<cs_xdef_array_context_t *>(def->context);
        assert(c->stride == 1);
        assert(cs_flag_test(c->value_location, cs_flag_primal_face));

        nsb->pressure_bc_val[i] = c->values[bf_id];
      } break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        switch (nsp->dof_reduction_mode) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_cw_eval_at_xyz_by_analytic(cm,
                                             1,
                                             cm->face[f].center,
                                             t_eval,
                                             def->context,
                                             nsb->pressure_bc_val + i);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_cw_eval_scalar_face_avg_by_analytic(
            cm, f, t_eval, def->context, def->qtype, nsb->pressure_bc_val + i);
          break;

        default:
          bft_error(__FILE__,
                    __LINE__,
                    0,
                    _(" %s: Invalid type of reduction.\n"
                      " Stop computing the Dirichlet value.\n"),
                    __func__);

        } /* switch on reduction */
        break;

      default:
        bft_error(__FILE__,
                  __LINE__,
                  0,
                  _(" %s: Invalid type of definition.\n"
                    " Stop computing the Dirichlet value.\n"),
                  __func__);
        break;

      } /* def->type */
    }
    else
      nsb->pressure_bc_val[i] = 0.;

  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mass flux playing the role of the advection field in
 *         the Navier-Stokes equations
 *         One considers the mass flux across primal faces which relies on the
 *         velocity vector defined on each face.
 *
 * \param[in]      nsp         set of parameters to define the NavSto system
 * \param[in]      quant       set of additional geometrical quantities
 * \param[in]      face_vel    velocity vectors for each face
 * \param[in, out] mass_flux   array of mass flux values to update (allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_mass_flux(const cs_navsto_param_t   *nsp,
                          const cs_cdo_quantities_t *quant,
                          const cs_real_t           *face_vel,
                          cs_real_t                 *mass_flux)
{
  if (mass_flux == nullptr)
    return;

  assert(face_vel != nullptr);
  assert(nsp->space_scheme == CS_SPACE_SCHEME_MACFB);
  assert(cs_property_is_uniform(nsp->mass_density));
  assert(nsp->mass_density->n_definitions == 1);

  const cs_real_t rho_val = nsp->mass_density->ref_value;

  /* Define the mass flux = rho * |f| * u_f */

#pragma omp parallel for if (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++) {
    mass_flux[f_id]
      = rho_val * cs_quant_get_face_surf(f_id, quant) * face_vel[f_id];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence in a cell of a vector-valued array defined at
 *         faces (values are defined both at interior and border faces).
 *         Variant based on the usage of \ref cs_cdo_quantities_t structure.
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

double
cs_macfb_navsto_cell_divergence(const cs_lnum_t            c_id,
                                const cs_cdo_quantities_t *quant,
                                const cs_adjacency_t      *c2f,
                                const cs_real_t           *f_vals)
{
  /* div(u) = 1/|c| * sum_f |f| * u_f * sign(f) * n.e^f */

  double div = 0.0;
  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id + 1]; f++) {

    const cs_lnum_t  f_id = c2f->ids[f];
    const cs_real_t *fq   = cs_quant_get_face_vector_area(f_id, quant);

    div += c2f->sgn[f] * f_vals[f_id] * fq[quant->face_axis[f_id]];

  } /* Loop on cell faces */

  div /= quant->cell_vol[c_id];

  return div;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the pressure at faces
 *
 * \param[in]       mesh      pointer to a cs_mesh_t structure
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]       ts        pointer to a cs_time_step_t structure
 * \param[in]       nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]       p_cell    value of the pressure inside each cell
 * \param[in, out]  p_face    value of the pressure at each face
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_compute_face_pressure(const cs_mesh_t           *mesh,
                                      const cs_cdo_connect_t    *connect,
                                      const cs_cdo_quantities_t *quant,
                                      const cs_time_step_t      *ts,
                                      const cs_navsto_param_t   *nsp,
                                      const cs_real_t           *p_cell,
                                      cs_real_t                 *p_face)
{
  cs_cdofb_navsto_compute_face_pressure(
    mesh, connect, quant, ts, nsp, p_cell, p_face);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the velocity at cells
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]       mass_flux scalar-valued mass flux for each face
 * \param[in, out]  u_cell    vector-value velocity velo at each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_compute_cell_velocity(const cs_cdo_connect_t    *connect,
                                      const cs_cdo_quantities_t *quant,
                                      const cs_real_t           *mass_flux,
                                      cs_real_t                 *u_cell)
{
  cs_reco_cell_vectors_by_face_dofs(connect->c2f, quant, mass_flux, u_cell);
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

void
cs_macfb_navsto_add_grad_div(short int       n_fc,
                             const cs_real_t zeta,
                             const cs_real_t div[],
                             cs_sdm_t       *mat)
{
  assert(n_fc == mat->n_rows && n_fc <= mat->n_cols);

  for (short int bi = 0; bi < n_fc; bi++) {

    const cs_real_t divi  = div[bi];
    const cs_real_t zt_di = zeta * divi;

    /* Begin with the diagonal block */

    const cs_real_t gd_coef_ii = zt_di * divi;

    mat->val[bi * mat->n_cols + bi] += gd_coef_ii;

    /* Continue with the extra-diag. blocks */

    for (short int bj = bi + 1; bj < n_fc; bj++) {

      const cs_real_t divj       = div[bj];
      const cs_real_t gd_coef_ij = zt_di * divj;

      /* Extra-diagonal: Use the symmetry of the grad-div */

      mat->val[bi * mat->n_cols + bj] += gd_coef_ij;
      mat->val[bj * mat->n_cols + bi] += gd_coef_ij;

    } /* Loop on column: bj */
  }   /* Loop on row: bi */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values
 *
 * \param[in]       nsp     pointer to a \ref cs_navsto_param_t structure
 * \param[in]       quant   pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]       ts      pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr      pointer to the pressure \ref cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_init_pressure(const cs_navsto_param_t   *nsp,
                              const cs_cdo_quantities_t *quant,
                              const cs_time_step_t      *ts,
                              cs_field_t                *pr)
{
  cs_cdofb_navsto_init_pressure(nsp, quant, ts, pr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values when the pressure is defined at
 *         faces
 *
 * \param[in]       nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]       connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]       ts        pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr_f      pointer to the pressure values at faces
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_init_face_pressure(const cs_navsto_param_t *nsp,
                                   const cs_cdo_connect_t  *connect,
                                   const cs_time_step_t    *ts,
                                   cs_real_t               *pr_f)
{
  cs_cdofb_navsto_init_face_pressure(nsp, connect, ts, pr_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a mean-value
 *         equal to the reference value
 *
 * \param[in]       nsp       pointer to a cs_navsto_param_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_rescale_pressure_to_ref(const cs_navsto_param_t   *nsp,
                                        const cs_cdo_quantities_t *quant,
                                        cs_real_t                  values[])
{
  cs_cdofb_navsto_rescale_pressure_to_ref(nsp, quant, values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a zero-mean
 *         average
 *
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_set_zero_mean_pressure(const cs_cdo_quantities_t *quant,
                                       cs_real_t                  values[])
{
  cs_cdofb_navsto_set_zero_mean_pressure(quant, values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform extra-operation related to MAC-Fb schemes when solving
 *         Navier-Stokes. Computation of the following quantities according to
 *         post-processing flags beeing activated.
 *         - The mass flux accross the boundaries.
 *         - The global mass in the computational domain
 *         - The norm of the velocity divergence
 *         - the cellwise mass flux balance
 *         - the kinetic energy
 *         - the velocity gradient
 *         - the pressure gradient
 *         - the vorticity
 *         - the helicity
 *         - the enstrophy
 *         - the stream function
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t struct.
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      quant         pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]      connect       pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      ts            pointer to a \ref cs_time_step_t struct.
 * \param[in,out]  time_plotter  pointer to a \ref cs_time_plot_t struct.
 * \param[in]      adv_field     pointer to a \ref cs_adv_field_t struct.
 * \param[in]      mass_flux     scalar-valued mass flux for each face
 * \param[in]      p_cell        scalar-valued pressure in each cell
 * \param[in]      u_face        vector-valued velocity on each face
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_extra_op(const cs_navsto_param_t   *nsp,
                         const cs_mesh_t           *mesh,
                         const cs_cdo_quantities_t *quant,
                         const cs_cdo_connect_t    *connect,
                         const cs_time_step_t      *ts,
                         cs_time_plot_t            *time_plotter,
                         const cs_adv_field_t      *adv_field,
                         const cs_real_t           *mass_flux,
                         const cs_real_t           *p_cell,
                         const cs_real_t           *u_face)
{
  CS_UNUSED(adv_field);

  /* 0. Compute full velocity at cells from normal velocity */

  cs_field_t *u_cell = cs_field_by_name("velocity");

  cs_macfb_navsto_compute_cell_velocity(connect, quant, mass_flux, u_cell->val);

  const cs_boundary_t *boundaries = nsp->boundaries;
  const cs_real_t     *bmass_flux = mass_flux + quant->n_i_faces;

  /* 1. Compute for each boundary the integrated mass flux to perform mass
   *    balance
   */

  bool *belong_to_default = nullptr;
  BFT_MALLOC(belong_to_default, quant->n_b_faces, bool);
  cs_array_bool_fill_true(quant->n_b_faces, belong_to_default);

  cs_real_t *boundary_fluxes = nullptr;
  BFT_MALLOC(boundary_fluxes, boundaries->n_boundaries + 1, cs_real_t);
  cs_array_real_fill_zero(boundaries->n_boundaries + 1, boundary_fluxes);

  for (int b_id = 0; b_id < boundaries->n_boundaries; b_id++) {

    const cs_zone_t *z = cs_boundary_zone_by_id(boundaries->zone_ids[b_id]);

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      const cs_lnum_t bf_id    = z->elt_ids[i];
      belong_to_default[bf_id] = false;
      boundary_fluxes[b_id] += bmass_flux[bf_id];
    }

  } /* Loop on domain boundaries */

  /* Update the flux through the default boundary */

  cs_lnum_t default_case_count = 0;
  for (cs_lnum_t i = 0; i < quant->n_b_faces; i++) {
    if (belong_to_default[i]) {
      default_case_count += 1;
      boundary_fluxes[boundaries->n_boundaries] += bmass_flux[i];
    }
  }

  /* Parallel synchronization if needed */

  cs_parall_sum(boundaries->n_boundaries + 1, CS_REAL_TYPE, boundary_fluxes);
  cs_parall_counter_max(&default_case_count, 1);

  /* Output result */

  cs_log_printf(CS_LOG_DEFAULT,
                "\n- Balance of the mass flux across the boundaries:\n");

  char descr[32];
  for (int b_id = 0; b_id < boundaries->n_boundaries; b_id++) {

    const cs_zone_t *z = cs_boundary_zone_by_id(boundaries->zone_ids[b_id]);

    cs_boundary_get_type_descr(boundaries, boundaries->types[b_id], 32, descr);

    cs_log_printf(CS_LOG_DEFAULT,
                  "b %-32s | %-32s |% -8.6e\n",
                  descr,
                  z->name,
                  boundary_fluxes[b_id]);

  } /* Loop on boundaries */

  /* Default boundary (if something to do) */

  if (default_case_count > 0) {

    cs_boundary_get_type_descr(boundaries, boundaries->default_type, 32, descr);
    cs_log_printf(CS_LOG_DEFAULT,
                  "b %-32s | %-32s |% -8.6e\n",
                  descr,
                  "default boundary",
                  boundary_fluxes[boundaries->n_boundaries]);
  }

  /* Free temporary buffers */

  BFT_FREE(belong_to_default);
  BFT_FREE(boundary_fluxes);

  /* Predefined post-processing */
  /* ========================== */

  /* There are five values if all flags are activated for the monitoring plot */

  int       n_cols      = 0;
  cs_real_t col_vals[5] = { 0, 0, 0, 0, 0 };

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE) {

    double      div_norm2 = 0.;
    cs_field_t *vel_div   = cs_field_by_name("velocity_divergence");
    assert(vel_div != nullptr);

    if (nsp->coupling != CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY)
      cs_field_current_to_previous(vel_div);

      /* Only the face velocity is used */

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)                      \
  reduction(+ : div_norm2)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double div_c
        = cs_macfb_navsto_cell_divergence(c_id, quant, connect->c2f, u_face);

      vel_div->val[c_id] = div_c;
      div_norm2 += quant->cell_vol[c_id] * div_c * div_c;

    } /* Loop on cells */

    cs_parall_sum(1, CS_DOUBLE, &div_norm2);
    col_vals[n_cols++] = sqrt(div_norm2);

  } /* Velocity divergence */

  if (nsp->post_flag & CS_NAVSTO_POST_MASS_DENSITY) {

    double      mass_integral = 0.;
    cs_field_t *rho           = cs_field_by_name("mass_density");
    assert(rho != nullptr);

    cs_field_current_to_previous(rho);

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)                      \
  reduction(+ : mass_integral)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double boussi_coef = 1;

      for (int i = 0; i < nsp->n_boussinesq_terms; i++) {

        cs_navsto_param_boussinesq_t *bp = nsp->boussinesq_param + i;
        boussi_coef += -bp->beta * (bp->var[c_id] - bp->var0);

      } /* Loop on Boussinesq terms */

      double rho_c   = nsp->mass_density->ref_value * boussi_coef;
      rho->val[c_id] = rho_c;
      mass_integral += quant->cell_vol[c_id] * rho_c;

    } /* Loop on cells */

    cs_parall_sum(1, CS_DOUBLE, &mass_integral);
    col_vals[n_cols++] = mass_integral;

  } /* Mass density */

  if (nsp->post_flag & CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE) {

    cs_field_t *mf_balance = cs_field_by_name("mass_flux_balance");
    assert(mf_balance != nullptr);

    cs_field_current_to_previous(mf_balance);

    const cs_adjacency_t *c2f = connect->c2f;

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double balance = 0.;
      for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id + 1]; j++) {
        balance += c2f->sgn[j] * mass_flux[c2f->ids[j]];
      }
      mf_balance->val[c_id] = balance;

    } /* Loop on cells */

  } /* Cell mass flux balance */

  if (nsp->post_flag & CS_NAVSTO_POST_PRESSURE_GRADIENT) {

    cs_field_t *pr_grd = cs_field_by_name("pressure_gradient");
    assert(pr_grd != nullptr);

    cs_field_current_to_previous(pr_grd);

    /* Compute a face pressure */

    cs_real_t *p_face = nullptr;
    BFT_MALLOC(p_face, quant->n_faces, cs_real_t);

    cs_macfb_navsto_compute_face_pressure(
      mesh, connect, quant, ts, nsp, p_cell, p_face);

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      cs_reco_grad_cell_from_fb_dofs(
        c_id, connect, quant, p_cell, p_face, pr_grd->val + 3 * c_id);

    BFT_FREE(p_face);

  } /* Pressure gradient */

  if (nsp->post_flag & CS_NAVSTO_POST_KINETIC_ENERGY) {

    double      k_integral     = 0.;
    cs_field_t *kinetic_energy = cs_field_by_name("kinetic_energy");
    assert(kinetic_energy != nullptr);

    cs_field_current_to_previous(kinetic_energy);

    if (cs_property_is_uniform(nsp->mass_density)) {

      /* This can be any cell but one assumes that there is at least one
      cell by
         MPI rank */

      const cs_real_t rho = cs_property_get_cell_value(0, /* cell_id */
                                                       ts->t_cur,
                                                       nsp->mass_density);

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)                      \
  reduction(+ : k_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        double kc = 0.5 * rho * cs_math_3_square_norm(u_cell->val + 3 * c_id);

        kinetic_energy->val[c_id] = kc;
        k_integral += quant->cell_vol[c_id] * kc;
      }
    }
    else { /* Mass density is not uniform in space */

#pragma omp parallel for if (quant->n_cells > CS_THR_MIN)                      \
  reduction(+ : k_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        cs_real_t rho_c
          = cs_property_get_cell_value(c_id, ts->t_cur, nsp->mass_density);
        double kc = 0.5 * rho_c * cs_math_3_square_norm(u_cell->val + 3 * c_id);

        kinetic_energy->val[c_id] = kc;
        k_integral += quant->cell_vol[c_id] * kc;
      }
    }

    cs_parall_sum(1, CS_DOUBLE, &k_integral); /* Sync. parallel computations */
    col_vals[n_cols++] = k_integral;

  } /* Kinetic energy */

  if (cs_glob_rank_id < 1 && time_plotter != nullptr)
    cs_time_plot_vals_write(
      time_plotter, ts->nt_cur, ts->t_cur, n_cols, col_vals);

  /* Stream function */
  /* --------------- */

  if (nsp->post_flag & CS_NAVSTO_POST_STREAM_FUNCTION) {

    cs_equation_t *eq = cs_equation_by_name(CS_NAVSTO_STREAM_EQNAME);
    assert(eq != nullptr);
    cs_equation_solve_steady_state(mesh, eq);

    cs_equation_param_t *eqp = cs_equation_get_param(eq);
    if (eqp->n_bc_defs == 0) {

      /* Since this is an equation solved with only homogeneous Neumann BCs, one
       * substracts the mean value to get a unique solution */

      cs_real_t mean_value;
      cs_equation_integrate_variable(connect, quant, eq, &mean_value);
      mean_value /= quant->vol_tot;

      cs_real_t *psi_v = cs_equation_get_vertex_values(eq, false);
      for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
        psi_v[i] -= mean_value;

    } /* If homogeneous Neumann everywhere */

  } /* Computation of the stream function is requested */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the normal velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         an algebraic technique.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_block_dirichlet_alge(short int                  f,
                              const cs_equation_param_t *eqp,
                              const cs_cell_mesh_t      *cm,
                              const cs_property_data_t  *pty,
                              cs_cell_builder_t         *cb,
                              cs_cell_sys_t             *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(cm);
  CS_UNUSED(pty);

  double   *x_dir  = cb->values;
  double   *ax_dir = cb->values + 1;
  cs_sdm_t *m      = csys->mat;
  int       n_dofs = csys->n_dofs;
  int       n_fc   = cm->n_fc;

  assert(n_fc == m->n_rows && n_dofs == m->n_cols);

  /* Build x_dir */

  bool is_non_homogeneous = false; /* Assume homogeneous by default */

  memset(cb->values, 0, 2 * sizeof(double));

  if (csys->dof_flag[f] & CS_CDO_BC_DIRICHLET) {
    x_dir[0]           = csys->dir_values[f];
    is_non_homogeneous = true;
  }

  if (is_non_homogeneous) {

    for (int bi = 0; bi < n_fc; bi++) {

      if (bi == f)
        continue;

      ax_dir[0] = m->val[bi * m->n_cols + f] * x_dir[0];
      csys->rhs[bi] -= ax_dir[0];
    }

  } /* Non-homogeneous Dirichlet BC */

  /* Set RHS to the Dirichlet value for the related face */

  if (f < n_fc) {
    csys->rhs[f] = x_dir[0];
  }

  /* Second pass: Replace the Dirichlet block by a diagonal block and fill
   * with zero the remaining row and column */

  for (int bi = 0; bi < n_fc; bi++) {

    if (bi == f) { /* bi == f */

      /* Reset term (I==F,J) = Row F */

      for (int bj = 0; bj < n_dofs; bj++) {
        m->val[bi * m->n_cols + bj] = 0.;
      }

      m->val[bi * m->n_cols + bi] = 1;
    }
    else { /* bi != f */

      /* Reset term (I,F) = Col F */

      m->val[bi * m->n_cols + f] = 0.;
    }

  } /* Block bi */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a penalization technique (with a large coefficient).
 *         One assumes that static condensation has been performed and that
 *         the velocity-block has size n_fc
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_block_dirichlet_pena(short int                  f,
                              const cs_equation_param_t *eqp,
                              const cs_cell_mesh_t      *cm,
                              const cs_property_data_t  *pty,
                              cs_cell_builder_t         *cb,
                              cs_cell_sys_t             *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  CS_UNUSED(f);
  CS_UNUSED(eqp);
  CS_UNUSED(cb);
  CS_UNUSED(cm);
  CS_UNUSED(pty);

  assert(csys != nullptr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size (n_fc + 1)
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_block_dirichlet_weak(short int                  fb,
                              const cs_equation_param_t *eqp,
                              const cs_cell_mesh_t      *cm,
                              const cs_property_data_t  *pty,
                              cs_cell_builder_t         *cb,
                              cs_cell_sys_t             *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  assert(cm != nullptr && cb != nullptr && csys != nullptr && pty != nullptr);
  assert(pty->is_iso == true);

  CS_UNUSED(fb);
  CS_UNUSED(eqp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (symmetrized Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size (n_fc + 1)
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_block_dirichlet_wsym(short int                  fb,
                              const cs_equation_param_t *eqp,
                              const cs_cell_mesh_t      *cm,
                              const cs_property_data_t  *pty,
                              cs_cell_builder_t         *cb,
                              cs_cell_sys_t             *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  assert(cm != nullptr && cb != nullptr && csys != nullptr && pty != nullptr);
  assert(cs_equation_param_has_diffusion(eqp));
  assert(pty->is_iso == true);

  CS_UNUSED(fb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a boundary defined as 'symmetry' (treated as a
 *         sliding BCs on the three velocity components.)
 *         A weak penalization technique (symmetrized Nitsche) is used.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_symmetry(short int                  fb,
                  const cs_equation_param_t *eqp,
                  const cs_cell_mesh_t      *cm,
                  const cs_property_data_t  *pty,
                  cs_cell_builder_t         *cb,
                  cs_cell_sys_t             *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  assert(cm != nullptr && cb != nullptr && csys != nullptr && pty != nullptr);
  assert(pty->is_iso == true); /* if not the case something else TODO ? */
  assert(cs_equation_param_has_diffusion(eqp));

  CS_UNUSED(fb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a wall BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_fixed_wall(short int                  fb,
                    const cs_equation_param_t *eqp,
                    const cs_cell_mesh_t      *cm,
                    const cs_property_data_t  *pty,
                    cs_cell_builder_t         *cb,
                    cs_cell_sys_t             *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  CS_UNUSED(cb);
  CS_UNUSED(pty);
  CS_UNUSED(fb);
  CS_UNUSED(eqp);

  assert(cm != nullptr && csys != nullptr); /* Sanity checks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one has to do one more non-linear iteration.
 *         Test if performed on the relative norm on the increment between
 *         two iterations
 *
 * \param[in]      nl_algo_type   type of non-linear algorithm
 * \param[in]      pre_iterate    previous state of the mass flux iterate
 * \param[in]      cur_iterate    current state of the mass flux iterate
 * \param[in, out] algo           pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_macfb_navsto_nl_algo_cvg(cs_param_nl_algo_t nl_algo_type,
                            const cs_real_t   *pre_iterate,
                            cs_real_t         *cur_iterate,
                            cs_iter_algo_t    *algo)
{
  return cs_cdofb_navsto_nl_algo_cvg(
    nl_algo_type, pre_iterate, cur_iterate, algo);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the function pointer computing the source term in the momentum
 *         equation related to the gravity effect (hydrostatic pressure or the
 *         Boussinesq approximation)
 *
 * \param[in]  nsp          set of parameters for the Navier-Stokes system
 * \param[out] p_func       way to compute the gravity effect
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_set_gravity_func(const cs_navsto_param_t   *nsp,
                                 cs_macfb_navsto_source_t **p_func)
{
  if (nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ) {

    switch (cs_macfb_navsto_boussinesq_type) {

    case CS_CDOFB_NAVSTO_BOUSSINESQ_FACE_DOF:
      *p_func = cs_macfb_navsto_boussinesq_at_face;
      break;

    case CS_CDOFB_NAVSTO_BOUSSINESQ_CELL_DOF:
      *p_func = cs_macfb_navsto_boussinesq_at_cell;
      break;

    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                "%s: Invalid type of algorithm to compute the Boussinesq"
                " approximation.\n",
                __func__);
    }
  }
  else if (nsp->model_flag & CS_NAVSTO_MODEL_GRAVITY_EFFECTS) {
    *p_func = cs_macfb_navsto_gravity_term;
  }

  else
    *p_func = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the gravity effects.
 *         Compute and add the source term to the local RHS.
 *         This is a special treatment since of face DoFs are involved.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_gravity_term(const cs_navsto_param_t         *nsp,
                             const cs_cell_mesh_t            *cm,
                             const cs_macfb_navsto_builder_t *nsb,
                             cs_cell_sys_t                   *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  assert(nsp->model_flag & CS_NAVSTO_MODEL_GRAVITY_EFFECTS);

  CS_UNUSED(nsb);
  CS_UNUSED(cm);
  CS_UNUSED(csys);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This is the standard case where the face DoFs are used for the
 *         constant part rho0 . g[] and only the cell DoFs are involved for the
 *         remaining part (the Boussinesq approximation).
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_boussinesq_at_cell(const cs_navsto_param_t         *nsp,
                                   const cs_cell_mesh_t            *cm,
                                   const cs_macfb_navsto_builder_t *nsb,
                                   cs_cell_sys_t                   *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  CS_UNUSED(nsb);
  CS_UNUSED(cm);
  CS_UNUSED(csys);
  assert(nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This way to compute the Boussinesq approximation applies only to
 *         face DoFs. This should enable to keep a stable (no velocity) in
 *         case of a stratified configuration.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_boussinesq_at_face(const cs_navsto_param_t         *nsp,
                                   const cs_cell_mesh_t            *cm,
                                   const cs_macfb_navsto_builder_t *nsb,
                                   cs_cell_sys_t                   *csys)
{
  bft_error(__FILE__, __LINE__, 0, _(" %s: Not implemented\n"), __func__);

  CS_UNUSED(nsb);
  CS_UNUSED(cm);
  CS_UNUSED(csys);
  assert(nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the source term for computing the stream function.
 *         This relies on the prototype associated to the generic function
 *         pointer \ref cs_dof_func_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         nullptr or pointer to a structure cast
 * on-the-fly \param[in, out] retval        result of the function. Must be
 * allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_stream_source_term(cs_lnum_t        n_elts,
                                   const cs_lnum_t *elt_ids,
                                   bool             dense_output,
                                   void            *input,
                                   cs_real_t       *retval)
{
  cs_cdofb_navsto_stream_source_term(
    n_elts, elt_ids, dense_output, input, retval);
}

/*----------------------------------------------------------------------------*/

#undef _dp3
END_C_DECLS
