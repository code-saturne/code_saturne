/*============================================================================
 * Helper functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

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
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_advection_field.h"
#include "cs_array.h"
#include "cs_evaluate.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_priv.c

  \brief Helper functions dedicated to groundwater flows when using CDO schemes
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */

#define _dp3 cs_math_3_dot_product

#define CS_GWF_PRIV_DBG 0

/*============================================================================
 * Local definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the portion of boundary fluxes crossing the area closing a
 *        dual cell. Optimized way.
 *
 * \param[in]      xv            coordinates of vertices
 * \param[in]      xf            face barycenter
 * \param[in]      f_u_normal    unit normal vector of the face
 * \param[in]      bf2v_idx      boundary face to vertices index
 * \param[in]      bf2v_ids      vertex ids of the boundary faces
 * \param[in]      cell_vector   cell-wise velocity field
 * \param[in, out] fluxes        computed fluxes at each vertex
 */
/*----------------------------------------------------------------------------*/

static void
_compute_v_weighted_boundary_fluxes(const cs_real_t *xv,
                                    const cs_real_t  xf[3],
                                    const cs_nreal_t face_unitv[3],
                                    const cs_lnum_t  bf2v_idx[2],
                                    const cs_lnum_t *bf2v_ids,
                                    const cs_real_t  cell_vector[3],
                                    cs_real_t        fluxes[])
{
  const cs_lnum_t *ids  = bf2v_ids + bf2v_idx[0];
  const int        n_vf = bf2v_idx[1] - bf2v_idx[0];

  for (cs_lnum_t v = 0; v < n_vf; v++)
    fluxes[v] = 0.; /* Init */

  /* Scaled by 1/2 since the part for a vertex v from the triangle tef is
     weighted by 1/2. This flux is computed without the face area which is
     added on-the-fly in the same time as the weight associated to each
     vertex */

  const cs_real_t unflx = 0.5 * cs_math_3_dot_product(face_unitv, cell_vector);

  int _v0, _v1;
  for (cs_lnum_t v = 0; v < n_vf; v++) {

    if (v < n_vf - 1)
      _v0 = v, _v1 = v + 1;
    else
      _v0 = n_vf - 1, _v1 = 0;

    const double tef
      = cs_math_surftri(xv + 3 * ids[_v0], xv + 3 * ids[_v1], xf);

    fluxes[_v0] += tef * unflx;
    fluxes[_v1] += tef * unflx;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the values of (potential) fields needed for the update of
 *        the Darcy velocity/fluxes.
 *
 * \param[in]  eq           pointer to an equation structure
 * \param[out] p_dof_vals   double pointer to the values (degrees of freedom)
 * \param[out] p_cell_vals  double pointer to the values (cell values)
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_get_value_pointers(const cs_equation_t *eq,
                          cs_real_t          **p_dof_vals,
                          cs_real_t          **p_cell_vals)
{
  cs_real_t *dof_vals = nullptr, *cell_vals = nullptr;

  switch (cs_equation_get_space_scheme(eq)) {
  case CS_SPACE_SCHEME_CDOVB:
    dof_vals = cs_equation_get_vertex_values(eq, false);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    dof_vals  = cs_equation_get_vertex_values(eq, false);
    cell_vals = cs_equation_get_cell_values(eq, false);
    break;

  case CS_SPACE_SCHEME_CDOFB:
    dof_vals  = cs_equation_get_face_values(eq, false);
    cell_vals = cs_equation_get_cell_values(eq, false);
    break;

  default:
    /* Do nothing */
    break;
  }

  /* Returns pointers */

  *p_dof_vals  = dof_vals;
  *p_cell_vals = cell_vals;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in]      loc_flag   flag to define where the flux is defined
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_darcy_flux_t *
cs_gwf_darcy_flux_create(cs_flag_t loc_flag)
{
  cs_gwf_darcy_flux_t *darcy = nullptr;

  BFT_MALLOC(darcy, 1, cs_gwf_darcy_flux_t);

  darcy->flux_location     = loc_flag;
  darcy->adv_field         = nullptr;
  darcy->flux_val          = nullptr;
  darcy->boundary_flux_val = nullptr;

  return darcy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] p_darcy   pointer of pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_free(cs_gwf_darcy_flux_t **p_darcy)
{
  if (p_darcy == nullptr)
    return;
  if (*p_darcy == nullptr)
    return;

  cs_gwf_darcy_flux_t *darcy = *p_darcy;

  /* flux_val and boundary_flux_val are allocated only if the related advection
     field is defined by array. At the definition step, the GWF module has kept
     the ownership of the lifecycle of the flux_val and boundary_flux_val
     arrays. In this case, the lifecycle is not managed by the definition and
     thus one has to free the arrays now. */

  BFT_FREE(darcy->boundary_flux_val);
  BFT_FREE(darcy->flux_val);

  BFT_FREE(darcy);
  *p_darcy = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] darcy   pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_log(cs_gwf_darcy_flux_t *darcy)
{
  if (darcy == nullptr)
    return;
  assert(darcy->adv_field != nullptr);

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Darcy name: %s with flux location: %s\n",
                darcy->adv_field->name,
                cs_flag_str_location(darcy->flux_location));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the advection field attached to a
 *        \ref cs_gwf_darcy_flux_t structure
 *        If the function pointer is set to nullptr, then an automatic settings
 *        is done.
 *
 * \param[in]      connect         pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq            pointer to a cs_cdo_quantities_t structure
 * \param[in]      space_scheme    space discretization using this structure
 * \param[in]      update_context  pointer to the context for the update step
 * \param[in]      update_func     pointer to an update function or nullptr
 * \param[in, out] darcy           pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_define(const cs_cdo_connect_t    *connect,
                         const cs_cdo_quantities_t *cdoq,
                         cs_param_space_scheme_t    space_scheme,
                         void                      *update_context,
                         cs_gwf_darcy_update_t     *update_func,
                         cs_gwf_darcy_flux_t       *darcy)
{
  if (darcy == nullptr)
    return;

  cs_adv_field_t *adv = darcy->adv_field;
  assert(adv != nullptr);

  /* Set the pointer to the input structure for the update step */

  darcy->update_input = update_context;
  darcy->update_func  = update_func;

  if (update_func == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Need a function pointer to apply the update.",
              __func__);

  /* Additional set-up */

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB: {
    const cs_adjacency_t *bf2v = connect->bf2v;

    /* Define the flux of the advection field at the boundary */

    cs_flag_t array_location = CS_FLAG_SCALAR | cs_flag_dual_closure_byf;
    size_t    array_size     = bf2v->idx[cdoq->n_b_faces];

    BFT_MALLOC(darcy->boundary_flux_val, array_size, cs_real_t);
    cs_array_real_fill_zero(array_size, darcy->boundary_flux_val);

    cs_xdef_t *bdy_def
      = cs_advection_field_def_boundary_flux_by_array(adv,
                                                      nullptr, /* all cells */
                                                      array_location,
                                                      darcy->boundary_flux_val,
                                                      false, /* not owner */
                                                      true); /* full length */

    cs_xdef_array_set_adjacency(bdy_def, bf2v);

    /* Define the advection field in the volume */

    if (cs_flag_test(darcy->flux_location, cs_flag_dual_face_byc)) {

      const cs_adjacency_t *c2e = connect->c2e;

      array_location = CS_FLAG_SCALAR | darcy->flux_location;
      array_size     = c2e->idx[cdoq->n_cells];
      BFT_MALLOC(darcy->flux_val, array_size, cs_real_t);
      cs_array_real_fill_zero(array_size, darcy->flux_val);

      /* Do not transfer the ownership (automatically on the full domain) */

      cs_xdef_t *adv_def = cs_advection_field_def_by_array(
        adv, array_location, darcy->flux_val, false);

      cs_xdef_array_set_adjacency(adv_def, c2e);

      /* Reset the type of advection field */

      if (adv->status & CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR)
        adv->status -= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;
      adv->status |= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;
    }
    else if (cs_flag_test(darcy->flux_location, cs_flag_primal_cell)) {

      cs_field_t *cell_adv_field
        = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
      assert(cell_adv_field != nullptr);

      cs_advection_field_def_by_field(adv, cell_adv_field);

      /* Reset the type of advection field */

      if (adv->status & CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX)
        adv->status -= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;
      adv->status |= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;
    }
    else
      bft_error(__FILE__,
                __LINE__,
                0,
                " %s: Invalid location for the definition of the Darcy flux.",
                __func__);

  } break;

  case CS_SPACE_SCHEME_CDOFB: /* TODO */
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: CDO-Fb space scheme not fully implemented.",
              __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    break;

  } /* Switch on the space scheme */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate the balance by zone (relying on the splitting arising from
 *         the boundary settings) for the advection field attached to a \ref
 *         cs_gwf_darcy_flux_t structure
 *
 * \param[in]       connect       pointer to a cs_cdo_connect_t structure
 * \param[in]       quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]       eqp           pointer to the set of equation parameters
 * \param[in, out]  darcy         pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_balance(const cs_cdo_connect_t    *connect,
                          const cs_cdo_quantities_t *quant,
                          const cs_equation_param_t *eqp,
                          cs_gwf_darcy_flux_t       *darcy)
{
  if (darcy == nullptr)
    return;

  const cs_lnum_t       n_b_faces = quant->n_b_faces;
  const cs_adv_field_t *adv       = darcy->adv_field;
  assert(adv != nullptr);

  /* Try to retrieve the boundary flux values */

  const cs_field_t *nflx
    = cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);
  cs_real_t *flux_val
    = (nflx == nullptr) ? darcy->boundary_flux_val : nflx->val;

  if (flux_val == nullptr && n_b_faces > 0) /* No value on which operates */
    return;

  const cs_adjacency_t *bf2v = connect->bf2v;

  /* Define the balance by zone */

  bool *is_counted = nullptr;
  BFT_MALLOC(is_counted, n_b_faces, bool);
#pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (int i = 0; i < n_b_faces; i++)
    is_counted[i] = false;

  /* n_bc_defs + 1 to take into account the default boundary condition */

  cs_real_t *balances = nullptr;
  BFT_MALLOC(balances, eqp->n_bc_defs + 1, cs_real_t);

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {

    const cs_xdef_t *def = eqp->bc_defs[ibc];
    const cs_zone_t *z   = cs_boundary_zone_by_id(def->z_id);

    balances[ibc] = 0;

    if (nflx == nullptr) { /* The definition of the boundary flux relies on the
                           bf2v adjacency */

#if defined(DEBUG) && !defined(NDEBUG)
      cs_xdef_t               *_def = adv->bdy_flux_defs[0];
      cs_xdef_array_context_t *cx   = (cs_xdef_array_context_t *)_def->context;

      assert(adv->n_bdy_flux_defs == 1 && _def->type == CS_XDEF_BY_ARRAY);
      assert(cs_flag_test(cx->value_location, cs_flag_dual_closure_byf)
             == true);
#endif

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t bf_id = z->elt_ids[i];
        is_counted[bf_id]     = true;
        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id + 1]; j++)
          balances[ibc] += flux_val[j];
      }
    }
    else {

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t bf_id = z->elt_ids[i];
        is_counted[bf_id]     = true;
        balances[ibc] += flux_val[bf_id];
      }

    } /* nflux is nullptr ? */

  } /* Loop on BC definitions */

  /* Manage the default boundary condition if there is at least one face not
     considered in the previous loop */

  bool display             = false;
  balances[eqp->n_bc_defs] = 0.;
  for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
    if (is_counted[bf_id] == false) {

      display = true;
      if (nflx == nullptr) {

        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id + 1]; j++)
          balances[eqp->n_bc_defs] += flux_val[j];
      }
      else
        balances[eqp->n_bc_defs] += flux_val[bf_id];

    } /* Not already counted */
  } /* Loop on boundary faces */

  int display_flag = display ? 1 : 0;

  /* Parallel synchronizations */

  cs_parall_max(1, CS_INT_TYPE, &display_flag);
  cs_parall_sum(eqp->n_bc_defs + 1, CS_REAL_TYPE, balances);

  /* Output into the default log file */

  cs_log_printf(CS_LOG_DEFAULT,
                "-b- Balance of %s across the boundary zones:\n",
                adv->name);

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {
    const cs_zone_t *z = cs_boundary_zone_by_id((eqp->bc_defs[ibc])->z_id);
    cs_log_printf(
      CS_LOG_DEFAULT, "-b- %-32s: % -5.3e\n", z->name, balances[ibc]);
  }

  if (display_flag > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "-b- %-32s: % -5.3e\n",
                  "Remaining part of the boundary",
                  balances[eqp->n_bc_defs]);

  /* Free memory */

  BFT_FREE(is_counted);
  BFT_FREE(balances);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the associated Darcy flux over the boundary of the domain for
 *        each vertex of a boundary face.  Case of a vertex-based
 *        discretization and single-phase flows in porous media (saturated or
 *        not).
 *
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      eq       pointer to the equation related to this Darcy flux
 * \param[in, out] adv      pointer to the Darcy advection field
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update_on_boundary(cs_real_t            t_eval,
                                     const cs_equation_t *eq,
                                     cs_adv_field_t      *adv)
{
  if (adv->n_bdy_flux_defs > 1
      || adv->bdy_flux_defs[0]->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_xdef_t               *def      = adv->bdy_flux_defs[0];
  cs_xdef_array_context_t *cx       = (cs_xdef_array_context_t *)def->context;
  cs_real_t               *nflx_val = cx->values;

  if (cs_flag_test(cx->value_location, cs_flag_dual_closure_byf) == false)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_equation_compute_boundary_diff_flux(t_eval, eq, nflx_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the associated Darcy flux over the boundary of the domain for
 *        each vertex of a boundary face without using an equation (i.e. there
 *        is no associated boundary condition).
 *        Case of a vertex-based discretization.
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      cell_vel   Darcy velocity in each cell
 * \param[in, out] adv        pointer to the Darcy advection field to update
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update_on_boundary_wo_eq(const cs_cdo_connect_t    *connect,
                                           const cs_cdo_quantities_t *cdoq,
                                           cs_real_t                 *cell_vel,
                                           cs_adv_field_t            *adv)
{
  if (adv->n_bdy_flux_defs > 1
      || adv->bdy_flux_defs[0]->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_xdef_t               *def      = adv->bdy_flux_defs[0];
  cs_xdef_array_context_t *cx       = (cs_xdef_array_context_t *)def->context;
  cs_real_t               *nflx_val = cx->values;

  if (cs_flag_test(cx->value_location, cs_flag_dual_closure_byf) == false)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  const cs_adjacency_t *bf2v     = connect->bf2v;
  const cs_adjacency_t *f2c      = connect->f2c;
  const cs_lnum_t      *bf2c_ids = f2c->ids + f2c->idx[cdoq->n_i_faces];

#pragma omp parallel for if (cdoq->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t bf_id = 0; bf_id < cdoq->n_b_faces; bf_id++) {

    cs_lnum_t *bf2v_idx = bf2v->idx + bf_id;
    _compute_v_weighted_boundary_fluxes(cdoq->vtx_coord,
                                        cdoq->b_face_center + 3 * bf_id,
                                        cdoq->b_face_u_normal[bf_id],
                                        bf2v_idx,
                                        bf2v->ids,
                                        cell_vel + 3 * bf2c_ids[bf_id],
                                        nflx_val + bf2v_idx[0]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update head values (pressure head or head values for laws)
 *        Up to now, this is only used for single-phase flows in porous media
 *        (saturated or not case).
 *
 * \param[in]      connect         pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq            pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards        pointer to the Richards equation
 * \param[in]      option_flag     calculation option related to the GWF module
 * \param[in, out] pressure_head   pressure head field
 * \param[in, out] head_in_law     values of the head used in law
 * \param[in]      cur2prev        true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update_head(const cs_cdo_connect_t    *connect,
                   const cs_cdo_quantities_t *cdoq,
                   const cs_equation_t       *richards,
                   cs_flag_t                  option_flag,
                   cs_field_t                *pressure_head,
                   cs_real_t                  head_in_law[],
                   bool                       cur2prev)
{
  cs_param_space_scheme_t r_scheme = cs_equation_get_space_scheme(richards);
  cs_field_t             *hydraulic_head = cs_equation_get_field(richards);

  if (option_flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE) {

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB: {
      assert(hydraulic_head->location_id
             == cs_mesh_location_get_id_by_name("vertices"));

      cs_real_t domain_integral = cs_evaluate_scal_domain_integral_by_array(
        cs_flag_primal_vtx, hydraulic_head->val);

      const cs_real_t mean_value = domain_integral / cdoq->vol_tot;

#pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++)
        hydraulic_head->val[i] -= mean_value;
    } break;

    case CS_SPACE_SCHEME_CDOFB: {
      assert(hydraulic_head->location_id
             == cs_mesh_location_get_id_by_name("cells"));

      cs_real_t domain_integral = cs_evaluate_scal_domain_integral_by_array(
        cs_flag_primal_cell, hydraulic_head->val);

      const cs_real_t mean_value = domain_integral / cdoq->vol_tot;
#pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++)
        hydraulic_head->val[i] -= mean_value;
    } break;

    default:
      break; /* Nothing to do */
    }

  } /* Rescale hydraulic_head */

  if (option_flag & CS_GWF_GRAVITATION) { /* Update the pressure head (and if
                                          needed head_in_law) */

    cs_physical_constants_t *phys = cs_get_glob_physical_constants();

    if (pressure_head == nullptr)
      bft_error(__FILE__,
                __LINE__,
                0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */

    if (cur2prev)
      cs_field_current_to_previous(pressure_head);

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:

#pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        const cs_real_t gpot  = _dp3(cdoq->vtx_coord + 3 * i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }

      /* Update head_in_law */

      if (head_in_law != nullptr)
        cs_reco_scalar_v2c_full(
          connect->c2v, cdoq, pressure_head->val, head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB: {
      assert(hydraulic_head->location_id
             == cs_mesh_location_get_id_by_name("vertices"));

#pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        const cs_real_t gpot  = _dp3(cdoq->vtx_coord + 3 * i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }

      /* Update head_in_law */

      if (head_in_law != nullptr) {

        const cs_real_t *hydraulic_head_cells
          = cs_equation_get_cell_values(richards, false); /* current values */

#pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
          const cs_real_t gpot
            = _dp3(cdoq->cell_centers + 3 * i, phys->gravity);
          head_in_law[i] = hydraulic_head_cells[i] - gpot;
        }

      } /* head_in_law */

    } break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:

#pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
        const cs_real_t gpot  = _dp3(cdoq->cell_centers + 3 * i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }
      break;

      /* Head in law points either to hydraulic_head->val or pressure_head->val
       */

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

    } /* Switch on space scheme */
  }
  else { /* No gravity effect is taken into account */

    if (head_in_law == nullptr)
      return;

    /* Update head_in_law */

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cs_reco_scalar_v2c_full(
        connect->c2v, cdoq, hydraulic_head->val, head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB: {
      const cs_real_t *hydraulic_head_cells
        = cs_equation_get_cell_values(richards, false); /* current values */

      cs_array_real_copy(cdoq->n_cells, hydraulic_head_cells, head_in_law);
    } break;

    default:
      break; /* Nothing to do for CDO-Fb schemes and HHO schemes */

    } /* Switch on the space scheme related to the Richards equation */

  } /* Gravity is activated or not */
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
