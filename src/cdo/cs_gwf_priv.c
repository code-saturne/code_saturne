/*============================================================================
 * Helper functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_advection_field.h"
#include "cs_parall.h"

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
 * \brief  Compute the associated boundary Darcy flux for each vertex of
 *         boundary faces.
 *         Case of a vertex-based discretization and single-phase flows in
 *         porous media (saturated or not).
 *
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      eq       pointer to the equation related to this Darcy flux
 * \param[in, out] adv      pointer to the Darcy advection field
 */
/*----------------------------------------------------------------------------*/

static void
_update_vb_darcy_flux_at_boundary(cs_real_t                t_eval,
                                  const cs_equation_t     *eq,
                                  cs_adv_field_t          *adv)
{
  if (adv->n_bdy_flux_defs > 1 ||
      adv->bdy_flux_defs[0]->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_xdef_t  *def = adv->bdy_flux_defs[0];
  cs_xdef_array_context_t  *actx = (cs_xdef_array_context_t *)def->context;
  cs_real_t  *nflx_val = actx->values;

  if (cs_flag_test(actx->loc, cs_flag_dual_closure_byf) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_equation_compute_boundary_diff_flux(t_eval, eq, nflx_val);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
cs_gwf_darcy_flux_create(cs_flag_t     loc_flag)
{
  cs_gwf_darcy_flux_t  *darcy = NULL;

  BFT_MALLOC(darcy, 1, cs_gwf_darcy_flux_t);

  darcy->flux_location = loc_flag;
  darcy->adv_field = NULL;
  darcy->flux_val = NULL;
  darcy->boundary_flux_val = NULL;

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
cs_gwf_darcy_flux_free(cs_gwf_darcy_flux_t  **p_darcy)
{
  if (p_darcy == NULL)
    return;
  if (*p_darcy == NULL)
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
  *p_darcy = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] darcy   pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_log(cs_gwf_darcy_flux_t  *darcy)
{
  if (darcy == NULL)
    return;
  assert(darcy->adv_field != NULL);

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Darcy name: %s with flux location: %s\n",
                darcy->adv_field->name,
                cs_flag_str_location(darcy->flux_location));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the definition of the advection field attached to a
 *         \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in]       connect       pointer to a cs_cdo_connect_t structure
 * \param[in]       quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]       space_scheme  space discretization using this structure
 * \param[in, out]  darcy         pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_define(const cs_cdo_connect_t       *connect,
                         const cs_cdo_quantities_t    *quant,
                         cs_param_space_scheme_t       space_scheme,
                         cs_gwf_darcy_flux_t          *darcy)
{
  if (darcy == NULL)
    return;

  cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != NULL);

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_adjacency_t  *bf2v = connect->bf2v;

      /* Define the flux of the advection field at the boundary */

      cs_flag_t  array_location = CS_FLAG_SCALAR | cs_flag_dual_closure_byf;
      size_t  array_size = bf2v->idx[quant->n_b_faces];

      BFT_MALLOC(darcy->boundary_flux_val, array_size, cs_real_t);
      memset(darcy->boundary_flux_val, 0, array_size*sizeof(cs_real_t));

      /* Do not transfer the ownership */

      cs_advection_field_def_boundary_flux_by_array(adv,
                                                    NULL,
                                                    array_location,
                                                    darcy->boundary_flux_val,
                                                    false,
                                                    bf2v->idx,
                                                    bf2v->ids);

      /* Define the advection field in the volume */

      if (cs_flag_test(darcy->flux_location, cs_flag_dual_face_byc)) {

        const cs_adjacency_t  *c2e = connect->c2e;

        array_location = CS_FLAG_SCALAR | darcy->flux_location;
        array_size = c2e->idx[quant->n_cells];
        BFT_MALLOC(darcy->flux_val, array_size, cs_real_t);
        memset(darcy->flux_val, 0, array_size*sizeof(cs_real_t));

        /* Do not transfer the ownership */

        cs_advection_field_def_by_array(adv,
                                        array_location,
                                        darcy->flux_val,
                                        false, /* transfer ownership */
                                        c2e->idx,
                                        c2e->ids);

        /* Reset the type of advection field */

        if (adv->status & CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR)
          adv->status -= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;
        adv->status |= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

      }
      else if (cs_flag_test(darcy->flux_location, cs_flag_primal_cell)) {

        cs_field_t  *cell_adv_field =
          cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
        assert(cell_adv_field != NULL);

        cs_advection_field_def_by_field(adv, cell_adv_field);

        /* Reset the type of advection field */

        if (adv->status & CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX)
          adv->status -= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;
        adv->status |= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid location for the definition of the Darcy flux.",
                  __func__);

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:   /* TODO */
    bft_error(__FILE__, __LINE__, 0,
              " %s: CDO-Fb space scheme not fully implemented.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    break;

  } /* Switch on the space scheme */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the advection field/arrays related to the Darcy flux.
 *
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      eq        pointer to the equation related to this Darcy flux
 * \param[in]      cur2prev  true or false
 * \param[in, out] darcy     pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update(const cs_real_t              t_eval,
                         const cs_equation_t         *eq,
                         bool                         cur2prev,
                         cs_gwf_darcy_flux_t         *darcy)
{
  cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != NULL);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);

  assert(vel != NULL); /* Sanity check */
  if (cur2prev)
    cs_field_current_to_previous(vel);

  /* Update arrays related to the Darcy flux:
   * Compute the new darcian flux and darcian velocity inside each cell */

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:

    /* Update the array of flux values associated to the advection field */

    if (cs_flag_test(darcy->flux_location, cs_flag_dual_face_byc)) {

      assert(darcy->flux_val != NULL);
      if (adv->definition->type != CS_XDEF_BY_ARRAY)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid definition of the advection field", __func__);

      cs_equation_compute_diffusive_flux(eq,
                                         darcy->flux_location,
                                         t_eval,
                                         darcy->flux_val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_PRIV_DBG > 2
      cs_dbg_darray_to_listing("DARCIAN_FLUX_DFbyC",
                               connect->c2e->idx[cdoq->n_cells],
                               darcy->flux_val, 8);
#endif

      /* Set the new values of the vector field at cell centers */

      cs_advection_field_in_cells(adv, t_eval, vel->val);

    }
    else if (cs_flag_test(darcy->flux_location, cs_flag_primal_cell))
      cs_equation_compute_diffusive_flux(eq,
                                         darcy->flux_location,
                                         t_eval,
                                         vel->val);

    /* Update the Darcy flux at the boundary */

    _update_vb_darcy_flux_at_boundary(t_eval, eq, adv);
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    bft_error(__FILE__, __LINE__, 0, " TODO.");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

  } /* End of switch */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_PRIV_DBG > 1
  cs_dbg_darray_to_listing("DARCIAN_FLUX_CELL", 3*cdoq->n_cells, vel->val, 3);
#endif

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);

  }
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
cs_gwf_darcy_flux_balance(const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_equation_param_t    *eqp,
                          cs_gwf_darcy_flux_t          *darcy)
{
  if (darcy == NULL)
    return;

  const cs_lnum_t  n_b_faces = quant->n_b_faces;
  const cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != NULL);

  /* Try to retrieve the boundary flux values */

  const cs_field_t  *nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);
  cs_real_t  *flux_val = (nflx == NULL) ? darcy->boundary_flux_val : nflx->val;

  if (flux_val == NULL && n_b_faces > 0) /* No value on which operates */
    return;

  const cs_adjacency_t  *bf2v = connect->bf2v;

  /* Define the balance by zone */

  bool  *is_counted = NULL;
  BFT_MALLOC(is_counted, n_b_faces, bool);
# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (int i = 0; i < n_b_faces; i++) is_counted[i] = false;

  /* n_bc_defs + 1 to take into account the default boundary condition */

  cs_real_t  *balances = NULL;
  BFT_MALLOC(balances, eqp->n_bc_defs + 1, cs_real_t);

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {

    const cs_xdef_t  *def = eqp->bc_defs[ibc];
    const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);

    balances[ibc] = 0;

    if (nflx == NULL) { /* The definition of the boundary flux relies on the
                           bf2v adjacency */

#if defined(DEBUG) && !defined(NDEBUG)
      cs_xdef_t  *_def = adv->bdy_flux_defs[0];
      cs_xdef_array_context_t  *actx = _def->context;

      assert(adv->n_bdy_flux_defs == 1 && _def->type == CS_XDEF_BY_ARRAY);
      assert(cs_flag_test(actx->loc, cs_flag_dual_closure_byf) == true);
#endif

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  bf_id = z->elt_ids[i];
        is_counted[bf_id] = true;
        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++)
          balances[ibc] += flux_val[j];
      }

    }
    else {

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  bf_id = z->elt_ids[i];
        is_counted[bf_id] = true;
        balances[ibc] += flux_val[bf_id];
      }

    } /* nflux is NULL ? */

  } /* Loop on BC definitions */

  /* Manage the default boundary condition if there is at least one face not
     considered in the previous loop */

  bool  display = false;
  balances[eqp->n_bc_defs] = 0.;
  for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
    if (is_counted[bf_id] == false) {

      display = true;
      if (nflx == NULL) {

        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++)
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
                "-b- Balance of %s across the boundary zones:\n", adv->name);

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {
    const cs_zone_t  *z = cs_boundary_zone_by_id((eqp->bc_defs[ibc])->z_id);
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-32s: % -5.3e\n",
                  z->name, balances[ibc]);
  }

  if (display_flag > 0)
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-32s: % -5.3e\n",
                  "Remaining part of the boundary", balances[eqp->n_bc_defs]);

  /* Free memory */

  BFT_FREE(is_counted);
  BFT_FREE(balances);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
