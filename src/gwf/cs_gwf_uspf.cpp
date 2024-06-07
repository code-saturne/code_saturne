/*============================================================================
 * Main functions dedicated to groundwater flows in case of single-phase flows
 * in an unsaturated porous media (shorter notation = USPF)
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
#include <bft_printf.h>

#include "cs_array.h"
#include "cs_equation_bc.h"
#include "cs_gwf_soil.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_post.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_uspf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_uspf.c

  \brief Main high-level functions dedicated to groundwater flows when using
         CDO schemes for single-phase flows in an unsaturated porous media.
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_USPF_DBG 0

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
 * \brief Update the advection field/arrays related to the Darcy flux
 *        Case of CDO-Vb schemes and a Darcy flux defined at dual faces for an
 *        unsaturated porous media
 *
 *        cell_values could be set to nullptr when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or nullptr
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_update_darcy_arrays(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq,
                          const cs_real_t             *dof_values,
                          const cs_real_t             *cell_values,
                          cs_real_t                    t_eval,
                          bool                         cur2prev,
                          cs_gwf_darcy_flux_t         *darcy)
{
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  cs_adv_field_t  *adv = darcy->adv_field;

  assert(darcy->flux_val != nullptr);
  assert(adv != nullptr);
  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_uspf_t  *mc = (cs_gwf_uspf_t *)darcy->update_input;
  cs_equation_t  *eq = mc->richards;

  /* Update the array of flux values associated to the advection field */

  assert(eq != nullptr);
  cs_equation_compute_diffusive_flux(eq,
                                     nullptr, /* eqp --> default */
                                     nullptr, /* diff_pty --> default */
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux.
   * This field is always defined when the definition relies on a flux. */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != nullptr);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at the boundary (take into account the BCs) */

  cs_gwf_darcy_flux_update_on_boundary(t_eval, eq, adv);

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx
      != nullptr) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the model context structure in case of
 *         single-phase flows in an unsaturated porous media
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_uspf_t *
cs_gwf_uspf_create(void)
{
  cs_gwf_uspf_t *mc = nullptr;

  BFT_MALLOC(mc, 1, cs_gwf_uspf_t);

  mc->permeability_field = nullptr;
  mc->moisture_field     = nullptr;
  mc->capacity_field     = nullptr;
  mc->pressure_head      = nullptr;
  mc->head_in_law        = nullptr;

  /* Create a new equation structure for Richards' equation */

  mc->richards = cs_equation_add("Richards",       /* equation name */
                                 "hydraulic_head", /* variable name */
                                 CS_EQUATION_TYPE_GROUNDWATER,
                                 1,
                                 CS_BC_SYMMETRY);

  /* Define the Darcy flux structure.

     Add an advection field related to the darcian flux stemming from the
     Richards equation. This advection field is steady since the head is
     steady-state in this model. The flux is defined at the dual faces (split
     by primal cells)
  */

  mc->darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF |
    CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  mc->darcy->adv_field = cs_advection_field_add("darcy_field", adv_status);

  /* Add several properties:
   * - the moisture content
   * - the soil capacity --> for the unsteady term
   * - the (full) permeability = rel_permeability * abs_permeability
   */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  mc->soil_capacity = cs_property_add("soil_capacity", CS_PROPERTY_ISO);

  mc->permeability = nullptr; /* Done when the type of property is set */

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the model context structure in case of unsaturated single-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_free(cs_gwf_uspf_t   **p_mc)
{
  if (p_mc == nullptr)
    return;
  if (*p_mc == nullptr)
    return;

  cs_gwf_uspf_t  *mc = *p_mc;

  cs_gwf_darcy_flux_free(&(mc->darcy));

  BFT_FREE(mc);
  *p_mc = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of unsaturated
 *        single-phase flows in porous media
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_log_setup(cs_gwf_uspf_t   *mc)
{
  if (mc == nullptr)
    return;

  cs_gwf_darcy_flux_log(mc->darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of an unsaturated single-phase flows model in porous media
 *
 * \param[in, out] mc          pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_init(cs_gwf_uspf_t          *mc,
                 cs_property_type_t      perm_type)
{
  if (mc == nullptr)
    return;

  cs_equation_t  *eq = mc->richards;

  if (eq == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  assert(eqp != nullptr);

  /* Add the property related to the diffusion term */

  mc->permeability = cs_property_add("permeability", perm_type);

  /* Add the diffusion term to the Richards equation by associating the
     full permeability to the diffusion property of the Richards eq. */

  cs_equation_add_diffusion(eqp, mc->permeability);

  /* Associate the soil capacity to the unsteady term of the Richards
     equation */

  assert(mc->soil_capacity != nullptr);
  cs_equation_add_time(eqp, mc->soil_capacity);

  /* Default treatment of the Dirichlet BCs (algebraic since this is a pure
     diffusion equation) */

  cs_equation_param_set(eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");

  /* The steady diffusion eq. is solved once at the beginning and then the
     hydraulic head is used to compute the advective flux for the tracer
     equations. A good accuracy is thus required knowing that strong
     heterogeneities may be taken into account. Moreover, the system is SPD by
     construction. One can thus use a AMG + FCG as the default solver. */

  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RTOL, "1e-8");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "fcg");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_MAX_ITER, "5000");

  /* The default value of the relative tolerance is 1e-6 (please refer to
     cs_param_sles_create() to know all the default settings related to linear
     algebra).

     The default value of the max. number of iteration is 10000. The AMG/FCG
     solver is more efficient but each iteration is more costly. Thus, one
     reduces this max. number of iterations if one encounters a convergence
     issue.
  */

  cs_equation_predefined_create_field(1, eq); /* Keep two states */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for single-phase flows in an unsaturated porous
 *        media. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in]      flag         optional settings for the module
 * \param[in]      post_flag    optional postprocessing request(s)
 * \param[in]      perm_dim     dimension of the permeability (scalar, etc.)
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_init_setup(cs_flag_t            flag,
                       cs_flag_t            post_flag,
                       int                  perm_dim,
                       cs_gwf_uspf_t       *mc)
{
  if (mc == nullptr)
    return;

  if (mc->richards == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);
  assert(cs_equation_is_steady(mc->richards) == false);

  cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);
  assert(eqp != nullptr);

  /* Set additional fields */

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Handle gravity effects */

  if (flag & CS_GWF_GRAVITATION) {

    switch (eqp->space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          v_loc_id,
                                          1,
                                          true); /* has_previous */
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          true); /* has_previous */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(mc->pressure_head, log_key, 1);
    cs_field_set_key_int(mc->pressure_head, post_key, 1);

  } /* Gravitation effects are activated */

  /* Field for the liquid saturation */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;
  mc->moisture_field = cs_field_create("liquid_saturation",
                                       pty_mask,
                                       c_loc_id,
                                       1,     /* dimension */
                                       true); /* has_previous */

  cs_field_set_key_int(mc->moisture_field, log_key, 1);
  if (post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(mc->moisture_field, post_key, 1);

  /* Field for the permeability */

  mc->permeability_field = cs_field_create("permeability",
                                           pty_mask,
                                           c_loc_id,
                                           perm_dim,
                                           true); /* has_previous */

  if (post_flag & CS_GWF_POST_PERMEABILITY) {
    cs_field_set_key_int(mc->permeability_field, log_key, 1);
    cs_field_set_key_int(mc->permeability_field, post_key, 1);
  }

  /* Field for the soil capacity */

  mc->capacity_field = cs_field_create("soil_capacity",
                                       pty_mask,
                                       c_loc_id,
                                       1,   /* dimension */
                                       true);

  cs_field_set_key_int(mc->capacity_field, log_key, 1);
  if (post_flag & CS_GWF_POST_SOIL_CAPACITY)
    cs_field_set_key_int(mc->capacity_field, post_key, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of single-phase flows in an unsaturated
 *        porous media
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      flag      optional settings for the module
 * \param[in, out] mc        pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_finalize_setup(const cs_cdo_connect_t         *connect,
                           const cs_cdo_quantities_t      *cdoq,
                           cs_flag_t                       flag,
                           cs_gwf_uspf_t                  *mc)
{
  const cs_field_t  *hydraulic_head = cs_equation_get_field(mc->richards);
  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(mc->richards);
  const cs_lnum_t  n_cells = connect->n_cells;

  /* Set the Darcian flux (in the volume and at the boundary) */

  cs_gwf_darcy_flux_define(connect, cdoq,
                           richards_scheme,
                           mc,                          /* context */
                           _uspf_update_darcy_arrays,  /* update function */
                           mc->darcy);

  /* Allocate a head array defined at cells and used to update the soil
     properties */

  switch (richards_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    BFT_MALLOC(mc->head_in_law, n_cells, cs_real_t);
#if defined(DEBUG) && !defined(NDEBUG)
    cs_array_real_fill_zero(n_cells, mc->head_in_law);
#endif
    break;

  case CS_SPACE_SCHEME_CDOFB:

    if (flag & CS_GWF_GRAVITATION)
      mc->head_in_law = mc->pressure_head->val;
    else
      mc->head_in_law = hydraulic_head->val;

    bft_error(__FILE__, __LINE__, 0,
              "%s: Fb space scheme not fully implemented.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme.", __func__);
    break;

  } /* Switch on Richards scheme */

  /* Define the permeability property using a field */

  cs_property_def_by_field(mc->permeability, mc->permeability_field);

  /* Define the moisture content (liquid saturation) property using a field */

  cs_property_def_by_field(mc->moisture_content, mc->moisture_field);

  /* Define the soil capacity property using a field */

  cs_property_def_by_field(mc->soil_capacity, mc->capacity_field);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update quantities in the case of single-phase flows in an unsaturated
 *        porous media
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_update(const cs_mesh_t                *mesh,
                   const cs_cdo_connect_t         *connect,
                   const cs_cdo_quantities_t      *cdoq,
                   const cs_time_step_t           *ts,
                   cs_flag_t                       update_flag,
                   cs_flag_t                       option_flag,
                   cs_gwf_uspf_t                  *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = cs_equation_get_time_eval(ts, mc->richards);
    cur2prev = true;

  }

  /* Update head */

  cs_gwf_update_head(connect, cdoq, mc->richards,
                     option_flag,
                     mc->pressure_head,
                     mc->head_in_law,
                     cur2prev);

  /* Update the advection field related to the groundwater flow module */

  cs_gwf_darcy_flux_t  *darcy = mc->darcy;
  cs_real_t            *dof_vals = nullptr, *cell_vals = nullptr;

  cs_gwf_get_value_pointers(mc->richards, &dof_vals, &cell_vals);

  darcy->update_func(connect,
                     cdoq,
                     dof_vals,
                     cell_vals,
                     time_eval,
                     cur2prev,
                     darcy);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_USPF_DBG > 2
  if (cs_flag_test(darcy->flux_location, cs_flag_dual_face_byc))
    cs_dbg_darray_to_listing("DARCIAN_FLUX_DFbyC",
                             connect->c2e->idx[cdoq->n_cells],
                             darcy->flux_val, 8);
  else if (cs_flag_test(darcy->flux_location, cs_flag_primal_cell)) {
    cs_field_t  *vel = cs_advection_field_get_field(darcy->adv_field,
                                                    CS_MESH_LOCATION_CELLS);
    cs_dbg_darray_to_listing("DARCIAN_FLUX_CELL",
                             3*cdoq->n_cells, vel->val, 3);
  }
#endif

  /* Update properties related to soils.
   * Handle the permeability, the moisture content and the soil capacity */

  if (cur2prev) {

    cs_field_current_to_previous(mc->permeability_field);
    cs_field_current_to_previous(mc->moisture_field);
    cs_field_current_to_previous(mc->capacity_field);

  }

  /* Update soil properties with the new head values */

  cs_gwf_soil_update(time_eval, mesh, connect, cdoq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new state for the groundwater flows module.
 *        Case of single-phase flows in an unstaturated porous media.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a cs_time_step_t structure
 * \param[in]      flag        optional metadata for the module
 * \param[in, out] mc          pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_compute(const cs_mesh_t               *mesh,
                    const cs_cdo_connect_t        *connect,
                    const cs_cdo_quantities_t     *cdoq,
                    const cs_time_step_t          *time_step,
                    cs_flag_t                      flag,
                    cs_gwf_uspf_t                 *mc)
{
  CS_NO_WARN_IF_UNUSED(flag);

  if (mc == nullptr)
    return;

  cs_equation_t  *richards = mc->richards;
  assert(richards != nullptr);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(!cs_equation_is_steady(richards));

  bool cur2prev = true;

  /* Build and solve the algebraic system related to the Richards equation. By
     default, a current to previous operation is performed */

  cs_equation_solve(cur2prev, mesh, richards);

    /* Update the variables related to the groundwater flow system */

  cs_gwf_uspf_update(mesh, connect, cdoq, time_step,
                     CS_FLAG_CURRENT_TO_PREVIOUS,
                     flag,
                     mc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of single phase flows in an unsaturated porous media
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      post_flag   requested quantities to be postprocessed
 * \param[in, out] mc          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_extra_op(const cs_cdo_connect_t         *connect,
                     const cs_cdo_quantities_t      *cdoq,
                     cs_flag_t                       post_flag,
                     cs_gwf_uspf_t                  *mc)
{
  assert(mc != nullptr);

  if (cs_flag_test(post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE))
    cs_gwf_darcy_flux_balance(connect,
                              cdoq,
                              cs_equation_get_param(mc->richards),
                              mc->darcy);

  if (cs_flag_test(post_flag, CS_GWF_POST_SOIL_STATE))
    cs_gwf_soil_update_soil_state(cdoq->n_cells, mc->moisture_field->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of single-phase flows in an unsaturated porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] mc           pointer to the model context structure
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_uspf_extra_post(int                        mesh_id,
                       cs_lnum_t                  n_cells,
                       const cs_lnum_t            cell_ids[],
                       cs_flag_t                  post_flag,
                       const cs_gwf_uspf_t        *mc,
                       const cs_time_step_t       *time_step)
{
  CS_NO_WARN_IF_UNUSED(n_cells);
  CS_NO_WARN_IF_UNUSED(cell_ids);
  CS_NO_WARN_IF_UNUSED(mc);

  if (mesh_id != CS_POST_MESH_VOLUME)
    return; /* Only postprocessings in the volume are defined */

  /* Note:
     - moisture fied (liquid saturation)
     - capacity field
     - permeability field
     are postprocessed if needed using the standard field process */

  if (post_flag & CS_GWF_POST_SOIL_STATE) {

    const int  *soil_state = cs_gwf_soil_get_soil_state();

    if (soil_state != nullptr)
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "soil_state",
                        1,
                        false, /* interlace */
                        true,  /* use_parent */
                        CS_POST_TYPE_int,
                        soil_state,
                        nullptr,
                        nullptr,
                        time_step);

  } /* Post-processing of the soil state */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
