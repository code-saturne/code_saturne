/*============================================================================
 * Main functions dedicated to groundwater flows in case of single-phase flows
 * in a saturated porous media (SSPF in short)
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
#include "cs_field.h"
#include "cs_gwf_priv.h"
#include "cs_gwf_soil.h"
#include "cs_log.h"
#include "cs_param_types.h"
#include "cs_post.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_sspf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_sspf.c

  \brief Main high-level functions dedicated to groundwater flows when using
         CDO schemes for single-phase flows in a saturated porous media.
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_SSPF_DBG 0

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
 *        Case of CDO-Vb schemes and a Darcy flux defined at dual faces for
 *        a saturated porous media.
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
_sspf_update_darcy_arrays(const cs_cdo_connect_t      *connect,
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

  cs_gwf_sspf_t  *mc = (cs_gwf_sspf_t *)darcy->update_input;
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
 * \brief Allocate and initialize the model context structure in case of
 *        single-phase flows in a saturated porous media
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_sspf_t *
cs_gwf_sspf_create(void)
{
  cs_gwf_sspf_t *mc = nullptr;

  BFT_MALLOC(mc, 1, cs_gwf_sspf_t);

  mc->pressure_head = nullptr;

  /* Create a new equation structure for Richards' equation */

  mc->richards = cs_equation_add("Richards",       /* equation name */
                                 "hydraulic_head", /* variable name */
                                 CS_EQUATION_TYPE_GROUNDWATER,
                                 1,
                                 CS_BC_SYMMETRY);

  /* Add a property related to the moisture content (It should be a constant
     in each soil). This is closely related to the soil porosity. */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  /* Define the Darcy flux structure.

     Add an advection field related to the darcian flux stemming from the
     Richards equation. This advection field is steady since the head is
     steady-state in this model. The flux is defined at the dual faces (split
     by primal cells)
  */

  mc->darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF |
    CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX |
    CS_ADVECTION_FIELD_STEADY;

  mc->darcy->adv_field = cs_advection_field_add("darcy_field", adv_status);

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the model context structure in case of a saturated single-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_free(cs_gwf_sspf_t   **p_mc)
{
  if (p_mc == nullptr)
    return;
  if ((*p_mc) == nullptr)
    return;

  cs_gwf_sspf_t  *mc = *p_mc;

  cs_gwf_darcy_flux_free(&(mc->darcy));

  BFT_FREE(mc);
  *p_mc = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of saturated single-phase
 *        flows in porous media
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_log_setup(cs_gwf_sspf_t   *mc)
{
  if (mc == nullptr)
    return;

  cs_gwf_darcy_flux_log(mc->darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a saturated single-phase flows model in porous media
 *
 * \param[in, out] mc          pointer to the model context structure
 * \param[in]      abs_perm    property struct. for the absolute permeability
 * \param[in]      flag        optional metadata
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_init(cs_gwf_sspf_t       *mc,
                 cs_property_t       *abs_perm,
                 cs_flag_t            flag)
{
  if (mc == nullptr)
    return;

  cs_equation_t  *eq = mc->richards;

  if (eq == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);
  assert(cs_equation_is_steady(eq) == true);

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  assert(eqp != nullptr);

  /* Add the diffusion term to the Richards equation by associating the
     absolute permeability to the diffusion property of the Richards eq. */

  cs_equation_add_diffusion(eqp, abs_perm);

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

  /* Add the variable field */

  if (flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    cs_equation_predefined_create_field(1, eq); /* Keep two states */
  else
    cs_equation_predefined_create_field(0, eq); /* Keep only one state */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for saturated single-phase flows in a porous
 *        media. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in]      flag         optional settings for the module
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_init_setup(cs_flag_t           flag,
                       cs_gwf_sspf_t      *mc)
{
  if (mc == nullptr)
    return;

  cs_equation_t  *eq = mc->richards;
  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  assert(eqp != nullptr);
  assert(cs_equation_is_steady(eq) == true);

  /* Set the "has_previous" variable */

  bool  has_previous = false;
  if (flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    has_previous = true;

  /* Add new fields if needed */

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
                                          has_previous);
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          has_previous);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(mc->pressure_head, log_key, 1);
    cs_field_set_key_int(mc->pressure_head, post_key, 1);

  } /* Gravitation effects are activated */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of single-phase flows in a saturated
 *        porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_finalize_setup(const cs_cdo_connect_t        *connect,
                           const cs_cdo_quantities_t     *cdoq,
                           cs_gwf_sspf_t                 *mc)
{
  if (mc == nullptr)
    return;

  /* Set the Darcian flux (in the volume and at the boundary) */

  cs_gwf_darcy_flux_define(connect,
                           cdoq,
                           cs_equation_get_space_scheme(mc->richards),
                           mc,                        /* context */
                           _sspf_update_darcy_arrays, /* update function */
                           mc->darcy);

  /* Set the soil porosity each soil from the moisture content */

  cs_gwf_soil_define_sspf_property(mc->moisture_content);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of single-phase flows in a
 *        saturated porous media
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  type of operation(s) to do
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_update(const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq,
                   const cs_time_step_t         *ts,
                   cs_flag_t                     update_flag,
                   cs_flag_t                     option_flag,
                   cs_gwf_sspf_t                *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = cs_equation_get_time_eval(ts, mc->richards);
    cur2prev = true;

  }

  /* Update the pressure head */

  cs_gwf_update_head(connect,
                     cdoq,
                     mc->richards,
                     option_flag,
                     mc->pressure_head,
                     nullptr, /* there is no head_in_law for saturated soils */
                     cur2prev);

  /* Update the advection field (the Darcy flux related to the groundwater flow
     module) */

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_SSPF_DBG > 2
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

  /* Other properties are constant in saturated soils. Therefore, there is no
     need to update them (e.g. permeability and moisture content) */

  /* Check if an update of some soil properties has to be done with the new
     head values for instance if a user-defined model has been set. */

  cs_gwf_soil_update(time_eval, mesh, connect, cdoq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the steady-state of the groundwater flows module in case of
 *        single-phase flows in a saturated porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      flag       type of additional treatment(s) to do
 * \param[in, out] mc         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_compute_steady_state(const cs_mesh_t             *mesh,
                                 const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *cdoq,
                                 const cs_time_step_t        *time_step,
                                 cs_flag_t                    flag,
                                 cs_gwf_sspf_t               *mc)
{
  if (mc == nullptr)
    return;

  cs_equation_t  *richards = mc->richards;
  assert(richards != nullptr);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  if (cs_equation_is_steady(richards) ||
      flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

    /* Build and solve the linear system related to the Richards equations */

    cs_equation_solve_steady_state(mesh, richards);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_sspf_update(mesh, connect, cdoq, time_step,
                       CS_FLAG_CURRENT_TO_PREVIOUS,
                       flag,
                       mc);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new hydraulic state for the groundwater flows module.
 *        Case of single-phase flows in a saturated porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      flag       type of additional treatment(s) to do
 * \param[in, out] mc         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_compute(const cs_mesh_t              *mesh,
                    const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    const cs_time_step_t         *time_step,
                    cs_flag_t                     flag,
                    cs_gwf_sspf_t                *mc)
{
  if (mc == nullptr)
    return;

  cs_equation_t  *richards = mc->richards;
  assert(richards != nullptr);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  bool cur2prev = true;

  if (!cs_equation_is_steady(richards)) {

    /* Build and solve the linear system related to the Richards equations. By
       default, a current to previous operation is performed. */

    cs_equation_solve(cur2prev, mesh, richards);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_sspf_update(mesh, connect, cdoq, time_step,
                       CS_FLAG_CURRENT_TO_PREVIOUS,
                       flag,
                       mc);

  }
  else {

    /* Richards is steady but one can force the resolution */

    if (flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

      /* Build and solve the linear system related to the Richards equations */

      cs_equation_solve_steady_state(mesh, richards);

      /* Update the variables related to the groundwater flow system */

      cs_gwf_sspf_update(mesh, connect, cdoq, time_step,
                         CS_FLAG_CURRENT_TO_PREVIOUS,
                         flag,
                         mc);

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of single phase flows in a saturated porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_extra_op(const cs_cdo_connect_t         *connect,
                     const cs_cdo_quantities_t      *cdoq,
                     cs_flag_t                       post_flag,
                     cs_gwf_sspf_t                  *mc)
{
  assert(mc != nullptr);

  if (cs_flag_test(post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  cs_gwf_darcy_flux_balance(connect, cdoq,
                            cs_equation_get_param(mc->richards),
                            mc->darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of single-phase flows in a saturated porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] abs_perm     property for the absolute permeability
 * \param[in] mc           pointer to the model context structure
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_sspf_extra_post(int                        mesh_id,
                       cs_lnum_t                  n_cells,
                       const cs_lnum_t            cell_ids[],
                       cs_flag_t                  post_flag,
                       const cs_property_t       *abs_perm,
                       const cs_gwf_sspf_t       *mc,
                       const cs_time_step_t      *time_step)
{
  if (mesh_id != CS_POST_MESH_VOLUME)
    return; /* Only postprocessings in the volume are defined */

  assert(mc != nullptr);

  if (post_flag & CS_GWF_POST_PERMEABILITY) {

    cs_real_t *permeability = nullptr;
    int  dim = cs_property_get_dim(abs_perm);
    int  post_dim = (dim == 1) ? 1 : 9;

    if (dim > 1) {

      BFT_MALLOC(permeability, post_dim*n_cells, cs_real_t);

      for (cs_lnum_t i = 0; i < n_cells; i++) {

        cs_real_t  tensor[3][3];

        cs_property_get_cell_tensor(cell_ids[i],
                                    time_step->t_cur,
                                    abs_perm,
                                    false, /* inversion */
                                    tensor);

        cs_real_t  *_cell_perm = permeability + post_dim*i;
        for (int ki = 0; ki < 3; ki++)
          for (int kj = 0; kj < 3; kj++)
            _cell_perm[3*ki+kj] = tensor[ki][kj];

      }

    }
    else {

      BFT_MALLOC(permeability, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        permeability[c_id] = cs_property_get_cell_value(cell_ids[c_id],
                                                        time_step->t_cur,
                                                        abs_perm);

    }

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "permeability",
                      post_dim,
                      true,  /* interlace */
                      false, /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      permeability,
                      nullptr,
                      nullptr,
                      time_step);

    BFT_FREE(permeability);

  } /* Postprocessing of the permeability */

  if (post_flag & CS_GWF_POST_LIQUID_SATURATION) {

    cs_real_t *liquid_saturation = nullptr;
    BFT_MALLOC(liquid_saturation, n_cells, cs_real_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      liquid_saturation[c_id] =
        cs_property_get_cell_value(cell_ids[c_id],
                                   time_step->t_cur,
                                   mc->moisture_content);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "liquid_saturation",
                      1,
                      false, /* interlaced */
                      false, /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      liquid_saturation,
                      nullptr,
                      nullptr,
                      time_step);

    BFT_FREE(liquid_saturation);

  } /* Postprocessing of the liquid saturation */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
