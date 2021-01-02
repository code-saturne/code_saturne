/*============================================================================
 * Routines to handle the resolution of the turbulence modelling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_equation.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_turbulence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_cdo_turbulence.c
 *
 *  \brief  Routines to handle the resoltion of the turbulence modelling
 *          within the CDO framework
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDO_TURBULENCE_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/* --------------------------------------------------------------------------
 * Context stucture for the k-epsilon turbulence modelling
 * -------------------------------------------------------------------------- */

typedef struct {

  /* Equations */

  cs_equation_t   *tke;
  cs_equation_t   *eps;

  /* Properties associated to the two equations */

  cs_property_t   *tke_diffusivity; /* mu_t/sigma_k */
  cs_property_t   *eps_diffusivity; /* mu_t/sigma_e */
  cs_property_t   *sigma_k;         /* TKE Schmidt  */
  cs_property_t   *sigma_eps;       /* epsilon Schmidt  */

  cs_property_t   *tke_reaction;    /* eps/tke by default + ... if needed */
  cs_property_t   *eps_reaction;    /* by default + ... if needed */

} cs_turb_context_k_eps_t;


/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for defining a quantity at known locations
 *         Here at cells with a function.  elt_ids is optional. If not NULL,
 *         the function works on a sub-list of elements. Moreover, it enables
 *         to fill retval with an indirection if dense_output is set to false
 *         Case of k-epsilon model
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        result of the function. Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_tke_source_term(cs_lnum_t            n_elts,
                 const cs_lnum_t     *elt_ids,
                 bool                 dense_output,
                 void                *input,
                 cs_real_t           *retval)
{
  /* TBD */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for defining a quantity at known locations
 *         Here at cells with a function.  elt_ids is optional. If not NULL,
 *         the function works on a sub-list of elements. Moreover, it enables
 *         to fill retval with an indirection if dense_output is set to false
 *         Case of k-epsilon model with linear production
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        result of the function. Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_tke_lin_source_term(cs_lnum_t            n_elts,
                     const cs_lnum_t     *elt_ids,
                     bool                 dense_output,
                     void                *input,
                     cs_real_t           *retval)
{
  /* TBD */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure storing the set of parameters for the
 *         turbulence modelling
 *
 * \return a pointer to a new allocated cs_turbulence_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_turbulence_param_t *
cs_turbulence_param_create(void)
{
  cs_turbulence_param_t  *tbp = NULL;

  BFT_MALLOC(tbp, 1, cs_turbulence_param_t);

  /* The following structures are shared with the Legacy part */

  tbp->model = cs_get_glob_turb_model();
  tbp->rans_param = cs_get_glob_turb_rans_model();
  tbp->les_param = cs_get_glob_turb_les_model();
  tbp->reference_values = cs_get_glob_turb_ref_values();

  return tbp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure managing the turbulence modelling
 *
 * \param[in] tbp       pointer to a \ref cs_turbulence_param_t structure
 *
 * \return a pointer to a new allocated cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

cs_turbulence_t *
cs_turbulence_create(cs_turbulence_param_t    *tbp)
{
  cs_turbulence_t  *turb = NULL;

  BFT_MALLOC(turb, 1, cs_turbulence_t);

  /* All the members of the following structures are shared with the Legacy
   * part. This structure is owned by cs_navsto_param_t
   */
  turb->param = tbp;

  /* Properties */
  turb->mu_tot = NULL;          /* Total viscosity */
  turb->mu_l = NULL;            /* Laminar viscosity */
  turb->mu_t = NULL;            /* Turbulent viscosity */

  turb->mu_tot_array = NULL;

  /* Fields */
  turb->mu_t_field = NULL;      /* Turbulent viscosity */
  turb->rij = NULL;             /* Reynolds stress tensor */

  /* Main structure (cast on-the-fly according to the turbulence model) */
  turb->context = NULL;

  /* Function pointers */
  turb->init_context = NULL;
  turb->free_context = NULL;
  turb->compute = NULL;
  turb->update = NULL;

  return turb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure managing the turbulence modelling
 *
 * \param[in, out]  p_turb_struct   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_free(cs_turbulence_t   **p_turb_struct)
{
  cs_turbulence_t  *turb = *p_turb_struct;

  /* Set of parameters (members are shared and freed elsewhere).
   * Properties, equations and fields are freed in an other part of the code
   */

  BFT_FREE(turb->mu_tot_array);

  if (turb->free_context != NULL)
    turb->context = turb->free_context(turb->context);

  assert(turb->context == NULL);
  BFT_FREE(turb);
  *p_turb_struct = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the structure managing the turbulence modelling
 *
 * \param[in, out]  turb   pointer to the structure to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_setup(cs_turbulence_t   *turb)
{
  if (turb == NULL)
    return;

  const cs_turbulence_param_t  *tbp = turb->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->iturb == CS_TURB_NONE)
    return; /* Nothing to do if there is a laminar flow */

  /* Set field metadata */
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const bool  has_previous = false;
  const int  field_post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;
  int  location_id = cs_mesh_location_get_id_by_name("cells");

  turb->mu_t_field = cs_field_find_or_create(CS_NAVSTO_TURB_VISCOSITY,
                                             field_mask,
                                             location_id,
                                             1, /* dimension */
                                             has_previous);

  /* Set default value for keys related to log and post-processing */
  cs_field_set_key_int(turb->mu_t_field, log_key, 1);
  cs_field_set_key_int(turb->mu_t_field, post_key, field_post_flag);

  /* Properties (shared) */
  turb->mu_tot = cs_property_by_name(CS_NAVSTO_TOTAL_VISCOSITY);
  turb->mu_l = cs_property_by_name(CS_NAVSTO_LAM_VISCOSITY);

  assert(turb->mu_l != NULL && turb->mu_tot != NULL);

  /* Add a mu_t property and define it with the associated field */
  turb->mu_t = cs_property_add(CS_NAVSTO_TURB_VISCOSITY, CS_PROPERTY_ISO);

  cs_property_def_by_field(turb->mu_t, turb->mu_t_field);

  /* Set function pointers and initialize the context structure */
  switch (model->iturb) {

  case CS_TURB_K_EPSILON:
  case CS_TURB_K_EPSILON_LIN_PROD:
    turb->init_context = cs_turb_init_k_eps_context;
    turb->free_context = cs_turb_free_k_eps_context;
    turb->compute = cs_turb_compute_k_eps;
    turb->update = NULL; /* TBD */

    turb->context = turb->init_context(model);
    break;

  case CS_TURB_NONE:
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid turbulence model with CDO schemes.\n"
              " Please check your settings.", __func__);
    break;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup of the turbulence modelling and especially the
 *         equations/properties and other related quantities
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] tbs        pointer to the turbulence main structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_finalize_setup(const cs_mesh_t            *mesh,
                             const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             const cs_time_step_t       *time_step,
                             cs_turbulence_t            *tbs)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_step);

  if (tbs == NULL)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->iturb == CS_TURB_NONE)
    return; /* Nothing to do */

  /* Define the property related to the total viscosity */
  BFT_MALLOC(tbs->mu_tot_array, quant->n_cells, cs_real_t);
  memset(tbs->mu_tot_array, 0, quant->n_cells*sizeof(cs_real_t));

  cs_property_def_by_array(tbs->mu_tot,
                           cs_flag_primal_cell,
                           tbs->mu_tot_array,
                           false, /* definition is owner ? */
                           NULL); /* no index */

  /* Last setup for each turbulence model */
  switch (model->iturb) {

  case CS_TURB_K_EPSILON:
    {
      /* Add a source term after having retrieve the equation param related to
         the tubulent kinetic energy equation */
      cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
      cs_equation_param_t  *tke_eqp = cs_equation_get_param(kec->tke);

      cs_equation_add_source_term_by_dof_func(tke_eqp,
                                              NULL, /* all cells */
                                              cs_flag_primal_cell,
                                              _tke_source_term,
                                              kec);
    }
    break;

  case CS_TURB_K_EPSILON_LIN_PROD:
    {
      /* Add a source term after having retrieve the equation param related to
         the tubulent kinetic energy equation */
      cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
      cs_equation_param_t  *tke_eqp = cs_equation_get_param(kec->tke);

      cs_equation_add_source_term_by_dof_func(tke_eqp,
                                              NULL, /* all cells */
                                              cs_flag_primal_cell,
                                              _tke_lin_source_term,
                                              kec);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid turbulence model with CDO schemes.\n"
              " Please check your settings.", __func__);
    break;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize quantities related to a turbulence model.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] tbs        pointer to the turbulence main structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_initialize(const cs_mesh_t            *mesh,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *quant,
                         const cs_time_step_t       *time_step,
                         cs_turbulence_t            *tbs)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  if (tbs == NULL)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->iturb == CS_TURB_NONE)
    return; /* Nothing to do */

  /* Initialize the total viscosity */
  const cs_real_t  *mut = tbs->mu_t_field->val;

  if (cs_property_is_uniform(tbs->mu_l)) {

    const cs_real_t  mul0 = cs_property_get_cell_value(0, time_step->t_cur,
                                                       tbs->mu_l);

    for (cs_lnum_t i = 0; i < quant->n_cells; i++)
      tbs->mu_tot_array[i] = mul0 + mut[i];

  }
  else {

    for (cs_lnum_t i = 0; i < quant->n_cells; i++) {

      const cs_real_t  mul = cs_property_get_cell_value(i, time_step->t_cur,
                                                        tbs->mu_l);

      tbs->mu_tot_array[i] = mul + mut[i];

    } /* Loop on cells */

  } /* laminar viscosity is uniform ? */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the context structure related to the
 *         k-epsilon turbulence model
 *
 * \param[in]  tbm         structure which defines the turbulence model
 *
 * \return a pointer to a new allocated context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_turb_init_k_eps_context(const cs_turb_model_t      *tbm)
{
  if (tbm == NULL)
    return NULL;

  /* Sanity checks */
  assert((tbm->iturb == CS_TURB_K_EPSILON) ||
         (tbm->iturb == CS_TURB_K_EPSILON_LIN_PROD));
  assert(tbm->type == CS_TURB_RANS);
  assert(tbm->order == CS_TURB_FIRST_ORDER);

  cs_turb_context_k_eps_t  *kec = NULL;

  BFT_MALLOC(kec, 1, cs_turb_context_k_eps_t);

  /* Add new equations for the tubulent kinetic energy (tke) and the dissipation
     (epsilon) */

  kec->tke = cs_equation_add("k", /* equation name */
                             "k", /* variable name */
                             CS_EQUATION_TYPE_NAVSTO,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  kec->eps = cs_equation_add("epsilon", /* equation name */
                             "epsilon", /* variable name */
                             CS_EQUATION_TYPE_NAVSTO,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  /* Add new related properties which will be associated to discretization terms
     in tke and epsilon */

  kec->tke_diffusivity = cs_property_add("k_diffusivity",
                                         CS_PROPERTY_ISO);

  kec->eps_diffusivity = cs_property_add("epsilon_diffusivity",
                                         CS_PROPERTY_ISO);

  /* Turbulent Schmidt coefficients : creation and set the reference value */

  kec->sigma_k = cs_property_add("k_turb_schmidt",
                                 CS_PROPERTY_ISO);
  cs_property_set_reference_value(kec->sigma_k, 1.0);

  kec->sigma_eps = cs_property_add("epsilon_turb_schmidt",
                                   CS_PROPERTY_ISO);
  cs_property_set_reference_value(kec->sigma_eps, 1.3);

  /* Reaction (implicit source terms) coefficients */

  kec->tke_reaction = cs_property_add("k_reaction",
                                                CS_PROPERTY_ISO);

  kec->eps_reaction = cs_property_add("epsilon_reaction",
                                                CS_PROPERTY_ISO);

  /* Retrieve the mass density */

  cs_property_t  *mass_density = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);

  /* Retrieve the advection field from Navier--Stokes (the mass flux) */

  cs_adv_field_t  *adv = cs_advection_field_by_name("mass_flux");
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Advection field not found.", __func__);

  /* Add terms to the TKE equation */

  cs_equation_param_t  *tke_eqp = cs_equation_get_param(kec->tke);

  cs_equation_add_time(tke_eqp, mass_density);
  cs_equation_add_diffusion(tke_eqp, kec->tke_diffusivity);
  cs_equation_add_reaction(tke_eqp, kec->tke_reaction);
  cs_equation_add_advection(tke_eqp, adv);

  /* Source term is defined elsewhere since it depends on the choice of the
   * sub-model */

  /* Add terms to the epsilon equation */

  cs_equation_param_t  *eps_eqp = cs_equation_get_param(kec->eps);

  cs_equation_add_time(eps_eqp, mass_density);
  cs_equation_add_diffusion(eps_eqp, kec->eps_diffusivity);
  cs_equation_add_reaction(tke_eqp, kec->tke_reaction);
  cs_equation_add_advection(tke_eqp, adv);

  return kec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to the k-epsilon turbulence model
 *
 * \param[in, out]  tbc   pointer to a structure context cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_turb_free_k_eps_context(void     *tbc)
{
  cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbc;

  if (kec == NULL)
    return kec;

  BFT_FREE(kec);

  return kec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the turbulence
 *         model. This means that all related equations are built and then
 *         solved.
 *
 * \param[in]      mesh      pointer to a \ref cs_mesh_t structure
 * \param[in]      tbp       pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] tbc       pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_k_eps(const cs_mesh_t              *mesh,
                      const cs_turbulence_param_t  *tpb,
                      void                         *tbc)
{
  if (tbc == NULL)
    return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
