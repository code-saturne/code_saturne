/*============================================================================
 * Functions to handle the resolution of the turbulence modelling
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

#include "cs_post.h"
#include "cs_turbulence_model.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_turbulence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_cdo_turbulence.c
 *
 *  \brief  Functions to handle the resolution of the turbulence modelling
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

  cs_xdef_t       *tke_reaction;    /* eps/tke by default + ... if needed */
  cs_xdef_t       *eps_reaction;    /* by default + ... if needed */

  cs_xdef_t       *tke_source_term; /* Production + buoyancy if needed for k */
  cs_xdef_t       *eps_source_term; /* Same for epsilon */

} cs_turb_context_k_eps_t;

/* --------------------------------------------------------------------------
 * Context structure for the k-epsilon turbulence modelling
 * -------------------------------------------------------------------------- */

typedef struct {

  /* High level structures */

  const cs_mesh_t            *mesh;
  const cs_cdo_connect_t     *connect;
  const cs_cdo_quantities_t  *quant;
  const cs_time_step_t       *time_step;

  /* Turbulence structure */

  cs_turbulence_t  *tbs;

  /* Velocities arrays */

  cs_real_t   *u_cell;
  cs_real_t   *u_face;

} cs_turb_source_term_t;


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
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_ke(const cs_mesh_t            *mesh,
            const cs_cdo_connect_t     *connect,
            const cs_cdo_quantities_t  *quant,
            const cs_time_step_t       *time_step,
            const cs_turbulence_t      *tbs,
            cs_real_t                  *tke_reaction,
            cs_real_t                  *eps_reaction,
            cs_real_t                  *tke_source_term,
            cs_real_t                  *eps_source_term)
{
  const cs_turbulence_param_t  *tbp = tbs->param;
  cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;

  const cs_real_t *u_cell = cs_equation_get_cell_values(tbs->mom_eq, false);
  const cs_real_t *u_face = cs_equation_get_face_values(tbs->mom_eq, false);

  cs_real_t *mu_t = tbs->mu_t_field->val;

  /* Get the evaluation of rho */

  int rho_stride = 0;
  cs_real_t *rho = NULL;

  const cs_real_t *tke_cell = cs_equation_get_cell_values(kec->tke, false);
  const cs_real_t *eps_cell = cs_equation_get_cell_values(kec->eps, false);

  /* Get mass density values in each cell */

  cs_property_iso_get_cell_values(time_step->t_cur, tbs->rho,
                                  &rho_stride, &rho);

  const cs_real_t d1s3 = 1./3.;
  const cs_real_t d2s3 = 2./3.;

  /* Production term */

  if (tbp->model->iturb == CS_TURB_K_EPSILON) {
#   pragma omp parallel for if (mesh->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

      /* Compute the velocity gradient */

      cs_real_t grd_uc[3][3];
      cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                        u_cell, u_face, (cs_real_t *)grd_uc);

      cs_real_t strain_sq = 2.
        * (  cs_math_pow2(  d2s3*grd_uc[0][0]
                          - d1s3*grd_uc[1][1]
                          - d1s3*grd_uc[2][2])
           + cs_math_pow2(- d1s3*grd_uc[0][0]
                          + d2s3*grd_uc[1][1]
                          - d1s3*grd_uc[2][2])
           + cs_math_pow2(- d1s3*grd_uc[0][0]
                          - d1s3*grd_uc[1][1]
                          + d2s3*grd_uc[2][2]))
        + cs_math_pow2(grd_uc[0][1] + grd_uc[1][0])
        + cs_math_pow2(grd_uc[0][2] + grd_uc[2][0])
        + cs_math_pow2(grd_uc[1][2] + grd_uc[2][1]);

      tke_source_term[c_id] = mu_t[c_id] * strain_sq;
    }
  }
  else if (tbp->model->iturb == CS_TURB_K_EPSILON_LIN_PROD) {
#   pragma omp parallel for if (mesh->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

      /* Compute the velocity gradient */

      cs_real_t grd_uc[3][3];
      cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                        u_cell, u_face, (cs_real_t *)grd_uc);

      cs_real_t strain_sq = 2.
        * (  cs_math_pow2(  d2s3*grd_uc[0][0]
                          - d1s3*grd_uc[1][1]
                          - d1s3*grd_uc[2][2])
           + cs_math_pow2(- d1s3*grd_uc[0][0]
                          + d2s3*grd_uc[1][1]
                          - d1s3*grd_uc[2][2])
           + cs_math_pow2(- d1s3*grd_uc[0][0]
                          - d1s3*grd_uc[1][1]
                          + d2s3*grd_uc[2][2]))
        + cs_math_pow2(grd_uc[0][1] + grd_uc[1][0])
        + cs_math_pow2(grd_uc[0][2] + grd_uc[2][0])
        + cs_math_pow2(grd_uc[1][2] + grd_uc[2][1]);

      cs_real_t strain = sqrt(strain_sq);
      const cs_real_t sqrcmu = sqrt(cs_turb_cmu);
      cs_real_t cmueta =
        fmin(cs_turb_cmu*tke_cell[c_id]/eps_cell[c_id] * strain, sqrcmu);

      tke_source_term[c_id] = rho[c_id*rho_stride] * cmueta * strain * tke_cell[c_id];

    }
  }

  /* Implicit dissipation term and epsilon source term */

# pragma omp parallel for if (mesh->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

    /* Inverse integral time scale */

    cs_real_t d_ttke = eps_cell[c_id] / tke_cell[c_id];

    /* Ce1 * eps/k * P */

    eps_source_term[c_id] = cs_turb_ce1 * d_ttke * tke_source_term[c_id];
    tke_reaction[c_id] = rho[c_id*rho_stride] * d_ttke;

    /* TODO Get Ce2 from curvature correction, to be merged with legacy */

    eps_reaction[c_id] = cs_turb_ce2 * rho[c_id*rho_stride] * d_ttke;

  }
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
  cs_turbulence_t  *tbs = NULL;

  BFT_MALLOC(tbs, 1, cs_turbulence_t);

  /* All the members of the following structures are shared with the Legacy
   * part. This structure is owned by cs_navsto_param_t
   */

  tbs->param = tbp;
  tbs->mom_eq = NULL;

  /* Properties */

  tbs->rho = NULL;             /* Mass density */
  tbs->mu_tot = NULL;          /* Total viscosity */
  tbs->mu_l = NULL;            /* Laminar viscosity */
  tbs->mu_t = NULL;            /* Turbulent viscosity */

  tbs->mu_tot_array = NULL;

  /* Fields */

  tbs->mu_t_field = NULL;      /* Turbulent viscosity */
  tbs->rij = NULL;             /* Reynolds stress tensor */

  /* Main structure (cast on-the-fly according to the turbulence model) */

  tbs->context = NULL;

  /* Function pointers */

  tbs->init_context = NULL;
  tbs->free_context = NULL;
  tbs->compute = NULL;
  tbs->compute_steady = NULL;
  tbs->update = NULL;

  return tbs;
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
  cs_turbulence_t  *tbs = *p_turb_struct;

  /* Set of parameters (members are shared and freed elsewhere).
   * Properties, equations and fields are freed in an other part of the code
   */

  BFT_FREE(tbs->mu_tot_array);

  if (tbs->free_context != NULL)
    tbs->context = tbs->free_context(tbs->context);

  assert(tbs->context == NULL);
  BFT_FREE(tbs);
  *p_turb_struct = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the structure managing the turbulence modelling
 *
 * \param[in, out]  tbs     pointer to the structure to initialize
 * \param[in]       mom_eq  pointer to the structure mom_eq
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_setup(cs_turbulence_t     *tbs,
                         cs_equation_t       *mom_eq)
{
  if (tbs == NULL)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->type == CS_TURB_NONE)
    return; /* Nothing to do if there is a laminar flow */

  tbs->mom_eq = mom_eq;

  /* Set field metadata */

  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const bool  has_previous = false;
  const int  field_post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;
  int  location_id = cs_mesh_location_get_id_by_name("cells");

  tbs->mu_t_field = cs_field_find_or_create(CS_NAVSTO_TURB_VISCOSITY,
                                            field_mask,
                                            location_id,
                                            1, /* dimension */
                                            has_previous);

  /* Set default value for keys related to log and post-processing */

  cs_field_set_key_int(tbs->mu_t_field, log_key, 1);
  cs_field_set_key_int(tbs->mu_t_field, post_key, field_post_flag);

  /* Properties (shared) */

  tbs->rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  tbs->mu_tot = cs_property_by_name(CS_NAVSTO_TOTAL_VISCOSITY);
  tbs->mu_l = cs_property_by_name(CS_NAVSTO_LAM_VISCOSITY);

  assert(tbs->rho != NULL && tbs->mu_l != NULL && tbs->mu_tot != NULL);

  /* Add a mu_t property and define it with the associated field */

  tbs->mu_t = cs_property_add(CS_NAVSTO_TURB_VISCOSITY, CS_PROPERTY_ISO);

  cs_property_def_by_field(tbs->mu_t, tbs->mu_t_field);

  /* Set function pointers and initialize the context structure */

  switch (model->iturb) {

  case CS_TURB_K_EPSILON:
  case CS_TURB_K_EPSILON_LIN_PROD:
    tbs->init_context = cs_turb_init_k_eps_context;
    tbs->free_context = cs_turb_free_k_eps_context;
    tbs->compute = cs_turb_compute_k_eps;
    tbs->update = cs_turb_update_k_eps;

    tbs->context = tbs->init_context(model);
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
  CS_UNUSED(time_step);

  if (tbs == NULL)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->type == CS_TURB_NONE)
    return; /* Nothing to do */

  /* Define the property related to the total viscosity */

  BFT_MALLOC(tbs->mu_tot_array, quant->n_cells, cs_real_t);
  memset(tbs->mu_tot_array, 0, quant->n_cells*sizeof(cs_real_t));

  cs_property_def_by_array(tbs->mu_tot,
                           cs_flag_primal_cell,
                           tbs->mu_tot_array,
                           false, /* definition is owner ? */
                           NULL, NULL); /* no index/ids */

  /* Last setup for each turbulence model */

  switch (model->iturb) {

  case CS_TURB_K_EPSILON:
  case CS_TURB_K_EPSILON_LIN_PROD:
    {
      /* Add a source term after having retrieved the equation param related to
         the turbulent kinetic energy equation */

      cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
      cs_equation_param_t  *tke_eqp = cs_equation_get_param(kec->tke);
      kec->tke_source_term =
        cs_equation_add_source_term_by_array(tke_eqp,
                                             NULL, /* all cells */
                                             cs_flag_primal_cell,
                                             NULL,
                                             false, /*is owner */
                                             NULL, NULL); /* no index/ids */

      kec->tke_reaction =
        cs_property_def_by_array(cs_property_by_name("k_reaction"),
                                 cs_flag_primal_cell,
                                 NULL,
                                 false, /* definition is owner ? */
                                 NULL, NULL); /* no index/ids */

      cs_equation_param_t  *eps_eqp = cs_equation_get_param(kec->eps);
      kec->eps_source_term =
        cs_equation_add_source_term_by_array(eps_eqp,
                                             NULL, /* all cells */
                                             cs_flag_primal_cell,
                                             NULL,
                                             false, /*is owner */
                                             NULL, NULL); /* no index/ids */

      kec->eps_reaction =
        cs_property_def_by_array(cs_property_by_name("eps_reaction"),
                                 cs_flag_primal_cell,
                                 NULL,
                                 false, /* definition is owner ? */
                                 NULL, NULL); /* no index/ids */

      cs_property_def_by_array(tbs->mu_tot,
                               cs_flag_primal_cell,
                               tbs->mu_tot_array,
                               false, /* definition is owner ? */
                               NULL, NULL); /* no index/ids */

      /* Initialize TKE */

      cs_turb_ref_values_t *t_ref= cs_get_glob_turb_ref_values();
      cs_real_t tke_ref = 1.5 * cs_math_pow2(0.02 * t_ref->uref);
      cs_equation_add_ic_by_value(tke_eqp,
                                  NULL,
                                  &tke_ref);

      /* Initialize epsilon */

      cs_real_t eps_ref = powf(tke_ref, 1.5) * cs_turb_cmu / t_ref->almax;
      cs_equation_add_ic_by_value(eps_eqp,
                                  NULL,
                                  &eps_ref);
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
 * \brief  Initialize the values of quantities related to a turbulence model.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] tbs        pointer to the turbulence main structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_values(const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *time_step,
                          cs_turbulence_t             *tbs)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  if (tbs == NULL)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->iturb == CS_TURB_NONE)
    return; /* Nothing to do */

  tbs->update(mesh, connect, quant, time_step, tbs);
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

  assert((tbm->iturb == CS_TURB_K_EPSILON) ||
         (tbm->iturb == CS_TURB_K_EPSILON_LIN_PROD));
  assert(tbm->type == CS_TURB_RANS);
  assert(tbm->order == CS_TURB_FIRST_ORDER);

  cs_turb_context_k_eps_t  *kec = NULL;

  BFT_MALLOC(kec, 1, cs_turb_context_k_eps_t);

  /* Add new equations for the turbulent kinetic energy (tke) and the
     dissipation (epsilon) */

  kec->tke = cs_equation_add("k", /* equation name */
                             "k", /* variable name */
                             CS_EQUATION_TYPE_NAVSTO, /* related to NS */
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  kec->eps = cs_equation_add("epsilon", /* equation name */
                             "epsilon", /* variable name */
                             CS_EQUATION_TYPE_NAVSTO,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  /* Add new related properties which will be associated to discretization
     terms in tke and epsilon */

  kec->tke_diffusivity = cs_property_add("k_diffusivity", CS_PROPERTY_ISO);

  kec->eps_diffusivity = cs_property_add("epsilon_diffusivity",
                                         CS_PROPERTY_ISO);

  /* Turbulent Schmidt coefficients : creation and set the reference value */

  kec->sigma_k = cs_property_add("k_turb_schmidt", CS_PROPERTY_ISO);
  cs_property_set_reference_value(kec->sigma_k, 1.0);

  kec->sigma_eps = cs_property_add("epsilon_turb_schmidt", CS_PROPERTY_ISO);
  cs_property_set_reference_value(kec->sigma_eps, 1.3);

  /* Reaction (implicit source terms) coefficients */

  cs_property_t *k_reaction
    = cs_property_add("k_reaction", CS_PROPERTY_ISO);
  cs_property_t *eps_reaction
    = cs_property_add("epsilon_reaction", CS_PROPERTY_ISO);

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
  cs_equation_add_reaction(tke_eqp, k_reaction);
  cs_equation_add_advection(tke_eqp, adv);

  /* Source term is defined elsewhere since it depends on the choice of the
   * sub-model */

  /* Add terms to the epsilon equation */

  cs_equation_param_t  *eps_eqp = cs_equation_get_param(kec->eps);

  cs_equation_add_time(eps_eqp, mass_density);
  cs_equation_add_diffusion(eps_eqp, kec->eps_diffusivity);
  cs_equation_add_reaction(eps_eqp, eps_reaction);
  cs_equation_add_advection(eps_eqp, adv);

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
 * \brief  Update for the current time step the new state for the turbulence
 *         model. This is used to update the turbulent viscosity.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in]      tbs        pointer to a \ref cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_update_k_eps(const cs_mesh_t              *mesh,
                     const cs_cdo_connect_t       *connect,
                     const cs_cdo_quantities_t    *quant,
                     const cs_time_step_t         *time_step,
                     const cs_turbulence_t        *tbs)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  if (tbs == NULL)
    return;

  cs_lnum_t n_cells = mesh->n_cells;
  CS_UNUSED(n_cells); /* avoid a compiler warning without OpenMP */

  cs_turb_context_k_eps_t  *kec =
    (cs_turb_context_k_eps_t *)tbs->context;

  /* Update turbulent viscosity field */

  cs_real_t *mu_t = tbs->mu_t_field->val;
  cs_real_t *k = cs_equation_get_cell_values(kec->tke, false);
  cs_real_t *eps = cs_equation_get_cell_values(kec->eps, false);

  /* Get the evaluation of rho */

  int rho_stride = 0;
  cs_real_t *rho = NULL;

  /* Get mass density values in each cell */

  cs_property_iso_get_cell_values(time_step->t_cur, tbs->rho,
                                  &rho_stride, &rho);

  /* Get laminar viscosity values in each cell */

  int mu_stride = 0;
  cs_real_t *mu_l = NULL;
  cs_property_iso_get_cell_values(time_step->t_cur, tbs->mu_l,
                                  &mu_stride, &mu_l);


  /* Compute mu_t in each cell and define mu_tot = mu_t + mu_l */

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {

    mu_t[cell_id] = cs_turb_cmu * rho[cell_id*rho_stride] *
      cs_math_pow2(k[cell_id]) / eps[cell_id];

    tbs->mu_tot_array[cell_id] = mu_t[cell_id] + mu_l[cell_id*mu_stride];

  }

  /* Free memory */

  BFT_FREE(rho);
  BFT_FREE(mu_l);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for the current time step the new state for the turbulence
 *         model. This means that all related equations are built and then
 *         solved.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in, out] tbs        pointer to turbulence structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_k_eps(const cs_mesh_t            *mesh,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      const cs_time_step_t       *time_step,
                      cs_turbulence_t            *tbs)
{
  if (tbs == NULL)
    return;

  /* Get k epsilon context */

  cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
  cs_equation_t *tke_eq = kec->tke;
  cs_equation_t *eps_eq = kec->eps;
  assert(kec != NULL);

  /* Prepare source term and reaction term */

  cs_real_t *tke_source_term = NULL, *eps_source_term = NULL;
  cs_real_t *tke_reaction = NULL, *eps_reaction = NULL;
  BFT_MALLOC(tke_source_term, mesh->n_cells, cs_real_t);
  BFT_MALLOC(eps_source_term, mesh->n_cells, cs_real_t);
  BFT_MALLOC(tke_reaction, mesh->n_cells, cs_real_t);
  BFT_MALLOC(eps_reaction, mesh->n_cells, cs_real_t);

  /* Set xdefs */

  cs_xdef_set_array(kec->tke_reaction,
                    false, /* is_owner */
                    tke_reaction);


  cs_xdef_set_array(kec->eps_reaction,
                    false, /* is_owner */
                    eps_reaction);

  cs_xdef_set_array(kec->tke_source_term,
                    false, /* is_owner */
                    tke_source_term);

  cs_xdef_set_array(kec->eps_source_term,
                    false, /* is_owner */
                    eps_source_term);

  _prepare_ke(mesh,
              connect,
              quant,
              time_step,
              tbs,
              tke_reaction,
              eps_reaction,
              tke_source_term,
              eps_source_term);

  /* Solve k */

  cs_equation_solve(true, /* cur2prev */
                    mesh,
                    tke_eq);

  /* Solve epsilon */

  cs_equation_solve(true, /* cur2prev */
                    mesh,
                    eps_eq);

  /* Free memory */

  BFT_FREE(tke_source_term);
  BFT_FREE(eps_source_term);
  BFT_FREE(tke_reaction);
  BFT_FREE(eps_reaction);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
