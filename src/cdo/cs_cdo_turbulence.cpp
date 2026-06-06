/*============================================================================
 * Functions to handle the resolution of the turbulence modelling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <cassert>
#include <cstring>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_mem.h"
#include "base/cs_post.h"
#include "base/cs_wall_functions.h"
#include "cdo/cs_reco.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdo_turbulence.h"

/*----------------------------------------------------------------------------*/

/*!
 *  \file cs_cdo_turbulence.cpp
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
 *         Here at cells with a function.  elt_ids is optional. If non-null,
 *         the function works on a sub-list of elements. Moreover, it enables
 *         to fill retval with an indirection if dense_output is set to false
 *         Case of k-epsilon model
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
  cs_real_t *rho        = nullptr;

  const cs_real_t *tke_cell = cs_equation_get_cell_values(kec->tke, false);
  const cs_real_t *eps_cell = cs_equation_get_cell_values(kec->eps, false);

  /* Get mass density values in each cell */

  cs_property_iso_get_cell_values(time_step->t_cur, tbs->rho,
                                  &rho_stride, &rho);

  const cs_real_t d1s3 = 1./3.;
  const cs_real_t d2s3 = 2./3.;

  /* Production term */

  if (tbp->model->model == CS_TURB_K_EPSILON) {
#   pragma omp parallel for if (mesh->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

      /* Compute the velocity gradient */

      cs_real_t grd_uc[3][3];
      cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                        u_cell, u_face, (cs_real_t *)grd_uc);

      cs_real_t strain_sq =
        2. * (cs::pow2(d2s3 * grd_uc[0][0] - d1s3 * grd_uc[1][1] -
                       d1s3 * grd_uc[2][2]) +
              cs::pow2(-d1s3 * grd_uc[0][0] + d2s3 * grd_uc[1][1] -
                       d1s3 * grd_uc[2][2]) +
              cs::pow2(-d1s3 * grd_uc[0][0] - d1s3 * grd_uc[1][1] +
                       d2s3 * grd_uc[2][2])) +
        cs::pow2(grd_uc[0][1] + grd_uc[1][0]) +
        cs::pow2(grd_uc[0][2] + grd_uc[2][0]) +
        cs::pow2(grd_uc[1][2] + grd_uc[2][1]);

      tke_source_term[c_id] = mu_t[c_id] * strain_sq;
    }
  }
  else if (tbp->model->model == CS_TURB_K_EPSILON_LIN_PROD) {
#   pragma omp parallel for if (mesh->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

      /* Compute the velocity gradient */

      cs_real_t grd_uc[3][3];
      cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                        u_cell, u_face, (cs_real_t *)grd_uc);

      cs_real_t strain_sq =
        2. * (cs::pow2(d2s3 * grd_uc[0][0] - d1s3 * grd_uc[1][1] -
                       d1s3 * grd_uc[2][2]) +
              cs::pow2(-d1s3 * grd_uc[0][0] + d2s3 * grd_uc[1][1] -
                       d1s3 * grd_uc[2][2]) +
              cs::pow2(-d1s3 * grd_uc[0][0] - d1s3 * grd_uc[1][1] +
                       d2s3 * grd_uc[2][2])) +
        cs::pow2(grd_uc[0][1] + grd_uc[1][0]) +
        cs::pow2(grd_uc[0][2] + grd_uc[2][0]) +
        cs::pow2(grd_uc[1][2] + grd_uc[2][1]);

      cs_real_t       strain = std::sqrt(strain_sq);
      const cs_real_t sqrcmu = std::sqrt(cs_turb_cmu);
      cs_real_t       cmueta =
        cs::min(cs_turb_cmu * tke_cell[c_id] / eps_cell[c_id] * strain, sqrcmu);

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function used to define the exchange coefficients for tangential and
 *        normal components. (1-velocity scales method)
 *
 * \param[in]      pena_bc_coeff  pena coefficient used as Dirichlet BC
 * \param[in]      nu             laminar kinematic viscosity
 * \param[in]      k              turbulent kinetic energy
 * \param[in]      hfc            distance from cell center to the wall
 * \param[in]      uct            norm of tangential components of cell velocity
 * \param[in]      uft            norm of tangential components of face velocity
 * \param[in, out] res            value of the resulting exchange coefficients
 */
/*----------------------------------------------------------------------------*/

static void
_wall_function_1scale_log(const double    pena_bc_coeff,
                          const double    nu,
                          const double    k,
                          const double    hfc,
                          const double    uct,
                          const double    uft,
                          cs_real_t      *res)
{
  CS_UNUSED(k);

  const double ypluli = cs_get_glob_wall_functions()->ypluli;

  double ustar, ustarwer, ustarmin, ustaro;
  double y_over_nu = hfc/nu;

  double re = uct*y_over_nu;

  const double eps = 1e-4;
  const int n_max_iter = 100;

  if (re <= ypluli * ypluli) {

    /* In the viscous sub-layer: U+ = y+ */
    res[0] = pena_bc_coeff;
  }
  else {
    /* In the logaritim laye */

    /* The initial value is Wener or the minimun ustar to ensure convergence */

    ustarwer = pow(cs::abs(uct) / cs_turb_apow / pow(y_over_nu, cs_turb_bpow),
                   cs_turb_dpow);
    ustarmin = std::exp(-cs_turb_cstlog * cs_turb_xkappa) / y_over_nu;
    ustaro = cs::max(ustarwer, ustarmin);
    ustar =
      (cs_turb_xkappa * uct + ustaro) /
      (std::log(y_over_nu * ustaro) + cs_turb_xkappa * cs_turb_cstlog + 1.);

    int iter = 0;
    while (cs::abs(ustar - ustaro) >= eps * ustaro && iter < n_max_iter) {
      ustaro = ustar;
      ustar =
        (cs_turb_xkappa * uct + ustaro) /
        (std::log(y_over_nu * ustaro) + cs_turb_xkappa * cs_turb_cstlog + 1.);
      iter ++;
    }

    if (iter >= n_max_iter) {
      bft_printf(_("WARNING: non-convergence in the computation\n"
                   "******** of the friction velocity\n\n"
                   "friction vel: %f \n" ), ustar);
    }

    double h_t = ustar*ustar/uft;
    res[0] = cs::max(h_t, 0.0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function used to define the exchange coefficients for tangential and
 *        normal components. (2-velocity scales method)
 *
 * \param[in]      pena_bc_coeff  pena coefficient used as Dirichlet BC
 * \param[in]      nu             laminar kinematic viscosity
 * \param[in]      k              turbulent kinetic energy
 * \param[in]      hfc            distance from cell center to the wall
 * \param[in]      uct            norm of tangential components of cell velocity
 * \param[in]      uft            norm of tangential components of face velocity
 * \param[in, out] res            value of the resulting exchange coefficients
 */
/*----------------------------------------------------------------------------*/

static void
_wall_function_2scales_log(const double    pena_bc_coeff,
                           const double    nu,
                           const double    k,
                           const double    hfc,
                           const double    uct,
                           const double    uft,
                           cs_real_t      *res)
{
  const double ypluli = cs_get_glob_wall_functions()->ypluli;
  const double re     = std::sqrt(k) * hfc / nu;
  const double g      = std::exp(-re / 11.);
  const double uk =
    std::sqrt((1. - g) * std::sqrt(cs_turb_cmu) * k + g * nu * uct / hfc);
  const double yplus = hfc * uk / nu;

  if (yplus > ypluli) { /* In the logarithm zone */
    double ustar = uct / (std::log(yplus) / cs_turb_xkappa + cs_turb_cstlog);
    // FIXME: I set res[0] = 0 if uft = 0.0 but I am not sure
    if (cs::abs(uft) > 0.0) {
      double h_t = uk * ustar / uft;

      res[0] = cs::max(h_t, 0.0);
    }
    else {
      res[0] = 0.0;
    }
  }
  else /* In the linear zone */
    res[0] = pena_bc_coeff;
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
  cs_turbulence_param_t *tbp = nullptr;

  CS_MALLOC(tbp, 1, cs_turbulence_param_t);

  /* The following structures are shared with the Legacy part */

  tbp->model = cs_get_glob_turb_model();
  tbp->rans_param = cs_get_glob_turb_rans_model();
  tbp->les_param = cs_get_glob_turb_les_model();
  tbp->reference_values = cs_get_glob_turb_ref_values();

  if (cs_param_cdo_has_cdo_and_fv()) {
    tbp->shared_from_legacy = true;
  }
  else {
    tbp->shared_from_legacy = false;
  }

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
  cs_turbulence_t *tbs = nullptr;

  CS_MALLOC(tbs, 1, cs_turbulence_t);

  /* All the members of the following structures are shared with the Legacy
   * part. This structure is owned by cs_navsto_param_t
   */

  tbs->param = tbp;
  tbs->mom_eq = nullptr;

  /* Properties */

  tbs->rho    = nullptr; /* Mass density */
  tbs->mu_tot = nullptr; /* Total viscosity */
  tbs->mu_l   = nullptr; /* Laminar viscosity */
  tbs->mu_t   = nullptr; /* Turbulent viscosity */

  tbs->mu_tot_array = nullptr;

  /* Fields */

  tbs->mu_t_field = nullptr; /* Turbulent viscosity */
  tbs->rij        = nullptr; /* Reynolds stress tensor */

  /* Main structure (cast on-the-fly according to the turbulence model) */

  tbs->context = nullptr;

  /* Function pointers */

  tbs->init_context   = nullptr;
  tbs->free_context   = nullptr;
  tbs->compute        = nullptr;
  tbs->compute_steady = nullptr;
  tbs->update         = nullptr;

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

  CS_FREE(tbs->mu_tot_array);

  if (tbs->free_context != nullptr)
    tbs->context = tbs->free_context(tbs->context);

  assert(tbs->context == nullptr);
  CS_FREE(tbs);
  *p_turb_struct = nullptr;
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
  if (tbs == nullptr)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->type == CS_TURB_NONE)
    return; /* Nothing to do if there is a laminar flow */

  tbs->mom_eq = mom_eq;

  /* Set field metadata */

  const int  log_key         = cs_field_key_id("log");
  const int  post_key        = cs_field_key_id("post_vis");
  const bool has_previous    = false;
  const int  field_post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
  int        location_id     = cs_mesh_location_get_id_by_name("cells");

  if (tbp->shared_from_legacy) {
    int field_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;

    tbs->mu_t_field = cs_field_find_or_create(CS_NAVSTO_TURB_VISCOSITY,
                                              field_mask,
                                              location_id,
                                              1, /* dimension */
                                              has_previous);

    /* Properties (shared) */

    tbs->rho    = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
    tbs->mu_tot = cs_property_by_name(CS_NAVSTO_TOTAL_VISCOSITY);
    tbs->mu_l   = cs_property_by_name(CS_NAVSTO_LAM_VISCOSITY);
    tbs->mu_t = cs_property_by_name(CS_NAVSTO_TURB_VISCOSITY);

    assert(tbs->rho != nullptr && tbs->mu_l != nullptr &&
           tbs->mu_tot != nullptr && tbs->mu_t != nullptr);
  }
  else {
    int field_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

    tbs->mu_t_field = cs_field_find_or_create(CS_NAVSTO_TURB_VISCOSITY,
                                              field_mask,
                                              location_id,
                                              1, /* dimension */
                                              has_previous);

    /* Set default value for keys related to log and post-processing */

    tbs->mu_t_field->set_key_int(log_key, 1);
    tbs->mu_t_field->set_key_int(post_key, field_post_flag);

    /* Properties (shared) */

    tbs->rho    = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
    tbs->mu_tot = cs_property_by_name(CS_NAVSTO_TOTAL_VISCOSITY);
    tbs->mu_l   = cs_property_by_name(CS_NAVSTO_LAM_VISCOSITY);

    assert(tbs->rho != nullptr && tbs->mu_l != nullptr &&
           tbs->mu_tot != nullptr);

    /* Add a mu_t property and define it with the associated field */

    tbs->mu_t = cs_property_add(CS_NAVSTO_TURB_VISCOSITY, CS_PROPERTY_ISO);

    cs_property_def_by_field(tbs->mu_t, tbs->mu_t_field);
  }

  /* Set function pointers and initialize the context structure */

  if (!tbp->shared_from_legacy) {
    switch (model->model) {
      case CS_TURB_K_EPSILON:
      case CS_TURB_K_EPSILON_LIN_PROD:
        tbs->init_context = cs_turb_init_k_eps_context;
        tbs->free_context = cs_turb_free_k_eps_context;
        tbs->compute      = cs_turb_compute_k_eps;
        tbs->update       = cs_turb_update_k_eps;

        tbs->context = tbs->init_context(model);
        break;

      case CS_TURB_NONE:
        break;

      default:
        bft_error(__FILE__,
                  __LINE__,
                  0,
                  " %s: Invalid turbulence model with CDO schemes.\n"
                  " Please check your settings.",
                  __func__);
        break;
    }
  }
  else {
    tbs->init_context = nullptr;
    tbs->free_context = nullptr;
    tbs->compute      = nullptr;
    tbs->update       = cs_turb_update_shared_legacy;
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

  if (tbs == nullptr)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->type == CS_TURB_NONE)
    return; /* Nothing to do */

  /* Define the property related to the total viscosity */

  CS_MALLOC(tbs->mu_tot_array, quant->n_cells, cs_real_t);
  std::memset(tbs->mu_tot_array, 0, quant->n_cells * sizeof(cs_real_t));

  cs_property_def_by_array(tbs->mu_tot,
                           nullptr, /* all cells */
                           cs_flag_primal_cell,
                           tbs->mu_tot_array,
                           false, /* definition is owner ? */
                           true); /* full length */

  /* Turbulence model is computed by FV */
  if (tbp->shared_from_legacy) {
    return;
  }

  /* Last setup for each turbulence model */

  switch (model->model) {

  case CS_TURB_K_EPSILON:
  case CS_TURB_K_EPSILON_LIN_PROD:
    {
      /* Add a source term after having retrieved the equation param related to
         the turbulent kinetic energy equation */

      cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
      cs_equation_param_t  *tke_eqp = cs_equation_get_param(kec->tke);
      kec->tke_source_term
        = cs_equation_add_source_term_by_array(tke_eqp,
                                               nullptr, /* all cells */
                                               cs_flag_primal_cell,
                                               nullptr,
                                               false, /* is owner */
                                               true); /* full length */

      kec->tke_reaction
        = cs_property_def_by_array(cs_property_by_name("k_reaction"),
                                   nullptr, /* all cells */
                                   cs_flag_primal_cell,
                                   nullptr,
                                   false, /* is owner */
                                   true); /* full length */

      cs_equation_param_t  *eps_eqp = cs_equation_get_param(kec->eps);
      kec->eps_source_term
        = cs_equation_add_source_term_by_array(eps_eqp,
                                               nullptr, /* all cells */
                                               cs_flag_primal_cell,
                                               nullptr,
                                               false, /* is owner */
                                               true); /* full length */

      kec->eps_reaction
        = cs_property_def_by_array(cs_property_by_name("eps_reaction"),
                                   nullptr, /* all cells */
                                   cs_flag_primal_cell,
                                   nullptr,
                                   false, /* is owner ? */
                                   true); /* full length */

      cs_property_def_by_array(tbs->mu_tot,
                               nullptr, /* all cells */
                               cs_flag_primal_cell,
                               tbs->mu_tot_array,
                               false, /* is owner */
                               true); /* full length */

      /* Initialize TKE */

      cs_turb_ref_values_t *t_ref= cs_get_glob_turb_ref_values();
      cs_real_t             tke_ref = 1.5 * cs::pow2(0.02 * t_ref->uref);
      cs_equation_add_ic_by_value(tke_eqp, nullptr, &tke_ref);

      /* Initialize epsilon */

      cs_real_t eps_ref = powf(tke_ref, 1.5) * cs_turb_cmu / t_ref->almax;
      cs_equation_add_ic_by_value(eps_eqp, nullptr, &eps_ref);
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
 * \brief  Indicate whether use Legacy solved turbulent viscosity
 * \param[in, out] tbs        pointer to the turbulence main structure
 * \param[in]      is_shared  boolean parameter to set
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_set_shared_from_fv(cs_turbulence_t *tbs, bool is_shared)
{
  cs_turbulence_param_t *tbp = tbs->param;
  tbp->shared_from_legacy    = is_shared;
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

  if (tbs == nullptr)
    return;

  const cs_turbulence_param_t  *tbp = tbs->param;
  const cs_turb_model_t  *model = tbp->model;

  if (model->model == CS_TURB_NONE)
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
  if (tbm == nullptr)
    return nullptr;

  assert((tbm->model == CS_TURB_K_EPSILON) ||
         (tbm->model == CS_TURB_K_EPSILON_LIN_PROD));
  assert(tbm->type == CS_TURB_RANS);
  assert(tbm->order == CS_TURB_FIRST_ORDER);

  cs_turb_context_k_eps_t *kec = nullptr;

  CS_MALLOC(kec, 1, cs_turb_context_k_eps_t);

  /* Add new equations for the turbulent kinetic energy (tke) and the
     dissipation (epsilon) */

  kec->tke = cs_equation_add("k",                     /* equation name */
                             "k",                     /* variable name */
                             CS_EQUATION_TYPE_NAVSTO, /* related to NS */
                             1,
                             CS_BC_HMG_NEUMANN);

  kec->eps = cs_equation_add("epsilon", /* equation name */
                             "epsilon", /* variable name */
                             CS_EQUATION_TYPE_NAVSTO,
                             1,
                             CS_BC_HMG_NEUMANN);

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

  cs_property_t *k_reaction = cs_property_add("k_reaction",
                                              CS_PROPERTY_ISO);
  cs_property_t *eps_reaction = cs_property_add("epsilon_reaction",
                                                CS_PROPERTY_ISO);

  /* Retrieve the mass density */

  cs_property_t  *mass_density = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);

  /* Retrieve the advection field from Navier--Stokes (the mass flux) */

  cs_adv_field_t  *adv = cs_advection_field_by_name("mass_flux");
  if (adv == nullptr)
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
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_turb_free_k_eps_context(void     *tbc)
{
  cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbc;

  if (kec == nullptr)
    return kec;

  CS_FREE(kec);

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

  if (tbs == nullptr)
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
  cs_real_t *rho        = nullptr;

  /* Get mass density values in each cell */

  cs_property_iso_get_cell_values(time_step->t_cur, tbs->rho,
                                  &rho_stride, &rho);

  /* Get laminar viscosity values in each cell */

  int mu_stride = 0;
  cs_real_t *mu_l      = nullptr;
  cs_property_iso_get_cell_values(time_step->t_cur, tbs->mu_l,
                                  &mu_stride, &mu_l);


  /* Compute mu_t in each cell and define mu_tot = mu_t + mu_l */

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    mu_t[cell_id] = cs_turb_cmu * rho[cell_id * rho_stride] *
                    cs::pow2(k[cell_id]) / eps[cell_id];

    tbs->mu_tot_array[cell_id] = mu_t[cell_id] + mu_l[cell_id*mu_stride];

  }

  /* Free memory */

  CS_FREE(rho);
  CS_FREE(mu_l);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update for the current time step the new state for the turbulence
 *         model. directly update the turbulent viscosity from Legacy.
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  structure managing the time stepping
 * \param[in]      tbs        pointer to a \ref cs_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_update_shared_legacy(const cs_mesh_t           *mesh,
                             const cs_cdo_connect_t    *connect,
                             const cs_cdo_quantities_t *quant,
                             const cs_time_step_t      *time_step,
                             const cs_turbulence_t     *tbs)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  if (tbs == nullptr)
    return;

  const cs_lnum_t n_cells = mesh->n_cells;
  cs_real_t      *mu_t    = tbs->mu_t_field->val;

  /* Get laminar viscosity values in each cell */

  int        mu_stride = 0;
  cs_real_t *mu_l      = nullptr;
  cs_property_iso_get_cell_values(time_step->t_cur,
                                  tbs->mu_l,
                                  &mu_stride,
                                  &mu_l);

  /* Compute mu_t in each cell and define mu_tot = mu_t + mu_l */

#pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    tbs->mu_tot_array[cell_id] = mu_t[cell_id] + mu_l[cell_id * mu_stride];
  }

  /* Free memory */
  CS_FREE(mu_l);
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
  if (tbs == nullptr)
    return;

  /* Get the k-epsilon context */

  cs_turb_context_k_eps_t  *kec = (cs_turb_context_k_eps_t *)tbs->context;
  cs_equation_t *tke_eq = kec->tke;
  cs_equation_t *eps_eq = kec->eps;
  assert(kec != nullptr);

  /* Prepare source term and reaction term */

  cs_real_t *tke_source_term = nullptr, *eps_source_term = nullptr;
  cs_real_t *tke_reaction = nullptr, *eps_reaction = nullptr;

  CS_MALLOC(tke_source_term, mesh->n_cells, cs_real_t);
  CS_MALLOC(eps_source_term, mesh->n_cells, cs_real_t);
  CS_MALLOC(tke_reaction, mesh->n_cells, cs_real_t);
  CS_MALLOC(eps_reaction, mesh->n_cells, cs_real_t);

  /* Set the array values for each cs_xdef_t structures */

  cs_xdef_array_set_values(kec->tke_reaction,
                           false, /* is_owner */
                           tke_reaction);

  cs_xdef_array_set_values(kec->eps_reaction,
                           false, /* is_owner */
                           eps_reaction);

  cs_xdef_array_set_values(kec->tke_source_term,
                           false, /* is_owner */
                           tke_source_term);

  cs_xdef_array_set_values(kec->eps_source_term,
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

  CS_FREE(tke_source_term);
  CS_FREE(eps_source_term);
  CS_FREE(tke_reaction);
  CS_FREE(eps_reaction);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function used to define the exchange coefficients for tangential and
 *        normal components.
 *
 * \param[in]      eqp  pointer to a cs_equation_param_t
 * \param[in]      nu   laminar kinematic viscosity
 * \param[in]      k    turbulent kinetic energy
 * \param[in]      hfc  distance from cell center to the wall
 * \param[in]      uct  norm of tangential components of cell velocity
 * \param[in]      uft  norm of tangential components of face velocity
 * \param[in, out] res  value of the resulting exchange coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_wall_bc_coeffs(const cs_equation_param_t  *eqp,
                               const double                nu,
                               const double                k,
                               const double                hfc,
                               const double                uct,
                               const double                uft,
                               cs_real_t                  *res)
{
  const cs_wall_functions_t *glob_wf = cs_get_glob_wall_functions();
  const cs_wall_f_type_t iwallf = glob_wf->iwallf;

  switch (iwallf) {
    case CS_WALL_F_DISABLED:
      res[0] = eqp->strong_pena_bc_coeff;
      break;
    case CS_WALL_F_1SCALE_LOG:
      _wall_function_1scale_log(eqp->strong_pena_bc_coeff,
                                nu,
                                k,
                                hfc,
                                uct,
                                uft,
                                res);
      break;
    case CS_WALL_F_SCALABLE_2SCALES_LOG:
    case CS_WALL_F_2SCALES_LOG:
      _wall_function_2scales_log(eqp->strong_pena_bc_coeff,
                                 nu,
                                 k,
                                 hfc,
                                 uct,
                                 uft,
                                 res);
      break;
    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid wall function type.\n"
                " Expected wall function types in CDO turbulence: \n"
                " CS_WALL_F_DISABLED, CS_WALL_F_1SCALE_LOG"
                " or CS_WALL_F_SCALABLE_2SCALES_LOG.", __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return true if a wall function is used for turbulence.
 *
 * \param[in] turbulence  pointer to a cs_turbulence_param_t
 *
 *\return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_turb_wall_functions_is_activated(const cs_turbulence_param_t *turbulence)
{
  if (turbulence == nullptr)
    return false;

  if (turbulence->model->model != CS_TURB_NONE) {
    const cs_wall_functions_t *glob_wf = cs_get_glob_wall_functions();
    const cs_wall_f_type_t     iwallf  = glob_wf->iwallf;
    if (iwallf != CS_WALL_F_UNSET && iwallf != CS_WALL_F_DISABLED) {
      return true;
    }
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute and add the turbulent part of the boundray stress
 *
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in] boundary_stress pointer to a \ref cs_field_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
_cs_turbulence_t::add_boundary_stress(const cs_mesh_t *m,
                                      cs_field_t      *boundary_stress) const
{
  /* This a copy of cs_solve_navier_stokes.cpp */

  /* -2/3 rho * grad(k) for Eddy viscosity models with k defined
   * Note: we do not take the gradient of (rho k), as this would make
   *       the handling of BC's more complex...
   *
   * It is not clear whether the extrapolation in time is useful.
   *
   * This explicit term is computed once, at the first iteration on
   * cs_solve_navier_stokes: it is saved in a field if it must be extrapolated
   * in time ; it goes into trava if we do not extrapolate or iterate on
   * cs_solve_navier_stokes. */

  // index of the iteration on Navier-Stokes
  const int iterns = 1;

  if ((cs_glob_turb_model->order == CS_TURB_FIRST_ORDER &&
       cs_glob_turb_model->type == CS_TURB_RANS && CS_F_(k) != nullptr) &&
      cs_glob_turb_rans_model->igrhok == 1 && iterns == 1) {
    cs_dispatch_context ctx;

    cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

    cs_lnum_t n_b_faces   = m->n_b_faces;
    cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

    const cs_lnum_t *b_face_cells = m->b_face_cells;

    const cs_rreal_3_t *restrict diipb           = mq->diipb;
    const cs_nreal_3_t *restrict b_face_u_normal = mq->b_face_u_normal;

    cs_real_3_t *b_stress = (cs_real_3_t *)boundary_stress->val;

    cs_real_3_t *grad_k = nullptr;
    CS_MALLOC_HD(grad_k, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    const cs_field_t *fk      = CS_F_(k);
    cs_real_t        *cvara_k = fk->val;

    const cs_real_t *crom = CS_F_(rho)->val;

    cs_field_gradient_scalar(fk, true, 1, grad_k);

    constexpr cs_real_t d2s3 = 2.0 / 3.0;

    /* Calculation of wall stresses (part 3/5), if requested */
    const cs_real_t *coefa_k = fk->bc_coeffs->a;
    const cs_real_t *coefb_k = fk->bc_coeffs->b;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE(cs_lnum_t f_id) {
      const cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t       xkb =
        cvara_k[c_id] + cs_math_3_dot_product(diipb[f_id], grad_k[c_id]);

      xkb = coefa_k[f_id] + coefb_k[f_id] * xkb;
      xkb = d2s3 * crom[c_id] * xkb;
      for (cs_lnum_t i = 0; i < 3; i++)
        b_stress[f_id][i] -= xkb * b_face_u_normal[f_id][i];
    });
    ctx.wait();

    CS_FREE(grad_k);
  }
};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the convergence of the pseudo-steady algorithm
 *        when the unsteady Navier-Stokes system with a CDO face-based scheme
 *        is used.
 *
 * \param[in]  quant         pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]  ts            pointer to a \ref cs_time_step_t structure
 * \param[in]  psp           pointer to a \ref cs_param_psteady_t struct.
 *
 * \return returns true if the pseudo-steady algorithm has converged else
 * false
 *
 */
/*----------------------------------------------------------------------------*/

bool
_cs_turbulence_t::check_convergence(const cs_cdo_quantities_t *quant,
                                    const cs_time_step_t      *ts,
                                    const cs_param_psteady_t  &psp) const
{
  bool cvg = false;

  if (cs_glob_turb_model->order == CS_TURB_FIRST_ORDER &&
      cs_glob_turb_model->type == CS_TURB_RANS) {
    cs_dispatch_context ctx;

    const cs_field_t *f_k = CS_F_(k);

    assert(f_k != nullptr);

    const cs_real_t *k_curr = f_k->val;
    const cs_real_t *k_prev = f_k->val_pre;

    if (k_curr == nullptr || k_prev == nullptr) {
      return false;
    }

    const auto cell_vol = quant->cell_vol;
    cs_lnum_t  n_cells  = quant->n_cells;

    const cs_real_t dt = ts->dt[0];

    struct cs_double_n<2>      rd;
    struct cs_reduce_sum_nr<2> reducer;

    ctx.parallel_for_reduce(
      n_cells,
      rd,
      reducer,
      [=] CS_F_HOST_DEVICE(cs_lnum_t c_id, cs_double_n<2> & res) {
        const cs_real_t diff = k_curr[c_id] - k_prev[c_id];

        res.r[0] = cell_vol[c_id] * diff * diff;
        res.r[1] = cell_vol[c_id] * k_curr[c_id] * k_curr[c_id];
      });

    ctx.wait();

    cs::parall::sum(rd.r[0], rd.r[1]);

    // Simulate temporal derivative
    const cs_real_t norm2_k_diff = std::sqrt(rd.r[0]) / dt,
                    norm2_k      = std::sqrt(rd.r[1]);

    if (norm2_k_diff < psp.rtol * cs::max(1.0, norm2_k)) {
      cvg = true;
    }

    if (norm2_k_diff < psp.atol) {
      cvg = true;
    }

    cs_log_printf(CS_LOG_DEFAULT,
                  "- turbulence k: residual %5.3g (with "
                  "tolerence relative %5.3g and absolute %5.3g)\n",
                  norm2_k_diff,
                  psp.rtol * cs::max(1.0, norm2_k),
                  psp.atol);
  }

  return cvg;
};

/*----------------------------------------------------------------------------*/
