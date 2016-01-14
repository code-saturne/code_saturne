/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_post.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_param.h"
#include "cs_cdo_toolbox.h"
#include "cs_reco.h"
#include "cs_hodge.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_groundwater.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Set of parameters related to a tracer equation */
typedef struct {

  int              eq_id;

  /* Product of the bulk density by the distribution coefficient */
  double           rho_kd;

  /* Longitudinal and transversal dispersivity */
  double           alpha_l;
  double           alpha_t;

  /* Water molecular diffusivity */
  double           wmd;

  /* First order decay coefficient (related to the reaction term) */
  double           reaction_rate;


} cs_gw_tracer_t;

/* Set of parameters describing a type of soil */
typedef struct {

  int          ml_id;    /* id associated to a mesh location structure
                            The support entities are cells */

  /* Physical modelling adopted for this soil */
  cs_groundwater_model_t  model;

  /* Parameters for predefined models */
  cs_gw_genuchten_t       genuchten_param; /* Van-Genuchten-Mualen law */
  cs_gw_tracy_t           tracy_param;     /* Tracy law */

  /* Main soil properties */
  double                  residual_moisture;       /* theta_r */
  double                  saturated_moisture;      /* theta_s */
  cs_get_t                saturated_permeability;  /* k_s [m.s^-1] */

} cs_gw_soil_t;

/* Set of parameters related to the groundwater module */
struct _groundwater_t {

  cs_flag_t                flag;       /* Compact information */
  int                      post_freq;  /* Frequency for post-processing */

  cs_lnum_t                n_cells;    /* number of cells (useful for accessing
                                          soil_id) */

  cs_groundwater_model_t   global_model;

  /* Gravity effect */
  bool                     with_gravitation;
  cs_real_3_t              gravity;
  cs_real_t               *gravity_source_term;

  /* Set of equations associated to this module */
  int                      richards_eq_id;
  int                      n_tracers;
  cs_gw_tracer_t          *tracer_param;

  /* Moisture content variable and attached quantities */
  cs_field_t              *moisture_content;

  /* Physical parameters related to each kind of soil considered.
     If n_soils > 1, soil_id array stores the id giving access to the soil
     parameters related to each cell of the mesh */
  int                      n_soils;
  cs_gw_soil_t            *soil_param;
  short int               *soil_id;        /* NULL or allocated to n_cells */

  /* Permeability is the diffusion property related to Richards equation but
     this property plays also a role in the diffusion of tracer equations */
  cs_property_t           *permeability;  /* shared with domain */

  /* Scan the c2e connectivity index to get the darcian flux related to
     each dual face when CDO vertex-based scheme is activated */
  cs_real_t               *darcian_flux;
  cs_adv_field_t          *adv_field;    /* shared with domain */

  /* Work buffer */
  cs_real_t  *work;

};

/* List of available keys for setting the groundwater module */
typedef enum {

  GWKEY_GRAVITATION,
  GWKEY_OUTPUT_MOISTURE,
  GWKEY_POST_FREQ,
  GWKEY_ERROR

} gwkey_t;

typedef enum {

  SOILKEY_SATURATED_MOISTURE,
  SOILKEY_RESIDUAL_MOISTURE,
  SOILKEY_TRACY_HS,
  SOILKEY_TRACY_HR,
  SOILKEY_ERROR

} soilkey_t;


/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_gw[] =
  " Stop execution. The structure related to the groundwater module is empty.\n"
  " Please check your settings.\n";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding groundwater key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_gwkey(gwkey_t  key)
{
  switch (key) {

  case GWKEY_GRAVITATION:
    return "gravity";
  case GWKEY_OUTPUT_MOISTURE:
    return "output_moisture";
  case GWKEY_POST_FREQ:
    return "post_freq";

  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of an groundwater key.
 *         If not found, return a key error.
 *
 * \param[in] keyname    name of the key
 *
 * \return a gwkey_t
 */
/*----------------------------------------------------------------------------*/

static gwkey_t
_get_gwkey(const char  *keyname)
{
  gwkey_t  key = GWKEY_ERROR;

  if (strcmp(keyname, "gravity") == 0)
    key = GWKEY_GRAVITATION;
  else if (strcmp(keyname, "output_moisture") == 0)
    key = GWKEY_OUTPUT_MOISTURE;
  else if (strcmp(keyname, "post_freq") == 0)
    key = GWKEY_POST_FREQ;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding soil key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_soilkey(soilkey_t  key)
{
  switch (key) {

  case SOILKEY_SATURATED_MOISTURE:
    return "saturated_moisture";
  case SOILKEY_RESIDUAL_MOISTURE:
    return "residual_moisture";
  case SOILKEY_TRACY_HS:
    return "tracy_hs";
  case SOILKEY_TRACY_HR:
    return "tracy_hr";

  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of a soil key.
 *         If not found, return a key error.
 *
 * \param[in] keyname    name of the key
 *
 * \return a soilkey_t
 */
/*----------------------------------------------------------------------------*/

static soilkey_t
_get_soilkey(const char  *keyname)
{
  soilkey_t  key = SOILKEY_ERROR;

  if (strcmp(keyname, "saturated_moisture") == 0)
    key = SOILKEY_SATURATED_MOISTURE;
  else if (strcmp(keyname, "residual_moisture") == 0)
    key = SOILKEY_RESIDUAL_MOISTURE;
  else if (strcmp(keyname, "tracy_hs") == 0)
    key = SOILKEY_TRACY_HS;
  else if (strcmp(keyname, "tracy_hr") == 0)
    key = SOILKEY_TRACY_HR;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_gw_soil_t structure (already allocated)
 *
 * \param[in]      ml_name         name of the mesh location
 * \param[in]      model_kw        keyword related to the modelling
 * \param[in, out] soil            pointer to a cs_gw_soil_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_soil(const char     *ml_name,
           const char     *model_kw,
           cs_gw_soil_t   *soil)
{
  if (soil == NULL)
    return;

  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  if (cs_mesh_location_get_type(ml_id) != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of mesh location for mesh location  %s.\n"
                " The expected type is CS_MESH_LOCATION_CELLS.\n"), ml_name);

  soil->ml_id = ml_id;

  /* Set the model associated to this soil */
  if (strcmp(model_kw, "saturated") == 0) {

    soil->model = CS_GROUNDWATER_MODEL_SATURATED;
    soil->saturated_moisture = 1.0;
    soil->residual_moisture = 0.0;

  }
  else if (strcmp(model_kw, "genutchten") == 0) {

    soil->model = CS_GROUNDWATER_MODEL_GENUCHTEN;

    /* Default initialization */
    soil->saturated_moisture = 0.75;
    soil->residual_moisture = 0.15;

    double  n = 1.56;
    soil->genuchten_param.n = n;
    soil->genuchten_param.m = 1 - 1/n;
    soil->genuchten_param.scale = 0.036;
    soil->genuchten_param.tortuosity = 0.5;

  }
  else if (strcmp(model_kw, "tracy") == 0) {

    soil->model = CS_GROUNDWATER_MODEL_TRACY;

    /* Default initialization */
    soil->saturated_moisture = 0.75;
    soil->residual_moisture = 0.15;
    soil->tracy_param.h_r = -100;
    soil->tracy_param.h_s = 0;

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible model for a soil in the groundwater module.\n"
              " Value given: %s\n"
              " Availaible models: saturated, genutchen, tracy", model_kw);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_gw_tracer_t structure
 *
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 * \param[in, out] tp              pointer to a cs_gw_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_tracer_param(int                 tracer_eq_id,
                  double              wmd,
                  double              alpha_l,
                  double              alpha_t,
                  double              bulk_density,
                  double              distrib_coef,
                  double              reaction_rate,
                  cs_gw_tracer_t     *tp)
{
  assert(tp != NULL); /* Sanity check */

  tp->eq_id = tracer_eq_id;

  tp->wmd = wmd;
  tp->rho_kd = bulk_density * distrib_coef;
  tp->reaction_rate = reaction_rate;
  tp->alpha_l = alpha_l;
  tp->alpha_t = alpha_t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the diffusion term of the
 *         tracer equation
 *
 * \param[in]      theta          value of the moisture content
 * \param[in]      v              value of the local velocity
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_tracer_diffusion_tensor(double          theta,
                             const double    v[],
                             const void     *tracer_struc,
                             cs_get_t       *result)
{
  const cs_gw_tracer_t  *tp = (const cs_gw_tracer_t  *)tracer_struc;

  const double  vxx = v[0]*v[0], vyy = v[1]*v[1], vzz = v[2]*v[2];
  const double  vxy = v[0]*v[1], vxz = v[0]*v[1], vyz = v[1]*v[2];
  const double  vnorm = sqrt(vxx + vyy + vzz);

  assert(vnorm > cs_get_zero_threshold());
  const double  onv = 1/vnorm;
  const double  delta_coef = (tp->alpha_l - tp->alpha_t)*onv;

  /* Extra diagonal terms */
  result->tens[0][1] = result->tens[1][0] = delta_coef * vxy;
  result->tens[0][2] = result->tens[2][0] = delta_coef * vxz;
  result->tens[1][2] = result->tens[2][1] = delta_coef * vyz;

  /* Diagonal terms */
  const double  diag_coef = tp->wmd * theta;

  result->tens[0][0] = tp->alpha_l*vxx + tp->alpha_t*vyy + tp->alpha_t*vzz;
  result->tens[1][1] = tp->alpha_t*vxx + tp->alpha_l*vyy + tp->alpha_t*vzz;
  result->tens[2][2] = tp->alpha_t*vxx + tp->alpha_t*vyy + tp->alpha_l*vzz;

  for (int k = 0; k < 3; k++) {
    result->tens[k][k] *= onv;
    result->tens[k][k] += diag_coef;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *
 * \param[in]      theta          value of the moisture content
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_tracer_time_coeff(double        theta,
                       const void   *tracer_struc,
                       cs_get_t     *result)
{
  const cs_gw_tracer_t  *tracer = (const cs_gw_tracer_t  *)tracer_struc;

  result->val = theta + tracer->rho_kd;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in reaction term for the simulation
 *         of tracer equations
 *
 * \param[in]      theta          value of the moisture content
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_tracer_reaction_coeff(double        theta,
                           const void   *tracer_struc,
                           cs_get_t     *result)
{
  const cs_gw_tracer_t  *tracer = (const cs_gw_tracer_t  *)tracer_struc;

  result->val = tracer->reaction_rate * (theta + tracer->rho_kd);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         van Genuchten-Mualen law
 *
 * \param[in]      h             value of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_permeability_by_genuchten_law(double        h,
                               const void   *soil_struc,
                               cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_genuchten_t  law = soil->genuchten_param;

  /* Up to now, only isotropic values are considered */
  double  isoval = soil->saturated_permeability.val;

  if (h < 0) {

    /* S_e(h) = [1 + |alpha*h|^n]^(-m) */
    const double  alpha_h = fabs(law.scale*h);
    const double  one_alpha_hn = 1 + pow(alpha_h, law.n);
    const double  se_pow_overm = 1/one_alpha_hn;
    const double  se_L = pow(one_alpha_hn, -law.m*law.tortuosity);
    const double  coef_base = 1 - pow(1 - se_pow_overm, law.m);

    isoval *= se_L * coef_base*coef_base;

  }

  /* Build the related tensor (permeability is always defined as a tensor) */
  for (int k = 0; k < 3; k++) {
    result->tens[k][k] = isoval;
    for (int l = k+1; l < 3; l++) {
      result->tens[k][l] = result->tens[l][k] = 0;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         Tracy law
 *
 * \param[in]      h             value of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_permeability_by_tracy_law(double        h,
                           const void   *soil_struc,
                           cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_tracy_t  law = soil->tracy_param;

  /* Up to now, only isotropic values are considered */
  const double  ks = soil->saturated_permeability.val;
  const double  isoval = ks * (h - law.h_r)/(law.h_s - law.h_r);

  /* Build the related tensor (in bottom laywer permeability is always defined
     as a tensor) */
  for (int k = 0; k < 3; k++) {
    result->tens[k][k] = isoval;
    for (int l = k+1; l < 3; l++) {
      result->tens[k][l] = result->tens[l][k] = 0;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the moisture content using the Tracy law
 *
 * \param[in]      h             value of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_moisture_by_tracy_law(double        h,
                       const void   *soil_struc,
                       cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_tracy_t  law = soil->tracy_param;
  const double  delta_theta = soil->saturated_moisture - soil->residual_moisture;

  double  k_r = (h - law.h_r)/(law.h_s - law.h_r);
  double  moisture = k_r * delta_theta + soil->residual_moisture;

  result->val = moisture;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the darcian flux playing the role of advection field in
 *         groundwater flows
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards    pointer to the Richards equation structure
 * \param[in, out] gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcian_flux(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *cdoq,
                     const cs_equation_t         *richards,
                     cs_groundwater_t            *gw)
{
  cs_lnum_t  i, _i, c_id;

  /* Sanity checks */
  if (richards == NULL || gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Groundwater module or Richards eq. is not allocated.");

  const cs_field_t  *h = cs_equation_get_field(richards);
  const cs_equation_param_t  *eqp = cs_equation_get_param(richards);
  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  bool  diff_tensor_uniform = cs_property_is_uniform(eqp->diffusion_property);
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_real_t  *loc_grdh = gw->work;
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect,
                                                  eqp->diffusion_hodge);

  /* Define the flux by cellwise contributions
     loc_flux = - loc_hodge * loc_gradient(h) */
  for (c_id = 0; c_id < cdoq->n_cells; c_id++) {

    if (c_id == 0 || diff_tensor_uniform == false) {
      cs_property_get_cell_tensor(c_id,
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);
      cs_hodge_builder_set_tensor(hb, (const cs_real_t (*)[3])diff_tensor);
    }

    /* Build a local discrete Hodge op. and return a local dense matrix */
    const cs_locmat_t  *hloc = cs_hodge_build_local(c_id, connect, cdoq, hb);

    for (i = c2e->idx[c_id], _i = 0; i < c2e->idx[c_id+1]; i++, _i++) {

      const cs_lnum_t  e_shft = 2*c2e->ids[i];
      const cs_lnum_t  v1_id = e2v->col_id[e_shft];
      const cs_lnum_t  v2_id = e2v->col_id[e_shft+1];
      const short int  sgn_v1 = e2v->sgn[e_shft];
      const short int  sgn_v2 = e2v->sgn[e_shft+1];

      loc_grdh[_i] = -(sgn_v1*h->val[v1_id] + sgn_v2*h->val[v2_id]);

    } // Loop on cell edges

    cs_locmat_matvec(hloc, loc_grdh,
                     gw->darcian_flux + c2e->idx[c_id]); // local flux

  } // Loop on cells

  /* Free builder */
  hb = cs_hodge_builder_free(hb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the moisture content from the value of the hydraulic head
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards    pointer to the Richards equation structure
 * \param[in, out] gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_moisture_content(const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *cdoq,
                         const cs_equation_t         *richards,
                         cs_groundwater_t            *gw)
{
  cs_lnum_t  i, c_id;
  double  val_xc;
  cs_get_t  get;

  /* Sanity checks */
  if (richards == NULL || gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Groundwater module or Richards eq. is not allocated.");

  const cs_field_t  *h = cs_equation_get_field(richards);
  cs_field_t  *moisture = gw->moisture_content;

  /* Sanity checks */
  if (moisture == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " The field related to the moisture content is not allocated.");

  /* Copy current field values to previous values */
  cs_field_current_to_previous(moisture);

  /* Moisture content is define in each cell while the hydraulic head is
     defined at vertices (CDO vertex-based is used) */
  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gw_soil_t  *soil = gw->soil_param + soil_id;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

    switch (soil->model) {

    case CS_GROUNDWATER_MODEL_TRACY:

      if (elt_ids == NULL) {

        /* Sanity checks */
        assert((cdoq->n_cells == n_elts[0]) && (gw->n_cells == cdoq->n_cells));

        for (c_id = 0; c_id <  gw->n_cells; c_id++) {

          /* Reconstruct (or interpolate) value at the current cell center */
          cs_reco_pv_at_cell_center(c_id,
                                    connect->c2v, cdoq, h->val,
                                    &val_xc);

          _moisture_by_tracy_law(val_xc, (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on cells

      }
      else {

        for (i = 0; i <  n_elts[0]; i++) {

          /* Reconstruct (or interpolate) value at the current cell center */
          c_id = elt_ids[i];
          cs_reco_pv_at_cell_center(c_id,
                                    connect->c2v, cdoq, h->val,
                                    &val_xc);

          _moisture_by_tracy_law(val_xc, (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on selected cells

      } // elt_ids != NULL
      break; // Tracy model

    case CS_GROUNDWATER_MODEL_SATURATED:

      if (elt_ids == NULL) {

        /* Sanity checks */
        assert((cdoq->n_cells == n_elts[0]) && (gw->n_cells == cdoq->n_cells));

        for (c_id = 0; c_id <  gw->n_cells; c_id++)
          moisture->val[c_id] = soil->saturated_moisture;

      }
      else
        for (i = 0; i <  n_elts[0]; i++)
          moisture->val[elt_ids[i]] = soil->saturated_moisture;

      break; // Saturated model

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of model for computing the moisture content."));
      break; // Nothing to do

    } /* Switch according to the kind of modelling */

  } /* Loop on soils */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \param[in]  n_cells    number of cells in the computational domain
 *
 * \return a pointer to a new allocated cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_groundwater_create(cs_lnum_t    n_cells)
{
  cs_groundwater_t  *gw = NULL;

  BFT_MALLOC(gw, 1, cs_groundwater_t);

  /* Default initialization */
  gw->flag = 0;
  gw->post_freq = -1;

  gw->n_soils = 0;
  gw->soil_param = NULL;
  gw->n_cells = n_cells;
  gw->soil_id = NULL;
  gw->global_model = CS_GROUNDWATER_N_MODELS;

  gw->with_gravitation = false;
  gw->gravity[0] = 0, gw->gravity[1] = 0, gw->gravity[2] = 0;
  gw->gravity_source_term = NULL;

  gw->richards_eq_id = -1;
  gw->n_tracers = 0;
  gw->tracer_param = NULL;

  gw->darcian_flux = NULL;
  gw->adv_field = NULL;
  gw->work = NULL;

  return gw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to groundwater flows
 *
 * \param[in, out]  gw     pointer to a cs_groundwater_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_groundwater_finalize(cs_groundwater_t   *gw)
{
  if (gw == NULL)
    return NULL;

  BFT_FREE(gw->soil_param);
  BFT_FREE(gw->tracer_param);
  BFT_FREE(gw->darcian_flux);
  BFT_FREE(gw->work);

  if (gw->with_gravitation)
    BFT_FREE(gw->gravity_source_term);

  if (gw->n_soils > 1)
    BFT_FREE(gw->soil_id);

  BFT_FREE(gw);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       keyname   name of key related to the member of adv to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_param(cs_groundwater_t    *gw,
                         const char          *keyname,
                         const char          *keyval)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gwkey_t  key = _get_gwkey(keyname);

  if (key == GWKEY_ERROR) {

    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < GWKEY_ERROR; i++) {
      bft_printf("%s ", _print_gwkey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the groundwater module.\n"
                " Please read listing for more details and modify your"
                " settings."));

  } /* Error message */

  switch(key) {

  case GWKEY_GRAVITATION:
    gw->with_gravitation = true;
    if (strcmp(keyval, "x") == 0)
      gw->gravity[0] = 1., gw->gravity[1] = gw->gravity[2] = 0.;
    else if (strcmp(keyval, "-x") == 0)
      gw->gravity[0] = -1., gw->gravity[1] = gw->gravity[2] = 0.;
    else if (strcmp(keyval, "y") == 0)
      gw->gravity[1] = 1., gw->gravity[0] = gw->gravity[2] = 0.;
    else if (strcmp(keyval, "-y") == 0)
      gw->gravity[1] = -1., gw->gravity[0] = gw->gravity[2] = 0.;
    else if (strcmp(keyval, "z") == 0)
      gw->gravity[2] = 1., gw->gravity[0] = gw->gravity[1] = 0.;
    else if (strcmp(keyval, "-z") == 0)
      gw->gravity[2] = -1., gw->gravity[0] = gw->gravity[1] = 0.;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid choice of gravitation axis: %s.\n"
                  " Available choices are 'x', 'y' and 'z'\n"
                  " Please check your settings."), keyval);
    break;

  case GWKEY_OUTPUT_MOISTURE:
    if (strcmp(keyval, "false")) // not "false"
      gw->flag |= CS_GROUNDWATER_POST_MOISTURE;
    break;

  case GWKEY_POST_FREQ:
    gw->post_freq = atoi(keyval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_groundwater_t structure
 *
 * \param[in]  gw     pointer to a cs_groundwater_t struct. to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_summary(const cs_groundwater_t   *gw)
{
  if (gw == NULL)
    return;

  bft_printf("\n");
  bft_printf("%s", lsepline);
  bft_printf("\tSummary of the groundwater module\n");
  bft_printf("%s", lsepline);

  if (gw->with_gravitation)
    bft_printf("  <GW/Gravitation> true -- Axis = [%.2f %.2f %.2f]\n",
               gw->gravity[0], gw->gravity[1], gw->gravity[2]);
  else
    bft_printf("  <GW/Gravitation> false\n");

  bft_printf("  <GW/Tracer> n_tracer_equations %d\n", gw->n_tracers);
  bft_printf("  <GW/Soils>  n_soils %d\n", gw->n_soils);

  const cs_property_t  *permeability = gw->permeability;

  for (int i = 0; i < gw->n_soils; i++) {

    const cs_gw_soil_t  soil = gw->soil_param[i];
    const char *ml_name = cs_mesh_location_get_name(soil.ml_id);
    const cs_get_t  sat_perm = soil.saturated_permeability;

    bft_printf("  <GW/Soil %s>", ml_name);
    bft_printf(" residual_moisture %5.3e", soil.residual_moisture);
    bft_printf(" saturated_moisture %5.3e\n", soil.saturated_moisture);
    bft_printf("  <GW/Soil %s>", ml_name);

    switch (cs_property_get_type(permeability)) {

    case CS_PROPERTY_ISO:
      bft_printf(" saturated_permeability (iso) %5.3e\n", sat_perm.val);
      break;

    case CS_PROPERTY_ORTHO:
      bft_printf(" saturated_permeability (ortho) %5.3e %5.3e %5.3e\n",
                 sat_perm.vect[0], sat_perm.vect[1], sat_perm.vect[2]);
      break;

    case CS_PROPERTY_ANISO:
      bft_printf(" saturated_permeability (aniso) %-5.3e %5.3e %5.3e\n"
                 "                                %-5.3e %5.3e %5.3e\n"
                 "                                %-5.3e %5.3e %5.3e\n",
                 sat_perm.tens[0][0], sat_perm.tens[0][1], sat_perm.tens[0][2],
                 sat_perm.tens[1][0], sat_perm.tens[1][1], sat_perm.tens[1][2],
                 sat_perm.tens[2][0], sat_perm.tens[2][1], sat_perm.tens[2][2]);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of property for %s."),
                cs_property_get_name(permeability));
      break;

    } // Switch on property type

    bft_printf("  <GW/Soil %s>", ml_name);
    switch (soil.model) {
    case CS_GROUNDWATER_MODEL_GENUCHTEN:
      bft_printf(" model VanGenuchten-Mualen\n");
      break;
    case CS_GROUNDWATER_MODEL_SATURATED:
      bft_printf(" model saturated\n");
      break;
    case CS_GROUNDWATER_MODEL_TRACY:
      bft_printf(" model Tracy\n");
      break;
    case CS_GROUNDWATER_MODEL_USER:
      bft_printf(" model User-defined\n");
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid model for a soil in the groundwater module.\n"
                " Please check your settings.");
    } // Switch model

    if (gw->n_soils > 1) {

      switch (gw->global_model) {
      case CS_GROUNDWATER_MODEL_COMPOSITE:
        bft_printf("  <GW/Global model> composite model\n");
        break;
      case CS_GROUNDWATER_MODEL_GENUCHTEN:
        bft_printf("  <GW/Global model> model VanGenuchten-Mualen\n");
        break;
      case CS_GROUNDWATER_MODEL_SATURATED:
        bft_printf("  <GW/Global model> model saturated\n");
        break;
      case CS_GROUNDWATER_MODEL_TRACY:
        bft_printf("  <GW/Global model> model Tracy\n");
        break;
      case CS_GROUNDWATER_MODEL_USER:
        bft_printf("  <GW/Global model> model User-defined\n");
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid model for groundwater module.\n"
                  " Please check your settings.");
    } // Switch model

    }

  } // Loop on soils

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]      connect          pointer to a cs_cdo_connect_t structure
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_groundwater_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_initialize(const cs_cdo_connect_t  *connect,
                          int                      richards_eq_id,
                          cs_property_t           *permeability,
                          cs_property_t           *soil_capacity,
                          cs_adv_field_t          *adv_field,
                          cs_groundwater_t        *gw)
{
  cs_equation_t  *eq = NULL;

  const cs_connect_index_t  *c2e = connect->c2e;

  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gw->richards_eq_id = richards_eq_id;

  /* Create a new equation structure for Richards' equation */
  eq = cs_equation_create("Richards",                   // equation name
                          "hydraulic_head",             // variable name
                          CS_EQUATION_TYPE_GROUNDWATER, // type of equation
                          CS_PARAM_VAR_SCAL,            // type of variable
                          CS_PARAM_BC_HMG_NEUMANN);     // default BC

  /* Associate soil_capacity to the unsteady term of the Richards eq. */
  if (soil_capacity != NULL)
    cs_equation_link(eq, "time", soil_capacity);

  /* Associate permeability to the diffusion property of the Richards eq. */
  gw->permeability = permeability;
  cs_equation_link(eq, "diffusion", permeability);

  /* Advection field induced by the hydraulic head */
  gw->adv_field = adv_field;

  /* Up to now Richards equation is only set with CS_SPACE_SCHEME_CDOVB */
  BFT_MALLOC(gw->darcian_flux, c2e->idx[connect->c_info->n_ent], cs_real_t);
  for (cs_lnum_t i = 0; i < c2e->idx[connect->c_info->n_ent]; i++)
    gw->darcian_flux[i] = 0;

  /* Work (temporary) buffer */
  BFT_MALLOC(gw->work, connect->n_max_ebyc, cs_real_t);

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new type of soil to consider in the groundwater module
 *
 * \param[in, out] gw         pointer to a cs_groundwater_t structure
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      model_kw   keyword related to the model used
 * \param[in]      ks         value(s) of the saturated permeability
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_add_soil_by_value(cs_groundwater_t   *gw,
                                 const char         *ml_name,
                                 const char         *model_kw,
                                 const char         *pty_val)
{
  if (gw == NULL)
    return;

  int  soil_id = gw->n_soils;

  gw->n_soils += 1;
  BFT_REALLOC(gw->soil_param, gw->n_soils, cs_gw_soil_t);

  cs_gw_soil_t  *soil = gw->soil_param + soil_id;

  _init_soil(ml_name, model_kw, soil);

  /* Set the saturated permeability */
  switch (cs_property_get_type(gw->permeability)) {

  case CS_PROPERTY_ISO:
    cs_param_set_get(CS_PARAM_VAR_SCAL, (const void *)pty_val,
                     &(soil->saturated_permeability));
    break;

  case CS_PROPERTY_ORTHO:
    cs_param_set_get(CS_PARAM_VAR_VECT, (const void *)pty_val,
                     &(soil->saturated_permeability));
    break;

  case CS_PROPERTY_ANISO:
    cs_param_set_get(CS_PARAM_VAR_TENS, (const void *)pty_val,
                     &(soil->saturated_permeability));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of property for %s."),
              cs_property_get_name(gw->permeability));
    break;

  } /* End of switch */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       ml_name   name of the mesh location associated to this soil
 * \param[in]       keyname   name of key related to the member of adv to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_soil_param(cs_groundwater_t    *gw,
                              const char          *ml_name,
                              const char          *keyname,
                              const char          *keyval)
{
  int  i;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int ml_id = -1;

  if (ml_name != NULL) {

    ml_id = cs_mesh_location_get_id_by_name(ml_name);
    if (ml_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location name %s.\n"
                  " This mesh location is not already defined.\n"), ml_name);

  }

  soilkey_t  key = _get_soilkey(keyname);

  if (key == SOILKEY_ERROR) {

    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (i = 0; i < SOILKEY_ERROR; i++) {
      bft_printf("%s ", _print_soilkey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the soil parameters in groundwater"
                " module.\n"
                " Please read the listing for more details and modify your"
                " settings."));

  } /* Error message */

  /* Set different available keys */
  switch(key) {

  case SOILKEY_SATURATED_MOISTURE:
    {
      double  theta_s = atof(keyval);

      if (ml_id == -1)
        for (i = 0; i < gw->n_soils; i++)
          gw->soil_param[i].saturated_moisture = theta_s;
      else
        for (i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].saturated_moisture = theta_s;

    }
    break;

  case SOILKEY_RESIDUAL_MOISTURE:
    {
      double  theta_r = atof(keyval);

      if (ml_id == -1)
        for (i = 0; i < gw->n_soils; i++)
          gw->soil_param[i].residual_moisture = theta_r;
      else
        for (i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].residual_moisture = theta_r;

    }
    break;

  case SOILKEY_TRACY_HS:
    {
      double  h_s = atof(keyval);

      if (ml_id == -1) {

        for (i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GROUNDWATER_MODEL_TRACY)
            gw->soil_param[i].tracy_param.h_s = h_s;

      }
      else {

        for (i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].tracy_param.h_s = h_s;

      }
    }
    break;

  case SOILKEY_TRACY_HR:
    {
      double  h_r = atof(keyval);

      if (ml_id == -1) {

        for (i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GROUNDWATER_MODEL_TRACY)
            gw->soil_param[i].tracy_param.h_r = h_r;

      }
      else {

        for (i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].tracy_param.h_r = h_r;

      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_groundwater_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      eqname          name of the equation
 * \param[in]      varname         name of the related variable
 * \param[in]      diff_pty        related property for the diffusion term
 * \param[in]      time_pty        related property for the time-dependent term
 * \param[in]      reac_pty        related property for the reaction term
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_add_tracer(cs_groundwater_t    *gw,
                          int                  tracer_eq_id,
                          const char          *eqname,
                          const char          *varname,
                          cs_property_t       *diff_pty,
                          cs_property_t       *time_pty,
                          cs_property_t       *reac_pty,
                          double               wmd,
                          double               alpha_l,
                          double               alpha_t,
                          double               bulk_density,
                          double               distrib_coef,
                          double               reaction_rate)
{
  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_equation_t  *eq = NULL;

  BFT_REALLOC(gw->tracer_param, gw->n_tracers + 1, cs_gw_tracer_t);

  cs_gw_tracer_t  *tp = gw->tracer_param + gw->n_tracers;

  _set_tracer_param(tracer_eq_id,
                    wmd,
                    alpha_l,
                    alpha_t,
                    bulk_density,
                    distrib_coef,
                    reaction_rate,
                    tp);

  gw->n_tracers += 1;

  eq = cs_equation_create(eqname,                       // equation name
                          varname,                      // variable name
                          CS_EQUATION_TYPE_GROUNDWATER, // type of equation
                          CS_PARAM_VAR_SCAL,            // type of variable
                          CS_PARAM_BC_HMG_NEUMANN);     // default BC

  /* Associate the property for the unsteady term */
  cs_equation_link(eq, "time", time_pty);
  cs_property_def_by_law(time_pty, _get_tracer_time_coeff);
  cs_property_set_struct(time_pty, (const void *)tp);

  /* Associate the advection field for the advection term */
  assert(gw->adv_field != NULL); /* Sanity check */
  cs_equation_link(eq, "advection", gw->adv_field);

  /* Associate the property for the diffusion term */
  if (diff_pty != NULL) {

    cs_equation_link(eq, "diffusion", diff_pty);
    cs_property_def_by_scavec_law(diff_pty, _get_tracer_diffusion_tensor);
    cs_property_set_struct(diff_pty, (const void *)tp);

  }

  /* Associate the property for the reaction term */
  if (reac_pty != NULL) {

    cs_equation_add_reaction(eq, "decay", "linear", reac_pty);
    cs_property_def_by_law(reac_pty, _get_tracer_reaction_coeff);
    cs_property_set_struct(reac_pty, (const void *)tp);

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the module dedicated to groundwater flows
 *
 * \param[in, out] equations    pointer to the array of cs_equation_t struct.
 * \param[in, out] gw           pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_automatic_settings(cs_equation_t      **equations,
                                  cs_groundwater_t    *gw)
{
  cs_lnum_t  i, j, c_id;
  cs_flag_t  flag;

  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->n_soils == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is activated but no soil is defined."));

  cs_property_t  *permeability = gw->permeability;
  cs_equation_t  *richards = equations[gw->richards_eq_id];
  cs_field_t  *hydraulic_head = cs_equation_get_field(richards);

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_space_scheme(richards) == CS_SPACE_SCHEME_CDOVB);

  /* Moisture content (defined in each cell) */
  bool has_previous = cs_equation_is_steady(richards) ? false:true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  int  location_id = cs_mesh_location_get_id_by_name(N_("cells"));

  gw->moisture_content = cs_field_create("moisture_content",
                                         field_mask,
                                         location_id,
                                         1,        // dimension
                                         true,     // interleave
                                         has_previous);
  cs_field_allocate_values(gw->moisture_content);

  if (gw->with_gravitation) { /* Gravitation effect */

    int  ml_id = hydraulic_head->location_id;
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
    const cs_lnum_t  n_vertices = n_elts[0];

    BFT_MALLOC(gw->gravity_source_term, n_vertices, cs_real_t);
    for (i = 0; i < n_vertices; i++)
      gw->gravity_source_term[i] = 0;

    flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;
    cs_equation_add_gravity_source_term(richards,
                                        ml_id,
                                        flag,
                                        gw->gravity_source_term);

  }

  /* Set the values for the permeability and the moisture content
     and if needed set also the value of the soil capacity */
  if (gw->n_soils == 1) {

    cs_gw_soil_t  *soil = gw->soil_param;

    gw->global_model = soil->model;

    switch (soil->model) {

    case CS_GROUNDWATER_MODEL_GENUCHTEN:
      cs_property_def_by_law(permeability, _permeability_by_genuchten_law);
      flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;
      cs_property_set_array(permeability, flag, hydraulic_head->val);
      cs_property_set_struct(permeability, (const void *)soil);

      /* Soil capacity settings (related to unsteady term) */
      if (has_previous)
        bft_error(__FILE__, __LINE__, 0, "Not implemented. To do.");

      break;

    case CS_GROUNDWATER_MODEL_TRACY:
      cs_property_def_by_law(permeability, _permeability_by_tracy_law);
      flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;
      cs_property_set_array(permeability, flag, hydraulic_head->val);
      cs_property_set_struct(permeability, (const void *)soil);

      /* Soil capacity settings (related to unsteady term) */
      if (has_previous) {

        cs_property_t  *capacity = cs_equation_get_time_property(richards);

        const cs_gw_tracy_t  law = soil->tracy_param;
        const double  delta_h = law.h_s - law.h_r;
        const double  delta_theta =
          soil->saturated_moisture - soil->residual_moisture;

        cs_property_iso_def_by_value(capacity, delta_theta/delta_h);

      } /* Set the soil capacity */
      break;

    case CS_GROUNDWATER_MODEL_SATURATED:
      switch (cs_property_get_type(permeability)) {

      case CS_PROPERTY_ISO:
        cs_property_iso_def_by_value(permeability,
                                     soil->saturated_permeability.val);
        break;

      case CS_PROPERTY_ORTHO:
        cs_property_ortho_def_by_value(permeability,
                                       soil->saturated_permeability.vect);
        break;

      case CS_PROPERTY_ANISO:
        cs_property_aniso_def_by_value(permeability,
                      (const double (*)[3])soil->saturated_permeability.tens);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
        break;

      } // End of switch on permeability type

      if (has_previous)
        bft_error(__FILE__, __LINE__, 0,
                  " Saturated model should lead to steady Richards equation.");

      break; // Switch on saturated model

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Incompatible model for groundwater flows.\n"
                " Availaible models: saturated, genutchen, tracy");
      // Nothing to do
      break;

    } // switch on modelling

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

    assert(elt_ids == NULL && gw->n_cells == n_elts[0]);
#endif

    /* Default initialization of the moisture content */
    const double  theta_s = soil->saturated_moisture;

    for (i = 0; i < gw->n_cells; i++)
      gw->moisture_content->val[i] = theta_s;

  }
  else { /* n_soils > 1 */

    cs_groundwater_model_t  model = gw->soil_param[0].model;

    BFT_MALLOC(gw->soil_id, gw->n_cells, short int);
    for (i = 0; i < gw->n_cells; i++)
      gw->soil_id[i] = -1; // default initialization (not set)

    gw->global_model = model;

    /* Is there a unique model ? */
    for (i = 0; i < gw->n_soils; i++) {

      const cs_gw_soil_t  *soil = gw->soil_param + i;
      const double  theta_s = soil->saturated_moisture;

      if (soil->model != model)
        gw->global_model = CS_GROUNDWATER_MODEL_COMPOSITE;

      const char  *ml_name = cs_mesh_location_get_name(soil->ml_id);
      const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);
      const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

      /* Sanity check */
      assert(elt_ids != NULL);

      /* Initialization of soil_id and moisture content */
      for (j = 0; j < n_elts[0]; j++) {

        c_id = elt_ids[j];
        gw->soil_id[c_id] = i;
        gw->moisture_content->val[c_id] = theta_s;

      } // Loop on selected cells

      switch (soil->model) {

      case CS_GROUNDWATER_MODEL_TRACY:

        cs_property_def_subdomain_by_law(permeability,
                                         ml_name,
                                         (const void *)soil,
                                         _permeability_by_tracy_law);
        flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;
        cs_property_set_array(permeability, flag, hydraulic_head->val);

        /* Soil capacity settings (related to unsteady term) */
        if (has_previous) {

          cs_property_t  *capacity = cs_equation_get_time_property(richards);

          const cs_gw_tracy_t  law = soil->tracy_param;
          const double  delta_h = law.h_s - law.h_r;
          const double  delta_theta =
            soil->saturated_moisture - soil->residual_moisture;

          cs_property_iso_def_subdomain_by_value(capacity,
                                                 ml_name,
                                                 delta_theta/delta_h);

        } /* Set the soil capacity */

        break; // Tracy model

      case CS_GROUNDWATER_MODEL_SATURATED:

        switch (cs_property_get_type(permeability)) {

        case CS_PROPERTY_ISO:
          cs_property_iso_def_subdomain_by_value(permeability, ml_name,
                                       soil->saturated_permeability.val);
          break;

        case CS_PROPERTY_ORTHO:
          cs_property_ortho_def_subdomain_by_value(permeability, ml_name,
                                        soil->saturated_permeability.vect);
          break;

        case CS_PROPERTY_ANISO:
          cs_property_aniso_def_subdomain_by_value(permeability, ml_name,
                   (const double (*)[3])soil->saturated_permeability.tens);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
          break;

        } // End of switch on permeability type

        if (has_previous)
          bft_error(__FILE__, __LINE__, 0,
                    " Saturated model should lead to steady Richards equation.");

        break; // Saturated model

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid type of model for groundwater module.."));
        break; // Nothing to do

      } /* Switch depending on the type of model */

    } // Loop on the different type of soils

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_lnum_t  n_unset_cells = 0;
    for (i = 0; i < gw->n_cells; i++)
      if (gw->soil_id[i] == -1) n_unset_cells++;

    if (n_unset_cells > 0)
      bft_error(__FILE__, __LINE__, 0,
                _(" %d cells have no related soil.\n"
                  " Please check your settings."), n_unset_cells);
#endif
  } /* n_soils > 1 */

  /* Define and then link the advection field to each tracer equations */
  flag = CS_PARAM_FLAG_FACE | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_BY_CELL;
  flag |= CS_PARAM_FLAG_SCAL;

  cs_advection_field_def_by_array(gw->adv_field, flag, gw->darcian_flux);

  for (i = 0; i < gw->n_tracers; i++) {

    cs_gw_tracer_t  tp = gw->tracer_param[i];
    cs_equation_t  *eq = equations[tp.eq_id];
    cs_flag_t  eq_flag = cs_equation_get_flag(eq);

    const char  *eqname = cs_equation_get_name(eq);

    /* Set time property */
    cs_property_t  *time_pty = cs_equation_get_time_property(eq);

    flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL;
    cs_property_set_array(time_pty, flag, gw->moisture_content->val);

    /* Set diffusion property */
    if (eq_flag & CS_EQUATION_DIFFUSION) {

      cs_property_t  *diff_pty = cs_equation_get_diffusion_property(eq);

      flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL;
      cs_property_set_array(diff_pty, flag, gw->moisture_content->val);

      flag = CS_PARAM_FLAG_FACE | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_BY_CELL;
      flag |= CS_PARAM_FLAG_SCAL;
      cs_property_set_second_array(diff_pty, flag, gw->darcian_flux);

    }

    /* Add reaction term if needed */
    if (eq_flag & CS_EQUATION_REACTION) {

      cs_property_t  *reac_pty = cs_equation_get_reaction_property(eq,
                                                                   "decay");

      flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL;
      cs_property_set_array(reac_pty, flag, gw->moisture_content->val);

    }

    /* TODO: add predefined source term */

  } // Loop on tracer equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] eqs        array of pointers to cs_equation_t structures
 * \param[in, out] gw         pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_compute(const cs_mesh_t              *mesh,
                       const cs_time_step_t         *time_step,
                       double                        dt_cur,
                       const cs_cdo_connect_t       *connect,
                       const cs_cdo_quantities_t    *cdoq,
                       cs_equation_t                *eqs[],
                       cs_groundwater_t             *gw)
{
  int  i;

  if (gw == NULL)
    return;

  const int  nt_cur = time_step->nt_cur;

  /* Solve the Richards equation */
  cs_equation_t  *richards = eqs[gw->richards_eq_id];

  /* Sanity check */
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  if (nt_cur == 0) {

    /* Initialize system before resolution for all equations
       - create system builder
       - initialize field according to initial conditions
       - initialize source term */
    cs_equation_init_system(mesh, connect, cdoq, time_step, richards);

    /* Build and solve the linear system related to the Richards equations */
    if (cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(time_step, richards);

      /* Compute the darcian flux */
      _update_darcian_flux(connect, cdoq, richards, gw);

      /* Update the moisture content */
      _update_moisture_content(connect, cdoq, richards, gw);

    }

    for (i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_param[i].eq_id];

      cs_equation_init_system(mesh, connect, cdoq, time_step, eq);

      if (cs_equation_is_steady(eq)) {

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(time_step, eq);

      } /* Solve this equation which is steady */

    } /* Loop on tracer equations */

  }
  else { /* nt_cur > 0 */

    /* Build and solve the linear system related to the Richards equations */
    if (!cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      if (cs_equation_needs_build(richards)) // unsteady ?
        cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(time_step, richards);

      /* Compute the darcian flux */
      _update_darcian_flux(connect, cdoq, richards, gw);

      /* Update the moisture content */
      _update_moisture_content(connect, cdoq, richards, gw);

    }

    for (i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_param[i].eq_id];

      cs_equation_init_system(mesh, connect, cdoq, time_step, eq);

      if (!cs_equation_is_steady(eq)) { // unsteady ?

        /* Define the algebraic system */
        if (cs_equation_needs_build(eq))
          cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(time_step, eq);

      } /* Solve this equation which is steady */

    } /* Loop on tracer equations */

  } /* nt_cur > 0 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined postprocessing for the groundwater module
 *
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 * \param[in]  gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_post(const cs_time_step_t      *time_step,
                    const cs_groundwater_t    *gw)
{
  if (gw == NULL)
    return;

  const int  nt_cur = time_step->nt_cur;

  /* Cases where a post-processing is not required */
  if (gw->post_freq == -1)
    return;
  if (nt_cur == 0) {
    if (gw->global_model != CS_GROUNDWATER_MODEL_SATURATED)
      return;
  }
  else { /* nt_cur > 0 */
    if (gw->global_model == CS_GROUNDWATER_MODEL_SATURATED)
      return;
    if (gw->post_freq == 0)
      return;
    if (nt_cur % gw->post_freq > 0)
      return;
  }

  if (gw->flag & CS_GROUNDWATER_POST_MOISTURE) {

    cs_field_t  *f = gw->moisture_content;

    cs_post_write_var(-1,              // id du maillage de post
                      f->name,
                      1,               // dim
                      true,            // interlace
                      true,            // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      f->val,          // values on cells
                      NULL,            // values at internal faces
                      NULL,            // values at border faces
                      time_step);      // time step structure
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
