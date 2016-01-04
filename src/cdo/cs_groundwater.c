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

/* Parameters defining the van Genuchten-Mualen law */
typedef struct {

  double  n;          // 1.25 < n < 6
  double  m;          // m = 1 - 1/n
  double  scale;      // scale parameter [m^-1]
  double  tortuosity; // tortuosity param. for saturated hydraulic conductivity

} cs_gw_genuchten_t;

typedef struct {

  double   h_r;
  double   h_s;

} cs_gw_tracy_t;

/* Set of parameters related to a tracer equation */
typedef struct {

  int              eq_id;

  cs_real_3_t      dispersivity;
  double           bulk_density;
  double           distrib_coef;
  double           reaction_rate;

} cs_gw_tracer_t;

/* set of parameters related to the groundwater module */
struct _groundwater_t {

  cs_groundwater_model_t   model;      /* Physical modelling */
  cs_flag_t                flag;       /* Compact information */
  int                      post_freq;  /* Frequency for post-processing */

  /* Physical parameters related to this module */
  double                   residual_moisture;      // theta_r
  double                   saturated_moisture;     // theta_s
  double                   saturated_permeability; // k_s [m.s^-1]

  /* Parameters for predefined models */
  cs_gw_genuchten_t        genuchten_param;
  cs_gw_tracy_t            tracy_param;

  /* Gravity effect */
  bool                     with_gravitation;
  cs_real_3_t              gravity;

  /* Set of equations associated to this module */
  int                      richards_eq_id;
  int                      n_tracers;
  cs_gw_tracer_t          *tracer_param;

  /* Moisture content variable and attached quantities */
  cs_field_t              *moisture_content;

  /* Scan the c2e connectivity index to get the darcian flux related to
     each dual face when CDO vertex-based scheme is activated */
  cs_real_t               *darcian_flux;
  cs_adv_field_t          *adv_field; /* shared with domain */

  /* Work buffer */
  cs_real_t  *work;

};

/* List of available keys for setting the groundwater module */
typedef enum {

  GWKEY_SATURATED_PERMEABILITY,
  GWKEY_MAX_MOISTURE,
  GWKEY_RESIDUAL_MOISTURE,
  GWKEY_TRACY_HS,
  GWKEY_TRACY_HR,
  GWKEY_POST_FREQ,
  GWKEY_OUTPUT_MOISTURE,
  GWKEY_ERROR

} gwkey_t;

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
 * \brief  Print the name of the corresponding advection key
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

  case GWKEY_SATURATED_PERMEABILITY:
    return "saturated_permeability";
  case GWKEY_MAX_MOISTURE:
    return "max_moisture";
  case GWKEY_RESIDUAL_MOISTURE:
    return "residual_moisture";
  case GWKEY_TRACY_HS:
    return "tracy_hs";
  case GWKEY_TRACY_HR:
    return "tracy_hr";
  case GWKEY_POST_FREQ:
    return "post_freq";
  case GWKEY_OUTPUT_MOISTURE:
    return "output_moisture";

  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of an advection key.
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

  if (strcmp(keyname, "saturated_permeability") == 0)
    key = GWKEY_SATURATED_PERMEABILITY;
  else if (strcmp(keyname, "max_moisture") == 0)
    key = GWKEY_MAX_MOISTURE;
  else if (strcmp(keyname, "residual_moisture") == 0)
    key = GWKEY_RESIDUAL_MOISTURE;
  else if (strcmp(keyname, "tracy_hs") == 0)
    key = GWKEY_TRACY_HS;
  else if (strcmp(keyname, "tracy_hr") == 0)
    key = GWKEY_TRACY_HR;
  else if (strcmp(keyname, "post_freq") == 0)
    key = GWKEY_POST_FREQ;
  else if (strcmp(keyname, "output_moisture") == 0)
    key = GWKEY_OUTPUT_MOISTURE;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_gw_tracer_t structure
 *
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      dispersivity    dispersivity for each axis (x, y, z]
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 * \param[in, out] tp              pointer to a cs_gw_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_tracer_param(int                 tracer_eq_id,
                  cs_real_3_t         dispersivity,
                  double              bulk_density,
                  double              distrib_coef,
                  double              reaction_rate,
                  cs_gw_tracer_t     *tp)
{
  assert(tp != NULL); /* Sanity check */

  tp->eq_id = tracer_eq_id;

  for (int k = 0; k < 3; k++)
    tp->dispersivity[k] = dispersivity[k];

  tp->bulk_density = bulk_density;
  tp->distrib_coef = distrib_coef;
  tp->reaction_rate = reaction_rate;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default parametrization of Van Genuchten-Mualen laws
 */
/*----------------------------------------------------------------------------*/

static cs_gw_genuchten_t
_set_default_genuchten_param(void)
{
  cs_gw_genuchten_t  law;

  law.n = 1.56;
  law.m = 1 - 1/law.n;
  law.scale = 0.036;
  law.tortuosity = 0.5;

  return law;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         van Genuchten-Mualen law
 *
 * \param[in]      h           value of the hydralic head
 * \param[in]      gw_struct   pointer to the groundwater structure
 * \param[in, out] result      pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_permeability_by_genuchten_law(double        h,
                               const void   *gw_struct,
                               cs_get_t     *result)
{
  const cs_groundwater_t  *gw = (const cs_groundwater_t *)gw_struct;
  const cs_gw_genuchten_t  law = gw->genuchten_param;

  double  isoval = gw->saturated_permeability;

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
 * \brief  Set the default parametrization of Tracy model
 */
/*----------------------------------------------------------------------------*/

static cs_gw_tracy_t
_set_default_tracy_param(void)
{
  cs_gw_tracy_t  law;

  law.h_r = -100;
  law.h_s = 0;

  return law;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         van Tracy law
 *
 * \param[in]      h           value of the hydralic head
 * \param[in]      gw_struct   pointer to the groundwater structure
 * \param[in, out] result      pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_permeability_by_tracy_law(double        h,
                           const void   *gw_struct,
                           cs_get_t     *result)
{
  const cs_groundwater_t  *gw = (const cs_groundwater_t *)gw_struct;
  const cs_gw_tracy_t  law = gw->tracy_param;

  const double  ks = gw->saturated_permeability;
  const double  isoval = ks * (h - law.h_r)/(law.h_s - law.h_r);

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
 * \brief  Define the moisture content using the Tracy law
 *
 * \param[in]      h           value of the hydralic head
 * \param[in]      gw_struct   pointer to the groundwater structure
 * \param[in, out] result      pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_moisture_by_tracy_law(double        h,
                       const void   *gw_struct,
                       cs_get_t     *result)
{
  const cs_groundwater_t  *gw = (const cs_groundwater_t *)gw_struct;
  const cs_gw_tracy_t  law = gw->tracy_param;
  const double  delta_theta = gw->saturated_moisture - gw->residual_moisture;

  double  k_r = (h - law.h_r)/(law.h_s - law.h_r);
  double  moisture = k_r * delta_theta + gw->residual_moisture;

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
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards    pointer to the Richards equation structure
 * \param[in, out] gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_moisture_content(const cs_cdo_quantities_t   *cdoq,
                         const cs_equation_t         *richards,
                         cs_groundwater_t            *gw)
{
  cs_lnum_t  i;
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

  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(moisture->location_id);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(moisture);

  /* Up to now groundwater module is discretized using vertex-based schemes */
  assert(cdoq->n_vertices == n_elts[0]);

  switch (gw->model) {

  case CS_GROUNDWATER_MODEL_TRACY:
    for (i = 0; i < cdoq->n_vertices; i++) {

      _moisture_by_tracy_law(h->val[i], (const void *)gw, &get);
      moisture->val[i] = get.val;

    } // Loop on vertices
    break;

  default: // include "saturated"
    break; // Nothing to do

  } /* Switch according to the kind of modelling */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \return a pointer to a new allocated cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_groundwater_create(void)
{
  cs_groundwater_t  *gw = NULL;

  BFT_MALLOC(gw, 1, cs_groundwater_t);

  /* Default initialization */
  gw->model = CS_GROUNDWATER_N_MODELS;
  gw->flag = 0;
  gw->post_freq = -1;

  gw->residual_moisture = 0.0;
  gw->saturated_moisture = 1.0;
  gw->saturated_permeability = 1.0;

  gw->genuchten_param = _set_default_genuchten_param();
  gw->tracy_param = _set_default_tracy_param();

  gw->with_gravitation = false;
  gw->gravity[0] = 0, gw->gravity[1] = 0, gw->gravity[2] = 0;

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

  BFT_FREE(gw->tracer_param);
  BFT_FREE(gw->darcian_flux);
  BFT_FREE(gw->work);

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

  case GWKEY_SATURATED_PERMEABILITY:
    gw->saturated_permeability = atof(keyval);
    break;
  case GWKEY_MAX_MOISTURE:
    gw->saturated_moisture = atof(keyval);
    break;
  case GWKEY_RESIDUAL_MOISTURE:
    gw->residual_moisture = atof(keyval);
    break;
  case GWKEY_TRACY_HS:
    gw->tracy_param.h_s = atof(keyval);
    break;
  case GWKEY_TRACY_HR:
    gw->tracy_param.h_r = atof(keyval);
    break;
  case GWKEY_POST_FREQ:
    gw->post_freq = atoi(keyval);
    break;
  case GWKEY_OUTPUT_MOISTURE:
    if (strcmp(keyval, "false")) // not "false"
      gw->flag |= CS_GROUNDWATER_POST_MOISTURE;
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

  bft_printf("  <GW/Tracer> n_tracer_equations %d\n", gw->n_tracers);
  bft_printf("  <GW/Parameters>");
  bft_printf(" residual_moisture: %5.3e", gw->residual_moisture);
  bft_printf(" saturated_moisture: %5.3e\n", gw->saturated_moisture);
  bft_printf("  <GW/Parameters>");
  bft_printf(" saturated_permeability: %5.3e\n", gw->saturated_permeability);

  if (gw->with_gravitation)
    bft_printf("  <GW/Gravitation> true\n");
  else
    bft_printf("  <GW/Gravitation> false\n");

  switch (gw->model) {
  case CS_GROUNDWATER_MODEL_SATURATED:
    bft_printf("  <GW/Model> saturated\n");
    break;
  case CS_GROUNDWATER_MODEL_GENUCHTEN:
    bft_printf("  <GW/Model> VanGenuchten-Mualen\n");
    break;
  case CS_GROUNDWATER_MODEL_TRACY:
    bft_printf("  <GW/Model> Tracy\n");
    break;
  case CS_GROUNDWATER_MODEL_USER:
    bft_printf("  <GW/Model> User-defined\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid model for groundwater module.\n"
              " Please check your settings.");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]      connect          pointer to a cs_cdo_connect_t structure
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in]      model            keyword related to the model used
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_groundwater_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_init(const cs_cdo_connect_t  *connect,
                    int                      richards_eq_id,
                    const char              *model,
                    cs_property_t           *permeability,
                    cs_property_t           *soil_capacity,
                    cs_adv_field_t          *adv_field,
                    cs_groundwater_t        *gw)
{
  cs_equation_t  *eq = NULL;

  const cs_connect_index_t  *c2e = connect->c2e;

  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Create a new equation structure for Richards' equation */
  gw->richards_eq_id = richards_eq_id;
  eq = cs_equation_create("Richards",                   // equation name
                          "hydraulic_head",             // variable name
                          CS_EQUATION_TYPE_GROUNDWATER, // type of equation
                          CS_PARAM_VAR_SCAL,            // type of variable
                          CS_PARAM_BC_HMG_NEUMANN);     // default BC

  /* Define and associate properties according to the type of modelling */
  if (strcmp(model, "saturated") == 0) {

    assert(soil_capacity == NULL); // Sanity check
    gw->model = CS_GROUNDWATER_MODEL_SATURATED;

    /* Default initialization */
    const char  identity_val[] = "1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0";
    cs_property_def_by_value(permeability, identity_val);

  }
  else {

    if (strcmp(model, "genutchten") == 0) {

      gw->model = CS_GROUNDWATER_MODEL_GENUCHTEN;
      cs_property_def_by_law(permeability, _permeability_by_genuchten_law);

    }
    else if (strcmp(model, "tracy") == 0) {

      gw->model = CS_GROUNDWATER_MODEL_TRACY;
      cs_property_def_by_law(permeability, _permeability_by_tracy_law);

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " Incompatible model for groundwater flows.\n"
                " Value given: %s\n"
                " Availaible models: saturated, genutchen, tracy", model);

    /* Associate soil_capacity to the unsteady term of the Richards eq. */
    cs_equation_link(eq, "time", soil_capacity);

  }

  /* Associate permeability to the diffusion property of the Richards eq. */
  cs_equation_link(eq, "diffusion", permeability);

  /* Advection field induced by the hydraulic head */
  gw->adv_field = adv_field;

  BFT_MALLOC(gw->darcian_flux, c2e->idx[connect->c_info->n_ent], cs_real_t);
  for (cs_lnum_t i = 0; i < c2e->idx[connect->c_info->n_ent]; i++)
    gw->darcian_flux[i] = 0;

  /* Field related to the moisture content */


  /* Work buffer */
  BFT_MALLOC(gw->work, connect->n_max_ebyc, cs_real_t);

  return eq;
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
 * \param[in]      diff_property   pointer to a cs_property_t struct.
 * \param[in]      dispersivity    dispersivity for each axis (x, y, z]
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
                          cs_property_t       *diff_property,
                          cs_real_3_t          dispersivity,
                          double               bulk_density,
                          double               distrib_coef,
                          double               reaction_rate)
{
  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_equation_t  *eq = NULL;

  BFT_REALLOC(gw->tracer_param, gw->n_tracers + 1, cs_gw_tracer_t);

  _set_tracer_param(tracer_eq_id,
                    dispersivity,
                    bulk_density,
                    distrib_coef,
                    reaction_rate,
                    gw->tracer_param + gw->n_tracers);

  gw->n_tracers += 1;

  eq = cs_equation_create(eqname,                       // equation name
                          varname,                      // variable name
                          CS_EQUATION_TYPE_GROUNDWATER, // type of equation
                          CS_PARAM_VAR_SCAL,            // type of variable
                          CS_PARAM_BC_HMG_NEUMANN);     // default BC

  /* Associate permeability to the diffusion property of the Richards eq. */
  cs_equation_link(eq, "diffusion", diff_property);

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
  cs_flag_t  flag;
  cs_equation_t  *richards = equations[gw->richards_eq_id];

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_space_scheme(richards) == CS_SPACE_SCHEME_CDOVB);

  /* Moisture content */
  bool has_previous = cs_equation_is_steady(richards) ? false:true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  int  location_id = cs_mesh_location_get_id_by_name(N_("vertices"));

  gw->moisture_content = cs_field_create("moisture_content",
                                         field_mask,
                                         location_id,
                                         1,     // dimension
                                         true,  // interleave
                                         has_previous);
  cs_field_allocate_values(gw->moisture_content);

  /* Default initialization of the moisture content */
  const cs_lnum_t *n_elts =
    cs_mesh_location_get_n_elts(gw->moisture_content->location_id);
  for (cs_lnum_t  i = 0; i < n_elts[0]; i++)
    gw->moisture_content->val[i] = gw->saturated_moisture;

  /* Permeability settings for the diffusion term */
  cs_property_t  *permeability = cs_equation_get_diffusion_property(richards);
  switch (gw->model) {
  case CS_GROUNDWATER_MODEL_GENUCHTEN:
  case CS_GROUNDWATER_MODEL_TRACY:
    {
      cs_field_t  *hydraulic_head = cs_equation_get_field(richards);

      flag = CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;
      cs_property_set_array(permeability, flag, hydraulic_head->val);
      cs_property_set_struct(permeability, (const void *)gw);

    }
    break;
  case CS_GROUNDWATER_MODEL_SATURATED:
    { /* Anisotropic by construction */
      double  tens[9] = {0.,0.,0., 0.,0.,0., 0.,0.,0.};

      tens[0] = tens[4] = tens[8] = gw->saturated_permeability;
      /* Modify the value of a property defined by value.
         Be careful when using this function (do not use it in a loop on cells
         for instance) */
      cs_property_set_value(permeability, tens);

    }
    break;

  default:
    // Nothing to do
    break;

  } // switch on modelling

  /* Soil capacity settings (related to unsteady term) */
  if (gw->model == CS_GROUNDWATER_MODEL_TRACY) {

    cs_property_t  *capacity = cs_equation_get_time_property(richards);

    const cs_gw_tracy_t  law = gw->tracy_param;
    const double  delta_theta = gw->saturated_moisture-gw->residual_moisture;
    const double  delta_h = law.h_s - law.h_r;

    char cval[16];
    sprintf(cval, "%10.8e", delta_theta/delta_h);
    cs_property_def_by_value(capacity, cval);

  }

  /* Define and then link the advection field to each tracer equations */
  flag = CS_PARAM_FLAG_FACE | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_BY_CELL;
  flag |= CS_PARAM_FLAG_SCAL;

  cs_advection_field_def_by_array(gw->adv_field, flag, gw->darcian_flux);

  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gw_tracer_t  tp = gw->tracer_param[i];
    cs_equation_t  *eq = equations[tp.eq_id];

    cs_equation_link(eq, "advection", gw->adv_field);

  }

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
      _update_moisture_content(cdoq, richards, gw);

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
      _update_moisture_content(cdoq, richards, gw);

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
    if (gw->model != CS_GROUNDWATER_MODEL_SATURATED)
      return;
  }
  else { /* nt_cur > 0 */
    if (gw->model == CS_GROUNDWATER_MODEL_SATURATED)
      return;
    if (gw->post_freq == 0)
      return;
    if (nt_cur % gw->post_freq > 0)
      return;
  }

  if (gw->flag & CS_GROUNDWATER_POST_MOISTURE) {

    cs_field_t  *f = gw->moisture_content;

    cs_post_write_vertex_var(-1,              // id du maillage de post
                             f->name,
                             1,               // dim
                             true,            // interlace
                             true,            // true = original mesh
                             CS_POST_TYPE_cs_real_t,
                             f->val,          // values on vertices
                             time_step);      // time step structure
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
