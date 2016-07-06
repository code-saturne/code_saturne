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

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_post.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_math.h"
#include "cs_param.h"
#include "cs_reco.h"
#include "cs_hodge.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_groundwater.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Tag dedicated to build a flag for the groundwater module */
/*   1: post the moisture content */
#define CS_GROUNDWATER_POST_MOISTURE  (1 <<  0)

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Set of parameters related to a tracer equation */
typedef struct {

  /* Bulk density times the distribution coefficient for each soil */
  double           rho_kd;

  /* Longitudinal and transversal dispersivity for each soil */
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

  /* Set of parameters for each tracer equation */
  cs_gw_tracer_t         *tracer_param;

  /* Physical modelling adopted for this soil */
  cs_groundwater_model_t  model;

  /* Parameters for predefined models */
  cs_gw_genuchten_t       genuchten_param; /* Van-Genuchten-Mualen law */
  cs_gw_tracy_t           tracy_param;     /* Tracy law */

  /* Main soil properties */
  double                  residual_moisture;       /* theta_r */
  double                  saturated_moisture;      /* theta_s */
  double                  delta_moisture;          /* theta_s - theta_r */
  cs_get_t                saturated_permeability;  /* k_s [m.s^-1] */

} cs_gw_soil_t;

/* Set of parameters related to the groundwater module */
struct _groundwater_t {

  cs_flag_t                flag;       /* Compact information */

  cs_groundwater_model_t   global_model;

  /* Physical parameters related to each kind of soil considered.
     If n_soils > 1, soil_id array stores the id giving access to the soil
     parameters related to each cell of the mesh */
  int            n_soils;
  int            n_max_soils;
  cs_gw_soil_t  *soil_param;
  cs_lnum_t      n_cells;     /* Number of cells (useful to get soil_id) */
  short int     *soil_id;     /* NULL or allocated to n_cells */

  /* Gravity effect */
  bool           with_gravitation;
  cs_real_3_t    gravity;

  /* Fields located at cells */
  cs_field_t    *hydraulic_head;  /* H (shared or owner according to the
                                     discretization space scheme) */
  cs_field_t    *pressure_head;   /* h = H - gravity_potential */

  /* Set of equations associated to this module */
  int            richards_eq_id;

  int            n_tracers;
  int            n_max_tracers;
  int           *tracer_eq_ids;

  /* Moisture content variable and attached quantities */
  cs_field_t          *moisture_content;

  /* Permeability is the diffusion property related to Richards equation but
     this property plays also a role in the diffusion of tracer equations */
  cs_property_t       *permeability;  /* shared with a cs_domain_t structure */

  /* Scan the c2e connectivity index to get the darcian flux related to
     each dual face when CDO vertex-based scheme is activated */
  cs_real_t           *darcian_flux;
  cs_adv_field_t      *adv_field;    /* shared with a cs_domain_t structure */

  /* Work buffer (allocated to n_max_ebyc) */
  cs_real_t  *work;

};

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
 * \brief  Set a cs_gw_tracer_t structure
 *
 * \param[in, out] tp              pointer to a cs_gw_tracer_t structure
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

static void
_set_tracer_param(cs_gw_tracer_t     *tp,
                  double              wmd,
                  double              alpha_l,
                  double              alpha_t,
                  double              bulk_density,
                  double              distrib_coef,
                  double              reaction_rate)
{
  assert(tp != NULL); /* Sanity check */

  tp->wmd = wmd;
  tp->rho_kd = bulk_density * distrib_coef;
  tp->reaction_rate = reaction_rate;
  tp->alpha_l = alpha_l;
  tp->alpha_t = alpha_t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_gw_soil_t structure (already allocated)
 *
 * \param[in]      ml_name         name of the mesh location
 * \param[in]      model_kw        keyword related to the modelling
 * \param[in]      n_tracers       number of related tracer eqs
 * \param[in, out] soil            pointer to a cs_gw_soil_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_soil(const char     *ml_name,
           const char     *model_kw,
           int             n_tracers,
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
  else if (strcmp(model_kw, "genuchten") == 0) {

    soil->model = CS_GROUNDWATER_MODEL_GENUCHTEN;

    /* Default initialization */
    soil->saturated_moisture = 0.75;
    soil->residual_moisture = 0.15;

    const double  n = 1.56;
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
              " Availaible models: saturated, genuchten, tracy", model_kw);

  soil->delta_moisture = soil->saturated_moisture - soil->residual_moisture;

  /* Set of parameters for each tracer which are related to this soil */
  BFT_MALLOC(soil->tracer_param, n_tracers, cs_gw_tracer_t);

  for (int i = 0; i < n_tracers; i++) /* default initialization */
    _set_tracer_param(soil->tracer_param + i,
                      0.0,  /* water molecular diffusivity */
                      0.0,  /* alpha_l */
                      0.0,  /* alpha_t */
                      1.0,  /* bulk density */
                      0.0,  /* Kd (distribution coef.) */
                      0.0); /* reaction rate */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the id in a cs_groundwater_t structure from the tracer
 *         equation id
 *
 * \param[in]   gw             pointer to a cs_groundwater_t structure
 * \param[in]   tracer_eq_id   tracer equation id
 *
 * \returns an id related to this tracer equation id
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_tracer_id(const cs_groundwater_t   *gw,
               int                       tracer_eq_id)
{
  int  tracer_id = -1;

  for (int id = 0; id < gw->n_tracers; id++) {
    if (gw->tracer_eq_ids[id] == tracer_eq_id) {
      tracer_id = id;
      break;
    }
  }

  if (tracer_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              " Stop setting a tracer equation. Its identification number has"
              " not been found in the groundwater flow module.\n"
              " Please check your settings.");

  return tracer_id;
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

static inline void
_get_tracer_diffusion_tensor(double          theta,
                             const double    v[],
                             const void     *tracer_struc,
                             cs_get_t       *result)
{
  const cs_gw_tracer_t  *tp = (const cs_gw_tracer_t *)tracer_struc;

  const double v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
  const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
  const double  coef1 = tp->wmd * theta + tp->alpha_t*vnorm;

  double  delta_coef = 0.;
  if (vnorm > cs_math_zero_threshold)
    delta_coef = (tp->alpha_l - tp->alpha_t)/vnorm;

  const double  dcv[3] = {delta_coef*v[0], delta_coef*v[1], delta_coef*v[2]};

  for (int ki = 0; ki < 3; ki++) {

    /* Diagonal terms */
    result->tens[ki][ki] = coef1 + delta_coef*v2[ki];

    /* Extra-diagonal terms */
    for (int kj = ki + 1; kj < 3; kj++)
      result->tens[ki][kj] = dcv[ki]*v[kj];
  }

  /* Diffusion tensor is symmetric by construction */
  result->tens[1][0] = result->tens[0][1];
  result->tens[2][0] = result->tens[0][2];
  result->tens[2][1] = result->tens[1][2];

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

static inline void
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

static inline void
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

static inline void
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
    const double  one_alpha_hn = 1 + pow(fabs(law.scale*h), law.n);
    const double  se = pow(one_alpha_hn, -law.m);
    const double  se_pow_L = pow(se, law.tortuosity);
    const double  se_pow_overm = pow(se, 1/law.m);
    const double  coef_base = 1 - pow(1 - se_pow_overm, law.m);

    isoval *= se_pow_L * coef_base*coef_base;

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
 * \brief  Define the moisture content using the Van Genuchten law
 *
 * \param[in]      h             value of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_moisture_by_genuchten_law(double        h,
                           const void   *soil_struc,
                           cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_genuchten_t  law = soil->genuchten_param;

  double  Se = 1; // dimensionless moisture

  if (h < 0) {
    double  tmp_coef = pow(fabs(law.scale * h), law.n);
    Se = pow(1 + tmp_coef, -law.m);
  }

  /* Return the computed moisture content */
  result->val = Se*soil->delta_moisture + soil->residual_moisture;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the moisture content using the Van Genuchten law
 *
 * \param[in]      h             value of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        pointer to a cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_capacity_by_genuchten_law(double        h,
                           const void   *soil_struc,
                           cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_genuchten_t  law = soil->genuchten_param;

  result->val = 0.; // default initialization for saturated soil

  if (h >= 0)
    return;

  const double  mult_coef = -law.n * law.m * soil->delta_moisture;
  const double  alpha_h_pow_n = pow(fabs(law.scale * h), law.n);
  const double  se_m1 = pow(1 + alpha_h_pow_n, -law.m-1);

  result->val = mult_coef * alpha_h_pow_n/h * se_m1;
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

static inline void
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

static inline void
_moisture_by_tracy_law(double        h,
                       const void   *soil_struc,
                       cs_get_t     *result)
{
  const cs_gw_soil_t  *soil = (const cs_gw_soil_t  *)soil_struc;
  const cs_gw_tracy_t  law = soil->tracy_param;

  double  k_r = (h - law.h_r)/(law.h_s - law.h_r);
  double  moisture = k_r * soil->delta_moisture + soil->residual_moisture;

  result->val = moisture;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure head from the value of the hydraulic head
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      richards    pointer to the Richards equation structure
 * \param[in, out] gw          pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_head(const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *cdoq,
             const cs_equation_t         *richards,
             cs_groundwater_t            *gw)
{
  /* Sanity checks */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0,  _err_empty_gw);
  if (richards == NULL)
    bft_error(__FILE__, __LINE__, 0," Richards eq. is not allocated.");

  cs_field_t  *h_head = gw->hydraulic_head;

  /* Sanity check */
  assert(h_head->location_id == cs_mesh_location_get_id_by_name(N_("cells")));

  switch (cs_equation_get_space_scheme(richards)) {
  case CS_SPACE_SCHEME_CDOVB:
    {
      /* Copy current field values to previous values */
      cs_field_current_to_previous(h_head);

      cs_field_t  *h_head_vtx = cs_equation_get_field(richards);

      /* Sanity check */
      assert(h_head_vtx->location_id ==
             cs_mesh_location_get_id_by_name(N_("vertices")));

      /* Reconstruct (or interpolate) values at the cell centers */
# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++)
        cs_reco_pv_at_cell_center(c_id, connect->c2v, cdoq, h_head_vtx->val,
                                  h_head->val + c_id);

    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *hh_vals = cs_equation_get_cell_values(richards);

      /* Copy current field values to previous values */
      cs_field_current_to_previous(h_head);

      /* Copy values at cells center into this field */
      memcpy(h_head->val, hh_vals, sizeof(cs_real_t)*cdoq->n_cells);
    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
    break; // Nothing to do (h_head is a pointer to richards field

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

  } // Switch on space scheme

  if (gw->with_gravitation) {

    cs_field_t  *p_head = gw->pressure_head;

    /* Sanity checks */
    if (p_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */
    cs_field_current_to_previous(p_head);

# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      const cs_real_t  gpot = cs_math_3_dot_product(cdoq->cell_centers + 3*c_id,
                                                    gw->gravity);

      p_head->val[c_id] = h_head->val[c_id] - gpot;

    }

  } // Gravitation effect is activated

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
  /* Sanity checks */
  if (richards == NULL || gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Groundwater module or Richards eq. is not allocated.");

  cs_field_t  *moisture = gw->moisture_content;

  /* Sanity checks */
  if (moisture == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " The field related to the moisture content is not allocated.");

  /* Copy current field values to previous values */
  cs_field_current_to_previous(moisture);

  /* Moisture content is define in each cell along with the hydraulic (or
     pressure head */

  const cs_field_t  *h;
  if (gw->with_gravitation)
    h = gw->pressure_head;
  else
    h = gw->hydraulic_head;

  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gw_soil_t  *soil = gw->soil_param + soil_id;
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

    switch (soil->model) {

    case CS_GROUNDWATER_MODEL_GENUCHTEN:

      if (elt_ids == NULL) {
# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

          cs_get_t  get;

          _moisture_by_genuchten_law(h->val[c_id], (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on cells

      }
      else {

        const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

          cs_get_t  get;
          const cs_lnum_t  c_id = elt_ids[i];

          _moisture_by_genuchten_law(h->val[c_id], (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on cells

      } // cell selection is a part of all cells
      break;

    case CS_GROUNDWATER_MODEL_TRACY:
      if (elt_ids == NULL) {

# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

          cs_get_t  get;

          _moisture_by_tracy_law(h->val[c_id], (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on cells

      }
      else {

        const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

          cs_get_t  get;
          const cs_lnum_t  c_id = elt_ids[i];

          _moisture_by_tracy_law(h->val[c_id], (const void *)soil, &get);
          moisture->val[c_id] = get.val;

        } // Loop on cells

      } // cell selection is a part of all cells
      break;

    case CS_GROUNDWATER_MODEL_SATURATED:
      if (elt_ids == NULL) {

# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id <  cdoq->n_cells; c_id++)
          moisture->val[c_id] = soil->saturated_moisture;

      }
      else {

        const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++)
          moisture->val[elt_ids[i]] = soil->saturated_moisture;

      } // cell selection is a part of all cells
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of model for estimating the moisture"
                  " content."));
      break; // Nothing to do

    } // Switch on the type of soil modelling

  } // Loop on soils

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the setting is correct
 *
 * \param[in]  gw     pointer to a cs_groundwater_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_check_settings(const cs_groundwater_t  *gw)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  for (int i = 0; i < gw->n_max_tracers; i++)
    if (gw->tracer_eq_ids[i] == -1)
      bft_error(__FILE__, __LINE__, 0,
                " At least one tracer equation has not been set.\n"
                " %d tracer equations have heen initialy requested.\n"
                " Please check your settings.", gw->n_max_tracers);


  if (gw->n_soils < gw->n_max_soils)
    bft_error(__FILE__, __LINE__, 0,
              " %d soils are defined but %d have been initially requested."
              " Please check your settings.", gw->n_soils, gw->n_max_soils);

  if (gw->n_soils > 1) {

    cs_lnum_t  n_unset_cells = 0;
    int  cell_ml_id = cs_mesh_location_get_id_by_name(N_("cells"));
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(cell_ml_id);

# pragma omp parallel reduction(+:n_unset_cells) if (n_elts[0] > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts[0]; i++)
      if (gw->soil_id[i] == gw->n_max_soils) n_unset_cells++;

    if (n_unset_cells > 0)
      bft_error(__FILE__, __LINE__, 0,
                " %d cells are not associated to any soil.\n"
                " Please check your settings.", n_unset_cells);

  } // n_soils > 1

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
  gw->flag = 0;

  gw->n_soils = 0;
  gw->n_max_soils = 0;
  gw->soil_param = NULL;
  gw->soil_id = NULL;

  gw->global_model = CS_GROUNDWATER_N_MODELS;

  gw->with_gravitation = false;
  gw->gravity[0] = 0, gw->gravity[1] = 0, gw->gravity[2] = 0;
  gw->hydraulic_head = NULL;
  gw->pressure_head = NULL;

  gw->richards_eq_id = -1;
  gw->n_tracers = 0;
  gw->n_max_tracers = 0;
  gw->tracer_eq_ids = NULL;

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

  BFT_FREE(gw->tracer_eq_ids);
  BFT_FREE(gw->darcian_flux);
  BFT_FREE(gw->work);

  for (int i = 0; i < gw->n_soils; i++) {
    cs_gw_soil_t *soil = gw->soil_param + i;
    BFT_FREE(soil->tracer_param);
  }

  if (gw->n_soils > 1)
    BFT_FREE(gw->soil_id);

  BFT_FREE(gw->soil_param);

  BFT_FREE(gw);

  /* Fields, advection fields and properties are freed elsewhere */

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of requested soils
 *
 * \param[in]  gw        pointer to a cs_groundwater_t structure
 *
 * \return the number of requested soils
 */
/*----------------------------------------------------------------------------*/

int
cs_groundwater_get_n_soils(const cs_groundwater_t    *gw)
{
  if (gw == NULL)
    return 0;

  return gw->n_max_soils;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_groundwater_t structure
 *
 * \param[in, out]  gw        pointer to a cs_groundwater_t structure
 * \param[in]       key       key related to the member of gw to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_param(cs_groundwater_t      *gw,
                         cs_groundwater_key_t   key,
                         const char            *keyval)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_GWKEY_GRAVITATION:
    gw->with_gravitation = true;
    gw->gravity[0] = gw->gravity[1] = gw->gravity[2] = 0.;
    if (strcmp(val, "x") == 0)       gw->gravity[0] =  1.;
    else if (strcmp(val, "-x") == 0) gw->gravity[0] = -1.;
    else if (strcmp(val, "y") == 0)  gw->gravity[1] =  1.;
    else if (strcmp(val, "-y") == 0) gw->gravity[1] = -1.;
    else if (strcmp(val, "z") == 0)  gw->gravity[2] =  1.;
    else if (strcmp(val, "-z") == 0) gw->gravity[2] = -1.;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid choice of gravitation axis: %s.\n"
                  " Available choices are 'x', 'y' and 'z'\n"
                  " Please check your settings."), val);
    break;

  case CS_GWKEY_OUTPUT_MOISTURE:
    if (strcmp(val, "false")) // not "false"
      gw->flag |= CS_GROUNDWATER_POST_MOISTURE;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Key not implemented yet."));

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

  /* Sanity checks */
  _check_settings(gw);

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

  } // Loop on soils

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

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]      connect          pointer to a cs_cdo_connect_t structure
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in]      n_soils          number of soils to consider
 * \param[in]      n_tracers        number of tracers to consider
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
                          int                      n_soils,
                          int                      n_tracer_eqs,
                          cs_property_t           *permeability,
                          cs_property_t           *soil_capacity,
                          cs_adv_field_t          *adv_field,
                          cs_groundwater_t        *gw)
{
  cs_equation_t  *eq = NULL;

  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_lnum_t  n_cells = connect->c_info->n_elts;

  /* Sanity check */
  assert(n_soils > 0);

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

  /* Up to now Richards equation is only set with vertex-based schemes */
  BFT_MALLOC(gw->darcian_flux, c2e->idx[n_cells], cs_real_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2e->idx[n_cells]; i++)
    gw->darcian_flux[i] = 0;

  /* Work (temporary) buffer */
  BFT_MALLOC(gw->work, connect->n_max_ebyc, cs_real_t);

  /* Quantities related to soils */
  gw->n_soils = 0;           /* No soil is set at the beginning */
  gw->n_max_soils = n_soils; /* Max. number of soils allocated */
  BFT_MALLOC(gw->soil_param, n_soils, cs_gw_soil_t);

  if (n_soils > 1) {
    BFT_MALLOC(gw->soil_id, n_cells, short int);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      gw->soil_id[i] = n_soils; /* Default value => not set */
  }

  gw->n_tracers = 0;
  gw->n_max_tracers = n_tracer_eqs;
  BFT_MALLOC(gw->tracer_eq_ids, n_tracer_eqs, int);
  for (int i = 0; i < n_tracer_eqs; i++)
    gw->tracer_eq_ids[i] = -1; /* Default initialization = not set */

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
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  soil_id = gw->n_soils;

  gw->n_soils += 1;
  if (gw->n_soils > gw->n_max_soils)
    bft_error(__FILE__, __LINE__, 0,
              " Maximum number of soils is reached.\n"
              " Stop adding a new soil by value (mesh location: %s).\n"
              " Please modify your settings.", ml_name);

  cs_gw_soil_t  *soil = gw->soil_param + soil_id;

  _init_soil(ml_name, model_kw, gw->n_max_tracers, soil);

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
 * \param[in]       key       key related to a member of the soil to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_soil_param(cs_groundwater_t          *gw,
                              const char                *ml_name,
                              cs_groundwater_soilkey_t   key,
                              const char                *keyval)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  ml_id = -1;
  if (ml_name != NULL) {
    ml_id = cs_mesh_location_get_id_by_name(ml_name);
    if (ml_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location name %s.\n"
                  " This mesh location is not already defined.\n"), ml_name);
  }

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  /* Set different available keys */
  switch(key) {

  case CS_SOILKEY_SAT_MOISTURE:
    {
      double  theta_s = atof(val);

      if (ml_id == -1) {

        for (int i = 0; i < gw->n_soils; i++) {
          gw->soil_param[i].saturated_moisture = theta_s;
          /* Update delta_moisture */
          gw->soil_param[i].delta_moisture =
            theta_s - gw->soil_param[i].residual_moisture;
        } // Loop on soils

      } // All soils are updated
      else {

        for (int i = 0; i < gw->n_soils; i++) {
          if (ml_id == gw->soil_param[i].ml_id) {
            gw->soil_param[i].saturated_moisture = theta_s;
            /* Update delta_moisture */
            gw->soil_param[i].delta_moisture =
              theta_s - gw->soil_param[i].residual_moisture;
          }
        } // Loop on soils

      }
    }
    break;

  case CS_SOILKEY_RES_MOISTURE:
    {
      double  theta_r = atof(val);

      if (ml_id == -1) {

        for (int i = 0; i < gw->n_soils; i++) {
          gw->soil_param[i].residual_moisture = theta_r;
          /* Update delta_moisture */
          gw->soil_param[i].delta_moisture =
            gw->soil_param[i].saturated_moisture - theta_r;
        } // Loop on soils

      }
      else {

        for (int i = 0; i < gw->n_soils; i++) {
          if (ml_id == gw->soil_param[i].ml_id) {
            gw->soil_param[i].residual_moisture = theta_r;
            /* Update delta_moisture */
            gw->soil_param[i].delta_moisture =
              gw->soil_param[i].saturated_moisture - theta_r;
          }
        } // Loop on soils

      }
    }
    break;

  case CS_SOILKEY_TRACY_SAT_H:
    {
      double  h_s = atof(val);

      if (ml_id == -1) {
        for (int i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GROUNDWATER_MODEL_TRACY)
            gw->soil_param[i].tracy_param.h_s = h_s;
      }
      else {
        for (int i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].tracy_param.h_s = h_s;
      }

    }
    break;

  case CS_SOILKEY_TRACY_RES_H:
    {
      double  h_r = atof(val);

      if (ml_id == -1) {

        for (int i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GROUNDWATER_MODEL_TRACY)
            gw->soil_param[i].tracy_param.h_r = h_r;

      }
      else {

        for (int i = 0; i < gw->n_soils; i++)
          if (ml_id == gw->soil_param[i].ml_id)
            gw->soil_param[i].tracy_param.h_r = h_r;

      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Key not implemented"));

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a specific unsteady advection/diffusion/reaction eq.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion and reaction parameters result from a physical modelling.
 *
 * \param[in, out] gw              pointer to a cs_groundwater_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      eqname          name of the equation
 * \param[in]      varname         name of the related variable
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_groundwater_add_tracer(cs_groundwater_t    *gw,
                          int                  tracer_eq_id,
                          const char          *eqname,
                          const char          *varname)
{
  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->n_soils != gw->n_max_soils)
    bft_error(__FILE__, __LINE__, 0,
              " Add a tracer but not all soils are defined (%d/%d).\n"
              " Stop adding a new tracer equation (%s).\n"
              " Please check your settings.",
              gw->n_soils, gw->n_max_soils, eqname);

  int  tracer_id = gw->n_tracers;
  cs_equation_t  *eq = NULL;

  eq = cs_equation_create(eqname,                       // equation name
                          varname,                      // variable name
                          CS_EQUATION_TYPE_GROUNDWATER, // type of equation
                          CS_PARAM_VAR_SCAL,            // type of variable
                          CS_PARAM_BC_HMG_NEUMANN);     // default BC

  gw->n_tracers += 1;
  if (gw->n_tracers > gw->n_max_tracers)
    bft_error(__FILE__, __LINE__, 0,
              _(" Maximum number of tracers is reached.\n"
                " Stop adding the tracer equation %s.\n"
                " Please modify your settings."), eqname);

  BFT_REALLOC(gw->tracer_eq_ids, gw->n_tracers, int);
  gw->tracer_eq_ids[tracer_id] = tracer_eq_id;

  /* Associate the advection field for the advection term */
  assert(gw->adv_field != NULL); /* Sanity check */
  cs_equation_link(eq, "advection", gw->adv_field);

  /* Set default option */
  cs_equation_set_param(eq, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
  cs_equation_set_param(eq, CS_EQKEY_ITSOL, "bicg");
  cs_equation_set_param(eq, CS_EQKEY_BC_ENFORCEMENT, "weak");

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
 * \param[in]      ml_name         name of the related mesh location
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_set_tracer_param(cs_groundwater_t    *gw,
                                int                  tracer_eq_id,
                                const char          *ml_name,
                                double               wmd,
                                double               alpha_l,
                                double               alpha_t,
                                double               bulk_density,
                                double               distrib_coef,
                                double               reaction_rate)
{
  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  tracer_id = _get_tracer_id(gw, tracer_eq_id);

  /* Look for the related soil */
  if (ml_name == NULL) { /* All soils have to be set for this tracer */

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gw_soil_t  *soil = gw->soil_param + soil_id;

      /* Set tracer parameters */
      _set_tracer_param(soil->tracer_param + tracer_id,
                        wmd,
                        alpha_l, alpha_t,
                        bulk_density, distrib_coef,
                        reaction_rate);

    } // Loop on soils

  }
  else { /* Set this tracer equation for a specific soil */

    int ml_id = cs_mesh_location_get_id_by_name(ml_name);
    if (ml_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _( " Stop setting a tracer equation."
                   " Invalid mesh location name %s.\n"
                   " This mesh location is not already defined.\n"), ml_name);

    int  soil_id = -1;
    for (int id = 0; id < gw->n_soils; id++) {
      if (gw->soil_param[id].ml_id == ml_id) {
        soil_id = id;
        break;
      }
    }

    if (soil_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Stop setting a tracer equation."
                  " No soil related to mesh location %s has been found.\n"
                  " Please check your settings."), ml_name);

    cs_gw_soil_t  *soil = gw->soil_param + soil_id;

    /* Set tracer parameters */
    _set_tracer_param(soil->tracer_param + tracer_id,
                      wmd,
                      alpha_l, alpha_t,
                      bulk_density, distrib_coef,
                      reaction_rate);

  } /* Set a specific couple (tracer,soil) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the Richards equation
 *
 * \param[in, out] gw        pointer to a cs_groundwater_t structure
 * \param[in, out] richards  pointer to the related cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_richards_setup(cs_groundwater_t    *gw,
                              cs_equation_t       *richards)
{
  /* Sanity checks */
  assert(richards != NULL);
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->n_soils == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is activated but no soil is defined."));

  const cs_space_scheme_t  cdo_scheme = cs_equation_get_space_scheme(richards);
  if (cdo_scheme == CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0,
              _(" Richards eq. is only available for vertex-based schemes."));

  cs_property_t  *permeability = gw->permeability;
  bool has_previous = cs_equation_is_steady(richards) ? false:true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  int  cell_loc_id = cs_mesh_location_get_id_by_name(N_("cells"));

  /* Create a moisture field attached to cells */
  gw->moisture_content = cs_field_create("moisture_content",
                                         field_mask,
                                         cell_loc_id,
                                         1,        // dimension
                                         has_previous);

  /* Allocate and initialize values */
  cs_field_allocate_values(gw->moisture_content);

  /* Create or shared a field attached to cell for the hydraulic head */
  if (cdo_scheme == CS_SPACE_SCHEME_CDOFB) // For a future usage
    gw->hydraulic_head = cs_equation_get_field(richards);

  else { // Create a new field

    /* Create a hydraulic head field attached to cells (field on vertices
       are already defined thanks to the Richards equation in this case) */
    gw->hydraulic_head = cs_field_create("hydraulic_head_cell",
                                         field_mask,
                                         cell_loc_id,
                                         1,        // dimension
                                         has_previous);

    /* Allocate and initialize values */
    cs_field_allocate_values(gw->hydraulic_head);

  }

  if (gw->with_gravitation) { /* Gravitation effect */

    gw->pressure_head = cs_field_create("pressure_head",
                                        field_mask,
                                        cell_loc_id,
                                        1,
                                        has_previous);

    /* Allocate and initialize values */
    cs_field_allocate_values(gw->pressure_head);

  }

  /* Set the values for the permeability and the moisture content
     and if needed set also the value of the soil capacity */

  gw->global_model = gw->soil_param[0].model;

  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gw_soil_t  *soil = gw->soil_param + soil_id;
    const double  theta_s = soil->saturated_moisture;

    /* Is there a unique model ? */
    if (soil->model != gw->global_model)
      gw->global_model = CS_GROUNDWATER_MODEL_COMPOSITE;

    const char  *ml_name = cs_mesh_location_get_name(soil->ml_id);
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

    /* Initialization of soil_id and moisture content */
    if (gw->n_soils > 1) {

      assert(elt_ids != NULL); /* sanity check */

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

        cs_lnum_t  c_id = elt_ids[i];

        gw->moisture_content->val[c_id] = theta_s;
        gw->soil_id[c_id] = soil_id;

      } // Loop on selected cells

    }
    else { /* n_soils == 1 => all cells are selected */

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts[0]; i++)
        gw->moisture_content->val[i] = theta_s;

    }

    switch (soil->model) {

    case CS_GROUNDWATER_MODEL_GENUCHTEN:
      {
        cs_desc_t  desc = {.location = CS_FLAG_SCAL | cs_cdo_primal_cell,
                           .state = CS_FLAG_STATE_POTENTIAL};

        /* Set permeability tensor behavior */
        cs_property_def_by_law(permeability,
                               ml_name,
                               (const void *)soil,
                               _permeability_by_genuchten_law);

        if (gw->with_gravitation)
          cs_property_set_array(permeability, desc, gw->pressure_head->val);
        else
          cs_property_set_array(permeability, desc, gw->hydraulic_head->val);

        /* Soil capacity settings (related to unsteady term) */
        if (has_previous) {

          cs_property_t *capacity = cs_equation_get_time_property(richards);

          cs_property_def_by_law(cs_equation_get_time_property(richards),
                                 ml_name,
                                 (const void *)soil,
                                 _capacity_by_genuchten_law);

          if (gw->with_gravitation)
            cs_property_set_array(capacity, desc, gw->pressure_head->val);
          else
            cs_property_set_array(capacity, desc, gw->hydraulic_head->val);

        }

      }
      break;

    case CS_GROUNDWATER_MODEL_TRACY:
      {
        cs_desc_t  desc = {.location = CS_FLAG_SCAL | cs_cdo_primal_cell,
                           .state = CS_FLAG_STATE_POTENTIAL};

        /* Set permeability tensor behavior */
        cs_property_def_by_law(permeability,
                               ml_name,
                               (const void *)soil,
                               _permeability_by_tracy_law);

        if (gw->with_gravitation)
          cs_property_set_array(permeability, desc, gw->pressure_head->val);
        else
          cs_property_set_array(permeability, desc, gw->hydraulic_head->val);

        /* Soil capacity settings (related to unsteady term) */
        if (has_previous) {

          cs_property_t  *capacity = cs_equation_get_time_property(richards);

          const cs_gw_tracy_t  law = soil->tracy_param;
          const double  dh = law.h_s - law.h_r;
          const double  dtheta = soil->delta_moisture;

          cs_property_iso_def_by_value(capacity, ml_name, dtheta/dh);

        } /* Set the soil capacity */

      }
      break; // Tracy model

    case CS_GROUNDWATER_MODEL_SATURATED:

      switch (cs_property_get_type(permeability)) {

      case CS_PROPERTY_ISO:
        cs_property_iso_def_by_value(permeability,
                                     ml_name,
                                     soil->saturated_permeability.val);
        break;

      case CS_PROPERTY_ORTHO:
        cs_property_ortho_def_by_value(permeability,
                                       ml_name,
                                       soil->saturated_permeability.vect);
        break;

      case CS_PROPERTY_ANISO:
        cs_property_aniso_def_by_value(permeability,
                                       ml_name,
                  (const double (*)[3])soil->saturated_permeability.tens);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
        break;

      } // End of switch on permeability type

      if (has_previous)
        bft_error(__FILE__, __LINE__, 0,
                  " Saturated model must yield a steady Richards equation.");

      break; // Saturated model

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of model for groundwater module.."));
      break; // Nothing to do

    } /* Switch depending on the type of model */

  } // Loop on the different type of soils

  { /* Define and then link the advection field to each tracer equations */
    cs_desc_t  desc = {.location = CS_FLAG_SCAL | cs_cdo_dual_face_byc,
                       .state = CS_FLAG_STATE_FLUX};

    cs_advection_field_def_by_array(gw->adv_field, desc, gw->darcian_flux);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a reaction term for a given tracer
 *
 * \param[in] gw         pointer to a cs_groundwater_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_groundwater_tracer_needs_reaction(const cs_groundwater_t    *gw,
                                     int                        eq_id)
{
  int  tracer_id = _get_tracer_id(gw, eq_id);
  bool  is_needed = false;

  /* Loop on soils to check in a reaction term is needed */
  for (int soil_id = 0; soil_id < gw->n_soils && is_needed == false; soil_id++)
    {
      cs_gw_soil_t  *soil = gw->soil_param + soil_id;

      if (soil->tracer_param[tracer_id].reaction_rate > 0) is_needed = true;
    }

  return is_needed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a diffusion term for a given tracer
 *
 * \param[in] gw         pointer to a cs_groundwater_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_groundwater_tracer_needs_diffusion(const cs_groundwater_t    *gw,
                                      int                        eq_id)
{
  int  tracer_id = _get_tracer_id(gw, eq_id);
  bool  is_needed = false;

  /* Loop on soils to check in a reaction term is needed */
  for (int soil_id = 0; soil_id < gw->n_soils && is_needed == false; soil_id++)
    {
      cs_gw_soil_t  *soil = gw->soil_param + soil_id;

      if (soil->tracer_param[tracer_id].alpha_t > 0) is_needed = true;
      if (soil->tracer_param[tracer_id].alpha_l > 0) is_needed = true;
      if (soil->tracer_param[tracer_id].wmd > 0) is_needed = true;
    }

  return is_needed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for a tracer equation
 *
 * \param[in]      tracer_eq_id  id of the equation related to this tracer
 * \param[in, out] eq            pointer to the related cs_equation_t structure
 * \param[in, out] gw            pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_tracer_setup(int                  tracer_eq_id,
                            cs_equation_t       *eq,
                            cs_groundwater_t    *gw)
{
  const cs_flag_t  eq_flag = cs_equation_get_flag(eq);
  const cs_desc_t  desc1 = {.location = CS_FLAG_SCAL | cs_cdo_primal_cell,
                            .state = CS_FLAG_STATE_DENSITY};

  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Set time property */
  const int  tracer_id = _get_tracer_id(gw, tracer_eq_id);
  cs_property_t  *time_pty = cs_equation_get_time_property(eq);

  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    cs_gw_soil_t  *soil = gw->soil_param + soil_id;
    cs_gw_tracer_t  *tp = soil->tracer_param + tracer_id;

    cs_property_def_by_law(time_pty,
                           cs_mesh_location_get_name(soil->ml_id),
                           (const void *)tp,
                           _get_tracer_time_coeff);

  } // Loop on soils

  cs_property_set_array(time_pty, desc1, gw->moisture_content->val);

  /* Add a diffusion property */
  if (eq_flag & CS_EQUATION_DIFFUSION) {

    cs_property_t  *diff_pty = cs_equation_get_diffusion_property(eq);

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gw_soil_t  *soil = gw->soil_param + soil_id;
      cs_gw_tracer_t  *tp = soil->tracer_param + tracer_id;

      cs_property_def_by_scavec_law(diff_pty,
                                    cs_mesh_location_get_name(soil->ml_id),
                                    (const void *)tp,
                                    _get_tracer_diffusion_tensor);

    } // Loop on soils

    cs_property_set_array(diff_pty, desc1, gw->moisture_content->val);

    cs_desc_t  desc2 = {.location = CS_FLAG_SCAL | cs_cdo_dual_face_byc,
                        .state = CS_FLAG_STATE_FLUX};
    cs_property_set_second_array(diff_pty, desc2, gw->darcian_flux);

  } /* Diffusion term has to be set */

  /* Add a reaction property */
  if (eq_flag & CS_EQUATION_REACTION) {

    cs_property_t  *reac_pty = cs_equation_get_reaction_property(eq, "decay");

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gw_soil_t  *soil = gw->soil_param + soil_id;
      cs_gw_tracer_t  *tp = soil->tracer_param + tracer_id;

      cs_property_def_by_law(reac_pty,
                             cs_mesh_location_get_name(soil->ml_id),
                             (const void *)tp,
                             _get_tracer_reaction_coeff);

    } // Loop on soils

    cs_property_set_array(reac_pty, desc1, gw->moisture_content->val);

  } /* Reaction term has to be set */

  if (eq_flag & CS_EQUATION_DIFFUSION)
    cs_equation_set_param(eq, CS_EQKEY_ADV_SCHEME, "sg");

  /* TODO: add predefined source term */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      do_logcvg  output information on convergence or not
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
                       bool                          do_logcvg,
                       cs_equation_t                *eqs[],
                       cs_groundwater_t             *gw)
{
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
    cs_equation_init_system(mesh, richards);

    /* Update hydraulic/pressure head */
    _update_head(connect, cdoq, richards, gw);

    /* Build and solve the linear system related to the Richards equations */
    if (cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(richards, do_logcvg);

      /* Update hydraulic/pressure head */
      _update_head(connect, cdoq, richards, gw);

      /* Compute the darcian flux */
      cs_equation_compute_diff_flux(richards, gw->darcian_flux);

      /* Update the moisture content */
      _update_moisture_content(cdoq, richards, gw);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_eq_ids[i]];

      cs_equation_init_system(mesh, eq);

      if (cs_equation_is_steady(eq)) {

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(eq, do_logcvg);

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
      cs_equation_solve(richards, do_logcvg);

      /* Update hydraulic/pressure head */
      _update_head(connect, cdoq, richards, gw);

      /* Compute the darcian flux */
      cs_equation_compute_diff_flux(richards, gw->darcian_flux);

      /* Update the moisture content */
      _update_moisture_content(cdoq, richards, gw);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_eq_ids[i]];

      if (!cs_equation_is_steady(eq)) { // unsteady ?

        /* Define the algebraic system */
        if (cs_equation_needs_build(eq))
          cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(eq, do_logcvg);

      } /* Solve this equation which is steady */

    } /* Loop on tracer equations */

  } /* nt_cur > 0 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_groundwater_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_list    list of cells (1 to n)
 * \param[in]      i_face_list  list of interior faces (1 to n)
 * \param[in]      b_face_list  list of boundary faces (1 to n)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_groundwater_extra_post(void                      *input,
                          int                        mesh_id,
                          int                        cat_id,
                          int                        ent_flag[5],
                          cs_lnum_t                  n_cells,
                          cs_lnum_t                  n_i_faces,
                          cs_lnum_t                  n_b_faces,
                          const cs_lnum_t            cell_list[],
                          const cs_lnum_t            i_face_list[],
                          const cs_lnum_t            b_face_list[],
                          const cs_time_step_t      *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_list);
  CS_UNUSED(i_face_list);
  CS_UNUSED(b_face_list);

  if (input == NULL)
    return;

  if (mesh_id != -1) /* Post-processing only on the generic volume mesh */
    return;

  const cs_groundwater_t  *gw = (const cs_groundwater_t *)input;

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

  if (gw->with_gravitation) { /* Post-process pressure head */
    cs_field_t  *f = gw->pressure_head;
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
