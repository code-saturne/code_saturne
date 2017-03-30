/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_cdo.h"
#include "cs_field.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_DBG 0

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Set of parameters related to a tracer equation */
/* ---------------------------------------------- */

typedef struct {

  double     rho_kd;        // Bulk density times the distribution coefficient
  double     alpha_l;       // Longitudinal dispersivity
  double     alpha_t;       // Transversal dispersivity
  double     wmd;           // Water molecular diffusivity
  double     reaction_rate; /* First order decay coefficient (related to the
                               reaction term) */

} cs_gwf_tracer_t;

/* Set of parameters describing a type of soil */
/* ------------------------------------------- */
typedef struct {

  int                       ml_id;    /* id related to a mesh location structure
                                         The support entities are cells */

  /* Set of parameters for each tracer equation */
  cs_gwf_tracer_t          *tracer_param;

  /* Physical modelling adopted for this soil */
  cs_gwf_hydraulic_model_t  model;

  /* Parameters for predefined hydraulic models */
  cs_gwf_genuchten_t        genuchten_param; /* Van-Genuchten-Mualen law */
  cs_gwf_tracy_t            tracy_param;     /* Tracy law */

  /* Main soil properties */
  double                    bulk_density;
  double                    residual_moisture;       /* theta_r */
  double                    saturated_moisture;      /* theta_s */
  double                    delta_moisture;          /* theta_s - theta_r */
  cs_get_t                  saturated_permeability;  /* k_s [m.s^-1] */

} cs_gwf_soil_t;

/* Set of parameters related to the groundwater module */
/* --------------------------------------------------- */
struct _gwf_t {

  cs_gwf_hydraulic_model_t   global_model;

  /* Physical parameters related to each kind of soil considered.
     If n_soils > 1, soil_id array stores the id giving access to the soil
     parameters related to each cell of the mesh */
  int              n_soils;
  int              n_max_soils;
  cs_gwf_soil_t   *soil_param;

  cs_lnum_t        n_cells;   // Number of cells (useful to get soil_id)
  short int       *soil_id;   // NULL if only one soil or allocated to n_cells

  /* Gravity effect */
  bool             with_gravitation;
  cs_real_3_t      gravity;

  /* Related fields */
  cs_field_t      *moisture_content; /* Always located at cells */
  cs_field_t      *pressure_head;    /* Allocated only if gravitation is on.
                                        Location depends on the discretization
                                        scheme used to solve Richards eq.
                                        pressure head is denoted by h
                                        hydraulic head (solved in Richards eq.)
                                        is denoted by H.
                                        h = H - gravity_potential */

  cs_real_t       *head_in_law;       /* Array used as an input in laws */

  /* Set of equations associated to this module */
  int              richards_eq_id;
  int              n_tracers;
  int              n_max_tracers;
  int             *tracer_eq_ids;

  /* Permeability is the diffusion property related to Richards equation but
     this property plays also a role in the diffusion of tracer equations */
  cs_property_t   *permeability;  /* shared with a cs_domain_t structure */

  /* Settings related to the advection field stemming from the darcian flux */
  cs_flag_t        flux_location; /* indicate where the array is defined */
  cs_real_t       *darcian_flux;  /* array defining the advection field */
  cs_adv_field_t  *adv_field;     /* shared with a cs_domain_t structure */

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
 * \brief  Set a cs_gwf_tracer_t structure
 *
 * \param[in, out] tp              pointer to a cs_gwf_tracer_t structure
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      bulk_density    value of the bulk density
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_tracer_param(cs_gwf_tracer_t    *tp,
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
 * \brief  Get the id in a cs_gwf_t structure from the tracer
 *         equation id
 *
 * \param[in]   gw             pointer to a cs_gwf_t structure
 * \param[in]   tracer_eq_id   tracer equation id
 *
 * \returns an id related to this tracer equation id
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_tracer_id(const cs_gwf_t   *gw,
               int               tracer_eq_id)
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
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts         number of elements to treat
 * \param[in]      elt_ids        list of element ids (NULL if no indirection)
 * \param[in]      theta          values of the moisture content
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         array storing the result
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_tracer_time_coeff(cs_lnum_t          n_elts,
                       const cs_lnum_t    elt_ids[],
                       const cs_real_t    theta[],
                       const void        *tracer_struc,
                       cs_real_t         *result)
{
  const cs_gwf_tracer_t  *tracer = (const cs_gwf_tracer_t  *)tracer_struc;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  id = elt_ids[i];
      result[id] = theta[id] + tracer->rho_kd;
    } // Loop on cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++)
      result[i] = theta[i] + tracer->rho_kd;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in reaction term for the simulation
 *         of tracer equations
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts         number of elements to treat
 * \param[in]      elt_ids        list of element ids (NULL if no indirection)
 * \param[in]      theta          values of the moisture content
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         array storing the result
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_tracer_reaction_coeff(cs_lnum_t          n_elts,
                           const cs_lnum_t    elt_ids[],
                           const cs_real_t    theta[],
                           const void        *tracer_struc,
                           cs_real_t         *result)
{
  const cs_gwf_tracer_t  *tracer = (const cs_gwf_tracer_t  *)tracer_struc;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  id = elt_ids[i];
      result[id] = tracer->reaction_rate * (theta[id] + tracer->rho_kd);
    } // Loop on cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++)
      result[i] = tracer->reaction_rate * (theta[i] + tracer->rho_kd);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the diffusion term of the
 *         tracer equation
 *         This function fits the generic prototype of cs_twovar_law_func_t
 *
 * \param[in]      n_elts         number of elements to treat
 * \param[in]      elt_ids        list of element ids (NULL if no indirection)
 * \param[in]      theta          values of the moisture content
 * \param[in]      velocity       values of the local velocity
 * \param[in]      tracer_struc   pointer to a soil structure
 * \param[in, out] result         array storing the result
 */
/*----------------------------------------------------------------------------*/

static void
_get_tracer_diffusion_tensor(cs_lnum_t          n_elts,
                             const cs_lnum_t    elt_ids[],
                             const cs_real_t    theta[],
                             const cs_real_t    velocity[],
                             const void        *tracer_struc,
                             cs_real_t         *result)
{
  const cs_gwf_tracer_t  *tp = (const cs_gwf_tracer_t *)tracer_struc;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i]; // cell_id
      const cs_real_t  *v = velocity + 3*id;
      const double  v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
      const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
      const double  coef1 = tp->wmd * theta[id] + tp->alpha_t*vnorm;

      double  delta = 0.;
      if (vnorm > cs_math_zero_threshold)
        delta = (tp->alpha_l - tp->alpha_t)/vnorm;

      const double  dcv[3] = {delta*v[0], delta*v[1], delta*v[2]};

      cs_real_t  *_r = result + 9*id;
      for (int ki = 0; ki < 3; ki++) {

        /* Diagonal terms */
        _r[3*ki+ki] = coef1 + delta*v2[ki];

        /* Extra-diagonal terms */
        for (int kj = ki + 1; kj < 3; kj++) {
          _r[3*ki+kj] = dcv[ki]*v[kj];
          _r[3*kj+ki] = _r[3*ki+kj]; /* tensor is symmetric by construction */
        }
      }

    } // Loop on cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  *v = velocity + 3*i;
      const double  v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
      const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
      const double  coef1 = tp->wmd * theta[i] + tp->alpha_t*vnorm;

      double  delta = 0.;
      if (vnorm > cs_math_zero_threshold)
        delta = (tp->alpha_l - tp->alpha_t)/vnorm;

      const double  dcv[3] = {delta*v[0], delta*v[1], delta*v[2]};

      cs_real_t  *_r = result + 9*i;
      for (int ki = 0; ki < 3; ki++) {

        /* Diagonal terms */
        _r[3*ki+ki] = coef1 + delta*v2[ki];

        /* Extra-diagonal terms */
        for (int kj = ki + 1; kj < 3; kj++) {
          _r[3*ki+kj] = dcv[ki]*v[kj];
          _r[3*kj+ki] = _r[3*ki+kj]; /* tensor is symmetric by construction */
        }
      }

    } // Loop on all cells
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         van Genuchten-Mualen law
 *         One assumes that pressure head is located at cells
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts         number of elements to treat
 * \param[in]      elt_ids        list of element ids (NULL if no indirection)
 * \param[in]      h_vals         values of the pressure head
 * \param[in]      soil_struc     pointer to a soil structure
 * \param[in, out] result         array storing the result
 */
/*----------------------------------------------------------------------------*/

static void
_genuchten_permeability_from_c_head(cs_lnum_t          n_elts,
                                    const cs_lnum_t    elt_ids[],
                                    const cs_real_t    h_vals[],
                                    const void        *soil_struc,
                                    cs_real_t         *result)
{
  const cs_gwf_soil_t  *soil = (const cs_gwf_soil_t  *)soil_struc;
  const cs_gwf_genuchten_t  law = soil->genuchten_param;

  /* Up to now, only isotropic values are considered */
  const  double  iso_satval = soil->saturated_permeability.val;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i]; // cell_id
      const cs_real_t  h = h_vals[id];
      double  isoval = iso_satval;

      if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */
        const double  one_alpha_hn = 1 + pow(fabs(law.scale*h), law.n);
        const double  se = pow(one_alpha_hn, -law.m);
        const double  se_pow_L = pow(se, law.tortuosity);
        const double  se_pow_overm = pow(se, 1/law.m);
        const double  coef_base = 1 - pow(1 - se_pow_overm, law.m);

        isoval *= se_pow_L * coef_base*coef_base;
      }

      /* Build the permeability which is always defined as a tensor */
      cs_real_t  *_r = result + 9*id;
      for (int k = 0; k < 3; k++) {
        _r[3*k+k] = isoval;
        for (int l = k+1; l < 3; l++)
          _r[3*k+l] = _r[3*l+k] = 0;
      }

    } // Loop on selected cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  h = h_vals[i];
      double  isoval = iso_satval;

      if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */
        const double  one_alpha_hn = 1 + pow(fabs(law.scale*h), law.n);
        const double  se = pow(one_alpha_hn, -law.m);
        const double  se_pow_L = pow(se, law.tortuosity);
        const double  se_pow_overm = pow(se, 1/law.m);
        const double  coef_base = 1 - pow(1 - se_pow_overm, law.m);

        isoval *= se_pow_L * coef_base*coef_base;
      }

      /* Build the permeability which is always defined as a tensor */
      cs_real_t  *_r = result + 9*i;
      for (int k = 0; k < 3; k++) {
        _r[3*k+k] = isoval;
        for (int l = k+1; l < 3; l++)
          _r[3*k+l] = _r[3*l+k] = 0;
      }

    } // Loop on all cells
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the moisture content using the Van Genuchten law
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts        number of elements to treat
 * \param[in]      elt_ids       list of element ids (NULL if no indirection)
 * \param[in]      h_vals        values of the pressure head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        array storing the result
 */
/*----------------------------------------------------------------------------*/

static void
_genuchten_moisture_from_c_head(cs_lnum_t          n_elts,
                                const cs_lnum_t    elt_ids[],
                                const cs_real_t    h_vals[],
                                const void        *soil_struc,
                                cs_real_t          result[])
{
  const cs_gwf_soil_t  *soil = (const cs_gwf_soil_t  *)soil_struc;
  const cs_gwf_genuchten_t  law = soil->genuchten_param;

  if (elt_ids != NULL) {

# pragma omp parallel for if (n_elts > CS_THR_MIN) default(none) \
  shared(n_elts, elt_ids, h_vals, soil, result)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i];
      const cs_real_t  h = h_vals[id];
      double  Se = 1; // dimensionless moisture

      if (h < 0) {
        const double  coef = pow(fabs(law.scale * h), law.n);
        Se = pow(1 + coef, -law.m);
      }

      /* Computed moisture content */
      result[id] = Se*soil->delta_moisture + soil->residual_moisture;

    } // Loop on selected cells

  }
  else {

# pragma omp parallel for if (n_elts > CS_THR_MIN) default(none)     \
  shared(n_elts, elt_ids, h_vals, soil, result)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  h = h_vals[i];
      double  Se = 1; // dimensionless moisture

      if (h < 0) {
        const double  coef = pow(fabs(law.scale * h), law.n);
        Se = pow(1 + coef, -law.m);
      }

      /* Computed moisture content */
      result[i] = Se*soil->delta_moisture + soil->residual_moisture;

    } // Loop on all cells

  } // elt_ids == NULL

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the moisture content using the Van Genuchten law
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts         number of elements to treat
 * \param[in]      elt_ids        list of element ids (NULL if no indirection)
 * \param[in]      h_vals         values of the pressure head
 * \param[in]      soil_struc     pointer to a soil structure
 * \param[in, out] result         array storing the result
 */
/*----------------------------------------------------------------------------*/

static void
_genuchten_capacity_from_c_head(cs_lnum_t          n_elts,
                                const cs_lnum_t    elt_ids[],
                                const cs_real_t    h_vals[],
                                const void        *soil_struc,
                                cs_real_t          result[])
{
  const cs_gwf_soil_t  *soil = (const cs_gwf_soil_t  *)soil_struc;
  const cs_gwf_genuchten_t  law = soil->genuchten_param;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i];
      const cs_real_t  h = h_vals[id];

      if (h >= 0)
        result[id] = 0.; // default initialization for saturated soil
      else {

        const double  mult_coef = -law.n * law.m * soil->delta_moisture;
        const double  alpha_h_pow_n = pow(fabs(law.scale * h), law.n);
        const double  se_m1 = pow(1 + alpha_h_pow_n, -law.m-1);

        result[id] = mult_coef * alpha_h_pow_n/h * se_m1;

      }

    } // Loop on selected cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  h = h_vals[i];

      if (h >= 0)
        result[i] = 0.; // default initialization for saturated soil
      else {

        const double  mult_coef = -law.n * law.m * soil->delta_moisture;
        const double  alpha_h_pow_n = pow(fabs(law.scale * h), law.n);
        const double  se_m1 = pow(1 + alpha_h_pow_n, -law.m-1);

        result[i] = mult_coef * alpha_h_pow_n/h * se_m1;

      }

    } // Loop on all cells

  } // elt_ids == NULL

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the permeability (or hydraulic conductivity) using the
 *         Tracy law
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts        number of elements to treat
 * \param[in]      elt_ids       list of element ids (NULL if no indirection)
 * \param[in]      h_vals        values of the hydralic head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        array storing the result
 */
/*----------------------------------------------------------------------------*/

static inline void
_tracy_permeability_from_c_head(cs_lnum_t          n_elts,
                                const cs_lnum_t    elt_ids[],
                                const cs_real_t    h_vals[],
                                const void        *soil_struc,
                                cs_real_t          result[])
{
  const cs_gwf_soil_t  *soil = (const cs_gwf_soil_t  *)soil_struc;
  const cs_gwf_tracy_t  law = soil->tracy_param;

  /* Up to now, only isotropic values are considered */
  const double  ks = soil->saturated_permeability.val;

  if (elt_ids != NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i];
      const cs_real_t  h = h_vals[id];
      const double  isoval = ks * (h - law.h_r)/(law.h_s - law.h_r);

      /* Permeability is always defined as a tensor */
      cs_real_t  *_r = result + 9*id;
      for (int k = 0; k < 3; k++) {
        _r[3*k+k] = isoval;
        for (int l = k+1; l < 3; l++)
          _r[3*k+l] = _r[3*l+k] = 0;
      }

    } // Loop on selected cells

  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  h = h_vals[i];
      const double  isoval = ks * (h - law.h_r)/(law.h_s - law.h_r);

      /* Permeability is always defined as a tensor */
      cs_real_t  *_r = result + 9*i;
      for (int k = 0; k < 3; k++) {
        _r[3*k+k] = isoval;
        for (int l = k+1; l < 3; l++)
          _r[3*k+l] = _r[3*l+k] = 0;
      }

    } // Loop on all cells

  } // elt_ids == NULL

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the moisture content using the Tracy law
 *         This function fits the generic prototype of cs_onevar_law_func_t
 *
 * \param[in]      n_elts        number of elements to treat
 * \param[in]      elt_ids       list of element ids (NULL if no indirection)
 * \param[in]      h_vals        values of the pressure head
 * \param[in]      soil_struc    pointer to a soil structure
 * \param[in, out] result        array storing the result
 */
/*----------------------------------------------------------------------------*/

static void
_tracy_moisture_from_c_head(cs_lnum_t         n_elts,
                            const cs_lnum_t   elt_ids[],
                            const cs_real_t   h_vals[],
                            const void       *soil_struc,
                            cs_real_t         result[])
{
  const cs_gwf_soil_t  *soil = (const cs_gwf_soil_t  *)soil_struc;
  const cs_gwf_tracy_t  law = soil->tracy_param;

  if (elt_ids != NULL) {

# pragma omp parallel for if (n_elts > CS_THR_MIN) default(none) \
  shared(n_elts, elt_ids, h_vals, soil, result)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  id = elt_ids[i];
      const cs_real_t  h = h_vals[id];
      const double  k_r = (h - law.h_r)/(law.h_s - law.h_r);

      result[id] = k_r * soil->delta_moisture + soil->residual_moisture;

    } // Loop on selected cells

  }
  else {

# pragma omp parallel for if (n_elts > CS_THR_MIN) default(none) \
  shared(n_elts, h_vals, soil, result)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_real_t  h = h_vals[i];
      const double  k_r = (h - law.h_r)/(law.h_s - law.h_r);

      result[i] = k_r * soil->delta_moisture + soil->residual_moisture;

    } // Loop on all cells

  } // elt_ids == NULL

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the groundwater system (pressure head, head in law, moisture
 *         content, darcian velocity)
 *
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      richards    pointer to the Richards equation structure
 * \param[in]      cur2prev    true or false
 * \param[in, out] gw          pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_system(const cs_cdo_quantities_t   *cdoq,
               const cs_cdo_connect_t      *connect,
               const cs_equation_t         *richards,
               bool                         cur2prev,
               cs_gwf_t                    *gw)
{
  /* Sanity checks */
  if (richards == NULL || gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Groundwater module or Richards eq. is not allocated.");

  const cs_field_t  *hydraulic_head = cs_equation_get_field(richards);
  cs_field_t  *pressure_head = gw->pressure_head;

  if (gw->with_gravitation) { /* Update the pressure head */

    /* Sanity checks */
    if (pressure_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */
    if (cur2prev)
      cs_field_current_to_previous(gw->pressure_head);

    switch (cs_equation_get_space_scheme(richards)) {

    case CS_SPACE_SCHEME_CDOVB:
      assert(hydraulic_head->location_id ==
             cs_mesh_location_get_id_by_name("vertices"));

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {

        const cs_real_t  gpot = cs_math_3_dot_product(cdoq->vtx_coord + 3*i,
                                                      gw->gravity);

        pressure_head->val[i] = hydraulic_head->val[i] - gpot;

      }

      /* Update head_in_law */
      cs_reco_pv_at_cell_centers(connect->c2v,
                                 cdoq,
                                 pressure_head->val,
                                 gw->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("vertices"));

#       pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {

          const cs_real_t  gpot = cs_math_3_dot_product(cdoq->vtx_coord + 3*i,
                                                        gw->gravity);

          pressure_head->val[i] = hydraulic_head->val[i] - gpot;

        }

        /* Update head_in_law */
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards);

#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {

          const cs_real_t  gpot =
            cs_math_3_dot_product(cdoq->cell_centers + 3*i, gw->gravity);

          gw->head_in_law[i] = hydraulic_head_cells[i] - gpot;

        }

      }
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO:
      assert(hydraulic_head->location_id ==
             cs_mesh_location_get_id_by_name("cells"));

#     pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {

        const cs_real_t  gpot = cs_math_3_dot_product(cdoq->cell_centers + 3*i,
                                                      gw->gravity);

        pressure_head->val[i] = hydraulic_head->val[i] - gpot;

      }
      break; // Nothing to do (h_head is a pointer to richards field)

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

    } // Switch on space scheme

  }
  else { // No gravity effect id taken into account

    /* Update head_in_law */
    switch(cs_equation_get_space_scheme(richards)) {

    case CS_SPACE_SCHEME_CDOVB:
      cs_reco_pv_at_cell_centers(connect->c2v,
                                 cdoq,
                                 hydraulic_head->val,
                                 gw->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards);

        memcpy(gw->head_in_law, hydraulic_head_cells,
               sizeof(cs_real_t)*cdoq->n_cells);
      }
      break;

    default:
      break; // Nothing to do

    } // Switch on the space scheme related to the Richards equation

  } /* Gravity is activated or not */

  /* Update the moisture content */
  cs_field_t  *moisture = gw->moisture_content;

  if (moisture == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " The field related to the moisture content is not allocated.");

  if (cur2prev)
    cs_field_current_to_previous(moisture);

  /* Moisture content is define in each cell along with the hydraulic (or
     pressure head */
  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);

    switch (soil->model) {

    case CS_GWF_HYDRAULIC_GENUCHTEN:
      _genuchten_moisture_from_c_head(n_elts[0], elt_ids, gw->head_in_law,
                                      (const void *)soil,
                                      moisture->val);
      break;

    case CS_GWF_HYDRAULIC_TRACY:
      _tracy_moisture_from_c_head(n_elts[0], elt_ids, gw->head_in_law,
                                  (const void *)soil,
                                  moisture->val);
      break;

    case CS_GWF_HYDRAULIC_SATURATED:
      if (elt_ids == NULL) {

#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id <  cdoq->n_cells; c_id++)
          moisture->val[c_id] = soil->saturated_moisture;

      }
      else {

#       pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  cs_dump_array_to_listing("MOISTURE_CONTENT",
                           cdoq->n_cells,
                           gw->moisture_content->val, 8);
#endif

  /* Update the advection field related to the groundwater flow module */
  cs_field_t  *vel = cs_advection_field_get_field(gw->adv_field,
                                                  CS_MESH_LOCATION_CELLS);

  if (cur2prev)
    cs_field_current_to_previous(vel);

  /* Compute the darcian flux and the darcian velocity inside each cell */
  switch (cs_equation_get_space_scheme(richards)) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:

    /* Update the array gw->darcian_flux associated to the advection field */
    cs_equation_compute_diff_flux_cellwise(richards,
                                           gw->flux_location,
                                           gw->darcian_flux);

    /* Set the new values */
    cs_advection_field_at_cells(gw->adv_field, vel->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
    if (cs_test_flag(location, cs_cdo_dual_face_byc))
      cs_dump_array_to_listing("DARCIAN_FLUX_DFbyC",
                               connect->c2e->idx[cdoq->n_cells],
                               gw->darcian_flux, 8);
    else if (cs_test_flag(location, cs_cdo_primal_cell))
      cs_dump_array_to_listing("DARCIAN_FLUX_CELL",
                               3*cdoq->n_cells,
                               gw->darcian_flux, 3);
#endif
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    bft_error(__FILE__, __LINE__, 0, " TODO.");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

  } // End of switch

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_gwf_soil_t structure (already allocated)
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      ml_name    name of the mesh location
 * \param[in]      model      hydraulic modelling for this soil
 * \param[in]      tetha_s    saturated moisture
 * \param[in]      rho        bulk density
 *
 * \return a pointer to an initialized cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_soil_t *
_init_soil(cs_gwf_t                   *gw,
           const char                 *ml_name,
           cs_gwf_hydraulic_model_t    model,
           double                      theta_s,
           double                      rho)
{
  if (gw == NULL)
    return  NULL;

  int  soil_id = gw->n_soils;

  gw->n_soils += 1;
  if (gw->n_soils > gw->n_max_soils)
    bft_error(__FILE__, __LINE__, 0,
              " Maximum number of soils is reached.\n"
              " Stop adding a new soil by value (mesh location: %s).\n"
              " Please modify your settings.", ml_name);

  cs_gwf_soil_t  *soil = gw->soil_param + soil_id;

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
  soil->model = model;
  soil->saturated_moisture = theta_s;
  soil->bulk_density = rho;
  soil->residual_moisture = 0.15; // Default value

  /* Additionnal default initialization */
  switch (model) {

  case CS_GWF_HYDRAULIC_SATURATED: // Nothing to do
    break;

  case CS_GWF_HYDRAULIC_GENUCHTEN:
    {
      const double  n = 1.56;
      soil->genuchten_param.n = n;
      soil->genuchten_param.m = 1 - 1/n;
      soil->genuchten_param.scale = 0.036;
      soil->genuchten_param.tortuosity = 0.5;
    }
    break;

  case CS_GWF_HYDRAULIC_TRACY:
    soil->tracy_param.h_r = -100;
    soil->tracy_param.h_s = 0;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible model for a soil in the groundwater module.\n"
              " Availaible models: saturated, genuchten, tracy");

  } // Switch on soil modeling

  soil->delta_moisture = soil->saturated_moisture - soil->residual_moisture;

  /* Set of parameters for each tracer which are related to this soil */
  BFT_MALLOC(soil->tracer_param, gw->n_max_tracers, cs_gwf_tracer_t);

  for (int i = 0; i < gw->n_max_tracers; i++) /* default initialization */
    _set_tracer_param(soil->tracer_param + i,
                      0.0,  /* water molecular diffusivity */
                      0.0,  /* alpha_l */
                      0.0,  /* alpha_t */
                      1.0,  /* bulk density (not useful here since Kd = 0 */
                      0.0,  /* Kd (distribution coef.) */
                      0.0); /* reaction rate */

  return soil;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \return a pointer to a new allocated cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_create(void)
{
  cs_gwf_t  *gw = NULL;

  BFT_MALLOC(gw, 1, cs_gwf_t);

  /* Default initialization */
  gw->n_soils = 0;
  gw->n_max_soils = 0;
  gw->soil_param = NULL;
  gw->soil_id = NULL;

  gw->global_model = CS_GWF_N_HYDRAULIC_MODELS;

  gw->with_gravitation = false;
  gw->gravity[0] = 0, gw->gravity[1] = 0, gw->gravity[2] = 0;
  gw->pressure_head = NULL;
  gw->head_in_law = NULL;

  gw->richards_eq_id = -1;
  gw->n_tracers = 0;
  gw->n_max_tracers = 0;
  gw->tracer_eq_ids = NULL;

  gw->flux_location = cs_cdo_primal_cell;
  gw->darcian_flux = NULL;
  gw->adv_field = NULL;

  return gw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to groundwater flows
 *
 * \param[in, out]  gw     pointer to a cs_gwf_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_finalize(cs_gwf_t   *gw)
{
  if (gw == NULL)
    return NULL;

  BFT_FREE(gw->tracer_eq_ids);
  BFT_FREE(gw->darcian_flux);
  if (gw->head_in_law != NULL)
    BFT_FREE(gw->head_in_law);

  for (int i = 0; i < gw->n_soils; i++) {
    cs_gwf_soil_t *soil = gw->soil_param + i;
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
 * \param[in]  gw        pointer to a cs_gwf_t structure
 *
 * \return the number of requested soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(const cs_gwf_t    *gw)
{
  if (gw == NULL)
    return 0;

  return gw->n_max_soils;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the gravity and set the gravitaty vector
 *
 * \param[in, out]  gw        pointer to a cs_gwf_t structure
 * \param[in]       gvec      values of the gravity vector
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_gravity_vector(cs_gwf_t              *gw,
                          const cs_real_3_t      gvec)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gw->with_gravitation = true;
  gw->gravity[0] = gvec[0];
  gw->gravity[1] = gvec[1];
  gw->gravity[2] = gvec[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Advanced setting: indicate where the darcian flux is stored
 *         cs_cdo_primal_cell is the default setting
 *         cs_cdo_dual_face_byc is a valid choice for vertex-based schemes
 *
 * \param[in, out]  gw              pointer to a cs_gwf_t structure
 * \param[in]       location_flag   where the flux is defined
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_darcian_flux_location(cs_gwf_t      *gw,
                                 cs_flag_t      location_flag)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gw->flux_location = location_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_gwf_t structure
 *
 * \param[in]  gw     pointer to a cs_gwf_t struct. to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_summary(const cs_gwf_t   *gw)
{
  if (gw == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the groundwater module\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  if (gw->with_gravitation)
    cs_log_printf(CS_LOG_SETUP,
                  "  <GW/Gravitation> true -- Axis = [%.2f %.2f %.2f]\n",
                  gw->gravity[0], gw->gravity[1], gw->gravity[2]);
  else
    cs_log_printf(CS_LOG_SETUP, "  <GW/Gravitation> false\n");

  cs_log_printf(CS_LOG_SETUP,
                "  <GW/Tracer> n_tracer_equations %d\n", gw->n_tracers);
  cs_log_printf(CS_LOG_SETUP, "  <GW/Soils>  n_soils %d\n", gw->n_soils);

  const cs_property_t  *permeability = gw->permeability;

  for (int i = 0; i < gw->n_soils; i++) {

    const cs_gwf_soil_t  soil = gw->soil_param[i];
    const char *ml_name = cs_mesh_location_get_name(soil.ml_id);
    const cs_get_t  sat_perm = soil.saturated_permeability;

    cs_log_printf(CS_LOG_SETUP, "  <GW/Soil %s>", ml_name);
    cs_log_printf(CS_LOG_SETUP,
                  " residual_moisture %5.3e", soil.residual_moisture);
    cs_log_printf(CS_LOG_SETUP,
                  " saturated_moisture %5.3e\n", soil.saturated_moisture);
    cs_log_printf(CS_LOG_SETUP, "  <GW/Soil %s>", ml_name);

    switch (cs_property_get_type(permeability)) {

    case CS_PROPERTY_ISO:
      cs_log_printf(CS_LOG_SETUP,
                    " saturated_permeability (iso) %5.3e\n", sat_perm.val);
      break;

    case CS_PROPERTY_ORTHO:
      cs_log_printf(CS_LOG_SETUP,
                    " saturated_permeability (ortho) %5.3e %5.3e %5.3e\n",
                    sat_perm.vect[0], sat_perm.vect[1], sat_perm.vect[2]);
      break;

    case CS_PROPERTY_ANISO:
      cs_log_printf(CS_LOG_SETUP,
                    " saturated_permeability (aniso) %-5.3e %5.3e %5.3e\n"
                    "                                %-5.3e %5.3e %5.3e\n"
                    "                                %-5.3e %5.3e %5.3e\n",
                    sat_perm.tens[0][0], sat_perm.tens[0][1],
                    sat_perm.tens[0][2], sat_perm.tens[1][0],
                    sat_perm.tens[1][1], sat_perm.tens[1][2],
                    sat_perm.tens[2][0], sat_perm.tens[2][1],
                    sat_perm.tens[2][2]);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of property for %s."),
                cs_property_get_name(permeability));
      break;

    } // Switch on property type

    cs_log_printf(CS_LOG_SETUP, "  <GW/Soil %s> Hydraulic model", ml_name);
    switch (soil.model) {
    case CS_GWF_HYDRAULIC_GENUCHTEN:
      cs_log_printf(CS_LOG_SETUP, " VanGenuchten-Mualen\n");
      break;
    case CS_GWF_HYDRAULIC_SATURATED:
      cs_log_printf(CS_LOG_SETUP, " saturated\n");
      break;
    case CS_GWF_HYDRAULIC_TRACY:
      cs_log_printf(CS_LOG_SETUP, " Tracy\n");
      break;
    case CS_GWF_HYDRAULIC_USER:
      cs_log_printf(CS_LOG_SETUP, " User-defined\n");
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid model for a soil in the groundwater module.\n"
                " Please check your settings.");
    } // Switch model

  } // Loop on soils

  if (gw->n_soils > 1) {
    const char *meta = "  <GW/Global Hydraulic Model>";
    switch (gw->global_model) {
    case CS_GWF_HYDRAULIC_COMPOSITE:
    cs_log_printf(CS_LOG_SETUP, "%s composite model\n", meta);
      break;
    case CS_GWF_HYDRAULIC_GENUCHTEN:
      cs_log_printf(CS_LOG_SETUP, "%s VanGenuchten-Mualen\n", meta);
      break;
    case CS_GWF_HYDRAULIC_SATURATED:
      cs_log_printf(CS_LOG_SETUP, "%s saturated\n", meta);
      break;
    case CS_GWF_HYDRAULIC_TRACY:
      cs_log_printf(CS_LOG_SETUP, "%s Tracy\n", meta);
      break;
    case CS_GWF_HYDRAULIC_USER:
      cs_log_printf(CS_LOG_SETUP, "%s User-defined\n", meta);
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
 * \param[in]      richards_eq_id   id related to the Richards equation
 * \param[in]      n_soils          number of soils to consider
 * \param[in]      n_tracers        number of tracers to consider
 * \param[in, out] permeability     pointer to a property structure
 * \param[in, out] soil_capacity    pointer to a property structure
 * \param[in, out] adv_field        pointer to a cs_adv_field_t structure
 * \param[in, out] gw               pointer to a cs_gwf_t structure
 *
 * \return a pointer to a new allocated equation structure (Richards eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_initialize(int                        richards_eq_id,
                  int                        n_soils,
                  int                        n_tracer_eqs,
                  cs_property_t             *permeability,
                  cs_property_t             *soil_capacity,
                  cs_adv_field_t            *adv_field,
                  cs_gwf_t                  *gw)
{
  cs_equation_t  *eq = NULL;

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
  gw->darcian_flux = NULL;

  /* Quantities related to soils */
  gw->n_soils = 0;           /* No soil is set at the beginning */
  gw->n_max_soils = n_soils; /* Max. number of soils allocated */

  BFT_MALLOC(gw->soil_param, n_soils, cs_gwf_soil_t);

  gw->n_tracers = 0;
  gw->n_max_tracers = n_tracer_eqs;
  BFT_MALLOC(gw->tracer_eq_ids, n_tracer_eqs, int);
  for (int i = 0; i < n_tracer_eqs; i++)
    gw->tracer_eq_ids[i] = -1; /* Default initialization = not set */

  /* Add default post-processing related to groundwater flow module */
  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post, gw);

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an isotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      k_s        value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_iso_soil_by_value(cs_gwf_t                   *gw,
                             cs_gwf_hydraulic_model_t    model,
                             const char                 *ml_name,
                             double                      k_s,
                             double                      theta_s,
                             double                      rho)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_gwf_soil_t  *soil = _init_soil(gw, ml_name, model, theta_s, rho);

  soil->saturated_permeability.val = k_s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an orthotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      ks         value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_ortho_soil_by_value(cs_gwf_t                  *gw,
                               cs_gwf_hydraulic_model_t   model,
                               const char                *ml_name,
                               cs_real_t                 *ks,
                               double                     theta_s,
                               double                     rho)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_gwf_soil_t  *soil = _init_soil(gw, ml_name, model, theta_s, rho);

  for (int k = 0; k < 3; k++)
    soil->saturated_permeability.vect[k] = ks[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new soil attached to an orthotropic permeability
 *
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      model      type of modeling for the hydraulic behavior
 * \param[in]      ml_name    name of the mesh location related to this soil
 * \param[in]      k_s        value of the saturated permeability
 * \param[in]      theta_s    saturated moisture
 * \param[in]      rho        bulk density
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_aniso_soil_by_value(cs_gwf_t                  *gw,
                               cs_gwf_hydraulic_model_t   model,
                               const char                *ml_name,
                               cs_real_t                 *ks,
                               double                     theta_s,
                               double                     rho)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_gwf_soil_t  *soil = _init_soil(gw, ml_name, model, theta_s, rho);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      soil->saturated_permeability.tens[i][j] = ks[3*i+j];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters related to a cs_gwf_t structure
 *
 * \param[in, out]  gw        pointer to a cs_gwf_t structure
 * \param[in]       ml_name   name of the mesh location associated to this soil
 * \param[in]       key       key related to a member of the soil to set
 * \param[in]       val       value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_soil_param(cs_gwf_t          *gw,
                      const char        *ml_name,
                      cs_gwf_soilkey_t   key,
                      const cs_real_t    val)
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

  /* Set different available keys */
  switch(key) {

  case CS_SOILKEY_SAT_MOISTURE:
    {
      const double  theta_s = val;

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
      const double  theta_r = val;

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
      const double  h_s = val;

      if (ml_id == -1) {
        for (int i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GWF_HYDRAULIC_TRACY)
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
      const double  h_r = val;

      if (ml_id == -1) {

        for (int i = 0; i < gw->n_soils; i++)
          if (gw->soil_param[i].model == CS_GWF_HYDRAULIC_TRACY)
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
 * \param[in, out] gw              pointer to a cs_gwf_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      eqname          name of the equation
 * \param[in]      varname         name of the related variable
 *
 * \return a pointer to a new allocated equation structure (Tracer eq.)
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_add_tracer(cs_gwf_t         *gw,
                  int               tracer_eq_id,
                  const char       *eqname,
                  const char       *varname)
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
 * \param[in, out] gw              pointer to a cs_gwf_t structure
 * \param[in]      tracer_eq_id    id related to the tracer equation
 * \param[in]      ml_name         name of the related mesh location
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_tracer_param(cs_gwf_t        *gw,
                        int              tracer_eq_id,
                        const char      *ml_name,
                        double           wmd,
                        double           alpha_l,
                        double           alpha_t,
                        double           distrib_coef,
                        double           reaction_rate)
{
  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  tracer_id = _get_tracer_id(gw, tracer_eq_id);

  /* Look for the related soil */
  if (ml_name == NULL) { /* All soils have to be set for this tracer */

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;

      /* Set tracer parameters */
      _set_tracer_param(soil->tracer_param + tracer_id,
                        wmd,
                        alpha_l, alpha_t,
                        soil->bulk_density,
                        distrib_coef,
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

    cs_gwf_soil_t  *soil = gw->soil_param + soil_id;

    /* Set tracer parameters */
    _set_tracer_param(soil->tracer_param + tracer_id,
                      wmd,
                      alpha_l, alpha_t,
                      soil->bulk_density,
                      distrib_coef,
                      reaction_rate);

  } /* Set a specific couple (tracer,soil) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the Richards equation
 *
 * \param[in, out] gw        pointer to a cs_gwf_t structure
 * \param[in, out] richards  pointer to the related cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_richards_setup(cs_gwf_t            *gw,
                      cs_equation_t       *richards)
{
  /* Sanity checks */
  assert(richards != NULL);
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->n_soils == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is activated but no soil is defined."));

  bool has_previous = cs_equation_is_steady(richards) ? false:true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");

  const cs_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(richards);

  /* Create a moisture field attached to cells */
  gw->moisture_content = cs_field_create("moisture_content",
                                         field_mask,
                                         c_loc_id,
                                         1,        // dimension
                                         has_previous);

  cs_field_set_key_int(gw->moisture_content, cs_field_key_id("log"), 1);
  cs_field_set_key_int(gw->moisture_content, cs_field_key_id("post_vis"), 1);

  switch (space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    {
      if (gw->with_gravitation) /* Add the pressure_head field */
        gw->pressure_head = cs_field_create("pressure_head",
                                            field_mask,
                                            v_loc_id,
                                            1,
                                            has_previous);

      /* Define and then link the advection field to each tracer equations */
      cs_desc_t  flux_desc =
        {.location = CS_FLAG_SCALAR | gw->flux_location,
         .state = CS_FLAG_STATE_FLUX};

      cs_advection_field_def_by_array(gw->adv_field, flux_desc);

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    if (gw->with_gravitation) /* Add the pressure_head field */
      gw->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          has_previous);

    /* TODO: Set the darcian flux */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");
  }

  if (gw->with_gravitation) { /* Set option to the pressure_head field */
    cs_field_set_key_int(gw->pressure_head, cs_field_key_id("log"), 1);
    cs_field_set_key_int(gw->pressure_head, cs_field_key_id("post_vis"), 1);
  }

  /* Set the values for the permeability and the moisture content
     and if needed set also the value of the soil capacity */
  cs_property_t  *permeability = gw->permeability;

  gw->global_model = gw->soil_param[0].model;

  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
    const char  *ml_name = cs_mesh_location_get_name(soil->ml_id);

    /* Is there a unique model ? */
    if (soil->model != gw->global_model)
      gw->global_model = CS_GWF_HYDRAULIC_COMPOSITE;

    switch (soil->model) {

    case CS_GWF_HYDRAULIC_GENUCHTEN:

      /* Set permeability tensor behavior */
      cs_property_def_by_onevar_law(permeability,
                                    ml_name,
                                    (const void *)soil,
                                    _genuchten_permeability_from_c_head);

      /* Soil capacity settings (related to unsteady term) */
      if (has_previous)
        cs_property_def_by_onevar_law(cs_equation_get_time_property(richards),
                                      ml_name,
                                      (const void *)soil,
                                      _genuchten_capacity_from_c_head);
      break;

    case CS_GWF_HYDRAULIC_TRACY:

      /* Set permeability tensor behavior */
      cs_property_def_by_onevar_law(permeability,
                                    ml_name,
                                    (const void *)soil,
                                    _tracy_permeability_from_c_head);

      /* Soil capacity settings (related to unsteady term) */
      if (has_previous) {

        cs_property_t  *capacity = cs_equation_get_time_property(richards);

        const cs_gwf_tracy_t  law = soil->tracy_param;
        const double  dh = law.h_s - law.h_r;
        const double  dtheta = soil->delta_moisture;

        cs_property_iso_def_by_value(capacity, ml_name, dtheta/dh);

      }
      break; // Tracy model

    case CS_GWF_HYDRAULIC_SATURATED:

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

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a reaction term for a given tracer
 *
 * \param[in] gw         pointer to a cs_gwf_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_tracer_needs_reaction(const cs_gwf_t    *gw,
                             int                eq_id)
{
  int  tracer_id = _get_tracer_id(gw, eq_id);
  bool  is_needed = false;

  /* Loop on soils to check in a reaction term is needed */
  for (int soil_id = 0; soil_id < gw->n_soils && is_needed == false; soil_id++)
    {
      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;

      if (soil->tracer_param[tracer_id].reaction_rate > 0) is_needed = true;
    }

  return is_needed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to add a diffusion term for a given tracer
 *
 * \param[in] gw         pointer to a cs_gwf_t structure
 * \param[in] eq_id      id of the equation related to this tracer
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_tracer_needs_diffusion(const cs_gwf_t      *gw,
                              int                  eq_id)
{
  int  tracer_id = _get_tracer_id(gw, eq_id);
  bool  is_needed = false;

  /* Loop on soils to check in a reaction term is needed */
  for (int soil_id = 0; soil_id < gw->n_soils && is_needed == false; soil_id++)
    {
      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;

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
 * \param[in, out] gw            pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_setup(int                tracer_eq_id,
                    cs_equation_t     *eq,
                    cs_gwf_t          *gw)
{
  const cs_flag_t  eq_flag = cs_equation_get_flag(eq);

  /* Sanity check */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Set time property */
  const int  tracer_id = _get_tracer_id(gw, tracer_eq_id);

  if (eq_flag & CS_EQUATION_UNSTEADY) {

    cs_property_t  *time_pty = cs_equation_get_time_property(eq);

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
      cs_gwf_tracer_t  *tp = soil->tracer_param + tracer_id;

      cs_property_def_by_onevar_law(time_pty,
                                    cs_mesh_location_get_name(soil->ml_id),
                                    (const void *)tp,
                                    _get_tracer_time_coeff);

    } // Loop on soils

  }

  /* Add a diffusion property */
  if (eq_flag & CS_EQUATION_DIFFUSION) {

    cs_property_t  *diff_pty = cs_equation_get_diffusion_property(eq);

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
      cs_gwf_tracer_t  *tp = soil->tracer_param + tracer_id;

      cs_property_def_by_twovar_law(diff_pty,
                                    cs_mesh_location_get_name(soil->ml_id),
                                    (const void *)tp,
                                    _get_tracer_diffusion_tensor);

    } // Loop on soils

  } /* Diffusion term has to be set */

  /* Add a reaction property */
  if (eq_flag & CS_EQUATION_REACTION) {

    cs_property_t  *reac_pty = cs_equation_get_reaction_property(eq, "decay");

    for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

      cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
      cs_gwf_tracer_t  *tp = soil->tracer_param + tracer_id;

      cs_property_def_by_onevar_law(reac_pty,
                                    cs_mesh_location_get_name(soil->ml_id),
                                    (const void *)tp,
                                    _get_tracer_reaction_coeff);

    } // Loop on soils

  } /* Reaction term has to be set */

  if (eq_flag & CS_EQUATION_DIFFUSION)
    cs_equation_set_param(eq, CS_EQKEY_ADV_SCHEME, "sg");

  /* TODO: add predefined source term */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last initialization step of the groundwater flow module
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      n_equations  number of equations in the list
 * \param[in, out] equations    pointer to a list of cs_equation_t structures
 * \param[in, out] gw           pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_final_initialization(const cs_cdo_connect_t    *connect,
                            int                        n_equations,
                            cs_equation_t            **equations,
                            cs_gwf_t                  *gw)
{
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Check if settings are correct */
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

  const cs_equation_t  *richards = equations[gw->richards_eq_id];
  const cs_field_t  *hydraulic_head = cs_equation_get_field(richards);
  const cs_space_scheme_t  ric_scheme = cs_equation_get_space_scheme(richards);

  if (ric_scheme == CS_SPACE_SCHEME_CDOFB ||
      ric_scheme == CS_SPACE_SCHEME_HHO)
    bft_error(__FILE__, __LINE__, 0,
              _(" Richards eq. is only available for vertex-based schemes."));

  const cs_lnum_t  n_cells = connect->n_cells;

  /* Up to now Richards equation is only set with vertex-based schemes
     TODO: Face-based schemes */
  switch (ric_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_connect_index_t  *c2e = connect->c2e;

      BFT_MALLOC(gw->head_in_law, n_cells, cs_real_t);

      /* Darcian flux settings */
      BFT_MALLOC(gw->darcian_flux, c2e->idx[n_cells], cs_real_t);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < c2e->idx[n_cells]; i++)
        gw->darcian_flux[i] = 0;

      cs_advection_field_set_array(gw->adv_field, gw->darcian_flux);
    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:

    BFT_MALLOC(gw->head_in_law, n_cells, cs_real_t);

    /* Darcian flux settings */
    BFT_MALLOC(gw->darcian_flux, 3*n_cells, cs_real_t);

#   pragma omp parallel for if (3*n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < 3*n_cells; i++)
      gw->darcian_flux[i] = 0;

    cs_advection_field_set_array(gw->adv_field, gw->darcian_flux);
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    if (gw->with_gravitation)
      gw->head_in_law = gw->pressure_head->val;
    else
      gw->head_in_law = hydraulic_head->val;

    /* Darcian flux settings */
    BFT_MALLOC(gw->darcian_flux, 3*n_cells, cs_real_t);

#   pragma omp parallel for if (3*n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < 3*n_cells; i++)
      gw->darcian_flux[i] = 0;
    cs_advection_field_set_array(gw->adv_field, gw->darcian_flux);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");
    break;
  }

  if (gw->n_soils > 1) { /* Default initialization of soil_id */

    BFT_MALLOC(gw->soil_id, n_cells, short int);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      gw->soil_id[i] = gw->n_max_soils; /* Default value => not set */

  }

  bool has_previous = cs_equation_is_steady(richards) ? false:true;
  cs_desc_t  head_desc = {.location = CS_FLAG_SCALAR | cs_cdo_primal_cell,
                          .state = CS_FLAG_STATE_POTENTIAL};

  /* Set the values for the permeability and the moisture content
     and if needed set also the value of the soil capacity */
  for (int soil_id = 0; soil_id < gw->n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = gw->soil_param + soil_id;
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(soil->ml_id);
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(soil->ml_id);

    const double  theta_s = soil->saturated_moisture;

    if (gw->n_soils == 1) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        gw->moisture_content->val[i] = theta_s;

    }
    else {

      if (elt_ids == NULL) {
        assert(cs_glob_n_ranks > 1); // Only possible in parallel

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          gw->moisture_content->val[c_id] = theta_s;
          gw->soil_id[c_id] = soil_id;
        }

      }
      else {

#       pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++) {
          cs_lnum_t  c_id = elt_ids[i];
          gw->moisture_content->val[c_id] = theta_s;
          gw->soil_id[c_id] = soil_id;
        }

      } // elt_ids != NULL

    } // n_soils > 1

    /* Set arrays if needed */
    switch (soil->model) {

    case CS_GWF_HYDRAULIC_GENUCHTEN:
      cs_property_set_array(gw->permeability, head_desc, gw->head_in_law);

      /* Soil capacity settings (related to unsteady term) */
      if (has_previous) {
        cs_property_set_array(cs_equation_get_time_property(richards),
                              head_desc,
                              gw->head_in_law);
      }
      break;

    case CS_GWF_HYDRAULIC_TRACY:
      cs_property_set_array(gw->permeability, head_desc, gw->head_in_law);
      break;

    case CS_GWF_HYDRAULIC_SATURATED:
      break; // Saturated model --> Nothing to do

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of model for groundwater module.."));
      break;

    } /* Switch depending on the type of model */

  } // Loop on soils

  if (gw->n_soils > 1) { /* Sanity check */

    cs_gnum_t  n_unset_cells = 0;

#   pragma omp parallel reduction(+:n_unset_cells) if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      if (gw->soil_id[i] == gw->n_max_soils) n_unset_cells++;

    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&n_unset_cells, 1);
    if (n_unset_cells > 0)
      bft_error(__FILE__, __LINE__, 0,
                " %lu cells are not associated to any soil.\n"
                " Please check your settings.", n_unset_cells);

  } // n_soils > 1

  /* Loop on tracer equations */
  /* ------------------------ */

  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    if (eq_id != gw->richards_eq_id) {
      if (cs_equation_get_type(equations[eq_id]) ==
          CS_EQUATION_TYPE_GROUNDWATER) {

        const cs_equation_t  *tr_eq = equations[eq_id];
        const cs_flag_t  eq_flag = cs_equation_get_flag(tr_eq);

        const cs_desc_t  c_desc =
          {.location = CS_FLAG_SCALAR | cs_cdo_primal_cell,
           .state = CS_FLAG_STATE_DENSITY};

        if (eq_flag & CS_EQUATION_UNSTEADY)
          cs_property_set_array(cs_equation_get_time_property(tr_eq),
                                c_desc,
                                gw->moisture_content->val);

        if (eq_flag & CS_EQUATION_DIFFUSION) {

          cs_property_t  *diff_pty = cs_equation_get_diffusion_property(tr_eq);

          cs_property_set_array(diff_pty, c_desc, gw->moisture_content->val);

          cs_desc_t  desc2 =
            {.location = CS_FLAG_SCALAR | gw->flux_location,
             .state = CS_FLAG_STATE_FLUX};
          cs_property_set_second_array(diff_pty, desc2, gw->darcian_flux);

        } // Diffusion part

        if (eq_flag & CS_EQUATION_REACTION)
          cs_property_set_array(cs_equation_get_reaction_property(tr_eq,
                                                                  "decay"),
                                c_desc,
                                gw->moisture_content->val);

      } // This equation is a tracer equation
    }

  } // Loop on equations

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
 * \param[in, out] eqs        array of pointers to cs_equation_t structures
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t                      *mesh,
                       const cs_time_step_t         *time_step,
                       double                        dt_cur,
                       const cs_cdo_connect_t       *connect,
                       const cs_cdo_quantities_t    *cdoq,
                       cs_equation_t                *eqs[],
                       cs_gwf_t                     *gw)
{
  if (gw == NULL)
    return;

  const int  nt_cur = time_step->nt_cur;

  /* Solve the Richards equation */
  cs_equation_t  *richards = eqs[gw->richards_eq_id];

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  if (nt_cur == 0) {

    /* Initialize system before resolution for all equations
       - create system builder
       - initialize field according to initial conditions
       - initialize source term */
    cs_equation_init_system(mesh, richards);

    /* Take into the initialization */
    _update_system(cdoq, connect, richards, false, gw);

    /* Build and solve the linear system related to the Richards equations */
    if (cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(richards);

      /* Update the variables related to the groundwater flow system */
      _update_system(cdoq, connect, richards, true, gw);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_eq_ids[i]];

      cs_equation_init_system(mesh, eq);

      if (cs_equation_is_steady(eq)) {

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(eq);

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
      cs_equation_solve(richards);

      /* Update the variables related to the groundwater flow system */
      _update_system(cdoq, connect, richards, true, gw);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_equation_t  *eq = eqs[gw->tracer_eq_ids[i]];

      if (!cs_equation_is_steady(eq)) { // unsteady ?

        /* Define the algebraic system */
        if (cs_equation_needs_build(eq))
          cs_equation_build_system(mesh, time_step, dt_cur, eq);

        /* Solve the algebraic system */
        cs_equation_solve(eq);

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
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post(void                      *input,
                  int                        mesh_id,
                  int                        cat_id,
                  int                        ent_flag[5],
                  cs_lnum_t                  n_cells,
                  cs_lnum_t                  n_i_faces,
                  cs_lnum_t                  n_b_faces,
                  const cs_lnum_t            cell_ids[],
                  const cs_lnum_t            i_face_ids[],
                  const cs_lnum_t            b_face_ids[],
                  const cs_time_step_t      *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == NULL)
    return;

  if (mesh_id != -1) /* Post-processing only on the generic volume mesh */
    return;

  const cs_gwf_t  *gw = (const cs_gwf_t *)input;

  if (gw->with_gravitation) { /* Post-process pressure head */

    cs_field_t  *f = gw->pressure_head;

    if (f->location_id == cs_mesh_location_get_id_by_name("cells"))
      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        f->name,
                        1,              // dim
                        true,           // interlace
                        true,           // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        f->val,         // values on cells
                        NULL,           // values at internal faces
                        NULL,           // values at border faces
                        time_step);     // time step structure

    else if (f->location_id == cs_mesh_location_get_id_by_name("vertices"))
      cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                               CS_POST_WRITER_ALL_ASSOCIATED,
                               f->name,
                               1,              // dim
                               false,          // interlace
                               true,           // true = original mesh
                               CS_POST_TYPE_cs_real_t,
                               f->val,         // values on vertices
                               time_step);     // time step management structure


  } // Post pressure head

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
