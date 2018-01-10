/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_field.h"
#include "cs_gwf_soil.h"
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

#include "cs_gwf_tracer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_TRACER_DBG 0

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_tracer[] =
  " Stop execution. The structure related to a tracer is empty.\n"
  " Please check your settings.\n";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_xdef_cell_eval_t
 *
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      elt_ids   list of element ids
 * \param[in]      compact   indirection for output (true or false)
 * \param[in]      mesh      pointer to a cs_mesh_t structure
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts        pointer to a cs_time_step_t structure
 * \param[in]      input     pointer to an input structure cast on-the_fly
 * \param[in, out] result    array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_tracer(cs_lnum_t                    n_elts,
                         const cs_lnum_t              elt_ids[],
                         bool                         compact,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts,
                         void                        *input,
                         cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);

  const cs_gwf_std_tracer_input_t  *law =
    (const cs_gwf_std_tracer_input_t *)input;

  /* Sanity check */
  assert(law != NULL);

  const cs_real_t  *theta = law->moisture_content->val;
  const short int  *c2s = cs_gwf_get_cell2soil();

  if (elt_ids != NULL && !compact) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      result[c_id] = theta[c_id] + law->rho_kd[c2s[c_id]];
    }

  }
  else if (elt_ids != NULL && compact) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      result[i] = theta[c_id] + law->rho_kd[c2s[c_id]];
    }

  }
  else {

    assert(elt_ids == NULL);
    for (cs_lnum_t i = 0; i < n_elts; i++)
      result[i] = theta[i] + law->rho_kd[c2s[i]];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      ts       pointer to a cs_time_step_t structure
 * \param[in]      input    pointer to an input structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_tracer_cw(const cs_cell_mesh_t       *cm,
                            const cs_time_step_t       *ts,
                            void                       *input,
                            cs_real_t                  *result)
{
  CS_UNUSED(ts);

  const cs_gwf_std_tracer_input_t  *law =
    (const cs_gwf_std_tracer_input_t *)input;
  const short int  *c2s = cs_gwf_get_cell2soil();

  /* Sanity check */
  assert(law != NULL);

  const cs_real_t  *theta = law->moisture_content->val;

  *result = theta[cm->c_id] + law->rho_kd[c2s[cm->c_id]];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the reaction term for the
 *         simulation of standard tracer equations.
 *         This function fits the generic prototype of cs_xdef_cell_eval_t
 *
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      elt_ids   list of element ids
 * \param[in]      compact   indirection for output (true or false)
 * \param[in]      mesh      pointer to a cs_mesh_t structure
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts        pointer to a cs_time_step_t structure
 * \param[in]      input     pointer to an input structure cast on-the_fly
 * \param[in, out] result    array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_tracer(cs_lnum_t                    n_elts,
                             const cs_lnum_t              elt_ids[],
                             bool                         compact,
                             const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_time_step_t        *ts,
                             void                        *input,
                             cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);

  const cs_gwf_std_tracer_input_t  *law =
    (const cs_gwf_std_tracer_input_t *)input;

  /* Sanity check */
  assert(law != NULL);

  const cs_real_t  *theta = law->moisture_content->val;
  const short int  *c2s = cs_gwf_get_cell2soil();

  if (elt_ids != NULL && !compact) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      const int s = c2s[c_id];
      result[c_id] = (theta[c_id] + law->rho_kd[s]) * law->reaction_rate[s];
    }

  }
  else if (elt_ids != NULL && compact) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      const int s = c2s[c_id];
      result[i] = (theta[c_id] + law->rho_kd[s]) * law->reaction_rate[s];
    }

  }
  else {

    assert(elt_ids == NULL);
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const int s = c2s[i];
      result[i] = (theta[i] + law->rho_kd[s]) * law->reaction_rate[s];
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the reaction term for the
 *         simulation of standard tracer equations.
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      ts       pointer to a cs_time_step_t structure
 * \param[in]      input    pointer to an input structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_tracer_cw(const cs_cell_mesh_t     *cm,
                                const cs_time_step_t     *ts,
                                void                     *input,
                                cs_real_t                *result)
{
  CS_UNUSED(ts);

  const cs_gwf_std_tracer_input_t  *law =
    (const cs_gwf_std_tracer_input_t *)input;

  /* Sanity check */
  assert(law != NULL);

  const short int  *c2s = cs_gwf_get_cell2soil();
  const int s = c2s[cm->c_id];
  const cs_real_t  *theta = law->moisture_content->val;

  *result = (theta[cm->c_id] + law->rho_kd[s]) * law->reaction_rate[s];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update physical properties for a standard tracer model.
 *         Only the diffusivity needs an update (reaction property and time
 *         property are defined by function).
 *         Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] input        pointer to a structure cast on-the-fly
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_diff_pty4std_tracer(void                        *input,
                            const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(ts);

  cs_gwf_tracer_t  *tracer = (cs_gwf_tracer_t *)input;
  if (tracer->diffusivity == NULL)
    return;

  cs_real_t  *values = tracer->diffusivity->val;
  cs_gwf_std_tracer_input_t  *law = (cs_gwf_std_tracer_input_t *)tracer->input;

  /* Sanity check */
  assert(law != NULL);
  assert(tracer->model == CS_GWF_TRACER_STANDARD);

  const cs_real_t  *theta = law->moisture_content->val;
  const cs_real_t  *velocity = law->darcy_velocity_field->val;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    const cs_volume_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);
    const double  wmd = law->wmd[soil_id];
    const double  at = law->alpha_t[soil_id];
    const double  al = law->alpha_l[soil_id];

    if (z->n_cells == quant->n_cells) { // No need to apply indirection

      for (cs_lnum_t c_id = 0; c_id < z->n_cells; c_id++) {

        const cs_real_t  *v = velocity + 3*c_id;
        const double  v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
        const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
        const double  coef1 = wmd * theta[c_id] + at*vnorm;

        double  delta = 0.;
        if (vnorm > cs_math_zero_threshold)
          delta = (al - at)/vnorm;

        const double  dcv[3] = {delta*v[0], delta*v[1], delta*v[2]};

        cs_real_t  *_r = values + 9*c_id;
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
    else {

      for (cs_lnum_t i = 0; i < z->n_cells; i++) {

        cs_lnum_t  c_id = z->cell_ids[i];
        const cs_real_t  *v = velocity + 3*c_id;
        const double  v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
        const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
        const double  coef1 = wmd * theta[c_id] + at*vnorm;

        double  delta = 0.;
        if (vnorm > cs_math_zero_threshold)
          delta = (al - at)/vnorm;

        const double  dcv[3] = {delta*v[0], delta*v[1], delta*v[2]};

        cs_real_t  *_r = values + 9*c_id;
        for (int ki = 0; ki < 3; ki++) {

          /* Diagonal terms */
          _r[3*ki+ki] = coef1 + delta*v2[ki];

          /* Extra-diagonal terms */
          for (int kj = ki + 1; kj < 3; kj++) {
            _r[3*ki+kj] = dcv[ki]*v[kj];
            _r[3*kj+ki] = _r[3*ki+kj]; /* tensor is symmetric by construction */
          }
        }

      } // Loop on cells attached to this soil

    } // Unique soil ?

  } // Loop on soils

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the input related to a standard tracer equation
 *         Rely on the generic prototype cs_gwf_tracer_free_input_t
 *
 * \param[in, out] input     pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

static void
_free_std_tracer(void      *input)
{
  cs_gwf_std_tracer_input_t  *sti = (cs_gwf_std_tracer_input_t *)input;

  if (sti == NULL)
    return;

  BFT_FREE(sti->rho_kd);
  BFT_FREE(sti->alpha_l);
  BFT_FREE(sti->alpha_t);
  BFT_FREE(sti->wmd);
  BFT_FREE(sti->reaction_rate);

  BFT_FREE(sti);

  /* moisture content and darcy velocity fields are freed thanks to another
     mechanism (They are only shared) */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new cs_gwf_tracer_t structure and initialize its members by
 *         default.
 *         Add a new equation related to the groundwater flow module.
 *         This equation is a specific transport equation.
 *         Tracer is advected thanks to the darcian velocity which is given
 *         by the resolution of the Richards equation.
 *         Diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in]   tracer_id   id number of the soil
 * \param[in]   eqname      name of the tracer equation
 * \param[in]   varname     name of the related variable
 * \param[in]   model       model related to this tracer
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_init(int                      tracer_id,
                   const char              *eq_name,
                   const char              *var_name,
                   cs_adv_field_t          *adv_field,
                   cs_gwf_tracer_model_t    model)
{
  cs_gwf_tracer_t  *tracer = NULL;

  BFT_MALLOC(tracer, 1, cs_gwf_tracer_t);

  tracer->id = tracer_id;
  tracer->eq = cs_equation_add(eq_name,
                               var_name,
                               CS_EQUATION_TYPE_GROUNDWATER,
                               1,
                               CS_PARAM_BC_HMG_NEUMANN);

  tracer->model = model;
  tracer->input = NULL;
  tracer->diffusivity = NULL;
  tracer->reaction_id = -1;

  tracer->update_properties = NULL;
  tracer->free_input = NULL;

  /* Add a new property related to the time-depedent term */
  char  *pty_name = NULL;
  int  len = strlen(eq_name) + strlen("_time") + 1;
  BFT_MALLOC(pty_name, len, char);
  sprintf(pty_name, "%s_time", eq_name);

  cs_property_t  *time_pty = cs_property_add(pty_name, CS_PROPERTY_ISO);

  BFT_FREE(pty_name);

  cs_equation_param_t  *tr_eqp = cs_equation_get_param(tracer->eq);

  cs_equation_add_time(tr_eqp,  time_pty);

  /* Associate the advection field for the advection term */
  assert(adv_field != NULL); /* Sanity check */
  cs_equation_add_advection(tr_eqp, adv_field);

  cs_equation_set_param(tr_eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
  cs_equation_set_param(tr_eqp, CS_EQKEY_ITSOL, "bicg");
  cs_equation_set_param(tr_eqp, CS_EQKEY_BC_ENFORCEMENT, "weak");
  cs_equation_set_param(tr_eqp, CS_EQKEY_ADV_SCHEME, "sg");

  const int  n_soils = cs_gwf_get_n_soils();

  switch (model) {
  case CS_GWF_TRACER_STANDARD:
    {
      cs_gwf_std_tracer_input_t  *input = NULL;

      BFT_MALLOC(input, 1, cs_gwf_std_tracer_input_t);

      BFT_MALLOC(input->rho_kd, n_soils, double);
      BFT_MALLOC(input->alpha_l, n_soils, double);
      BFT_MALLOC(input->alpha_t, n_soils, double);
      BFT_MALLOC(input->wmd, n_soils, double);
      BFT_MALLOC(input->reaction_rate, n_soils, double);

      input->darcy_velocity_field = NULL;
      input->moisture_content = NULL;

      tracer->input = input;

      tracer->update_properties = _update_diff_pty4std_tracer;
      tracer->free_input = _free_std_tracer;
    }
    break;

  case CS_GWF_TRACER_USER:
    break; // All is done during the finalization of the setup

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid model of tracer.");
  }

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_gwf_tracer_t structure
 *
 * \param[in, out]  tracer   pointer to a cs_gwf_tracer_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_free(cs_gwf_tracer_t     *tracer)
{
  if (tracer == NULL)
    return tracer;

  if (tracer->free_input != NULL)
    tracer->free_input(tracer->input);

  /* Tracer equation is freed with all equations in the same time */

  BFT_FREE(tracer);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a tracer for a specified soil when the tracer is attached to
 *         the default model
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or NULL if all
 *                                 soils are selected)
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 * \param[in]      reaction_rate   value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_standard_tracer(cs_gwf_tracer_t   *tracer,
                           const char        *soil_name,
                           double             wmd,
                           double             alpha_l,
                           double             alpha_t,
                           double             distrib_coef,
                           double             reaction_rate)
{
  /* Sanity check */
  if (tracer == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_tracer));
  if (tracer->model != CS_GWF_TRACER_STANDARD)
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible model of tracer.\n"
              " Expect a CS_GWF_TRACER_STANDARD tracer model.\n"
              " Please check your settings.");

  cs_gwf_std_tracer_input_t  *sti = (cs_gwf_std_tracer_input_t *)tracer->input;

  /* Look for the related soil */
  if (soil_name == NULL) { /* All soils have to be set for this tracer */

    const int n_soils = cs_gwf_get_n_soils();
    for (int soil_id = 0; soil_id < n_soils; soil_id++) {

      cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
      assert(soil != NULL);
      cs_real_t  bulk_density = cs_gwf_soil_get_bulk_density(soil);

      sti->rho_kd[soil_id] = bulk_density * distrib_coef;
      sti->alpha_l[soil_id] = alpha_l;
      sti->alpha_t[soil_id] = alpha_t;
      sti->wmd[soil_id] = wmd;
      sti->reaction_rate[soil_id] = reaction_rate;

    } // Loop on soils

  }
  else { /* Set this tracer equation for a specific soil */

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_name(soil_name);
    if (soil == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " Soil %s not found among the predefined soils.\n"
                " Please check your settings.", soil_name);
    cs_real_t  bulk_density = cs_gwf_soil_get_bulk_density(soil);

    sti->rho_kd[soil->id] = bulk_density * distrib_coef;
    sti->alpha_l[soil->id] = alpha_l;
    sti->alpha_t[soil->id] = alpha_t;
    sti->wmd[soil->id] = wmd;
    sti->reaction_rate[soil->id] = reaction_rate;

  } /* Set a specific couple (tracer, soil) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add terms to the algebraic system related to a tracer equation
 *         according to the settings.
 *         Case of the standar tracer modelling
 *         Rely on the generic function: cs_gwf_tracer_add_terms_t
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_standard_add_terms(cs_gwf_tracer_t     *tracer)
{
  /* Sanity check */
  if (tracer == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " At least one tracer equation has not been set.\n"
              " Please check your settings.");
  if (tracer->model != CS_GWF_TRACER_STANDARD)
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible model of tracer.\n"
              " Expect a CS_GWF_TRACER_STANDARD tracer model.\n"
              " Please check your settings.");

  cs_gwf_std_tracer_input_t  *param =
    (cs_gwf_std_tracer_input_t  *)tracer->input;
  cs_equation_param_t  *eqp = cs_equation_get_param(tracer->eq);

  const int n_soils = cs_gwf_get_n_soils();
  const double  thd = 100*DBL_MIN; // threshold to avoid a wrong activation
  const char *eq_name = cs_equation_get_name(tracer->eq);

  bool  do_diffusion = false, do_reaction = false;

  /* Loop on soils to check in a reaction term is needed */
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    if (fabs(param->alpha_t[soil_id]) > thd) do_diffusion = true;
    if (fabs(param->alpha_l[soil_id]) > thd) do_diffusion = true;
    if (param->wmd[soil_id] > thd) do_diffusion = true;
    if (fabs(param->reaction_rate[soil_id]) > thd) do_reaction = true;

  }

  int  max_len = 0;
  char  *pty_name = NULL;

  if (do_diffusion) { /* Add a new diffusion property for this equation */

    int  len = strlen(eq_name) + strlen("_diffusivity") + 1;
    if (len > max_len) {
      max_len = len;
      BFT_REALLOC(pty_name, len, char);
    }
    sprintf(pty_name, "%s_diffusivity", eq_name);

    cs_property_t *diff_pty = cs_property_add(pty_name, CS_PROPERTY_ANISO);

    cs_equation_add_diffusion(eqp, diff_pty);

    /* Create a new field related to this property */
    const int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
    const bool  pty_has_previous = false; // No storage of different snapshots
    const int  field_dim = 9; // anisotropic
    const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");

    tracer->diffusivity = cs_field_create(pty_name,
                                          pty_mask,
                                          c_loc_id,
                                          field_dim,
                                          pty_has_previous);

    cs_field_set_key_int(tracer->diffusivity, cs_field_key_id("log"), 1);

  } // diffusion

  if (do_reaction) { /* Add a new reaction property for this equation */

    int  len = strlen(eq_name) + strlen("_reaction") + 1;
    if (len > max_len) {
      max_len = len;
      BFT_REALLOC(pty_name, len, char);
    }
    sprintf(pty_name, "%s_reaction", eq_name);

    cs_property_t *r_pty = cs_property_add(pty_name, CS_PROPERTY_ISO);

    tracer->reaction_id = cs_equation_add_reaction(eqp, r_pty);

  } // reaction

  BFT_FREE(pty_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_standard_setup(const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             cs_gwf_tracer_t             *tracer)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  /* Sanity check */
  if (tracer == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " At least one tracer equation has not been set.\n"
              " Please check your settings.");
  if (tracer->model != CS_GWF_TRACER_STANDARD)
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible model of tracer.\n"
              " Expect a CS_GWF_TRACER_STANDARD tracer model.\n"
              " Please check your settings.");

  const int n_soils = cs_gwf_get_n_soils();
  const cs_flag_t  eq_flag = cs_equation_get_flag(tracer->eq);

  cs_gwf_std_tracer_input_t *sti = (cs_gwf_std_tracer_input_t *)tracer->input;

  /* Set additional (pre-defined) fields */
  sti->darcy_velocity_field = cs_field_by_name("darcian_flux_cells");
  sti->moisture_content = cs_field_by_name("moisture_content");

  cs_property_t  *pty = cs_equation_get_time_property(tracer->eq);
  assert(pty != NULL);

  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const cs_volume_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_property_def_by_func(pty,
                            z->name,
                            (void *)tracer->input,
                            _get_time_pty4std_tracer,
                            _get_time_pty4std_tracer_cw);

  } // Loop on soils

  if (eq_flag & CS_EQUATION_DIFFUSION) { /* Setup the diffusion property */

    assert(tracer->diffusivity != NULL &&
           tracer->diffusivity->val != NULL); // Should be done previously

    cs_property_t  *diff_pty = cs_equation_get_diffusion_property(tracer->eq);

    cs_property_def_by_field(diff_pty, tracer->diffusivity);

  } // diffusion

  if (eq_flag & CS_EQUATION_REACTION) { /* Setup the reaction property */

    for (int soil_id = 0; soil_id < n_soils; soil_id++) {

      const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
      const cs_volume_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      cs_property_t  *r_pty =
        cs_equation_get_reaction_property(tracer->eq, tracer->reaction_id);

      cs_property_def_by_func(r_pty,
                              z->name,
                              (void *)tracer->input,
                              _get_reaction_pty4std_tracer,
                              _get_reaction_pty4std_tracer_cw);

    } // Loop on soils

  } // reaction

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
