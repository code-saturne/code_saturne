/*============================================================================
 * Main functions dedicated to groundwater flows
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

#include "cs_array.h"
#include "cs_field.h"
#include "cs_gwf_soil.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_reco.h"
#include "cs_sles_it.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_tracer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_TRACER_DBG  0

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_tracer[] =
  " Stop execution. The structure related to a tracer is empty.\n"
  " Please check your settings.\n";

/* Store shared structures (tracers and decay chains) */

static int  _n_tracers = 0;
static cs_gwf_tracer_t **_tracers   = nullptr;

static int  _n_decay_chains = 0;
static cs_gwf_tracer_decay_chain_t **_decay_chains   = nullptr;

/* liquid saturation also called the moisture content, denoted by \theta (-) no
   unit. This array is shared across all tracers. It may be allocated inside
   this file or shared according to the type of hydraulic model. */

static cs_real_t *cs_shared_liquid_saturation = nullptr;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations.
 *         Case of a fully saturated model.
 *         This function fits the generic prototype of cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      context       nullptr or pointer to a structure cast
 * on-the_fly \param[in, out] result        array storing the result (must be
 * allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_sat_tracer(cs_lnum_t                    n_elts,
                             const cs_lnum_t              elt_ids[],
                             bool                         dense_output,
                             const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             cs_real_t                    t_eval,
                             void                        *context,
                             cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_lnum_t  c_id    = (elt_ids == nullptr) ? i : elt_ids[i];
    const cs_lnum_t  id = dense_output ? i : c_id;
    const short int  soil_id = c2s[c_id];
    const cs_real_t  saturated_moisture =
      cs_gwf_soil_get_saturated_moisture(soil_id);

    result[id] = saturated_moisture + tc->rho_kd[soil_id];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations.
 *         Case of a fully saturated model.
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      context  pointer to an context structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_sat_tracer_cw(const cs_cell_mesh_t    *cm,
                                cs_real_t                t_eval,
                                void                    *context,
                                cs_real_t               *result)
{
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();
  const short int  soil_id = c2s[cm->c_id];
  const cs_real_t  saturated_moisture =
    cs_gwf_soil_get_saturated_moisture(soil_id);

  *result = saturated_moisture + tc->rho_kd[soil_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      context       nullptr or pointer to a structure cast
 * on-the_fly \param[in, out] result        array storing the result (must be
 * allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_tracer(cs_lnum_t                    n_elts,
                         const cs_lnum_t              elt_ids[],
                         bool                         dense_output,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         cs_real_t                    t_eval,
                         void                        *context,
                         cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const cs_real_t  *theta = cs_shared_liquid_saturation;
  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  if (elt_ids == nullptr)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      result[i] = theta[i] + tc->rho_kd[c2s[i]];

  else { /* Loop on a selection of cells */

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      const cs_lnum_t  id = dense_output ? i : c_id;
      result[id] = theta[c_id] + tc->rho_kd[c2s[c_id]];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      context  pointer to an context structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_time_pty4std_tracer_cw(const cs_cell_mesh_t    *cm,
                            cs_real_t                t_eval,
                            void                    *context,
                            cs_real_t               *result)
{
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  const short int  *c2s = cs_gwf_soil_get_cell2soil();
  assert(tc != nullptr);

  *result = cs_shared_liquid_saturation[cm->c_id] + tc->rho_kd[c2s[cm->c_id]];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the reaction term for the
 *         simulation of standard tracer equations.
 *         Case of a fully saturated model
 *         This function fits the generic prototype of cs_xdef_cell_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  indirection for output (true or false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      context       nullptr or pointer to a structure cast
 * on-the_fly \param[in, out] result        array storing the result (must be
 * allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_sat_tracer(cs_lnum_t                    n_elts,
                                 const cs_lnum_t              elt_ids[],
                                 bool                         dense_output,
                                 const cs_mesh_t             *mesh,
                                 const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *quant,
                                 cs_real_t                    t_eval,
                                 void                        *context,
                                 cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_lnum_t  c_id               = (elt_ids == nullptr) ? i : elt_ids[i];
    const cs_lnum_t  id = (dense_output) ? i : c_id;
    const short int  s = c2s[c_id];
    const cs_real_t  saturated_moisture = cs_gwf_soil_get_saturated_moisture(s);

    result[id] = (saturated_moisture + tc->rho_kd[s]) * tc->decay_coef;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the reaction term for the
 *         simulation of standard tracer equations.
 *         Case of a fully saturated model.
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      context  nullptr or pointer to a structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_sat_tracer_cw(const cs_cell_mesh_t     *cm,
                                    cs_real_t                 t_eval,
                                    void                     *context,
                                    cs_real_t                *result)
{
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();
  const int s = c2s[cm->c_id];
  const cs_real_t  saturated_moisture = cs_gwf_soil_get_saturated_moisture(s);

  *result = (saturated_moisture + tc->rho_kd[s]) * tc->decay_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in the reaction term for the
 *         simulation of standard tracer equations.
 *         This function fits the generic prototype of cs_xdef_cell_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  indirection for output (true or false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      context       nullptr or pointer to a structure cast
 * on-the_fly \param[in, out] result        array storing the result (must be
 * allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_tracer(cs_lnum_t                    n_elts,
                             const cs_lnum_t              elt_ids[],
                             bool                         dense_output,
                             const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             cs_real_t                    t_eval,
                             void                        *context,
                             cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const cs_real_t  *theta = cs_shared_liquid_saturation;
  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  if (elt_ids == nullptr) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const int s = c2s[i]; /* soil_id */
      result[i]   = (theta[i] + tc->rho_kd[s]) * tc->decay_coef;
    }
  }
  else {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t c_id = elt_ids[i];
      const int       s    = c2s[c_id];
      const cs_lnum_t id   = (dense_output) ? i : c_id;

      result[id] = (theta[c_id] + tc->rho_kd[s]) * tc->decay_coef;
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
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      context  nullptr or pointer to a structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_reaction_pty4std_tracer_cw(const cs_cell_mesh_t *cm,
                                cs_real_t             t_eval,
                                void                 *context,
                                cs_real_t            *result)
{
  CS_UNUSED(t_eval);

  const cs_gwf_tracer_default_context_t *tc
    = (const cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int *c2s = cs_gwf_soil_get_cell2soil();
  const int        s   = c2s[cm->c_id];

  *result
    = (cs_shared_liquid_saturation[cm->c_id] + tc->rho_kd[s]) * tc->decay_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the diffusion property for a (non-user) tracer model.
 *         Only the diffusivity is updated (reaction property and time
 *         property are defined and updated thanks to a function).
 *         Case of a diffusity defined by a value.
 *         Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_diff_value(cs_gwf_tracer_t           *tracer,
                   void                      *context,
                   const cs_time_step_t      *ts,
                   const cs_mesh_t           *mesh,
                   const cs_cdo_connect_t    *connect,
                   const cs_cdo_quantities_t *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(quant);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->diffusivity == nullptr)
    return;

  cs_real_t                       *values = tracer->diffusivity->val;
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;
  assert(tc != nullptr);

  const int n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    cs_gwf_soil_t *soil = cs_gwf_soil_by_id(soil_id);

    const cs_zone_t *z   = cs_volume_zone_by_id(soil->zone_id);
    const double     wmd = tc->wmd[soil_id];

    assert(fabs(tc->alpha_t[soil_id]) < DBL_MIN
           && fabs(tc->alpha_l[soil_id]) < DBL_MIN);

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];

      values[c_id] = wmd;

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update physical properties for a (non-user) tracer model.
 *         Only the diffusivity is updated (reaction property and time
 *         property are defined by function).
 *         Case of an unsaturated model and diffusity defined by a tensor.
 *         Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_diff_tensor(cs_gwf_tracer_t             *tracer,
                    void                        *context,
                    const cs_time_step_t        *ts,
                    const cs_mesh_t             *mesh,
                    const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(quant);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->diffusivity == nullptr)
    return;

  cs_real_t  *values = tracer->diffusivity->val;
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;
  assert(tc != nullptr);

  const cs_real_t  *velocity = tc->darcy_velocity_field->val;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);
    const double  wmd = tc->wmd[soil_id];
    const double  at = tc->alpha_t[soil_id];
    const double  al = tc->alpha_l[soil_id];

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t   c_id  = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
      const cs_real_t  *v = velocity + 3*c_id;
      const double  v2[3] = {v[0]*v[0], v[1]*v[1], v[2]*v[2]};
      const double  vnorm = sqrt(v2[0] + v2[1] + v2[2]);
      const double  coef1 = wmd + at*vnorm;

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

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update physical properties for a non-user tracer model.
 *        Case of a tracer with the precipitation/dissolution modelling and a
 *        vertex-based scheme.
 *        Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_precipitation_vb(cs_gwf_tracer_t             *tracer,
                         void                        *context,
                         const cs_time_step_t        *ts,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant)
{
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(context);

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  assert(tc != nullptr);
  assert(tc->precip_mass != nullptr);
  assert(cs_shared_liquid_saturation != nullptr);

  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_real_t  *pvol_vc = quant->pvol_vc;
  const cs_real_t  *theta = cs_shared_liquid_saturation;
  const double  dt = ts->dt[0]; /* current time step */

  /* Retrieve the current values of the concentration of tracer in the liquid
     phase */

  cs_real_t  *c_l = cs_equation_get_vertex_values(tracer->equation, false);

  /* c_pcp = concentration in precipitate [mol/kg]
   * m_pcp = current number of moles stored as a precipitate [mol]
   * m_l_vc = estimated number of moles of tracer in the liquid phase after the
   *          resolution of the transport equation associated to the tracer.
   *          This quantity is stored on the pvol_vc (= |c cap dcell(v)|) submesh
   */

  cs_real_t  *c_pcp = tc->precip_field->val;
  cs_real_t  *m_pcp = tc->precip_mass;
  cs_real_t  *m_l_vc = nullptr;

  BFT_MALLOC(m_l_vc, c2v->idx[n_cells], cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t  theta_c = theta[c_id];
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      m_l_vc[j] = theta_c * pvol_vc[j] * c_l[c2v->ids[j]];
  }

  /* 2) Update c_l and m_pcp
   *    Update the value of concentration in precipitate in each cell */
  /*    ------------------------------------------------------------- */

  const int  n_soils = cs_gwf_get_n_soils();
  const double  lambda = tc->decay_coef;

  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    if (z->n_elts < 1)
      continue;

    const double  rho = tc->rho_bulk[soil->id];
    const cs_real_t  c_star = tc->conc_l_star[soil->id];

    for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on cells */

      const cs_lnum_t  c_id    = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
      const cs_real_t  theta_c = theta[c_id];

      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

        const cs_lnum_t  v_id = c2v->ids[j];

        double  delta_m = 0.;

        /* Radioactive decay
         * -----------------
         * A part of the stock disappears thanks to the radioactive decay
         * Simply solve:
         *
         * d/dt m_pcp = -lambda*m_pcp
         *
         * Integration between t^n and t^(n+1) yields
         *   m_pcp[j](t^(n+1)) = m_pcp[j](t^n) * exp(-lambda*dt)
         */

        m_pcp[j] *= exp(-lambda*dt);

        /* Precipitation/dissolution
         * ------------------------- */

        if (c_l[v_id] > c_star) { /* Precipitation */

          delta_m = theta_c*(c_l[v_id] - c_star)*pvol_vc[j];
          m_pcp[j] += delta_m;
          m_l_vc[j] -= delta_m;

        }
        else { /* c_l[v_id] <= c_star: dissolution ? */

          if (m_pcp[j] > 0) { /* Dissolution */

            /* Estimate the concentration in the liquid phase if all the mass
               of precipitate is disolved. Threshold given by c_star */

            double  c_l_max = c_l[v_id] + m_pcp[j]/(theta_c*pvol_vc[j]);

            c_l_max = CS_MIN(c_star, c_l_max);
            delta_m = theta_c*(c_l_max - c_l[v_id])*pvol_vc[j];
            m_pcp[j] -= delta_m;
            m_l_vc[j] += delta_m;

          }

        } /* Dissolution ? */

      } /* Loop on cell vertices */

      /* The concentration of precipitate in each cell [mol/kg] can be now
         computed to define the average concentration in precipitate in each
         cell (for a post-processing usage) */

      c_pcp[c_id] = 0;
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        c_pcp[c_id] += m_pcp[j];
      c_pcp[c_id] /= rho*quant->cell_vol[c_id]; /* Cp = mass(c)/(rho_c*|c|) */

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */

  /* Update C_l in each dual cell.
   * At the end of the first step, C_l stores the mass of tracer (in moles) in
   * each dual cell. Then, the mass is divided by the dual porous volume after
   * a parallel reduction
   */

  cs_array_real_fill_zero(n_vertices, c_l);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      c_l[c2v->ids[j]] += m_l_vc[j];

  /* Parallel synchronization (in case of dissolution) */

  if (connect->vtx_ifs != nullptr)
    cs_interface_set_sum(connect->vtx_ifs,
                         n_vertices,
                         1, false, /* stride, interlace (not useful here) */
                         CS_REAL_TYPE,
                         c_l);

  /* C_l is equal to the number of moles in the liquid phase divided by the
     volume of the liquid phase (the porous volume) in the dual cell volume */

  const double  *dpv = cs_gwf_soil_get_dual_porous_volume();
  assert(dpv != nullptr);
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    c_l[i] /= dpv[i];

  /* Free temporary buffer */

  BFT_FREE(m_l_vc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add and initialize quantities related to the precipitation model
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_precipitation(const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_gwf_tracer_t             *tracer)
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(tracer->equation);
  const cs_lnum_t  n_cells = quant->n_cells;

  /* Allocate and initialize the array storing the quantity of precipitate */

  cs_lnum_t  a_size = 0;

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    a_size = c2v->idx[n_cells];
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    a_size = c2v->idx[n_cells] + n_cells;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme.", __func__);

  }

  BFT_MALLOC(tc->precip_mass, a_size, cs_real_t);
  cs_array_real_fill_zero(a_size, tc->precip_mass);

  if (space_scheme == CS_SPACE_SCHEME_CDOVCB ||
      space_scheme == CS_SPACE_SCHEME_CDOVB) {

    /* Compute the dual volume weigthed by the saturation associated to each
     * cell. This quantity is useful to compute the concentration in the liquid
     * phase from the knowdledge of the mass of tracer inside a dual volume.
     *
     * C_l(v) = mass(dual(v)) / (SUM_c\in C_v \theta_c |c \cap dual(v)|)
     *
     * The vertex-based quantity SUM_c\in C_v \theta_c |c \cap dual(v)| is
     * denoted by dual porous volume
     */

    if (cs_gwf_soil_get_dual_porous_volume() == nullptr)
      cs_gwf_soil_build_dual_porous_volume(quant, connect);

  } /* space scheme with DoFs at vertices */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the integral of the molar concentration of a tracer field
 *        over a given set of cells. This integral turns out to be exact for
 *        linear functions. A parallel operation (a sum reduction) is performed
 *        inside this function.
 *        Case of a fully saturated hydraulic model with precipitation effects.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq        equation related to a tracer
 * \param[in]      z         pointer to a volume zone structure
 * \param[in]      context   pointer to a context structure for a tracer
 * \param[in, out] results   resulting array of values
 */
/*----------------------------------------------------------------------------*/

static void
_integrate_sat_precip_tracer(const cs_cdo_connect_t         *connect,
                             const cs_cdo_quantities_t      *cdoq,
                             const cs_equation_t            *eq,
                             const cs_zone_t                *z,
                             void                           *context,
                             double                          results[])
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
        const short int  s = c2s[c_id];
        const cs_real_t  sat_moisture = cs_gwf_soil_get_saturated_moisture(s);

        double  _integral = 0., precip = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _integral += cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];
          precip += tc->precip_mass[j];
        }

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipitate state */

        results[0] += (sat_moisture + tc->rho_kd[s]) * _integral;
        results[1] += precip;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVB */

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
        const short int  s = c2s[c_id];
        const cs_real_t  sat_moisture = cs_gwf_soil_get_saturated_moisture(s);

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
        */

        double  _integral = 0.25*cdoq->cell_vol[c_id]*c_vals[c_id];
        double  precip = tc->precip_mass[c_id];

        const cs_real_t  *precip_mass_shifted = tc->precip_mass + cdoq->n_cells;

        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _integral += 0.75 * cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];
          precip += precip_mass_shifted[j];
        }

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipitate state */

        results[0] += (sat_moisture + tc->rho_kd[s]) * _integral;
        results[1] += precip;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVCB */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme", __func__);
    break;

  } /* End of switch on the space scheme */

  /* Parallel synchronization (sum over all ranks) */

  cs_parall_sum(2, CS_DOUBLE, results);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the integral of the molar concentration of a tracer field
 *        over a given set of cells. This integral turns out to be exact for
 *        linear functions.
 *        Case of a fully saturated hydraulic model.
 *        Parallel synchronized is done inside this function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq        equation related to a tracer
 * \param[in]      z         pointer to a volume zone structure
 * \param[in]      context   pointer to a context structure for a tracer
 * \param[in, out] results   resulting array of values
 */
/*----------------------------------------------------------------------------*/

static void
_integrate_sat_tracer(const cs_cdo_connect_t           *connect,
                      const cs_cdo_quantities_t        *cdoq,
                      const cs_equation_t              *eq,
                      const cs_zone_t                  *z,
                      void                             *context,
                      double                            results[])
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
        const short int  s = c2s[c_id];
        const cs_real_t  sat_moisture = cs_gwf_soil_get_saturated_moisture(s);

        double  _integral = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          _integral += cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipate state */

        results[0] += (sat_moisture + tc->rho_kd[s]) * _integral;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVB */

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];
        const short int  s = c2s[c_id];
        const cs_real_t  sat_moisture = cs_gwf_soil_get_saturated_moisture(s);

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
        */

        double  _integral = 0.25*cdoq->cell_vol[c_id]*c_vals[c_id];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          _integral += 0.75 * cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipate state */

        results[0] += (sat_moisture + tc->rho_kd[s]) * _integral;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVCB */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme", __func__);
    break;

  } /* End of switch on the space scheme */

  /* Parallel synchronization (sum over all ranks) */

  cs_parall_sum(1, CS_DOUBLE, results);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the integral of the molar concentration of a tracer field
 *        over a given set of cells. This integral turns out to be exact for
 *        linear functions.
 *        Case of a single-phase unsaturated hydraulic model without
 *        precipitation.
 *        Parallel synchronized is done inside this function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq        equation related to a tracer
 * \param[in]      z         pointer to a volume zone structure
 * \param[in]      context   pointer to a context structure for a tracer
 * \param[in, out] results   resulting array of values
 */
/*----------------------------------------------------------------------------*/

static void
_integrate_tracer(const cs_cdo_connect_t          *connect,
                  const cs_cdo_quantities_t       *cdoq,
                  const cs_equation_t             *eq,
                  const cs_zone_t                 *z,
                  void                            *context,
                  double                           results[])
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();
  const cs_real_t  *moisture_val = cs_shared_liquid_saturation;

  if (moisture_val == nullptr)
    bft_error(__FILE__, __LINE__, 0, " %s: \"moisture_content\" not defined",
              __func__);

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];

        double  _integral = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          _integral += cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipate state */

        results[0] += (moisture_val[c_id] + tc->rho_kd[c2s[c_id]]) * _integral;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVB */

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
        */

        double  _integral = 0.25*cdoq->cell_vol[c_id]*c_vals[c_id];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          _integral += 0.75 * cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipate state */

        results[0] += (moisture_val[c_id] + tc->rho_kd[c2s[c_id]]) * _integral;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVCB */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme", __func__);
    break;

  } /* End of switch on the space discretization */

  /* Parallel synchronization (sum over all ranks) */

  cs_parall_sum(1, CS_DOUBLE, results);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the integral of the molar concentration of a tracer field
 *        over a given set of cells. This integral turns out to be exact for
 *        linear functions.
 *        Case of a single-phase unsaturated hydraulic model with precipitation.
 *        Parallel synchronized is done inside this function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq        equation related to a tracer
 * \param[in]      z         pointer to a volume zone structure
 * \param[in]      context   pointer to a context structure for a tracer
 * \param[in, out] results   resulting array of values
 */
/*----------------------------------------------------------------------------*/

static void
_integrate_precip_tracer(const cs_cdo_connect_t          *connect,
                         const cs_cdo_quantities_t       *cdoq,
                         const cs_equation_t             *eq,
                         const cs_zone_t                 *z,
                         void                            *context,
                         double                           results[])
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)context;
  assert(tc != nullptr);

  const short int  *c2s = cs_gwf_soil_get_cell2soil();
  const cs_real_t  *moisture_val = cs_shared_liquid_saturation;

  if (moisture_val == nullptr)
    bft_error(__FILE__, __LINE__, 0, " %s: \"moisture_content\" not defined",
              __func__);

  switch (cs_equation_get_space_scheme(eq)) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];

        double  _integral = 0., precip = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _integral += cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];
          precip += tc->precip_mass[j];
        }

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipitate state */

        results[0] += (moisture_val[c_id] + tc->rho_kd[c2s[c_id]]) * _integral;
        results[1] += precip;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVB */

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(eq, false);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(eq, false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t c_id = (z->elt_ids == nullptr) ? i : z->elt_ids[i];

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
        */

        double  _integral = 0.25*cdoq->cell_vol[c_id]*c_vals[c_id];
        double  precip = tc->precip_mass[c_id];

        const cs_real_t  *precip_mass_shifted = tc->precip_mass + cdoq->n_cells;

        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _integral += 0.75 * cdoq->pvol_vc[j] * v_vals[c2v->ids[j]];
          precip += precip_mass_shifted[j];
        }

        /* results[0] => quantity of tracer in the liquid phase
         * results[1] => quantity of tracer in the precipitate state */

        results[0] += (moisture_val[c_id] + tc->rho_kd[c2s[c_id]]) * _integral;
        results[1] += precip;

      } /* Loop on the selected cells */
    }
    break; /* CS_SPACE_SCHEME_CDOVCB */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme", __func__);
    break;

  } /* End of switch on the space discretization */

  /* Parallel synchronization (sum over all ranks) */

  cs_parall_sum(1, CS_DOUBLE, results);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context related to a standard tracer equation
 *         Rely on the generic prototype cs_gwf_tracer_free_context_t
 *
 * \param[in, out] tracer     pointer to a structure cs_gwf_tracer_t
 */
/*----------------------------------------------------------------------------*/

static void
_free_default_tracer_context(cs_gwf_tracer_t   *tracer)
{
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  if (tc == nullptr)
    return;

  BFT_FREE(tc->rho_bulk);
  BFT_FREE(tc->kd0);
  BFT_FREE(tc->rho_kd);
  BFT_FREE(tc->alpha_l);
  BFT_FREE(tc->alpha_t);
  BFT_FREE(tc->wmd);

  /* Sorption phenomena */

  if (tracer->model & CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS ||
      tracer->model & CS_GWF_TRACER_SORPTION_EK_5_PARAMETERS) {

    BFT_FREE(tc->k0_plus);
    BFT_FREE(tc->k0_minus);
    BFT_FREE(tc->conc_site2);

  }

  /* Precipitation phenomena */

  if (tracer->model & CS_GWF_TRACER_PRECIPITATION) {

    BFT_FREE(tc->conc_l_star);
    BFT_FREE(tc->precip_mass);

  }

  BFT_FREE(tc);
  tracer->context = nullptr;

  /* All fields are freed thanks to another mechanism */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context by default for a "standard" tracer
 *
 * \param[in, out]  tracer   pointer to a cs_gwf_tracer_t structure
 * \param[in]       lambda   value of the first order decay coefficient
 */
/*----------------------------------------------------------------------------*/

static void
_create_default_tracer_context(cs_gwf_tracer_t    *tracer,
                               double              lambda)
{
  if (tracer == nullptr)
    return;

  if ((tracer->model & CS_GWF_TRACER_USER) != 0) /* user-defined ? */
    return;

  cs_gwf_tracer_default_context_t *context = nullptr;

  BFT_MALLOC(context, 1, cs_gwf_tracer_default_context_t);

  context->decay_coef = lambda;

  /* One handles a standard tracer */

  const int  n_soils = cs_gwf_get_n_soils();

  BFT_MALLOC(context->rho_bulk, n_soils, double);
  BFT_MALLOC(context->kd0, n_soils, double);
  BFT_MALLOC(context->rho_kd, n_soils, double);
  BFT_MALLOC(context->alpha_l, n_soils, double);
  BFT_MALLOC(context->alpha_t, n_soils, double);
  BFT_MALLOC(context->wmd, n_soils, double);

  context->darcy_velocity_field = nullptr;

  /* Sorption members */

  context->k0_plus    = nullptr;
  context->k0_minus   = nullptr;
  context->conc_site2 = nullptr;

  if (tracer->model & CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS) {

    BFT_MALLOC(context->k0_minus, n_soils, double);
    BFT_MALLOC(context->k0_plus, n_soils, double);

  }

  /* Precipitation members */

  context->conc_l_star  = nullptr;
  context->precip_mass  = nullptr;
  context->precip_field = nullptr;

  if (tracer->model & CS_GWF_TRACER_PRECIPITATION)
    BFT_MALLOC(context->conc_l_star, n_soils, double);

  /* Set additional function pointers */

  tracer->update_precipitation = nullptr;

  switch (tracer->hydraulic_model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    if (tracer->model & CS_GWF_TRACER_PRECIPITATION) {
      tracer->integrate = _integrate_sat_precip_tracer;
      tracer->update_precipitation = _update_precipitation_vb;
    }
    else
      tracer->integrate = _integrate_sat_tracer;
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    if (tracer->model & CS_GWF_TRACER_PRECIPITATION) {

      tracer->update_precipitation = _update_precipitation_vb;
      tracer->integrate = _integrate_precip_tracer;

    }
    else
      tracer->integrate = _integrate_tracer;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Precipitation model not implemented in this case.\n",
              __func__);

  } /* Switch on hydraulic model */

  /* Common to all default tracers */

  tracer->context = context;
  tracer->update_diff_pty = nullptr;
  tracer->free_context = _free_default_tracer_context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_gwf_tracer_t structure and initialize its members by
 *        default
 *
 * \param[in]   tr_model    model related to this tracer
 * \param[in]   gwf_model   main model for the GWF module
 * \param[in]   eq_name     name of the tracer equation
 * \param[in]   var_name    name of the related variable
 * \param[in]   adv_field   pointer to a cs_adv_field_t structure
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_tracer_t *
_create_tracer(cs_gwf_tracer_model_t    tr_model,
               cs_gwf_model_type_t      gwf_model,
               const char              *eq_name,
               const char              *var_name,
               cs_adv_field_t          *adv_field)
{
  cs_gwf_tracer_t *tracer = nullptr;

  BFT_MALLOC(tracer, 1, cs_gwf_tracer_t);

  tracer->equation = cs_equation_add(eq_name,
                                     var_name,
                                     CS_EQUATION_TYPE_GROUNDWATER,
                                     1, /* scalar-valued equation */
                                     CS_BC_SYMMETRY);

  tracer->model = tr_model;
  tracer->hydraulic_model = gwf_model;
  tracer->diffusivity       = nullptr;
  tracer->reaction_id = -1;
  tracer->chain_id = -1;          /* Not in a chain by default */
  tracer->chain_position_id = -1; /* Not in a chain */

  /* Add a new property related to the time-depedent term */

  char *pty_name = nullptr;
  int  len = strlen(eq_name) + strlen("_time") + 1;
  BFT_MALLOC(pty_name, len, char);
  sprintf(pty_name, "%s_time", eq_name);

  cs_property_t  *time_pty = cs_property_add(pty_name, CS_PROPERTY_ISO);

  BFT_FREE(pty_name);

  cs_equation_param_t  *tr_eqp = cs_equation_get_param(tracer->equation);

  cs_equation_add_time(tr_eqp,  time_pty);

  /* Associate the advection field for the advection term */

  assert(adv_field != nullptr); /* Sanity check */
  cs_equation_add_advection(tr_eqp, adv_field);

  cs_equation_param_set(tr_eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");

  /* Space discretization */

  cs_equation_param_set(tr_eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
  cs_equation_param_set(tr_eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");

  /* Default treatment of the Dirichlet BCs (weak since there is an advection
     term) */

  cs_equation_param_set(tr_eqp, CS_EQKEY_BC_ENFORCEMENT, "weak");

  /* Default advection scheme: centered scheme with 25/100 of upwinding */

  cs_equation_param_set(tr_eqp, CS_EQKEY_ADV_FORMULATION, "conservative");
  cs_equation_param_set(tr_eqp, CS_EQKEY_ADV_SCHEME, "mix_centered_upwind");
  cs_equation_param_set(tr_eqp, CS_EQKEY_ADV_UPWIND_PORTION, "0.25");

  /* Linear algebra */

  cs_equation_param_set(tr_eqp, CS_EQKEY_SOLVER, "bicgs");
  cs_equation_param_set(tr_eqp, CS_EQKEY_PRECOND, "jacobi");
  cs_equation_param_set(tr_eqp, CS_EQKEY_SOLVER_RTOL, "1e-8");

  /* Function pointers */

  tracer->update_diff_pty       = nullptr;
  tracer->update_precipitation  = nullptr;
  tracer->update_decay_chain_st = nullptr;

  tracer->context        = nullptr;
  tracer->free_context   = nullptr;
  tracer->init_setup     = nullptr;
  tracer->finalize_setup = nullptr;

  tracer->integrate = nullptr;

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the source term taking into account the decay chain for the
 *        current tracer Case of CDO-Vb schemes with saturated soils and mole
 *        unit.
 *        Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vb_sat_decay_chain_molar_st(cs_gwf_tracer_t             *tracer,
                             void                        *context,
                             const cs_time_step_t        *ts,
                             const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->chain_id < 0)
    return;
  if (tracer->chain_position_id < 1)
    return; /* Nothing to do for the common ancestor */

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the decay chain structure */

  cs_gwf_tracer_decay_chain_t  *tdc =
    cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
  assert(tdc != nullptr);

  /* Definition associated to this tracer in the decay chain */

  cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];
  double  *st_values = cs_xdef_array_get_values(st_def);

  cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

  /* Retrieve information on the parent tracer */

  cs_gwf_tracer_t  *tr_parent = tdc->tracers[tracer->chain_position_id-1];
  assert(tr_parent != nullptr);
  const cs_field_t  *tr_field_parent =
    cs_equation_get_field(tr_parent->equation);
  const cs_real_t  *parent_vals = tr_field_parent->val;

  cs_gwf_tracer_default_context_t *tc_parent
    = (cs_gwf_tracer_default_context_t *)tr_parent->context;
  assert(tc_parent != nullptr);
  const double  lamb_parent = tc_parent->decay_coef;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    /* When the soil is saturated, the moisture content is equal to the
     * porosity
     */

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const double  coef =
      lamb_parent * (soil->porosity + tc_parent->rho_kd[soil_id]);
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

#if defined(DEBUG) && !defined(NDEBUG)
    if (z->n_elts > 0)
      assert(z->elt_ids != nullptr);
#endif

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t  c_id = z->elt_ids[i];
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        st_values[j] += parent_vals[c2v->ids[j]] * coef;

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the source term taking into account the decay chain for the
 *        current tracer Case of CDO-Vb schemes with unsaturated soils and mole
 *        unit.
 *        Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vb_decay_chain_molar_st(cs_gwf_tracer_t             *tracer,
                         void                        *context,
                         const cs_time_step_t        *ts,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->chain_id < 0)
    return;
  if (tracer->chain_position_id < 1)
    return; /* Nothing to do for the common ancestor */

  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_real_t  *theta = cs_shared_liquid_saturation;
  assert(theta != nullptr);

  /* Retrieve the decay chain structure */

  cs_gwf_tracer_decay_chain_t  *tdc =
    cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
  assert(tdc != nullptr);

  /* Definition associated to this tracer in the decay chain */

  cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];
  double  *st_values = cs_xdef_array_get_values(st_def);

  cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

  /* Retrieve information on the parent tracer */

  cs_gwf_tracer_t  *tr_parent = tdc->tracers[tracer->chain_position_id-1];
  assert(tr_parent != nullptr);
  const cs_field_t  *tr_field_parent =
    cs_equation_get_field(tr_parent->equation);
  const cs_real_t  *parent_vals = tr_field_parent->val;

  cs_gwf_tracer_default_context_t *tc_parent
    = (cs_gwf_tracer_default_context_t *)tr_parent->context;
  assert(tc_parent != nullptr);
  const double  lamb_parent = tc_parent->decay_coef;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const double  rhokd_parent = tc_parent->rho_kd[soil_id];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

#if defined(DEBUG) && !defined(NDEBUG)
    if (z->n_elts > 0)
      assert(z->elt_ids != nullptr);
#endif

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t  c_id = z->elt_ids[i];
      const double  coef = lamb_parent * (theta[c_id] + rhokd_parent);

      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        st_values[j] += parent_vals[c2v->ids[j]] * coef;

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the source term taking into account the decay chain for the
 *        current tracer. Case of CDO-Vb schemes with saturated soils and
 *        Becquerel unit.
 *        Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vb_sat_decay_chain_becqu_st(cs_gwf_tracer_t             *tracer,
                             void                        *context,
                             const cs_time_step_t        *ts,
                             const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->chain_id < 0)
    return;
  if (tracer->chain_position_id < 1)
    return; /* Nothing to do for the common ancestor */

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the decay chain structure */

  cs_gwf_tracer_decay_chain_t  *tdc =
    cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
  assert(tdc != nullptr);

  /* Definition associated to this tracer in the decay chain */

  cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];
  double  *st_values = cs_xdef_array_get_values(st_def);

  cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

  /* Retrieve information on the parent tracer */

  cs_gwf_tracer_t  *tr_parent = tdc->tracers[tracer->chain_position_id-1];
  assert(tr_parent != nullptr);
  const cs_field_t  *tr_field_parent =
    cs_equation_get_field(tr_parent->equation);
  const cs_real_t  *parent_vals = tr_field_parent->val;

  cs_gwf_tracer_default_context_t *tc_parent
    = (cs_gwf_tracer_default_context_t *)tr_parent->context;
  assert(tc_parent != nullptr);
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;
  assert(tc != nullptr);
  const double  lamb = tc->decay_coef;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    /* When the soil is saturated, the moisture content is equal to the
     * porosity
     */

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const double  coef = lamb * (soil->porosity + tc_parent->rho_kd[soil_id]);
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

#if defined(DEBUG) && !defined(NDEBUG)
    if (z->n_elts > 0)
      assert(z->elt_ids != nullptr);
#endif

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t  c_id = z->elt_ids[i];
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        st_values[j] += parent_vals[c2v->ids[j]] * coef;

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the source term taking into account the decay chain for the
 *        current tracer. Case of CDO-Vb schemes with unsaturated soils and
 *        Becquerel unit.
 *        Generic function relying on the prototype cs_gwf_tracer_update_t
 *
 * \param[in, out] tracer     pointer to a cs_gwf_tracer_structure
 * \param[in, out] context    nullptr or pointer to a structure cast on-the-fly
 * \param[in]      ts         pointer to a cs_time_step_t structure
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vb_decay_chain_becqu_st(cs_gwf_tracer_t             *tracer,
                         void                        *context,
                         const cs_time_step_t        *ts,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant)
{
  /* Parameters not used since it relies on a generic function pointer */

  CS_NO_WARN_IF_UNUSED(context);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(ts);

  assert(tracer != nullptr);
  if (tracer->chain_id < 0)
    return;
  if (tracer->chain_position_id < 1)
    return; /* Nothing to do for the common ancestor */

  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_real_t  *theta = cs_shared_liquid_saturation;

  /* Retrieve the decay chain structure */

  cs_gwf_tracer_decay_chain_t  *tdc =
    cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
  assert(tdc != nullptr);

  /* Definition associated to this tracer in the decay chain */

  cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];
  double  *st_values = cs_xdef_array_get_values(st_def);

  cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

  /* Retrieve information on the parent tracer */

  cs_gwf_tracer_t  *tr_parent = tdc->tracers[tracer->chain_position_id-1];
  assert(tr_parent != nullptr);
  const cs_field_t  *tr_field_parent =
    cs_equation_get_field(tr_parent->equation);
  const cs_real_t  *parent_vals = tr_field_parent->val;

  cs_gwf_tracer_default_context_t *tc_parent
    = (cs_gwf_tracer_default_context_t *)tr_parent->context;
  assert(tc_parent != nullptr);
  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;
  assert(tc != nullptr);
  const double  lamb = tc->decay_coef;

  const int  n_soils = cs_gwf_get_n_soils();
  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    /* When the soil is saturated, the moisture content is equal to the
     * porosity */

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const double  rhokd_parent = tc_parent->rho_kd[soil_id];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

#if defined(DEBUG) && !defined(NDEBUG)
    if (z->n_elts > 0)
      assert(z->elt_ids != nullptr);
#endif

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {

      const cs_lnum_t  c_id = z->elt_ids[i];
      const double  coef = lamb * (theta[c_id] + rhokd_parent);

      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        st_values[j] += parent_vals[c2v->ids[j]] * coef;

    } /* Loop on cells attached to this soil */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the members related to a tracer belonging to a decay chain
 *
 * \param[in, out] tracer          pointer to the tracer to update
 * \param[in]      chain_id        -1 or id of the associated decay chain
 * \param[in]      chain_position  -1 or id in the chain (position)
 */
/*----------------------------------------------------------------------------*/

static void
_set_decay_chain_members(cs_gwf_tracer_t    *tracer,
                         int                 chain_id,
                         int                 chain_position)
{
  if (tracer == nullptr)
    return;

  assert(chain_position > -1 && chain_id > -1);

  tracer->chain_id = chain_id;
  tracer->chain_position_id = chain_position;

  cs_gwf_tracer_decay_chain_t  *tdc = cs_gwf_tracer_decay_chain_by_id(chain_id);
  assert(tdc != nullptr);

  switch (tdc->unit) {

  case CS_GWF_TRACER_UNIT_BECQUEREL:
    if (tracer->hydraulic_model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
      tracer->update_decay_chain_st = _vb_sat_decay_chain_becqu_st;
    else
      tracer->update_decay_chain_st = _vb_decay_chain_becqu_st;
    break;

  case CS_GWF_TRACER_UNIT_MOLE:
    if (tracer->hydraulic_model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
      tracer->update_decay_chain_st = _vb_sat_decay_chain_molar_st;
    else
      tracer->update_decay_chain_st = _vb_decay_chain_molar_st;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Undefined unit for a decay chain.",
              __func__);
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display the main features related to each tracer
 */
/*----------------------------------------------------------------------------*/

static void
_log_decay_chains(void)
{
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Number of decay chains: %d\n", _n_decay_chains);

  if (_n_decay_chains == 0)
    return;
  assert(_decay_chains != nullptr);

  for (int i = 0; i < _n_decay_chains; i++) {

    cs_gwf_tracer_decay_chain_t  *tdc = _decay_chains[i];
    assert(tdc != nullptr);

    cs_log_printf(CS_LOG_SETUP, "\n  * GWF | Decay chain: %s\n", tdc->name);

    if (tdc->unit == CS_GWF_TRACER_UNIT_BECQUEREL)
      cs_log_printf(CS_LOG_SETUP, "  ** %s | Unit: Becquerel\n", tdc->name);
    else if (tdc->unit == CS_GWF_TRACER_UNIT_MOLE)
      cs_log_printf(CS_LOG_SETUP, "  ** %s | Unit: mole\n", tdc->name);

    cs_log_printf(CS_LOG_SETUP, "  ** %s | n_tracers: %d\n", tdc->name,
                  tdc->n_tracers);

    if (tdc->n_tracers > 0) {

      cs_log_printf(CS_LOG_SETUP, "  ** %s | Chain description\n", tdc->name);
      cs_log_printf(CS_LOG_SETUP, "  ** %s |", tdc->name);

      const cs_gwf_tracer_t  *tr0 = tdc->tracers[0];
      cs_log_printf(CS_LOG_SETUP, "%s", cs_equation_get_name(tr0->equation));

      for (int j = 1; j < tdc->n_tracers; j++) {

        const cs_gwf_tracer_t  *tr = tdc->tracers[j];
        cs_log_printf(CS_LOG_SETUP, " --> %s",
                      cs_equation_get_name(tr->equation));

      }

      cs_log_printf(CS_LOG_SETUP, "\n");

    }

  } /* Loop on decay chains */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all decay chain structures used to manage several linked tracers
 */
/*----------------------------------------------------------------------------*/

static void
_free_all_decay_chains(void)
{
  if (_n_decay_chains < 1)
    return;
  assert(_decay_chains != nullptr);

  for (int i = 0; i < _n_decay_chains; i++) {

    cs_gwf_tracer_decay_chain_t  *tdc = _decay_chains[i];

    if (tdc == nullptr)
      continue;

    BFT_FREE(tdc->name);
    BFT_FREE(tdc->tracers);
    BFT_FREE(tdc->st_defs);

    BFT_FREE(tdc);
    _decay_chains[i] = nullptr;
  }

  BFT_FREE(_decay_chains);
  _n_decay_chains = 0;
  _decay_chains   = nullptr;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the cs_gwf_tracer_t structure associated to
 *         the name given as parameter
 *
 * \param[in]  eq_name    name of the tracer equation
 *
 * \return the pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_by_name(const char   *eq_name)
{
  if (_n_tracers == 0)
    return nullptr;

  if (eq_name == nullptr)
    return nullptr;

  if (_tracers == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_tracer));

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = _tracers[i];

    const char *name_to_cmp = cs_equation_get_name(tracer->equation);
    if (strcmp(eq_name, name_to_cmp) == 0)
      return tracer;

  } /* Loop on tracer equations */

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_gwf_tracer_t structure and initialize its members.
 *        This creation of a new tracer is fully done in the case of a default
 *        tracer. Additional settings has to be done in the case of a
 *        user-defined tracer.
 *
 *        Add a new equation related to the groundwater flow module. This
 *        equation is a specific transport equation. The tracer is advected
 *        thanks to the darcian velocity (in the liquid phase) which is given
 *        by the resolution of the Richards equation. Diffusion and reaction
 *        coefficients result from a physical modelling.
 *
 * \param[in] tr_model        model related to this tracer
 * \param[in] gwf_model       main model for the GWF module
 * \param[in] eq_name         name of the tracer equation
 * \param[in] var_name        name of the related variable
 * \param[in] adv_field       pointer to a cs_adv_field_t structure
 * \param[in] lambda          value of the first order decay coefficient
 * \param[in] chain_position  -1 if not used or the id in the chain position
 * \param[in] chain_id        -1 or id of the associated decay chain
 * \param[in] init_setup      function pointer (predefined prototype)
 * \param[in] finalize_setup  function pointer (predefined prototype)
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_add(cs_gwf_tracer_model_t            tr_model,
                  cs_gwf_model_type_t              gwf_model,
                  const char                      *eq_name,
                  const char                      *var_name,
                  cs_adv_field_t                  *adv_field,
                  double                           lambda,
                  int                              chain_position,
                  int                              chain_id,
                  cs_gwf_tracer_init_setup_t      *init_setup,
                  cs_gwf_tracer_finalize_setup_t  *finalize_setup)
{
  int  tr_id = _n_tracers;

  cs_gwf_tracer_t  *tracer = _create_tracer(tr_model,
                                            gwf_model,
                                            eq_name,
                                            var_name,
                                            adv_field);

  assert(tracer != nullptr);

  tracer->init_setup = init_setup;
  tracer->finalize_setup = finalize_setup;

  /* If this tracer is embedded inside a decay chain */

  if (chain_position > -1 && chain_id > -1)
    _set_decay_chain_members(tracer, chain_id, chain_position);

  /* If this is not a user-defined tracer, initialize the default context */

  if ((tracer->model & CS_GWF_TRACER_USER) == 0)
    _create_default_tracer_context(tracer, lambda);

  /* The default value of the breakdown threshold may be too high when dealing
     with a tracer equation since the concentration of radionuclides are (very)
     small in general */

  cs_sles_it_set_breakdown_threshold(1e-36);

  /* Update the array storing all tracers */

  _n_tracers += 1;
  BFT_REALLOC(_tracers, _n_tracers, cs_gwf_tracer_t *);
  _tracers[tr_id] = tracer;

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all tracer structures and all decay chains
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_free_all(void)
{
  /* Destroy all decay chains */

  _free_all_decay_chains();

  /* Free tracer structures */

  if (_n_tracers == 0)
    return;
  assert(_tracers != nullptr);

  /* One assumes that all tracers share the same hydraulic model */

  cs_gwf_tracer_t  *tracer = _tracers[0];

  if (tracer->hydraulic_model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    BFT_FREE(cs_shared_liquid_saturation);
  cs_shared_liquid_saturation = nullptr; /* unset the pointer in all cases */

  for (int i = 0; i < _n_tracers; i++) {

    tracer = _tracers[i];
    if (tracer == nullptr)
      continue;

    if (tracer->free_context != nullptr)
      tracer->free_context(tracer);

    /* Tracer equation is freed with all equations at the same time */

    BFT_FREE(tracer);
    _tracers[i] = nullptr;

  } /* Loop on tracers */

  _n_tracers = 0;
  BFT_FREE(_tracers);
  _tracers = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the max. value of the theta parameter associated to a time
 *        scheme. Loop on all tracer equations.
 *
 * \return the computed value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tracer_get_time_theta_max(void)
{
  cs_real_t  theta = -1;

  if (_n_tracers == 0)
    return theta;

  assert(_tracers != nullptr);

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    theta = fmax(theta, cs_equation_get_theta_time_val(tracer->equation));

  } /* Loop on tracers */

  return theta;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the main parameters corresponding to a default modelling of a
 *        tracer transport equation for a specified soil
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or nullptr if all
 *                                 soils are selected)
 * \param[in]      wmd             value of the water molecular diffusivity
 * \param[in]      alpha_l         value of the longitudinal dispersivity
 * \param[in]      alpha_t         value of the transversal dispersivity
 * \param[in]      distrib_coef    value of the distribution coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_set_soil_param(cs_gwf_tracer_t   *tracer,
                             const char        *soil_name,
                             double             wmd,
                             double             alpha_l,
                             double             alpha_t,
                             double             distrib_coef)
{
  if (tracer == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_tracer));

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  /* Look for the related soil */

  if (soil_name == nullptr) { /* All soils have to be set for this tracer */

    const int n_soils = cs_gwf_get_n_soils();
    for (int soil_id = 0; soil_id < n_soils; soil_id++) {

      cs_gwf_soil_t *soil = cs_gwf_soil_by_id(soil_id);
      assert(soil != nullptr);

      tc->rho_bulk[soil_id] = soil->bulk_density;
      tc->kd0[soil_id]      = distrib_coef;
      tc->rho_kd[soil_id]   = soil->bulk_density * distrib_coef;
      tc->alpha_l[soil_id]  = alpha_l;
      tc->alpha_t[soil_id]  = alpha_t;
      tc->wmd[soil_id]      = wmd;

    } /* Loop on all soils */
  }
  else { /* Set this tracer equation for a specific soil */

    cs_gwf_soil_t *soil = cs_gwf_soil_by_name(soil_name);
    if (soil == nullptr)
      bft_error(__FILE__,
                __LINE__,
                0,
                " Soil %s not found among the predefined soils.\n"
                " Please check your settings.",
                soil_name);

    tc->rho_bulk[soil->id] = soil->bulk_density;
    tc->kd0[soil->id]      = distrib_coef;
    tc->rho_kd[soil->id]   = soil->bulk_density * distrib_coef;
    tc->alpha_l[soil->id]  = alpha_l;
    tc->alpha_t[soil->id]  = alpha_t;
    tc->wmd[soil->id]      = wmd;

  } /* Set a specific couple (tracer, soil) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a specified soil set the parameters corresponding to a
 *         precipitation modelling of a tracer transport
 *
 * \param[in, out] tracer          pointer to a cs_gwf_tracer_t structure
 * \param[in]      soil_name       name of the related soil (or nullptr if all
 *                                 soils are selected)
 * \param[in]      conc_l_star     value of the saturated concentration in the
 *                                 liquid phase
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_set_precip_param(cs_gwf_tracer_t *tracer,
                               const char      *soil_name,
                               double           conc_l_star)
{
  if (tracer == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_tracer));

  if ((tracer->model & CS_GWF_TRACER_PRECIPITATION) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Precipitation model has not been activated for this"
              " tracer", __func__);

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  /* Look for the related soil */

  if (soil_name == nullptr) { /* All soils have to be set for this tracer */

    const int n_soils = cs_gwf_get_n_soils();
    for (int soil_id = 0; soil_id < n_soils; soil_id++)
      tc->conc_l_star[soil_id] = conc_l_star;
  }
  else { /* Set this tracer equation for a specific soil */

    cs_gwf_soil_t *soil = cs_gwf_soil_by_name(soil_name);
    if (soil == nullptr)
      bft_error(__FILE__,
                __LINE__,
                0,
                " Soil %s not found among the predefined soils.\n"
                " Please check your settings.",
                soil_name);

    tc->conc_l_star[soil->id] = conc_l_star;

  } /* Set a specific couple (tracer, soil) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup step for tracer equations. Soils and equation
 *        parameters are defined at this stage.
 *        Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_init_setup(void)
{
  if (_n_tracers == 0)
    return;
  assert(_tracers != nullptr);

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    cs_equation_t *eq = tracer->equation;

    if (tracer->init_setup != nullptr)
      tracer->init_setup(tracer);

    /* Add the variable field (One assumes an unsteady behavior) */

    cs_equation_predefined_create_field(1, eq); /* Keep two states */

  } /* Loop on tracers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the tracer setup
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_finalize_setup(const cs_cdo_connect_t    *connect,
                             const cs_cdo_quantities_t *quant)
{
  if (_n_tracers == 0)
    return;
  assert(_tracers != nullptr);

  /* One assumes that all tracers share the same hydraulic model */

  cs_gwf_tracer_t *tracer = _tracers[0];

  /* Set the liquid saturation */

  switch (tracer->hydraulic_model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE: {
    const char     mc_pty_name[] = "moisture_content";
    cs_property_t *mc            = cs_property_by_name(mc_pty_name);
    if (mc == nullptr)
      bft_error(__FILE__,
                __LINE__,
                0,
                "%s: Expected property \"%s\" is not defined.\n",
                __func__,
                mc_pty_name);

    BFT_MALLOC(cs_shared_liquid_saturation, quant->n_cells, cs_real_t);

    /* For a saturated model there is no time evolution of the liquid
       saturation so that one can evaluate the moisture content (i.e. the
       liquid saturation) once and for all */

    cs_property_eval_at_cells(0, mc, cs_shared_liquid_saturation);
  } break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE: {
    cs_field_t *f = cs_field_by_name("liquid_saturation");
    assert(f != nullptr);
    cs_shared_liquid_saturation = f->val;
  } break;

  default:
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Invalid type of hydraulic model.\n",
              __func__);

  } /* End of switch */

  if (cs_shared_liquid_saturation == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Liquid saturation/moisture content is not set.\n",
              __func__);

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    cs_equation_param_t *eqp = cs_equation_get_param(tracer->equation);
    assert(eqp != nullptr);

    if (tracer->finalize_setup != nullptr)
      tracer->finalize_setup(connect, quant, eqp->adv_field, tracer);

  } /* Loop on tracers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the diffusion property related to each tracer equation
 *        The update strategy depends on the soil/tracer features and also
 *        on the hydraulic model.
 *
 * \param[in] ts        pointer to a cs_time_step_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 * \param[in] quant     pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_update_diff_pty(const cs_time_step_t      *ts,
                              const cs_mesh_t           *mesh,
                              const cs_cdo_connect_t    *connect,
                              const cs_cdo_quantities_t *quant)
{
  if (_n_tracers == 0)
    return;

  assert(_tracers != nullptr);

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    if (tracer->update_diff_pty != nullptr)
      tracer->update_diff_pty(tracer, nullptr, ts, mesh, connect, quant);

  } /* Loop on tracers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display the main features related to each tracer
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_log_all(void)
{
  /* Log details about decay chains */

  _log_decay_chains();

  /* Log details about tracer equations */

  cs_log_printf(
    CS_LOG_SETUP, "  * GWF | Number of tracer equations: %d\n", _n_tracers);

  if (_n_tracers == 0)
    return;

  assert(_tracers != nullptr);

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    cs_equation_t *eq     = tracer->equation;
    cs_field_t    *f      = cs_equation_get_field(eq);
    const char    *eqname = cs_equation_get_name(eq);

    cs_log_printf(CS_LOG_SETUP, "\n  ** %s | Variable: %s\n", eqname, f->name);

    if (tracer->model & CS_GWF_TRACER_USER)
      cs_log_printf(CS_LOG_SETUP, "  ** %s | User-defined model\n", eqname);

    else {

      cs_log_printf(CS_LOG_SETUP, "  ** %s | Default model\n", eqname);

      if (tracer->chain_id > -1 && tracer->chain_position_id > -1) {

        cs_gwf_tracer_decay_chain_t *tdc
          = cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
        assert(tdc != nullptr);

        cs_log_printf(CS_LOG_SETUP,
                      "  ** %s | Belongs to a decay chain \"%s\" at position"
                      " %d\n",
                      eqname,
                      tdc->name,
                      tracer->chain_position_id);
      }

      if (tracer->model & CS_GWF_TRACER_PRECIPITATION)
        cs_log_printf(
          CS_LOG_SETUP, "  ** %s | Add precipitation effects\n", eqname);

      if (tracer->model & CS_GWF_TRACER_SORPTION_EK_3_PARAMETERS)
        cs_log_printf(CS_LOG_SETUP,
                      "  ** %s | Use an EK model with 3 parameters\n",
                      eqname);
      else if (tracer->model & CS_GWF_TRACER_SORPTION_EK_5_PARAMETERS)
        cs_log_printf(CS_LOG_SETUP,
                      "  ** %s | Use an EK model with 5 parameters\n",
                      eqname);
    }

  } /* Loop on tracers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state for all tracer equations.
 *         Nothing is done if all equations are unsteady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_compute_steady_all(const cs_mesh_t           *mesh,
                                 const cs_time_step_t      *time_step,
                                 const cs_cdo_connect_t    *connect,
                                 const cs_cdo_quantities_t *cdoq)
{
  if (_n_tracers == 0)
    return;

  assert(_tracers != nullptr);

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    cs_equation_t *eq = tracer->equation;

    if (cs_equation_is_steady(eq)) {

      /* Solve the algebraic system */

      cs_equation_solve_steady_state(mesh, eq);

      if (tracer->update_precipitation != nullptr)
        tracer->update_precipitation(
          tracer, nullptr, time_step, mesh, connect, cdoq);

    } /* Solve this equation which is steady */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for all tracer equations.
 *         Nothing is done if all equations are steady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_compute_all(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq)
{
  if (_n_tracers == 0)
    return;

  assert(_tracers != nullptr);

  bool cur2prev = true;

  /* Loop on tracer equations */

  for (int i = 0; i < _n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = _tracers[i];

    if (tracer == nullptr)
      continue;

    cs_equation_t  *eq = tracer->equation;

    if (!cs_equation_is_steady(eq)) {

      /* Solve the algebraic system. By default, a current to previous operation
         is performed */

      cs_equation_solve(cur2prev, mesh, eq);

      if (tracer->update_precipitation != nullptr)
        tracer->update_precipitation(
          tracer, nullptr, time_step, mesh, connect, cdoq);

      if (tracer->update_decay_chain_st != nullptr)
        tracer->update_decay_chain_st(
          tracer, nullptr, time_step, mesh, connect, cdoq);

    } /* Solve this equation which is unsteady */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add terms to the algebraic system related to a tracer equation
 *         according to the settings.
 *         Case of the default tracer modelling
 *         Rely on the generic function: cs_gwf_tracer_add_terms_t
 *
 * \param[in, out] tracer       pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_default_init_setup(cs_gwf_tracer_t     *tracer)
{
  if (tracer == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " At least one tracer equation has not been set.\n"
              " Please check your settings.");

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;
  cs_equation_param_t  *eqp = cs_equation_get_param(tracer->equation);

  const int  n_soils = cs_gwf_get_n_soils();
  const double  thd = 100*DBL_MIN; /* threshold to avoid a wrong activation */
  const char *eq_name = cs_equation_get_name(tracer->equation);

  int  max_len = 0;
  char *name    = nullptr;

  const int  log_key = cs_field_key_id("log");
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  post_key = cs_field_key_id("post_vis");

  /* Add a diffusion term ? */
  /* ---------------------- */

  bool  do_diffusion = false;
  cs_property_type_t  diff_pty_type = CS_PROPERTY_ISO;

  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    if (fabs(tc->alpha_t[soil_id]) > thd)
      do_diffusion = true, diff_pty_type = CS_PROPERTY_ANISO;

    if (fabs(tc->alpha_l[soil_id]) > thd)
      do_diffusion = true, diff_pty_type = CS_PROPERTY_ANISO;

    if (tc->wmd[soil_id] > thd) do_diffusion = true;

  }

  if (do_diffusion) { /* Add a new diffusion property for this equation */

    int  len = strlen(eq_name) + strlen("_diffusivity") + 1;
    if (len > max_len) {
      max_len = len;
      BFT_REALLOC(name, len, char);
    }
    sprintf(name, "%s_diffusivity", eq_name);

    cs_property_t  *diff_pty = cs_property_add(name, diff_pty_type);

    cs_equation_add_diffusion(eqp, diff_pty);

    /* Create a new field related to this property */

    const int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
    const bool  pty_has_previous = false; /* no previous snapshot */

    int  field_dim = 9;         /* anisotropic case */
    if (cs_property_is_isotropic(diff_pty))
      field_dim = 1;            /* isotropic case */

    tracer->diffusivity = cs_field_create(name,
                                          pty_mask,
                                          c_loc_id,
                                          field_dim,
                                          pty_has_previous);

    cs_field_set_key_int(tracer->diffusivity, cs_field_key_id("log"), 1);

  } /* Has diffusion */

  /* Add a reaction term ? */
  /* --------------------- */

  bool  do_reaction = false;
  if (fabs(tc->decay_coef) > thd) do_reaction = true;

  if (do_reaction) { /* Add a new reaction property for this equation */

    int  len = strlen(eq_name) + strlen("_reaction") + 1;
    if (len > max_len) {
      max_len = len;
      BFT_REALLOC(name, len, char);
    }
    sprintf(name, "%s_reaction", eq_name);

    cs_property_t  *r_pty = cs_property_add(name, CS_PROPERTY_ISO);

    tracer->reaction_id = cs_equation_add_reaction(eqp, r_pty);

  } /* Has reaction */

  if (tracer->model & CS_GWF_TRACER_PRECIPITATION) {

    bool has_previous = false;  /* Not useful up to now */
    int  len = strlen(eq_name) + strlen("_precip") + 1;
    if (len > max_len) {
      max_len = len;
      BFT_REALLOC(name, len, char);
    }
    sprintf(name, "%s_precip", eq_name);

    tc->precip_field = cs_field_create(name,
                                       CS_FIELD_INTENSIVE | CS_FIELD_CDO,
                                       c_loc_id,
                                       1,
                                       has_previous);

    cs_field_set_key_int(tc->precip_field, log_key, 1);
    cs_field_set_key_int(tc->precip_field, post_key, 1);

  }

  if (tracer->chain_position_id > 0) {

    const cs_param_space_scheme_t  space_scheme =
      cs_equation_get_space_scheme(tracer->equation);

    if (space_scheme != CS_SPACE_SCHEME_CDOVB)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Case not handle yet.\n"
                "  Specify a CDO-Vb scheme for the equation \"%s\"\n",
                __func__, eq_name);

    /* The common ancestor has no source term. Add a source term for all
     * children. At this step, mesh has not been preprocessed. A second stage
     * will be done after this one in *_finalize_setup()
     *
     * A CDO-Vb scheme is assumed. The source term is defined on the
     * intersection between the primal and dual cell so that one can manage the
     * properties defined on cells and the variable defined at vertices.
     */

    cs_gwf_tracer_decay_chain_t  *tdc =
      cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
    assert(tdc != nullptr);

    tdc->st_defs[tracer->chain_position_id]
      = cs_equation_add_source_term_by_array(eqp,
                                             nullptr, /* all soils */
                                             cs_flag_dual_cell_byc,
                                             nullptr,
                                             false, /* ownership */
                                             true); /* full length */
  }

  BFT_FREE(name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      adv           pointer to an advection field structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_sat_finalize_setup(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_adv_field_t       *adv,
                                 cs_gwf_tracer_t            *tracer)
{
  if (tracer == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " At least one tracer equation has not been set.\n"
              " Please check your settings.");

  const int  n_soils = cs_gwf_get_n_soils();
  const cs_flag_t  eq_flag = cs_equation_get_flag(tracer->equation);

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  /* Set additional (predefined) fields */

  tc->darcy_velocity_field =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);

  /* We assume that the unsteady term is always activated */

  cs_property_t  *pty = cs_equation_get_time_property(tracer->equation);
  assert(pty != nullptr);

  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_property_def_by_func(pty,
                            z->name,
                            (void *)tracer->context,
                            _get_time_pty4std_sat_tracer,
                            _get_time_pty4std_sat_tracer_cw);

  } /* Loop on soils */

  if (eq_flag & CS_EQUATION_DIFFUSION) { /* Setup the diffusion property */

    assert(tracer->diffusivity != nullptr
           && tracer->diffusivity->val
                != nullptr); /* Should be done previously */

    cs_property_t  *diff_pty =
      cs_equation_get_diffusion_property(tracer->equation);

    if (cs_property_is_isotropic(diff_pty)) {

      tracer->update_diff_pty = nullptr; /* No need. Value is constant */

      for (int soil_id = 0; soil_id < n_soils; soil_id++) {

        cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
        const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

        cs_property_def_iso_by_value(diff_pty, z->name, tc->wmd[soil_id]);

      } /* Loop on soils */

      /* Store the value of the diffusivity inside the field values */

      cs_property_eval_at_cells(0, diff_pty, tracer->diffusivity->val);

    }
    else {

      tracer->update_diff_pty = _update_diff_tensor;
      cs_property_def_by_field(diff_pty, tracer->diffusivity);

    }

  } /* Has diffusion */

  if (eq_flag & CS_EQUATION_REACTION) { /* Setup the reaction property */

    for (int soil_id = 0; soil_id < n_soils; soil_id++) {

      const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
      const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      cs_property_t  *r_pty =
        cs_equation_get_reaction_property(tracer->equation,
                                          tracer->reaction_id);

      if (r_pty != nullptr) /* The default reaction property is defined */
        cs_property_def_by_func(r_pty,
                                z->name,
                                (void *)tracer->context,
                                _get_reaction_pty4std_sat_tracer,
                                _get_reaction_pty4std_sat_tracer_cw);

    } /* Loop on soils */

  } /* Has reaction */

  /* Precipitation modelling */

  if (tracer->model & CS_GWF_TRACER_PRECIPITATION)
    _add_precipitation(connect, quant, tracer);

  /* Handle the source term in the decay chain */

  if (tracer->chain_position_id > 0) {

    /* Second step: Associate the array to the source term */

    const cs_adjacency_t  *c2v = connect->c2v;

    double *st_values = nullptr;
    BFT_MALLOC(st_values, c2v->idx[quant->n_cells], double);
    cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

    cs_gwf_tracer_decay_chain_t  *tdc =
      cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
    assert(tdc != nullptr);

    cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];

    cs_xdef_array_set_values(st_def,
                             true,     /* transfer ownership */
                             st_values);

    cs_xdef_array_set_adjacency(st_def, c2v);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a standard tracer equation in case of
 *         an unsaturated flow model
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      adv           pointer to an advection field structure
 * \param[in, out] tracer        pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_unsat_finalize_setup(const cs_cdo_connect_t      *connect,
                                   const cs_cdo_quantities_t   *quant,
                                   const cs_adv_field_t        *adv,
                                   cs_gwf_tracer_t             *tracer)
{
  if (tracer == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " At least one tracer equation has not been set.\n"
              " Please check your settings.");

  const int  n_soils = cs_gwf_get_n_soils();
  const cs_flag_t  eq_flag = cs_equation_get_flag(tracer->equation);

  cs_gwf_tracer_default_context_t *tc
    = (cs_gwf_tracer_default_context_t *)tracer->context;

  /* Set additional (pre-defined) fields */

  tc->darcy_velocity_field =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);

  /* We assume that the unsteady term is always activated */

  cs_property_t  *pty = cs_equation_get_time_property(tracer->equation);
  assert(pty != nullptr);

  for (int soil_id = 0; soil_id < n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_property_def_by_func(pty,
                            z->name,
                            (void *)tracer->context,
                            _get_time_pty4std_tracer,
                            _get_time_pty4std_tracer_cw);

  } /* Loop on soils */

  if (eq_flag & CS_EQUATION_DIFFUSION) { /* Setup the diffusion property */

    assert(tracer->diffusivity != nullptr
           && tracer->diffusivity->val
                != nullptr); /* Should be done previously */

    cs_property_t  *diff_pty =
      cs_equation_get_diffusion_property(tracer->equation);

    cs_property_def_by_field(diff_pty, tracer->diffusivity);

    if (cs_property_is_isotropic(diff_pty))
      tracer->update_diff_pty = _update_diff_value;
    else
      tracer->update_diff_pty = _update_diff_tensor;

  } /* Has diffusion */

  if (eq_flag & CS_EQUATION_REACTION) { /* Setup the reaction property */

    for (int soil_id = 0; soil_id < n_soils; soil_id++) {

      const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);
      const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      cs_property_t  *r_pty =
        cs_equation_get_reaction_property(tracer->equation,
                                          tracer->reaction_id);

      if (r_pty != nullptr) /* The default reaction property is defined */
        cs_property_def_by_func(r_pty,
                                z->name,
                                (void *)tracer->context,
                                _get_reaction_pty4std_tracer,
                                _get_reaction_pty4std_tracer_cw);

    } /* Loop on soils */

  } /* Has reaction */

  /* Precipitation modelling */

  if (tracer->model & CS_GWF_TRACER_PRECIPITATION)
    _add_precipitation(connect, quant, tracer);

  /* Handle the source term in the decay chain */

  if (tracer->chain_position_id > 0) {

    /* Second step: Associate the array to the source term */

    const cs_adjacency_t  *c2v = connect->c2v;

    double *st_values = nullptr;
    BFT_MALLOC(st_values, c2v->idx[quant->n_cells], double);
    cs_array_real_fill_zero(c2v->idx[quant->n_cells], st_values);

    cs_gwf_tracer_decay_chain_t  *tdc =
      cs_gwf_tracer_decay_chain_by_id(tracer->chain_id);
    assert(tdc != nullptr);

    cs_xdef_t  *st_def = tdc->st_defs[tracer->chain_position_id];

    cs_xdef_array_set_values(st_def,
                             true,     /* transfer ownership */
                             st_values);

    cs_xdef_array_set_adjacency(st_def, c2v);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the integral of the tracer concentration over a given set of
 *        cells. This integral gives the number of moles of tracer inside the
 *        related volume. Moreover, it turns out to be exact for linear
 *        functions. A parallel operation (a sum reduction) is performed
 *        inside this function.
 *
 * \param[in]    connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]    cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]    tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]    z_name    name of the volumic zone where the integral is done
 *                         (if nullptr or "" all cells are considered)
 *
 * \return the value of the integral (number of moles in the zone)
 *         parallel synchronization is done
 */
/*----------------------------------------------------------------------------*/

double
cs_gwf_tracer_integrate(const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_gwf_tracer_t      *tracer,
                        const char                 *z_name)
{
  if (tracer == nullptr)
    return 0;

  const int  z_id = cs_volume_zone_id_by_name(z_name);
  const cs_zone_t  *zone = cs_volume_zone_by_id(z_id);

  if (tracer->model & CS_GWF_TRACER_USER)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of tracer.\n",
              __func__);

  double  results[2] = {0, 0};

  if (tracer->integrate == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: No function is set to integrate the tracer \"%s\"\n",
              __func__,
              (tracer->equation == nullptr)
                ? ""
                : cs_equation_get_name(tracer->equation));

  tracer->integrate(connect, cdoq, tracer->equation, zone, tracer->context,
                    results);

  return results[0] + results[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the quantity of tracer using the integral of the tracer
 *        concentration over a given set of cells. Two terms are computed: one
 *        for the quantity of moles inside the liquid phase and another one for
 *        the quantity of tracer inside the precipitation state (in moles). A
 *        parallel operation (a sum reduction) is performed inside this
 *        function.
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]      z_name    name of the volume zone where the integral is
 *                           done (if nullptr or "" all cells are considered)
 * \param[in, out] results   array of values. [0]= the quantity of moles
 *                           in the liquid phase, [1]= the quantity of
 *                           moles inside the precipitation state
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tracer_integrate_by_terms(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *cdoq,
                                 const cs_gwf_tracer_t      *tracer,
                                 const char                 *z_name,
                                 double                      results[])
{
  if (results == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid result array.", __func__);

  results[0] = 0., results[1] = 0.;

  if (tracer == nullptr)
    return;

  if (tracer->model & CS_GWF_TRACER_USER)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of tracer.\n",
              __func__);

  const int  z_id = cs_volume_zone_id_by_name(z_name);
  const cs_zone_t  *zone = cs_volume_zone_by_id(z_id);

  if (tracer->integrate == nullptr)
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: No function is set to integrate the tracer \"%s\"\n",
              __func__,
              (tracer->equation == nullptr)
                ? ""
                : cs_equation_get_name(tracer->equation));

  tracer->integrate(connect, cdoq, tracer->equation, zone, tracer->context,
                    results);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a decay chain structure to manage several linked tracers
 *
 * \param[in] n_tracers    number of tracers equations
 * \param[in] chain_name   name of the decay chain
 * \param[in] unit         type of unit used in the tracer equations
 *
 * \return a pointer to the new cs_gwf_tracer_decay_chain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_create_decay_chain(int                      n_tracers,
                                 const char              *chain_name,
                                 cs_gwf_tracer_unit_t     unit)
{
  if (n_tracers < 1)
    return nullptr;

  int  tdc_id = _n_decay_chains;

  cs_gwf_tracer_decay_chain_t *tdc = nullptr;

  BFT_MALLOC(tdc, 1, cs_gwf_tracer_decay_chain_t);

  tdc->n_tracers = n_tracers;
  tdc->unit = unit;

  size_t  len = strlen(chain_name);
  BFT_MALLOC(tdc->name, len + 1, char);
  strncpy(tdc->name, chain_name, len + 1); /* Last character is '\0' */

  tdc->id = tdc_id;

  BFT_MALLOC(tdc->tracers, n_tracers, cs_gwf_tracer_t *);
  BFT_MALLOC(tdc->st_defs, n_tracers, cs_xdef_t *);

  for (int i = 0; i < n_tracers; i++) {
    tdc->tracers[i] = nullptr;
    tdc->st_defs[i] = nullptr;
  }

  /* Update the array storing all decay chains */

  _n_decay_chains += 1;
  BFT_REALLOC(_decay_chains, _n_decay_chains, cs_gwf_tracer_decay_chain_t *);
  _decay_chains[tdc_id] = tdc;

  return tdc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the decay chain structure associated to the given id
 *        If not found, it returns the nullptr pointer.
 *
 * \param[in] id   id of the decay chain to retrieve
 *
 * \return a pointer to a new cs_gwf_tracer_decay_chain_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_decay_chain_by_id(int        id)
{
  if (_n_decay_chains < 1)
    return nullptr;
  if (id < 0)
    return nullptr;

  assert(_decay_chains != nullptr);

  return  _decay_chains[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the decay chain structure associated to the name given as
 *        parameter. If not found, it returns the nullptr pointer.
 *
 * \param[in] chain_name   name of the decay chain
 *
 * \return a pointer to a new cs_gwf_tracer_decay_chain_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_tracer_decay_chain_by_name(const char      *chain_name)
{
  if (_n_decay_chains < 1)
    return nullptr;
  if (chain_name == nullptr)
    return nullptr;

  cs_gwf_tracer_decay_chain_t *tdc = nullptr;

  size_t  len = strlen(chain_name);
  for (int i = 0; i < _n_decay_chains; i++) {

    cs_gwf_tracer_decay_chain_t  *_c = _decay_chains[i];
    if (_c == nullptr)
      continue;

    if (strlen(_c->name) == len)
      if (strcmp(_c->name, chain_name) == 0)
        return _c;

  }

  return tdc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the tracer structure for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a nullptr
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_gwf_tracer_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_decay_chain_get_tracer(cs_gwf_tracer_decay_chain_t  *tdc,
                                     int                           id)
{
  if (tdc == nullptr)
    return nullptr;

  if (id < 0 || id >= tdc->n_tracers)
    return nullptr;

  cs_gwf_tracer_t  *tracer = tdc->tracers[id];

  assert(tracer != nullptr); /* Sanity checks */
  assert(tracer->chain_position_id == id);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the equation structure for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a nullptr
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_equation_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_gwf_tracer_decay_chain_get_equation(cs_gwf_tracer_decay_chain_t  *tdc,
                                       int                           id)
{
  if (tdc == nullptr)
    return nullptr;

  if (id < 0 || id >= tdc->n_tracers)
    return nullptr;

  cs_gwf_tracer_t  *tracer = tdc->tracers[id];
  assert(tracer != nullptr);
  assert(tracer->chain_position_id == id);

  return tracer->equation;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the equation parameters for the tracer at the position "id"
 *        in the decay chain structure. If "id" is not valid, then a nullptr
 *        pointer is returned.
 *
 * \param[in] tdc   pointer to a decay chain structure
 * \param[in] id    position of the tracer in the decay chain
 *
 * \return a pointer to a cs_equation_param_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_gwf_tracer_decay_chain_get_equation_param(cs_gwf_tracer_decay_chain_t  *tdc,
                                             int                           id)
{
  if (tdc == nullptr)
    return nullptr;

  if (id < 0 || id >= tdc->n_tracers)
    return nullptr;

  cs_gwf_tracer_t  *tracer = tdc->tracers[id];
  assert(tracer != nullptr);
  assert(tracer->chain_position_id == id);

  return cs_equation_get_param(tracer->equation);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
