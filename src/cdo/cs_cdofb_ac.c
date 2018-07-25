/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it with an artificial compressibility algorithm
 *============================================================================*/

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_bc.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_navsto.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_post.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_ac.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_ac.c
 *
 * \brief Build an algebraic CDO face-based system for the Navier-Stokes
 * equations and solved it with an artificial compressibility algorithm
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_ac_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         vector-valued unknowns
 */

typedef struct {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_ac_t (owned by \ref cs_navsto_system_t)
   *  containing the settings related to an artificial compressibility (AC)
   *  algorithm or vector penalty projection (VPP)
   */

  cs_navsto_ac_t   *coupling_context;

  /*!
   * @name Main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   */

  /*! \var velocity
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the velocity
   */

  cs_field_t  *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t  *pressure;

  /*!
   * @}
   * @name Arrays storing face unknowns
   * @{
   */

  /*! \var face_velocity
   *  Degrees of freedom for the velocity at faces
   */

  cs_real_t   *face_velocity;

  /*!
   * @}
   * @name Parameters of the algorithm
   * Easy access to useful features and parameters of the algorithm
   * @{
   */

  /*! \var is_zeta_uniform
   *  Bool telling if the auxiliary parameter zeta is uniform. Not always
   *  necessary: zeta is tipically used in Artificial Compressibility algos
   */

  bool is_zeta_uniform;

  /*!
   * @}
   * @name Performance monitoring
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   * @{
   */

  /*! \var timer
   *  Cumulated elapsed time for building and solving the Navier--Stokes system
   */
  cs_timer_counter_t  timer;

  /*! @} */

} cs_cdofb_ac_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_AC_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_scal_ms;
static const cs_matrix_structure_t  *cs_shared_vect_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_init_common(const cs_cdo_quantities_t     *quant,
                        const cs_cdo_connect_t        *connect,
                        const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /*
    Matrix structure related to the algebraic system for scalar-valued equation
  */
  cs_shared_scal_ms = cs_cdofb_scaleq_matrix_structure();

  /*
    Matrix structure related to the algebraic system for vector-valued equation
  */
  cs_shared_vect_ms = cs_cdofb_vecteq_matrix_structure();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_cdofb_ac_t structure
 *
 * \param[in] nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a \ref cs_navsto_ac_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_ac_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_init_scheme_context(const cs_navsto_param_t     *nsp,
                                void                        *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  /* Cast the coupling context (CC) */
  cs_navsto_ac_t  *cc = (cs_navsto_ac_t  *)nsc_input;

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_ac_t  *sc = NULL;

  BFT_MALLOC(sc, 1, cs_cdofb_ac_t);

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");

  /* Only one vector equation */
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_context = (cs_cdofb_vecteq_t *)mom_eq->scheme_context;

  sc->face_velocity = mom_context->face_values;

  sc->is_zeta_uniform = cs_property_is_uniform(cc->zeta);

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_ac_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_ac_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t  *)scheme_context;

  if (sc == NULL)
    return sc;

  /* Free temporary buffers */
  if (sc->face_velocity != NULL) BFT_FREE(sc->face_velocity);

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the velocity values
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_init_velocity(const cs_navsto_param_t     *nsp,
                          void                        *scheme_context)
{
  CS_UNUSED(nsp);
  CS_UNUSED(scheme_context);

  /* Nothing to do. All is already done during the initialization of the
     momentum equation */
  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values
 *
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_init_pressure(const cs_navsto_param_t     *nsp,
                          void                        *scheme_context)
{
  /* Sanity checks */
  assert(nsp != NULL && scheme_context != NULL);

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;

  /* Initial conditions for the pressure */
  if (nsp->n_pressure_ic_defs > 0) {

    assert(nsp->pressure_ic_defs != NULL);
    assert(sc != NULL);

    const cs_time_step_t *ts = cs_shared_time_step;
    cs_field_t *pr  = sc->pressure;
    cs_real_t  *values = pr->val;

    const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;
    const cs_flag_t   dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;
    const cs_real_t  t_cur = ts->t_cur;

    for (int def_id = 0; def_id < nsp->n_pressure_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      cs_xdef_t  *def = nsp->pressure_ic_defs[def_id];

      /* Initialize face-based array */
      switch (def->type) {

        /* Evaluating the integrals: the averages will be taken care of at the
         * end when ensuring zero-mean valuedness */
      case CS_XDEF_BY_VALUE:
        cs_evaluate_density_by_value(dof_flag, def, values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          /* Forcing BARY so that it is equivalent to DeRham (JB?)*/
          cs_xdef_set_quadrature(def, CS_QUADRATURE_BARY);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          /* Restoring the original */
          cs_xdef_set_quadrature(def, nsp->qtype);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Incompatible reduction for the field %s.\n"),
                    __func__, pr->name);

        }  /* Switch on possible reduction types */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Incompatible way to initialize the field %s.\n"),
                  __func__, pr->name);
        break;

      }  /* Switch on possible type of definition */

    }  /* Loop on definitions */

  /* We should ensure that the mean of the pressure is zero. Thus we compute
   * it and subtract it from every value.
   * NOTES:
   *  - It could be useful to stored this average somewhere
   *  - The procedure is not optimized (we can avoid setting the average if
   *    it's a value), but it is the only way to allow multiple definitions
   *    and definitions that do not cover all the domain. Moreover, we need
   *    information (e.g. cs_cdo_quant) which we do not know here */
    cs_cdofb_navsto_ensure_zero_mean_and_avg(values, 1);

  }  /* Not the default */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         a Ac-Lagrangian Augmented approach.
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in]      dt_cur          current value of the time step
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_ac_compute(const cs_mesh_t              *mesh,
                    const cs_navsto_param_t      *nsp,
                    double                        dt_cur,
                    void                         *scheme_context)
{
  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */
  CS_UNUSED(dt_cur);
  CS_UNUSED(mesh);
  CS_UNUSED(nsp);
  CS_UNUSED(sc);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values of the velocity on the faces
 *
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 *
 * \return a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_ac_get_face_velocity(void    *scheme_context)
{
  if (scheme_context == NULL)
    return NULL;

  cs_cdofb_ac_t  *sc = (cs_cdofb_ac_t *)scheme_context;

  return sc->face_velocity;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
