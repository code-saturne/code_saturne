/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
 *============================================================================*/

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
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_post.h"
#include "cs_source_term.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_navsto.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_NAVSTO_DBG      0
#define CS_CDOFB_NAVSTO_MODULO  10

/*! \struct cs_cdofb_navsto_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *  vector-valued unknows
 */

typedef struct {

  /*!
   * @}
   * @name Context structures for each related equations
   * @{
   */

  /* \var vecteq_context
   * Structure storing how to build the vector equation
   */

  cs_cdofb_vecteq_t   *vecteq_context;

  /* \var scaleq_context
   * Structure storing how to build the vector equation (optional).
   * It depends on the coupling algo.
   */

  cs_cdofb_scaleq_t   *scaleq_context;

  /*!
   * @}
   * @name Arrays storing face unknowns
   * @{
   *
   */

  /* \var face_velocity
   * Degrees of freedom for the velocity at faces
   */

  cs_real_t  *face_velocity;

  /* \var face_pressure
   * Degrees of freedom for the pressure at faces. Not always allocated.
   * It depends on the type of algorithm used to couple the Navier-Stokes
   * system.
   */

  cs_real_t  *face_pressure;

  /*!
   * @}
   * @name Performance monitoring
   * @{
   *
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   */

  cs_timer_counter_t  timer; /*!< Cumulated elapsed time for building and
                              *   solving the Navier--Stokes system */

  /*! @} */

} cs_cdofb_navsto_t;

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_assembler_t  *cs_shared_scal_ma;
static const cs_matrix_structure_t  *cs_shared_scal_ms;
static const cs_matrix_assembler_t  *cs_shared_vect_ma;
static const cs_matrix_structure_t  *cs_shared_vect_ms;

static cs_cdofb_navsto_t  *cs_cdofb_navsto_context = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdofb_navsto_t structure by default
 *
 * \param[in] nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_navsto_t strcuture
 */
/*----------------------------------------------------------------------------*/

static cs_cdofb_navsto_t *
_create_navsto_context(const cs_navsto_param_t  *nsp)
{
  cs_cdofb_navsto_t  *nssc = NULL;

  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n",
              __func__);

  BFT_MALLOC(nssc, 1, cs_cdofb_navsto_t);

  nssc->vecteq_context = NULL;
  nssc->scaleq_context = NULL;

  nssc->face_velocity = NULL;
  nssc->face_pressure = NULL;

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(nssc->timer);

  return nssc;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members for CDO face-based
 *         schemes
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  sms         pointer to a cs_matrix_structure_t structure (scalar)
 * \param[in]  vms         pointer to a cs_matrix_structure_t structure (vector)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_common(const cs_cdo_quantities_t     *quant,
                            const cs_cdo_connect_t        *connect,
                            const cs_time_step_t          *time_step,
                            const cs_matrix_structure_t   *sms,
                            const cs_matrix_structure_t   *vms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /*
    Matrix structure related to the algebraic system for scalar-valued equation
  */
  cs_shared_scal_ms = sms;

  /*
    Matrix structure related to the algebraic system for vector-valued equation
  */
  cs_shared_vect_ms = vms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of a
 *         Uzawa-Augmented Lagrangian approach
 *
 * \param[in] nsp        pointer to a cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_uzawa_context(const cs_navsto_param_t     *nsp,
                                   const void                  *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  /* Navier-Navsto scheme context (NSSC) */
  cs_cdofb_navsto_t  *nssc = _create_navsto_context(nsp);

  const cs_navsto_coupling_uzawa_t  *nsc =
    (const cs_navsto_coupling_uzawa_t  *)nsc_input;

  cs_cdofb_navsto_context = nssc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of an
 *         Artificial Compressibility approach
 *
 * \param[in] nsp    pointer to a cs_navsto_param_t structure
 * \param[in] nsc    pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_ac_context(const cs_navsto_param_t   *nsp,
                                const void                *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  /* Navier-Navsto scheme context (NSSC) */
  cs_cdofb_navsto_t  *nssc = _create_navsto_context(nsp);

  const cs_navsto_coupling_ac_t  *nsc =
    (const cs_navsto_coupling_ac_t *)nsc_input;

  cs_cdofb_navsto_context = nssc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_navsto_t structure storing in the case of an
 *         incremental Projection approach
 *
 * \param[in] nsp        pointer to a cs_navsto_param_t structure
 * \param[in] nsc_input  pointer to a cs_navsto_coupling_uzawa_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_proj_context(const cs_navsto_param_t    *nsp,
                                  const void                 *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);

  /* Navier-Navsto scheme context (NSSC) */
  cs_cdofb_navsto_t  *nssc = _create_navsto_context(nsp);

  const cs_navsto_coupling_projection_t  *nsc =
    (const cs_navsto_coupling_projection_t  *)nsc_input;

  cs_cdofb_navsto_context = nssc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_navsto_t structure
 *
 * \param[in]      nsp        pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_free_context(const cs_navsto_param_t      *nsp)
{
  cs_cdofb_navsto_t  *nssc = cs_cdofb_navsto_context;

  if (nssc == NULL)
    return;

  /* Free temporary buffers */

  BFT_FREE(nssc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         a Uzawa-Lagrangian Augmented approach
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_uzawa_compute(const cs_mesh_t              *mesh,
                              double                        dt_cur,
                              const cs_navsto_param_t      *nsp,
                              void                         *nsc_input)
{
  CS_UNUSED(dt_cur);

  cs_cdofb_navsto_t  *nssc = cs_cdofb_navsto_context;
  cs_navsto_coupling_uzawa_t  *nscc = (cs_navsto_coupling_uzawa_t  *)nsc_input;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(nssc->timer), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         an Artificial Compressibility approach.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_ac_compute(const cs_mesh_t              *mesh,
                           double                        dt_cur,
                           const cs_navsto_param_t      *nsp,
                           void                         *nsc_input)
{
  CS_UNUSED(dt_cur);

  cs_cdofb_navsto_t  *nssc = cs_cdofb_navsto_context;
  cs_navsto_coupling_ac_t  *nscc = (cs_navsto_coupling_ac_t *)nsc_input;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(nssc->timer), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the Navier-Stokes system with a CDO face-based scheme using
 *         an incremental correction-projection approach
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      dt_cur      current value of the time step
 * \param[in]      nsp         pointer to a cs_navsto_param_t structure
 * \param[in, out] nsc_input   Navier-Stokes coupling context: pointer to a
 *                             structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_proj_compute(const cs_mesh_t              *mesh,
                             double                        dt_cur,
                             const cs_navsto_param_t      *nsp,
                             void                         *nsc_input)
{
  CS_UNUSED(dt_cur);

  cs_cdofb_navsto_t  *nssc = cs_cdofb_navsto_context;
  cs_navsto_coupling_projection_t  *nscc =
    (cs_navsto_coupling_projection_t  *)nsc_input;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(nssc->timer), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdofb_navsto_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_update_field(const cs_real_t              *solu,
                             const cs_real_t              *rhs,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *data,
                             cs_real_t                    *field_val)
{
  CS_UNUSED(rhs);

  cs_cdofb_navsto_t  *eqc = (cs_cdofb_navsto_t *)data;
  cs_timer_t  t0 = cs_timer_time();


  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_navsto_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname); // avoid a compilation warning
  CS_UNUSED(eqp);

  char *postlabel = NULL;
  cs_timer_t  t0 = cs_timer_time();


  // TODO

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
