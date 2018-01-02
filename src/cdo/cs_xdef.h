#ifndef __CS_XDEF_H__
#define __CS_XDEF_H__

/*============================================================================
 * Routines to handle extended definitions of quantities (cs_xdef_t structures)
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_param.h"
#include "cs_quadrature.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_XDEF_BY_ANALYTIC_FUNCTION,
  CS_XDEF_BY_ARRAY,
  CS_XDEF_BY_FIELD,
  CS_XDEF_BY_FUNCTION,
  CS_XDEF_BY_QOV,
  CS_XDEF_BY_TIME_FUNCTION,
  CS_XDEF_BY_VALUE,
  CS_N_XDEF_TYPES

} cs_xdef_type_t;

typedef enum {

  CS_XDEF_SUPPORT_TIME,      /* support for time step description */
  CS_XDEF_SUPPORT_BOUNDARY,  /* zones attached to boundary faces */
  CS_XDEF_SUPPORT_VOLUME,    /* zones attached to cells */
  CS_N_XDEF_SUPPORTS

} cs_xdef_support_t;

typedef struct {

  int                    dim;      /* dimension of the values attached to
                                      this description */
  cs_xdef_type_t         type;     /* type of definition */

  int                    z_id;     /* zone id related to this definition */
  cs_xdef_support_t      support;  /* support for this definition */

  cs_flag_t              state;    /* steady, uniform, cellwise... */
  cs_flag_t              meta;     /* Metadata about the descitption. These
                                      metadata may vary according to the object
                                      on which the description applies */
  cs_quadrature_type_t   qtype;    /* type of quadrature to use for evaluating
                                      the description */

  void                  *input;    /* pointer to metadat to complete
                                      the description */

} cs_xdef_t;

/* Input structure when an array is used for the definition */
typedef struct {

  int               stride;  /* array stride */
  cs_flag_t         loc;     /* flag to know where are defined array values */
  cs_real_t        *values;  /* array values */
  const cs_lnum_t  *index;   /* optional index for accessing to the values */

} cs_xdef_array_input_t;

/* Input structure when an analytic function is used for the definition */
typedef struct {

  void                *input;  // NULL or pointer to a structure cast on-the-fly
  cs_analytic_func_t  *func;   // function to call

} cs_xdef_analytic_input_t;

/* Input structure when an analytic function is used for the definition */
typedef struct {

  void                *input;  // NULL or pointer to a structure cast on-the-fly
  cs_timestep_func_t  *func;   // function to call

} cs_xdef_timestep_input_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on volumic
 *         elements
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_volume_create(cs_xdef_type_t    type,
                      int               dim,
                      int               z_id,
                      cs_flag_t         state,
                      cs_flag_t         meta,
                      void             *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure based on boundary
 *         elements
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_boundary_create(cs_xdef_type_t    type,
                        int               dim,
                        int               z_id,
                        cs_flag_t         state,
                        cs_flag_t         meta,
                        void             *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_xdef_t structure for setting the
 *         time step
 *
 * \param[in]  type       type of definition
 * \param[in]  dim        dimension of the values to define
 * \param[in]  z_id       volume zone id
 * \param[in]  state      flag to know if this uniform, cellwise, steady...
 * \param[in]  meta       metadata associated to this description
 * \param[in]  input      pointer to a structure storing the parameters (cast
 *                        on-the-fly according to the type of definition)
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_timestep_create(cs_xdef_type_t             type,
                        cs_flag_t                  state,
                        cs_flag_t                  meta,
                        void                      *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_xdef_t structure
 *
 * \param[in, out] d    pointer to a cs_xdef_t structure
 *
 * \return NULL
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_xdef_free(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in, out]  d       pointer to a cs_xdef_t structure
 * \param[in]       qtype   type of quadrature
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_set_quadrature(cs_xdef_t              *d,
                       cs_quadrature_type_t    qtype);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the type of quadrature to use for evaluating the given
 *         description
 *
 * \param[in]  d       pointer to a cs_xdef_t structure
 *
 * \return the type of quadrature
 */
/*----------------------------------------------------------------------------*/

cs_quadrature_type_t
cs_xdef_get_quadrature(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_xdef_type_t
cs_xdef_get_type(const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the flag dedicated to the state
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 *
 * \return the value of the flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_xdef_get_state_flag(const cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a a cs_xdef_t structure
 *
 * \param[in] d    pointer to a cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_log(cs_xdef_t     *d);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_H__ */
