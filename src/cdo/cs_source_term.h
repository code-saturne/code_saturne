#ifndef __CS_SOURCE_TERM_H__
#define __CS_SOURCE_TERM_H__

/*============================================================================
 * Functions and structures to deal with source term computation
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Types of source terms */
typedef enum {

  CS_SOURCE_TERM_GRAVITY,  // specific treatment
  CS_SOURCE_TERM_HEADLOSS, // specific treatment (not implemented yet)
  CS_SOURCE_TERM_MASS,     // specific treatment (not implemented yet)
  CS_SOURCE_TERM_USER,     // user-defined
  CS_N_SOURCE_TERM_TYPES

} cs_source_term_type_t;

/* Indicate where the reduction of source terms has to be done */
typedef enum {

  CS_SOURCE_TERM_REDUC_DUAL,   // reduction on dual entities
  CS_SOURCE_TERM_REDUC_PRIM,   // reduction on primal entities
  CS_N_SOURCE_TERM_REDUCTIONS

} cs_source_term_reduction_t;

typedef struct _cs_source_term_t cs_source_term_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                   const cs_cdo_connect_t       *connect,
                                   const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a cs_source_term_t structure
 *
 * \param[in] st_name     name of the related source term
 * \param[in] ml_id       id of the related mesh location
 * \param[in] st_type     type of source term to create
 * \param[in] red_type    type of reduction to apply
 * \param[in] var_type    type of variables (scalar, vector, tensor...)
 *
 * \return a pointer to a new allocated source term structure
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_source_term_create(const char                  *name,
                      int                          ml_id,
                      cs_source_term_type_t        st_type,
                      cs_source_term_reduction_t   red_type,
                      cs_param_var_type_t          var_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_source_term_free(cs_source_term_t   *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for computing the source term
 *
 * \param[in, out]  st          pointer to a cs_source_term_t structure
 * \param[in]       quad_type   type of quadrature to use
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_quadrature(cs_source_term_t  *st,
                              cs_quadra_type_t   quad_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set where the reduction has to be applied forn defining source term
 *
 * \param[in, out]  st        pointer to a cs_source_term_t structure
 * \param[in]       red_type  type of reduction to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_reduction(cs_source_term_t             *st,
                             cs_source_term_reduction_t    red_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the name related to a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return the name of the source term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_source_term_get_name(const cs_source_term_t   *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the type of reduction applied to a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return the type of reduction
 */
/*----------------------------------------------------------------------------*/

cs_source_term_reduction_t
cs_source_term_get_reduction(const cs_source_term_t   *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the content of a cs_source_term_t structure
 *
 * \param[in] eqname  name of the related equation
 * \param[in] st      pointer to a cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_summary(const char               *eqname,
                       const cs_source_term_t   *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic way to define the value of a cs_source_term_t structure
 *
 * \param[in, out]  pty     pointer to a cs_source_term_t structure
 * \param[in]       val     accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_value(cs_source_term_t    *st,
                            const char          *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an analytic function
 *
 * \param[in, out]  st      pointer to a cs_source_term_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_analytic(cs_source_term_t      *st,
                               cs_analytic_func_t    *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an array of values
 *
 * \param[in, out]  st       pointer to a cs_source_term_t structure
 * \param[in]       desc     description of the main feature of this array
 * \param[in]       array    pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_array(cs_source_term_t    *st,
                            cs_desc_t            desc,
                            cs_real_t           *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term
 *
 * \param[in]      dof_desc   description of the associated DoF
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute(cs_desc_t                     dof_desc,
                       const cs_source_term_t       *source,
                       double                       *p_values[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOURCE_TERM_H__ */
