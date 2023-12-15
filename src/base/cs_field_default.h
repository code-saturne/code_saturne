#ifndef __CS_FIELD_DEFAULT_H__
#define __CS_FIELD_DEFAULT_H__

/*============================================================================
 * Field utility functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_defs.h"
#include "cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a field shared between CDO and legacy schemes. This field is
 *         related to a general solved variable, with default options.
 *
 * \param[in]  name          field name
 * \param[in]  label         field default label, or empty
 * \param[in]  location_id   id of associated location
 * \param[in]  dim           field dimension
 * \param[in]  has_previous  no if lower than 1
 *
 * \return  newly defined field id
 */
/*----------------------------------------------------------------------------*/

int
cs_variable_cdo_field_create(const char  *name,
                             const char  *label,
                             int          location_id,
                             int          dim,
                             int          has_previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add field defining a general solved variable, with default options.
 *
 * \param[in]  name         field name
 * \param[in]  label        field default label, or empty
 * \param[in]  location_id  id of associated location
 * \param[in]  dim          field dimension
 *
 * \return  newly defined field id
 */
/*----------------------------------------------------------------------------*/

int
cs_variable_field_create(const char  *name,
                         const char  *label,
                         int          location_id,
                         int          dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Access a field's equation parameters.
 *
 * If the equation parameters were never initialized, they will be initialized
 * based on the current defaults.
 *
 * If the field does not have associated equation paremeters (i.e. is not
 * a variable field or is a CDO field (which is referenced by but does not
 * directly reference equations), NULL is returned.
 *
 * \param[in, out]  f  pointer to associated field
 *
 * \return  pointer to field's equation parameters, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_field_get_equation_param(cs_field_t  *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Access a field's equation parameters for read only.
 *
 * If the equation parameters were never initialized, the current default
 * parameters will be returned instead.
 *
 * If the field does not have associated equation parameters (i.e. is not
 * a variable field or is a CDO field (which is referenced by but does not
 * directly reference equations), NULL is returned.
 *
 * \param[in]  f  pointer to associated field
 *
 * \return  const-qualified pointer to field's equation parameters, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_equation_param_t *
cs_field_get_equation_param_const(const cs_field_t  *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief For a given field, returns field defined as its variance, if present.
 *
 * \param[in]  f  field
 *
 * \return  pointer to matching variance (variable) field, or NULL.
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_field_get_variance(const cs_field_t  *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and map boundary condition coefficients for all
 *        variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_build_bc_codes_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deallocate and unmap boundary condition coefficients for all
 *        variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_free_bc_codes_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize all field BC coefficients.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_map_and_init_bcs(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_DEFAULT_H__ */
