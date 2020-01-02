#ifndef __CS_FIELD_DEFAULT_H__
#define __CS_FIELD_DEFAULT_H__

/*============================================================================
 * Field utility functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

END_C_DECLS

#endif /* __CS_FIELD_DEFAULT_H__ */
