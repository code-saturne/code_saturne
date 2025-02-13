#ifndef __CS_LES_FILTER_H__
#define __CS_LES_FILTER_H__

/*============================================================================
 * Filters for dynamic models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute filters for dynamic models for a scalar.
 *
 * This function deals with the standard or extended neighborhood.
 *
 * parameters:
 *   val     <->  array of values to filter
 *   f_val   -->  array of filtered values
 *----------------------------------------------------------------------------*/

void
cs_les_filter_scalar(const cs_real_t  val[],
                     cs_real_t        f_val[]);

END_C_DECLS

#ifdef __cplusplus

/*----------------------------------------------------------------------------
 * Compute filters for dynamic models for vectors.
 *
 * This function deals with the standard or extended neighborhood.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   val     <->  array of values to filter
 *   f_val   -->  array of filtered values
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_les_filter_strided(const cs_real_t  val[][stride],
                      cs_real_t        f_val[][stride]);

/*----------------------------------------------------------------------------*/

#endif /* cplusplus */

#endif /* __CS_LES_FILTER_H__ */
