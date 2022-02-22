#ifndef __CS_LES_FILTER_H__
#define __CS_LES_FILTER_H__

/*============================================================================
 * Filters for dynamic models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Compute filters for dynamic models.
 *
 * This function deals with the standard or extended neighborhood.
 *
 * parameters:
 *   stride  <--  stride of array to filter
 *   val     <->  array of values to filter
 *   f_val   -->  array of filtered values
 *----------------------------------------------------------------------------*/

void
cs_les_filter(int        stride,
              cs_real_t  val[],
              cs_real_t  f_val[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LES_FILTER_H__ */
