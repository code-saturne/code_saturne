#ifndef __CS_XDEF_EVAL_AT_ZONE__
#define __CS_XDEF_EVAL_AT_ZONE__

/*============================================================================
 * Mass source terms computation.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
cs_xdef_eval_at_zone(const cs_xdef_t   *def,
                     const char        *caller_description,
                     cs_real_t          t_eval,
                     bool               dense,
                     cs_real_t         *values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_XDEF_EVAL_AT_ZONE__ */
