#ifndef __CS_RUNAWAY_H__
#define __CS_RUNAWAY_H__

/*============================================================================
 * Runaway (diverging) computation detection.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that defined field bounds are not exceeded.
 *
 * \return 0 if no bounds are exceeded, 1 if bounds are exceeded.
 */
/*----------------------------------------------------------------------------*/

int
cs_runaway_check(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define maximum value for a field, beyon which computation is aborted.
 *
 * Currently, only one field is handled, so calling this multiple times
 * replaced the previous setting. Using a negative field id removes
 * this check.
 */
/*----------------------------------------------------------------------------*/

void
cs_runaway_check_define_field_max(int        f_id,
                                  cs_real_t  max_allowed);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that defined field bounds are not exceeded.
 */
/*----------------------------------------------------------------------------*/

void
cs_runaway_check_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RUNAWAY_H__ */
