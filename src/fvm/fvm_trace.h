#ifndef __FVM_TRACE_H__
#define __FVM_TRACE_H__

/*============================================================================
 * Tracing utility functions for profiling and debugging
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Print memory usage status.
 */

#define FVM_TRACE_MEM_STATUS \
  {char str[256]; sprintf(str, "%s:%d", __FILE__, __LINE__); \
   fvm_trace_mem_status(str);}

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print memory usage status.
 *
 * If no associated description is given, a call number will be issued.
 *
 * parameters:
 *   descr <-- optional description, or NULL
 *----------------------------------------------------------------------------*/

void
fvm_trace_mem_status(const char  *descr);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TRACE_H__ */
