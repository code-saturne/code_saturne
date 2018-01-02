#ifndef __ECS_MEM_USAGE_H__
#define __ECS_MEM_USAGE_H__

/*============================================================================
 * Base memory usage information (System and Library dependent)
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

#include "ecs_def.h"

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Initialize memory usage count depending on system.
 *
 * This functions checks if it has already been called, so
 * it is safe to call more than once (though it is not
 * thread-safe). Only the first call is effective.
 */

void
ecs_mem_usage_init(void);

/*
 * End memory usage count depending on system.
 */

void
ecs_mem_usage_end(void);

/*
 * Indicates if ecs_mem_usage_...() functions are initialized.
 *
 * returns:
 *   1 if ecs_mem_usage_init has been called, 0 otherwise.
 */

int
ecs_mem_usage_initialized(void);

/*
 * Return current process memory use (in kB) depending on OS.
 */

size_t
ecs_mem_usage_pr_size(void);

/*
 * Return maximum process memory use (in kB) depending on OS.
 *
 * The returned value is the maximum returned by ecs_mem_usage_pr_size()
 * during the program's lifetime. With memory allocations which return
 * memory to the system (such as the GNU glibc on Linux systems),
 * this value will be correct only if allocation is tracked. This should
 * be the case if malloc hooks are used with the glibc allocation
 * functions (ECS library's default configuration/installation option),
 * but may give results lower than the true maximum in other cases.
 */

size_t
ecs_mem_usage_max_pr_size(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __ECS_MEM_USAGE_H__ */
