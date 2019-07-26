#ifndef __BFT_MEM_USAGE_H__
#define __BFT_MEM_USAGE_H__

/*============================================================================
 * Base memory usage information (System and Library dependent)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
bft_mem_usage_init(void);

/*
 * End memory usage count depending on system.
 */

void
bft_mem_usage_end(void);

/*
 * Indicates if bft_mem_usage_...() functions are initialized.
 *
 * returns:
 *   1 if bft_mem_usage_init has been called, 0 otherwise.
 */

int
bft_mem_usage_initialized(void);

/*
 * Return current process memory use (in kB) depending on OS.
 */

size_t
bft_mem_usage_pr_size(void);

/*
 * Return maximum process memory use (in kB) depending on OS.
 *
 * The returned value is the maximum memory used during the program's
 * lifetime.
 *
 * returns:
 *   maximum measured program size, or 0 if not available
 */

size_t
bft_mem_usage_max_pr_size(void);

/*
 * Return maximum process virtual memory use (in kB) depending on OS.
 *
 * returns:
 *   maximum measured virtual memory usage, or 0 if not available
 */

size_t
bft_mem_usage_max_vm_size(void);

/*
 * Return shared library memory use (in kB) depending on OS.
 *
 * returns:
 *   maximum measured shared library memory usage, or 0 if not available
 */

size_t
bft_mem_usage_shared_lib_size(void);

/*
 * Return counter to number of calls to malloc, realloc, and free.
 *
 * This function returns zeroes when the appropriate instrumentation
 * is not available.
 */

void
bft_mem_usage_n_calls(size_t count[3]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __BFT_MEM_USAGE_H__ */
