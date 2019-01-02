#ifndef __BFT_BACKTRACE_H__
#define __BFT_BACKTRACE_H__

/*============================================================================
 * Obtaining a stack backtrace
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/* BFT backtrace descriptor */

typedef struct _bft_backtrace_t bft_backtrace_t;

/* Pointers for backtrace print functions */

typedef void (bft_backtrace_print_t) (int  start_depth);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Build a backtrace description structure.
 *
 * returns:
 *   pointer to bft_backtrace_t backtrace descriptor (NULL in case of
 *   error, or if backtracing is unavailable on this architecture).
 */

bft_backtrace_t *
bft_backtrace_create(void);

/*
 * Free a backtrace description structure.
 *
 * parameters:
 *   bt: <-> pointer to backtrace description structure.
 *
 * returns:
 *   NULL pointer.
 */

bft_backtrace_t *
bft_backtrace_destroy(bft_backtrace_t  *bt);

/*!
 * Demangle a backtrace description structure (for C++).
 *
 * parameters:
 *   bt: <-> pointer to backtrace description structure.
 */

void
bft_backtrace_demangle(bft_backtrace_t  *bt);

/*
 * Return the total depth of a backtrace.
 *
 * parameters:
 *   bt: <-- pointer to backtrace description structure.
 *
 * returns:
 *   total backtrace depth.
 */

int
bft_backtrace_size(const bft_backtrace_t  *bt);

/*
 * Return file name associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bft_backtrace_size(bt)).
 *
 * returns:
 *   file name at the given depth, or NULL.
 */

const char *
bft_backtrace_file(bft_backtrace_t  *bt,
                   int               depth);

/*
 * Return function name associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bft_backtrace_size(bt)).
 *
 * returns:
 *   function name at the given depth, or NULL.
 */

const char *
bft_backtrace_function(bft_backtrace_t  *bt,
                       int               depth);

/*
 * Return address associated with a backtrace at a given depth.
 *
 * parameters:
 *   bt:    <-- pointer to backtrace description structure.
 *   depth: <-- index in backtrace structure (< bft_backtrace_size(bt)).
 *
 * returns:
 *   address at the given depth, or NULL.
 */

const char *
bft_backtrace_address(bft_backtrace_t  *bt,
                      int               depth);

/*
 * Print a backtrace.
 *
 * parameters:
 *   start_depth: <-- depth of backtrace at which to start printing
 *                    (0 for all, including backtrace print function)
 */

void
bft_backtrace_print(int  start_depth);

/*
 * Returns backtrace print function.
 *
 * returns:
 *   pointer to the backtrace print function.
 */

bft_backtrace_print_t *
bft_backtrace_print_get(void);

/*
 * Sets a backtrace print function.
 *
 * parameters:
 *   fct: <-- pointer to a bft_backtrace_print_t type function.
 */

void
bft_backtrace_print_set(bft_backtrace_print_t  *const fct);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __BFT_BACKTRACE_H__ */
