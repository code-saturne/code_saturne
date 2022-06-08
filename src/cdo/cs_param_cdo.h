#ifndef __CS_PARAM_CDO_H__
#define __CS_PARAM_CDO_H__

/*============================================================================
 * Manage the definition/setting of a computation
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Specifications for OpenMP loops */

#define CS_CDO_OMP_CHUNK_SIZE     128
#define CS_CDO_OMP_SCHEDULE       schedule(static, CS_CDO_OMP_CHUNK_SIZE)
#define CS_CDO_OMP_SYNC_SECTIONS  0 /* > 0 --> critical sections
                                       otherwise atomic sections */

/* Avoid issues with assert in some OpenMp contructs using gcc 9 */

#if defined(HAVE_OPENMP) && defined(__GNUC__)
  #if __GNUC__ == 9
    #define CS_CDO_OMP_ASSERT(e)
  #else
    #define CS_CDO_OMP_ASSERT(e)  assert(e)
  #endif
#else
  #define CS_CDO_OMP_ASSERT(e)  assert(e)
#endif

/* Size of the buffer used to collect global ids for rows and columns
   when assembling the values in the global matrix from the local cellwise
   matrices */

#define CS_CDO_ASSEMBLE_BUF_SIZE  200

/* The following limitation only results from an optimization in the size of
   the bit mask (can be changed if needed by changing the definition of
   the type cs_mask_t)
   Here is the max. number of reaction terms allowed in an equation */

#define CS_CDO_N_MAX_REACTIONS  8

#define CS_ALL_FACES   0        /* All faces: interior + border */
#define CS_BND_FACES   1        /* Boundary faces */
#define CS_INT_FACES   2        /* Interior faces */

/* Number of DoFs on faces and cells according to the polynomial space */

#define CS_N_DOFS_FACE_0TH  1
#define CS_N_DOFS_FACE_1ST  3
#define CS_N_DOFS_FACE_2ND  6

#define CS_N_DOFS_CELL_0TH  1
#define CS_N_DOFS_CELL_1ST  4
#define CS_N_DOFS_CELL_2ND  10

/*============================================================================
 * Type definitions
 *============================================================================*/

/* OpenMP STRATEGY FOR THE ASSEMBLY STEP */
/* ===================================== */

typedef enum {

  CS_PARAM_ASSEMBLE_OMP_ATOMIC,
  CS_PARAM_ASSEMBLE_OMP_CRITICAL,
  CS_PARAM_ASSEMBLE_OMP_N_STRATEGIES

} cs_param_assemble_omp_strategy_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_CDO_H__ */
