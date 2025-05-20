#ifndef __CS_PARAM_CDO_H__
#define __CS_PARAM_CDO_H__

/*============================================================================
 * High-level metadata related to CDO/HHO schemes
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Same mechanism as the one used in PETSC_DEFAULT */

#define CS_CDO_KEEP_DEFAULT  -2

/* Specifications for OpenMP loops/sections */

#define CS_CDO_OMP_CHUNK_SIZE   128
#define CS_CDO_OMP_SCHEDULE     schedule(static, CS_CDO_OMP_CHUNK_SIZE)
#define CS_CDO_OMP_SYNC_MODE    0 /* > 0 --> critical sections otherwise
                                   *         atomic sections is used */

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

/* The following limitation only results from an optimization in the size of
   the bit mask (can be changed if needed by changing the definition of
   the type cs_mask_t)
   Here is the max. number of reaction terms allowed in an equation */

#define CS_CDO_N_MAX_REACTIONS  8

#define CS_ALL_FACES   0        /* All faces: interior + border */
#define CS_BND_FACES   1        /* Boundary faces */
#define CS_INT_FACES   2        /* Interior faces */

/* HHO specific part:
 *
 * Number of DoFs on faces and cells according to the polynomial space
 */

#define CS_N_DOFS_FACE_0TH  1
#define CS_N_DOFS_FACE_1ST  3
#define CS_N_DOFS_FACE_2ND  6

#define CS_N_DOFS_CELL_0TH  1
#define CS_N_DOFS_CELL_1ST  4
#define CS_N_DOFS_CELL_2ND  10

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_PARAM_CDO_MODE_OFF     = -1,  /* CDO schemes are not used */
  CS_PARAM_CDO_MODE_WITH_FV =  1,  /* CDO and legacy FV schemes are used */
  CS_PARAM_CDO_MODE_NS_WITH_FV = 2,   /* CDO schemes activated for NSE
                                         and legacy FV schemes for
                                         other equations */
  CS_PARAM_CDO_MODE_ONLY = 3          /* Only CDO schemes are used */

} cs_param_cdo_mode_t;

/*============================================================================
 * Global variables
 *============================================================================*/

extern cs_param_cdo_mode_t  cs_glob_param_cdo_mode;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the global variable storing the mode of activation to apply to
 *        CDO/HHO schemes. Deprecated way to set the CDO mode.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_cdo_mode_set(cs_param_cdo_mode_t   mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the mode of activation for the CDO/HHO schemes.
 *
 * \return the mode of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

cs_param_cdo_mode_t
cs_param_cdo_mode_get(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print a welcome message indicating what is the current CDO status
 */
/*----------------------------------------------------------------------------*/

void
cs_param_cdo_log(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print generic parameters used with CDO/HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_param_cdo_setup_log(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tell if FVM is used to solve main equations
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_cdo_has_fv_main(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tell if FVM and CDO are used
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_cdo_has_cdo_and_fv(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tell if CDO is only used
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_cdo_has_cdo_only(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tell if FVM is oly used
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_cdo_has_fv_only(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_CDO_H__ */
