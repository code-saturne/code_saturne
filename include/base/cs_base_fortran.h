/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_BASE_FORTRAN_H__
#define __CS_BASE_FORTRAN_H__

/*============================================================================
 * Initializtion and handling of Fortran-related mechanisms
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This is necessary for good consistency of messages output from C or
 * from Fortran, and to handle parallel and serial logging options.
 *
 * parameters:
 *   r0_log_flag <-- redirection for rank 0 log;
 *                   0: not redirected; 1: redirected to "listing" file
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   0: not redirected; 1: redirected to "listing_n*" file;
 *                   2: redirected to "/dev/null" (suppressed)
 *----------------------------------------------------------------------------*/

void
cs_base_fortran_bft_printf_set(int r0_log_flag,
                               int rn_log_flag);

END_C_DECLS

#endif /* __CS_BASE_FORTRAN_H__ */
