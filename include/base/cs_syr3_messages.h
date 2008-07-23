/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_SYR3_MESSAGES_H__
#define __CS_SYR3_MESSAGES_H__

/*============================================================================
 * Manage messages for Syrthes coupling: sending, receiving and interpolation
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_syr3_coupling.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if Syrthes coupling continues or if we must finalize communications.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTSYR (IMSFIN)
 * *****************
 *
 * INTEGER          IMSFIN      : <-- : Indicates "end" message
 * INTEGER          NTMABS      : <-> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstsyr, TSTSYR)
(
 cs_int_t *const imsfin,
 cs_int_t *const ntmabs,
 cs_int_t *const ntcabs
);

/*----------------------------------------------------------------------------
 * Synchronize new time step message.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ITDSYR (NTCABS, NTMABS)
 * *****************
 *
 * INTEGER          NTMABS      : --> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(itdsyr, ITDSYR)
(
 cs_int_t   *const ntcabs,
 cs_int_t   *const ntmabs
);

/*----------------------------------------------------------------------------
 * Receive coupling variables from Syrthes
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI (NUMSYR, NOMRUB, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of Syrthes coupling
 * INTEGER          NBRENT      : --> : Number of elements
 * DOUBLE PRECISION TWALL       : <-- : Wall temerature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *const numsyr,
 cs_int_t   *const nbrent,
 cs_real_t  *const twall
);

/*----------------------------------------------------------------------------
 * Send coupling variables to Syrthes
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO (NUMSYR, NOMRUB, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of Syrthes coupling
 * INTEGER          NBRENT      : --> : Number of elements
 * REAL             TFLUID      : --> : Fluid temperature
 * REAL             HWALL       : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *const numsyr,
 cs_int_t   *const nbrent,
 cs_real_t  *const tfluid,
 cs_real_t  *const hwall
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SYR3_MESSAGES_H__ */
