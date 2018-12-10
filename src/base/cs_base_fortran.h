#ifndef __CS_BASE_FORTRAN_H__
#define __CS_BASE_FORTRAN_H__

/*============================================================================
 * Initializtion and handling of Fortran-related mechanisms
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

/*----------------------------------------------------------------------------
 * Call exit routine from Fortran code
 *
 * Fortran interface:
 *
 * subroutine csexit (status)
 * *****************
 *
 * integer          status      : <-- : 0 for success, 1+ for error
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *status
);

/*----------------------------------------------------------------------------
 * Elapsed time since execution start
 *
 * Fortran interface:
 *
 * subroutine dmtmps (tw)
 * *****************
 *
 * double precision tw          : <-- : elapsed time
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *tw
);

/*----------------------------------------------------------------------------
 * Create a directory, or check it exists.
 *
 * Fortran interface
 *
 * subroutine csmkdr (dirnam, dirlen)
 * *****************
 *
 * character*       dirnam      : --> : Directory name
 * integer          dirlen      : --> : Directory name length
 *----------------------------------------------------------------------------*/

void CS_PROCF (csmkdr, CSMKDR)
(
 const char       *dirnam,
 const cs_int_t   *dirlen
);

/*----------------------------------------------------------------------------
 * Compute the gamma function of x.
 *
 * Fortran interface
 *
 * subroutine csgamma (x, g)
 * ******************
 *
 * double precision  x     : <-- : input value
 * double precision  gamma : --> : output value
 *----------------------------------------------------------------------------*/

void CS_PROCF (csgamma, CSGAMMA)
(
 const cs_real_t  *x,
 cs_real_t        *gamma
 );

/*----------------------------------------------------------------------------
 * Compute the erf function of x.
 *
 * Fortran interface
 *
 * subroutine cserf (x, g)
 * ******************
 *
 * double precision  x     : <-- : input value
 * double precision  ferf  : --> : output value
 *----------------------------------------------------------------------------*/

void CS_PROCF (cserf, CSERF)
(
 const cs_real_t  *x,
 cs_real_t        *ferf
 );

/*----------------------------------------------------------------------------
 * Get log name file information.
 *
 * When log file output is suppressed, it returns the name of the
 * bit buck file ("/dev/null")
 *
 * Fortran interface
 *
 * subroutine cslogname (len, name)
 * ********************
 *
 * integer          len         : <-- : maximum string length
 * character*       name        : --> : Fortran string
 *----------------------------------------------------------------------------*/

void CS_PROCF (cslogname, CSLOGNAME)
(
 const cs_int_t   *len,
 char             *dir
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Get package data path information.
 *
 * The aim of this function is to aviod issues with Fortran array bounds
 * checking when compilers such as icc 11 consider a character array from C
 * as an array of 1-character length strings.
 *
 * Fortran interface
 *
 * subroutine csdatadir (len, dir)
 * ********************
 *
 * integer          len         : <-- : maximum string length
 * character*       dir         : --> : Fortran string
 *----------------------------------------------------------------------------*/

void CS_PROCF (csdatadir, CSDATADIR)
(
 const cs_int_t   *len,
 char             *dir
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This variant is designed to allow switching from C to Fortran output,
 * whithout disabling regular C stdout output when switched to Fortran.
 *
 * This allows redirecting or suppressing logging for different ranks.
 *
 * parameters:
 *   log_name    <-- base file name for log, or NULL for stdout
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   false:  to "/dev/null" (suppressed)
 *                   true: redirected to <log_name>_n*.log" file;
 *----------------------------------------------------------------------------*/

void
cs_base_fortran_bft_printf_set(const char  *log_name,
                               bool         rn_log_flag);

/*----------------------------------------------------------------------------
 * Switch bft_printf() mechanism to C output.
 *
 * This function may only be called after cs_base_fortran_bft_printf_set()
 *----------------------------------------------------------------------------*/

void
cs_base_fortran_bft_printf_to_c(void);

/*----------------------------------------------------------------------------
 * Switch bft_printf() mechanism to Fortran output.
 *
 * This function may only be called after cs_base_fortran_bft_printf_set()
 *----------------------------------------------------------------------------*/

void
cs_base_fortran_bft_printf_to_f(void);

/*----------------------------------------------------------------------------
 * Wrappers to cs_user_extra_operations
 *----------------------------------------------------------------------------*/

void
cs_user_extra_operations_initialize_wrapper(void);

void
cs_user_extra_operations_wrapper(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_FORTRAN_H__ */
