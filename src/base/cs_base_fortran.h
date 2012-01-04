#ifndef __CS_BASE_FORTRAN_H__
#define __CS_BASE_FORTRAN_H__

/*============================================================================
 * Initializtion and handling of Fortran-related mechanisms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * Create a directory, or check it exists.
 *
 * Fortran interface
 *
 * SUBROUTINE CSMKDR (DIRNAM, DIRLEN)
 * *****************
 *
 * CHARACTER*       DIRNAM      : --> : Directory name
 * INTEGER          DIRLEN      : --> : Directory name length
 *----------------------------------------------------------------------------*/

void CS_PROCF (csmkdr, CSMKDR)
(
 const char       *dirnam,
 const cs_int_t   *dirlen
);

/*----------------------------------------------------------------------------
 * Copy a Fortan string buffer to a C string buffer
 *
 * The aim of this function is to aviod issues with Fortran array bounds
 * checking when compilers such as icc 11 consider a character array from C
 * as an array of 1-character length strings.
 *
 * Fortran interface
 *
 * SUBROUTINE CSSF2C (LEN, CSTR, FSTR)
 * *****************
 *
 * INTEGER          LEN         : --> : String length
 * CHARACTER*       FSTR        : --> : Fortran string
 * CHARACTER*       CSTR        : <-- : C string
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssf2c, CSSF2C)
(
 const cs_int_t   *len,
 const char       *fstr,
 char             *cstr
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

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
