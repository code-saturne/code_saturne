#ifndef __CS_CDO_H__
#define __CS_CDO_H__

/*============================================================================
 * General functions or variables for the INNOV module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Static inline function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(DEBUG) && !defined(NDEBUG)
/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a cs_sdm_t structure which is defined by block
 *          Print into the file f if given otherwise open a new file named
 *          fname if given otherwise print into the standard output
 *          The usage of threshold allows one to compare more easier matrices
 *          without taking into account numerical roundoff.
 *
 * \param[in]  fp         pointer to a file structure or NULL
 * \param[in]  fname      filename or NULL
 * \param[in]  thd        threshold (below this value --> set 0)
 * \param[in]  n_elts     size of the array
 * \param[in]  array      list of values to dump
 * \param[in]  n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_array_fprintf(FILE             *fp,
                     const char       *fname,
                     cs_real_t         thd,
                     cs_lnum_t         n_elts,
                     const cs_real_t   array[],
                     int               n_cols);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In debug mode, dump an array of double into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_darray_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_real_t    array[],
                         int                n_cols);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In debug mode, dump an array of integer into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_iarray_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_lnum_t    array[],
                         int                n_cols);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_H__ */
