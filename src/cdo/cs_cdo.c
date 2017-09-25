/*============================================================================
 * General functions or variables for the CDO module
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Activation of the CDO/HHO module */
int  cs_cdo_activation_mode = CS_CDO_OFF;

/* Separation lines: long, medium, short */
const char lsepline[80] =
  "# =======================================================================\n";
const char msepline[60] =
  "# =========================================\n";
const char ssepline[40] =
  "# =================\n";

/* Default locations */
const cs_flag_t  cs_cdo_primal_vtx  = CS_FLAG_PRIMAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_cdo_primal_face = CS_FLAG_PRIMAL | CS_FLAG_FACE;
const cs_flag_t  cs_cdo_primal_cell = CS_FLAG_PRIMAL | CS_FLAG_CELL;
const cs_flag_t  cs_cdo_dual_vtx  = CS_FLAG_DUAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_cdo_dual_face = CS_FLAG_DUAL | CS_FLAG_FACE;
const cs_flag_t  cs_cdo_dual_cell = CS_FLAG_DUAL | CS_FLAG_CELL;
const cs_flag_t  cs_cdo_dual_face_byc =
  CS_FLAG_DUAL | CS_FLAG_FACE | CS_FLAG_BY_CELL;

/*=============================================================================
 * Local static variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(DEBUG) && !defined(NDEBUG)
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
cs_dump_array_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_real_t    array[],
                         int                n_cols)
{
  cs_log_printf(CS_LOG_DEFAULT, "\nDUMP>> %s\n", header);

  if (n_cols < 1) n_cols = 1;
  int  n_rows = size/n_cols;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    for (cs_lnum_t j = i*n_cols; j < (i+1)*n_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6.4e", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  if (n_rows*n_cols < size) {
    for (cs_lnum_t j = n_rows*n_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6.4e", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

}

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
cs_dump_integer_to_listing(const char        *header,
                           const cs_lnum_t    size,
                           const cs_lnum_t    array[],
                           int                n_cols)
{
  cs_log_printf(CS_LOG_DEFAULT, "\nDUMP>> %s\n", header);

  if (n_cols < 1) n_cols = 1;
  int  n_rows = size/n_cols;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    for (cs_lnum_t j = i*n_cols; j < (i+1)*n_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6d", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

  if (n_rows*n_cols < size) {
    for (cs_lnum_t j = n_rows*n_cols; j < size; j++)
      cs_log_printf(CS_LOG_DEFAULT, " (%04d) % 6d", j, array[j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

}
#endif  /* Only in debug mode */

/*----------------------------------------------------------------------------*/

END_C_DECLS
