#ifndef __MEI_MATH_UTIL_H__
#define __MEI_MATH_UTIL_H__

/*!
 * \file mei_math_util.h
 *
 * \brief Provides mathemathical functions facilities
 */

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
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <math.h>

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Structure associated to a single user data set for 1D interpolation
 */
/*----------------------------------------------------------------------------*/

typedef struct {

  int       ncols;           /*!< number of columns */
  int       nrows;           /*!< number of lines */
  double  **values;          /*!< values for the data set */
  char     *name;            /*!< name of the data file */
  char     *description;     /*!< user description */

} mei_user_data_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return an interpolated value.
 *
 * param [in] filename char
 * param [in] c1 double
 * param [in] c2 double
 * param [in] x variable to interpolate
 *
 * return interpolated value
 *----------------------------------------------------------------------------*/

double
mei_interp1d(const char  *filename,
             int          c1,
             int          c2,
             double       x);

/*-----------------------------------------------------------------------------*/
/*
 * Destroy all user data set for 1D interpolation.
 */
/*-----------------------------------------------------------------------------*/

void mei_data_free(void);

/*-----------------------------------------------------------------------------*/

#endif
