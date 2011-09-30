/*!
 * \file mei_math_util.c
 *
 * \brief Provides mathemathical functions facilities
 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "mei_math_util.h"

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Return the max value from two doubles.
 *
 * \param [in] x1 double
 * \param [in] x2 double
 * \return max value
 */

double
mei_max(double x1, double x2)
{
    if (x1 < x2)
        return x2;
    else
        return x1;
}

/*!
 * \brief Return the min value from two doubles.
 *
 * \param [in] x1 double
 * \param [in] x2 double
 * \return min value
 */

double
mei_min(double x1, double x2)
{
    if (x1 < x2)
        return x1;
    else
        return x2;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */


