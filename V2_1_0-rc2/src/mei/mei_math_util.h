#ifndef __MEI_MATH_UTIL_H__
#define __MEI_MATH_UTIL_H__

/*!
 * \file mei_math_util.h
 *
 * \brief Provides mathemathical functions facilities
 */

/*
  This file is part of the "Mathematical Expression Interpreter" library.

  Copyright (C) 2008-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*----------------------------------------------------------------------------
 * Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/

#include <math.h>

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the max value from two doubles.
 *
 * parameters:
 *   x1       <-- double
 *   x2       <-- double
 *----------------------------------------------------------------------------*/

double
mei_max(double x1, double x2);

/*----------------------------------------------------------------------------
 * Return the min value from two doubles.
 *
 * parameters:
 *   x1       <-- double
 *   x2       <-- double
 *----------------------------------------------------------------------------*/

double
mei_min(double x1, double x2);


#endif
