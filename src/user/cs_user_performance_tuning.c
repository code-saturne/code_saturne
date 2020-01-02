/*============================================================================
 * Definition of advanced options relative to parallelism.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_performance_tuning.c
 *
 * \brief Definition of advanced options relative to parallelism.
 *
 * See \subpage cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced mesh numbering options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_numbering(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced partitioning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_partition(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define parallel IO settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parallel_io(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse matrix tuning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_matrix_tuning(void)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
