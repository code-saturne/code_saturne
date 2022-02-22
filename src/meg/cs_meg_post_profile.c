/*============================================================================
 * Mathematical Expression Generator function stubs for output control.
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is used to define profile coordinates.
 *
 * \param[in]       name          name of matching profile
 * \param[in]       n_coords      number of point coordinates
 * \param[in, out]  coords        point coordinates
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_post_profiles
void
cs_meg_post_profiles(const char       *name,
                     int               n_coords,
                     cs_real_t         coords[][3])
{
  /* Avoid compilation warnings */
  CS_UNUSED(name);
  CS_UNUSED(n_coords);
  CS_UNUSED(coords);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
