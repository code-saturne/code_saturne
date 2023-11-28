/*----------------------------------------------------------------------------*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * \file cs_meg_immersed_boundaries_inout.c
 *
 * \brief This function is used to indicate whether a given point is within or
 *        outside a given solid
 *
 * \param[in, out] ipenal       indicator for cut cells algorithm
 * \param[in]      object_name  name of the solid object
 * \param[in]      xyz          pointer containing the point coordinates
 * \param[in]      t            time value
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_immersed_boundaries_inout
void
cs_meg_immersed_boundaries_inout(int         *ipenal,
                                 const char  *object_name,
                                 cs_real_3_t  xyz,
                                 cs_real_t    t)
{
  CS_NO_WARN_IF_UNUSED(ipenal);
  CS_NO_WARN_IF_UNUSED(object_name);
  CS_NO_WARN_IF_UNUSED(xyz);
  CS_NO_WARN_IF_UNUSED(t);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
