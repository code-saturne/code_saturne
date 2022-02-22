/*============================================================================
 * Mathematical Expression Generator function stubs for FSI internal coupling.
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
 * \brief This function is used to query FSI internal coupling structure values
 *        for a given boundary and structure.
 *
 * \param[in]      object_type   name of object type
 * \param[in]      name          name of matching boundary
 * \param[in]      fluid_f       array of fluid forces on the object
 * \param[in,out]  val[]         matrix or vector coefficients
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_fsi_struct
void
cs_meg_fsi_struct(const char       *object_type,
                  const char       *name,
                  const cs_real_t   fluid_f[],
                  cs_real_t         val[])
{
  CS_UNUSED(object_type);
  CS_UNUSED(name);
  CS_UNUSED(fluid_f);
  CS_UNUSED(val);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
