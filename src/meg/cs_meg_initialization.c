/*----------------------------------------------------------------------------*/

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
 * \file cs_meg_volume_initialization.c
 *
 * \brief This function is used for initalization of fields over a
 * given volume zone.
 *
 * The caller is responsible for freeing the associated array.
 *
 * \param[in, out]  field_name  associated field name
 * \param[in]       vz          pointer to associated volume zone
 *
 * \return  pointer to allocated initialization values.
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_meg_initialization(const char       *field_name,
                      const cs_zone_t  *vz)
{
  CS_UNUSED(field_name);
  CS_UNUSED(vz);

  return NULL; /* avoid a compilation warning */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
