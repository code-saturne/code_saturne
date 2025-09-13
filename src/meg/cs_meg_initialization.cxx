/*----------------------------------------------------------------------------*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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
 * \brief This function is used for the initalization of fields over a given
 *        volume zone. The mathematical expression is defined in the GUI.
 *
 * The caller is responsible for freeing the associated array.
 *
 * \param[in]  zone_name    name of a volume zone
 * \param[in]  n_elts       number of elements related to the zone
 * \param[in]  elt_ids      list of element ids related to the zone
 * \param[in]  xyz          list of coordinates related to the zone
 * \param[in]  field_name   associated variable field name
 * \param[out] retvals      array of computed values
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_initialization
void
cs_meg_initialization([[maybe_unused]] const char       *zone_name,
                      [[maybe_unused]] const cs_lnum_t   n_elts,
                      [[maybe_unused]] const cs_lnum_t  *elt_ids,
                      [[maybe_unused]] const cs_real_t   xyz[][3],
                      [[maybe_unused]] const char       *field_name,
                      [[maybe_unused]] cs_real_t        *retvals)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
