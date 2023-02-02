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
 * \file cs_meg_source_terms.c
 *
 * \brief This function is used to compute source terms over a volume zone
 *
 * The caller is responsible for freeing the returned array.
 *
 * \param[in] zone          pointer to cs_volume_zone_t
 * \param[in] field_name    variable field name
 * \param[in] source_type   source term type
 *
 * \returns a cs_real_t pointer containing the computed values
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_source_terms
cs_real_t *
cs_meg_source_terms(const char       *zone_name,
                    const cs_lnum_t   n_elts,
                    const cs_lnum_t  *elt_ids,
                    const cs_real_t   xyz[][3],
                    const char       *name,
                    const char       *source_type)
{
  CS_UNUSED(elt_ids);
  CS_UNUSED(n_elts);
  CS_UNUSED(name);
  CS_UNUSED(source_type);
  CS_UNUSED(xyz);
  CS_UNUSED(zone_name);

  return NULL; /* avoid a compilation warning */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
