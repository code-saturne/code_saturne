
/* VERS */

/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * User function for geometry and mesh modification.
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Modifiy geometry and mesh.
 *
 * The mesh structure is described in cs_mesh.h
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Example: modify vertex coordinates */
  /*------------------------------------*/

  /* Divide coordinates by 1000 (millimetres to metres).
   *
   * Warning:
   *
   *   This is incompatible with pre-processed periodicity,
   *   as the periodicity transformation is not updated.
   *
   *   With periodicity, using a coordinate transformation matrix
   *   in cs_user_mesh_input is preferred. */

  {
    cs_lnum_t  vtx_id;
    const double  coo_mult = 1. / 1000.;

    for (vtx_id = 0; vtx_id < mesh->n_vertices; vtx_id++) {
      mesh->vtx_coord[vtx_id*3]     *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 1] *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 2] *= coo_mult;
    }

    /* Set mesh modification flag it i should be saved for future re-use */

    mesh->modified = 1;
  }
}
