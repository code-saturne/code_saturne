/*============================================================================
 * Insert thin walls into the mesh.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh_boundary.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_thinwall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 * The created faces share vertices.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * \deprecated This function is replaced by \ref cs_mesh_boundary_insert.
 *
 * \param[in]       mesh           pointer to mesh structure to modify
 * \param[in, out]  face_list      list of selected (interior) faces (0 to n-1)
 * \param[in]       face_list_size number of selected (interior) faces
 */
/*----------------------------------------------------------------------------*/

void
cs_create_thinwall(cs_mesh_t  *mesh,
                   cs_lnum_t  *face_list,
                   cs_lnum_t   face_list_size)
{
  cs_mesh_boundary_insert_with_shared_vertices(mesh,
                                               face_list_size,
                                               face_list);
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
