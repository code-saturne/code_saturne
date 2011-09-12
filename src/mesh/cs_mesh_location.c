/*============================================================================
 * Mesh locations management.
 *============================================================================*/

/*
  This file is part of the Code_Saturne Kernel, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1998-2011 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Kernel is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Kernel is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Kernel; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_parall.h"

#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

struct _cs_mesh_location_t {

  char                      name[32];     /* Name */

  const cs_mesh_t          *mesh;         /* Pointer to associated mesh */
  cs_mesh_location_type_t   type;         /* Location type */

  cs_lnum_t                 n_elts[3];    /* Number of associated elements:
                                             0: local,
                                             1: with standard ghost elements,
                                             2: with extended ghost elements */

  cs_lnum_t                *elt_list;     /* List of associated elements,
                                             (1 to n numbering) if non
                                             trivial (i.e. a subset) */

};

#endif

/*============================================================================
 * Static global variables
 *============================================================================*/

static int                  _n_mesh_locations_max = 0;
static int                  _n_mesh_locations = 0;
static cs_mesh_location_t  *_mesh_location = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get a a pointer to a mesh location by its id.
 *
 * parameters:
 *   id            <-- id of mesh location
 *
 * returns:
 *   pointer to associated mesh location
 *----------------------------------------------------------------------------*/

static const cs_mesh_location_t *
_const_mesh_location_by_id(int id)
{
  const cs_mesh_location_t *retval = NULL;

  if (id < 0 || id > _n_mesh_locations)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested mesh location\n"
                "%d is not defined.\n"), id);
  else
    retval = _mesh_location + id;

  return retval;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of mesh locations defined.
 *
 * \return  number of mesh locations defined
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_n_locations(void)
{
  return _n_mesh_locations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize mesh location API.
 *
 * By default, 4 mesh locations are built, matching the 4 first values of
 * the cs_mesh_location_type_t enum: CS_MESH_LOCCATION_CELLS for the cells
 * of the (default) global mesh, CS_MESH_LOCATION_INTERIOR_FACES and
 * CS_MESH_LOCATION_BOUNDARY_FACES for its faces, and
 * CS_MESH_LOCATION_VERTICES for its vertices.
 *
 * Locations should thus be built once the global mesh is complete, and
 * its halo structures completed.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_initialize(void)
{
  const cs_mesh_t *m = cs_glob_mesh;

  assert(cs_glob_mesh != NULL);
  assert(m->halo != NULL || (cs_glob_n_ranks == 1 && m->periodicity == NULL));

  cs_mesh_location_define(m,
                          "cells",
                          CS_MESH_LOCATION_CELLS,
                          m->n_cells,
                          NULL);
  cs_mesh_location_define(m,
                          "interior_faces",
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          m->n_i_faces,
                          NULL);
  cs_mesh_location_define(m,
                          "boundary_faces",
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          m->n_b_faces,
                          NULL);
  cs_mesh_location_define(m,
                          "vertices",
                          CS_MESH_LOCATION_VERTICES,
                          m->n_vertices,
                          NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize mesh location API.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_finalize(void)
{
  int  i;

  for (i = 0; i < _n_mesh_locations; i++) {
    cs_mesh_location_t  *ml = _mesh_location;
    if (ml->elt_list != NULL)
      BFT_FREE(ml->elt_list);
  }

  _n_mesh_locations = 0;
  _n_mesh_locations_max = 0;

  BFT_FREE(_mesh_location);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location.
 *
 * If a list of associated elements is given (defining a subset of a main
 * location), its ownership is transferred to the mesh location.
 *
 * \param[in]  mesh      pointer to associated mesh structure
 * \param[in]  name      name of location to define
 * \param[in]  type      type of location to define
 * \param[in]  n_elts    local number of associated elements
 * \param[in]  elt_list  list of associated elements (1 to n numbering,
 *                       size n_elts) if the location is a subset of a
 *                       main location, or NULL
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_define(const cs_mesh_t          *mesh,
                        const char               *name,
                        cs_mesh_location_type_t   type,
                        cs_lnum_t                 n_elts,
                        cs_lnum_t                *elt_list)
{
  /* local variables */

  int    i;
  cs_mesh_location_t  *ml;

  int id = _n_mesh_locations;

  /* Allocate new locations if necessary */

  if (_n_mesh_locations >= _n_mesh_locations_max) {
    if (_n_mesh_locations_max == 0)
      _n_mesh_locations_max = 4;
    else
      _n_mesh_locations_max *= 2;
    BFT_REALLOC(_mesh_location,
                _n_mesh_locations_max,
                cs_mesh_location_t);
  }

  _n_mesh_locations++;

  /* Define mesh location */

  ml = _mesh_location + id;

  ml->mesh = mesh;

  strncpy(ml->name, name, 31);
  ml->name[31] = '\0';

  ml->type = type;

  for (i = 0; i < 3; i++)
    ml->n_elts[i] = n_elts;

  if (   type == CS_MESH_LOCATION_CELLS
      && n_elts == mesh->n_cells && elt_list == NULL) {
    ml->n_elts[0] = mesh->n_cells;
    if (mesh->halo != NULL) {
      ml->n_elts[1] = mesh->halo->n_elts[0];
      ml->n_elts[2] = mesh->halo->n_elts[1];
    }
  }

  ml->elt_list = elt_list;

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's name.
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to mesh location name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_mesh_location_get_name(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's type.
 *
 * \param[in]  id  id of mesh location
 *
 * \return  mesh location type
 */
/*----------------------------------------------------------------------------*/

cs_mesh_location_type_t
cs_mesh_location_get_type(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's number of elements.
 *
 * A pointer to a array of 3 values is returned:
 *   0: local number of elements
 *   1: with standard ghost elements (if applicable)
 *   2: with extended ghost elements (if applicable)
 *
 * \param[in]  id  id of mesh location
 *
 * \return  array of numbers of elements.
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_n_elts(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->n_elts;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's elements list, if present.
 *
 * A list of elements is defined if the location is a subset of a main
 * location type.
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to elements list.
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_elt_list(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->elt_list;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
