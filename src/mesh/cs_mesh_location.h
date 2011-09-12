#ifndef __CS_MESH_LOCATION_H__
#define __CS_MESH_LOCATION_H__

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Mesh location types */

typedef enum {

  CS_MESH_LOCATION_CELLS,
  CS_MESH_LOCATION_INTERIOR_FACES,
  CS_MESH_LOCATION_BOUNDARY_FACES,
  CS_MESH_LOCATION_VERTICES,
  CS_MESH_LOCATION_PARTICLES,
  CS_MESH_LOCATION_OTHER,
  CS_MESH_LOCATION_NONE

} cs_mesh_location_type_t;

/* Opaque mesh location object */

typedef struct _cs_mesh_location_t cs_mesh_location_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return number of mesh locations defined.
 *----------------------------------------------------------------------------*/

int
cs_mesh_location_n_locations(void);

/*----------------------------------------------------------------------------
 * Initialize mesh location API.
 *
 * By default, 4 mesh locations are built, matching the 4 first values of
 * the cs_mesh_location_type_t enum: CS_MESH_LOCCATION_CELLS for the cells
 * of the (default) global mesh, CS_MESH_LOCATION_INTERIOR_FACES and
 * CS_MESH_LOCATION_BOUNDARY_FACES for its faces, and
 * CS_MESH_LOCATION_VERTICES for its vertices.
 *
 * Locations should thus be built once the global mesh is complete, and
 * its halo structures completed.
 *----------------------------------------------------------------------------*/

void
cs_mesh_location_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize mesh location API.
 *----------------------------------------------------------------------------*/

void
cs_mesh_location_finalize(void);

/*----------------------------------------------------------------------------
 * Define a new mesh location.
 *
 * If a list of associated elements is given (defining a subset of a main
 * location), its ownership is transferred to the mesh location.
 *
 * parameters:
 *   mesh     <-- pointer to associated mesh structure
 *   name     <-- name of location to define
 *   type     <-- type of location to define
 *   n_elts   <-- local number of associated elements
 *   elt_list <-- list of associated elements (1 to n numbering,
 *                size n_elts) if the location is a subset of a
 *                main location, or NULL
 *
 * returns:
 *   id of newly defined created mesh location
 *----------------------------------------------------------------------------*/

int
cs_mesh_location_define(const cs_mesh_t          *mesh,
                        const char               *name,
                        cs_mesh_location_type_t   type,
                        cs_lnum_t                 n_elts,
                        cs_lnum_t                *elt_list);

/*----------------------------------------------------------------------------
 * Get a mesh location's name.
 *
 * parameters:
 *   id <-- id of mesh location
 *
 * returns:
 *   pointer to mesh location name
 *----------------------------------------------------------------------------*/

const char *
cs_mesh_location_get_name(int id);

/*----------------------------------------------------------------------------
 * Get a mesh location's type.
 *
 * parameters:
 *   id <-- id of mesh location
 *
 * returns:
 *    mesh location type
 *----------------------------------------------------------------------------*/

cs_mesh_location_type_t
cs_mesh_location_get_type(int id);

/*----------------------------------------------------------------------------
 * Get a mesh location's number of elements.
 *
 * A pointer to a array of 3 values is returned:
 *   0: local number of elements
 *   1: with standard ghost elements (if applicable)
 *   2: with extended ghost elements (if applicable)
 *
 * parameters:
 *   id <-- id of mesh location
 *
 * returns:
 *   array of numbers of elements.
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_n_elts(int id);

/*----------------------------------------------------------------------------
 * Get a mesh location's elements list, if present.
 *
 * A list of elements is defined if the location is a subset of a main
 * location type.
 *
 * parameters:
 *   id <-- id of mesh location
 *
 * returns:
 *   pointer to elements list.
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_elt_list(int id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_LOCATION_H__ */
