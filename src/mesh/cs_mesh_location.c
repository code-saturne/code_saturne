/*============================================================================
 * Mesh locations management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

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

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for cs_mesh_location.h) */

/*!
 * \typedef cs_mesh_location_select_t
 *
 * \brief Function pointer to mesh location elements selection definition.
 *
 * If non-empty and not containing all elements, a list of elements
 * of the parent mesh belonging to the location should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * This list's lifecycle is then managed by the mesh location object.
 *
 * \param [in]   m            pointer to associated mesh structure.
 * \param [in]   location_id  id of associated location.
 * \param [out]  n_elts       number of selected elements
 * \param [out]  elt_list     list of selected elements.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

struct _cs_mesh_location_t {

  char                        name[32];     /* Name */

  const cs_mesh_t            *mesh;         /* Pointer to associated mesh */

  cs_mesh_location_type_t     type;         /* Location type */

  char                       *select_str;   /* String */
  cs_mesh_location_select_t  *select_fp;    /* Function pointer */

  cs_lnum_t                 n_elts[3];    /* Number of associated elements:
                                             0: local,
                                             1: with standard ghost elements,
                                             2: with extended ghost elements */

  cs_lnum_t                *elt_list;     /* List of associated elements,
                                             (0 to n-1 numbering) if non
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

/*----------------------------------------------------------------------------
 * Define a new mesh location.
 *
 * If a list of associated elements is given (defining a subset of a main
 * location), its ownership is transferred to the mesh location.
 *
 * parameters:
 *   name      name of location to define
 *   type      type of location to define
 *
 * returns:
 *   id of newly defined created mesh location
 *----------------------------------------------------------------------------*/

static int
_mesh_location_define(const char               *name,
                      cs_mesh_location_type_t   type)
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

  ml->mesh = NULL;

  strncpy(ml->name, name, 31);
  ml->name[31] = '\0';

  ml->type = type;

  ml->select_str = NULL;
  ml->select_fp = NULL;

  for (i = 0; i < 3; i++)
    ml->n_elts[i] = 0;

  ml->elt_list = NULL;

  return id;
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
 * By default, 5 mesh locations are built, matching the 5 first values of
 * the cs_mesh_location_type_t enum: CS_MESH_LOCATION_NONE for global
 * values, CS_MESH_LOCCATION_CELLS for the cellsof the (default) global mesh,
 * CS_MESH_LOCATION_INTERIOR_FACES and CS_MESH_LOCATION_BOUNDARY_FACES for
 * its faces, and CS_MESH_LOCATION_VERTICES for its vertices.
 *
 * Locations should then be built once the global mesh is complete, and
 * its halo structures completed.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_initialize(void)
{
  cs_mesh_location_define(N_("global"),
                          CS_MESH_LOCATION_NONE,
                          NULL);
  cs_mesh_location_define(N_("cells"),
                          CS_MESH_LOCATION_CELLS,
                          NULL);
  cs_mesh_location_define(N_("interior_faces"),
                          CS_MESH_LOCATION_INTERIOR_FACES,
                          NULL);
  cs_mesh_location_define(N_("boundary_faces"),
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          NULL);
  cs_mesh_location_define(N_("vertices"),
                          CS_MESH_LOCATION_VERTICES,
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
    BFT_FREE(ml->elt_list);
    BFT_FREE(ml->select_str);
  }

  _n_mesh_locations = 0;
  _n_mesh_locations_max = 0;

  BFT_FREE(_mesh_location);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate mesh locations with a mesh.
 *
 * If mesh_id is negative, all defined mesh locations are associated
 * (which is useful for the common case where only one mesh is present).
 * If mesh_id is non-negative, only the location with the matching
 * id is associated (which may be useful when multiple meshes are defined).
 *
 * The number of elements are computed based on the underlying mesh,
 * and element lists are built for mesh subset locations.
 *
 * \param[in]  mesh  pointer to associated mesh structure
 * \param[in]  id    id of mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_build(cs_mesh_t  *mesh,
                       int         id)
{
  int  i;
  int id_start = 0, id_end = _n_mesh_locations;

  assert(mesh != NULL);
  assert(   mesh->halo != NULL
         || (cs_glob_n_ranks == 1 && mesh->periodicity == NULL));

  if (id >= 0) {
    id_start = id;
    if (id < _n_mesh_locations)
      id_end = id + 1;
  }

  for (i = id_start; i < id_end; i++) {

    int n_elts_max = 0;
    fvm_selector_t *selector = NULL;
    cs_mesh_location_t  *ml = _mesh_location + i;

    ml->mesh = mesh;

    if (ml->elt_list != NULL)
      BFT_FREE(ml->elt_list);

    switch(ml->type) {
    case CS_MESH_LOCATION_CELLS:
      selector = mesh->select_cells;
      n_elts_max = ml->mesh->n_cells;
      break;
    case CS_MESH_LOCATION_INTERIOR_FACES:
      selector = mesh->select_i_faces;
      n_elts_max = ml->mesh->n_i_faces;
      break;
    case CS_MESH_LOCATION_BOUNDARY_FACES:
      selector = mesh->select_b_faces;
      n_elts_max = ml->mesh->n_b_faces;
      break;
    case CS_MESH_LOCATION_VERTICES:
      n_elts_max = mesh->n_vertices;
      break;
    default:
      break;
    }

    if (ml->select_str != NULL) {
      if (selector != NULL) {
        BFT_MALLOC(ml->elt_list, n_elts_max, cs_lnum_t);
        int c_id = fvm_selector_get_list(selector,
                                         ml->select_str,
                                         ml->n_elts,
                                         ml->elt_list);
        if (ml->n_elts[0] == n_elts_max && ml->elt_list != NULL)
          BFT_FREE(ml->elt_list);
        else
          BFT_REALLOC(ml->elt_list, ml->n_elts[0], cs_lnum_t);
        if (fvm_selector_n_missing(selector, c_id) > 0) {
          const char *missing
            = fvm_selector_get_missing(selector, c_id, 0);
          cs_base_warn(__FILE__, __LINE__);
          bft_printf(_("The group \"%s\" in the selection criteria:\n"
                       "\"%s\"\n"
                       " does not correspond to any boundary face.\n"),
                     missing, ml->select_str);
        }
      }
      else
        bft_error
          (__FILE__, __LINE__, 0,
           _("A selection criteria is given but no associated selector\n"
             "is available for mesh location %d of type %d."),
           i, (int)ml->type);
    }
    else if (ml->select_fp != NULL)
      ml->select_fp(ml->mesh,
                    i,
                    ml->n_elts,
                    &(ml->elt_list));
    else
      ml->n_elts[0] = n_elts_max;

    ml->n_elts[1] = ml->n_elts[0];
    ml->n_elts[2] = ml->n_elts[0];

    if (ml->type == CS_MESH_LOCATION_CELLS && ml->n_elts[0] == mesh->n_cells) {
      if (mesh->halo != NULL) {
        assert(mesh->halo->n_local_elts == ml->n_elts[0]);
        ml->n_elts[1] += mesh->halo->n_elts[0];
        ml->n_elts[2] += mesh->halo->n_elts[1];
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location.
 *
 * So as to define a subset of mesh entities of a given type, an optional
 * selection criteria may be given.
 *
 * \param[in]  name      name of location to define
 * \param[in]  type      type of location to define
 * \param[in]  criteria  selection criteria for associated elements,
 *                       or NULL
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_define(const char               *name,
                        cs_mesh_location_type_t   type,
                        const char               *criteria)
{
  int m_id = _mesh_location_define(name, type);

  cs_mesh_location_t *ml = _mesh_location + m_id;

  if (criteria != NULL) {
    BFT_MALLOC(ml->select_str, strlen(criteria) + 1, char);
    strcpy(ml->select_str, criteria);
  }

  return m_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location with a associated selection function.
 *
 * So as to define a subset of mesh entities of a given type, a pointer
 * to a selection function may be given.
 *
 * This requires more programming but allows finer control than selection
 * criteria, as the function has access to the complete mesh structure.
 *
 * \param[in]  name  name of location to define
 * \param[in]  type  type of location to define
 * \param[in]  func  pointer to selection function for associated elements,
 *                   or NULL
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_define_by_func(const char                 *name,
                                cs_mesh_location_type_t     type,
                                cs_mesh_location_select_t  *func)
{
  int m_id = _mesh_location_define(name, type);

  cs_mesh_location_t *ml = _mesh_location + m_id;

  ml->select_fp = func;

  return m_id;
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
 * \return  pointer to elements list (0 to n-1 numbering).
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
