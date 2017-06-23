/*============================================================================
 * Mesh locations management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

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

struct _cs_mesh_location_t {

  char                        name[32];     /* Name */

  const cs_mesh_t            *mesh;         /* Pointer to associated mesh */

  cs_mesh_location_type_t     type;         /* Location type */

  char                       *select_str;   /* String */
  cs_mesh_location_select_t  *select_fp;    /* Function pointer */
  void                       *select_input; /* Optional input for the function pointer */
  int                         n_sub_ids;    /* Number of mesh location ids
                                               used to build this location */
  int                        *sub_ids;      /* List of mesh location ids */
  bool                        complement;   /* Take the complement ? */
  bool                        explicit_ids; /* Need explicit ids ? */

  cs_lnum_t                   n_elts[3];    /* Number of associated elements:
                                               0: local,
                                               1: with standard ghost elements,
                                               2: with extended ghost elements
                                            */

  cs_lnum_t                  *elt_list;     /* List of associated elements,
                                               (0 to n-1 numbering) if non
                                               trivial (i.e. a subset) */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names associated with location types */

const char  *cs_mesh_location_type_name[] = {N_("none"),
                                             N_("cells"),
                                             N_("interior faces"),
                                             N_("boundary faces"),
                                             N_("vertices"),
                                             N_("faces"),
                                             N_("edges"),
                                             N_("particles"),
                                             N_("other")};

/* Location definitions */

static int                  _n_mesh_locations_max = 0;
static int                  _n_mesh_locations = 0;
static cs_mesh_location_t  *_mesh_location = NULL;

/* Shared explicit ids */

static cs_lnum_t            _explicit_ids_size = 0;
static cs_lnum_t           *_explicit_ids = NULL;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Get a pointer to a mesh location by its id.
 *
 * \param[in]  id         id of mesh location
 *
 * \return  a pointer to the associated mesh location
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_location_t *
_mesh_location_by_id(int  id)
{
  const cs_mesh_location_t  *retval = NULL;

  if (id < 0 || id > _n_mesh_locations)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested mesh location\n"
                "%d is not defined.\n"), id);
  else
    retval = _mesh_location + id;

  return retval;
}

/*----------------------------------------------------------------------------
 * \brief Get a const pointer to a mesh location by its id.
 *
 * \param[in]  id         id of mesh location
 *
 * \return  a pointer to the associated mesh location
 */
/*----------------------------------------------------------------------------*/

static const cs_mesh_location_t *
_const_mesh_location_by_id(int  id)
{
  const cs_mesh_location_t  *retval = NULL;

  if (id < 0 || id > _n_mesh_locations)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested mesh location\n"
                "%d is not defined.\n"), id);
  else
    retval = _mesh_location + id;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the related location id from the location name
 *
 * \param[in]  ref_name    name of the location to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

static int
_find_id_by_name(const char  *ref_name)
{
  int  ml_id = -1, reflen = strlen(ref_name);

  for (int i = 0; i < _n_mesh_locations; i++) {

    cs_mesh_location_t  *ml = _mesh_location + i;
    int len = strlen(ml->name);

    if (reflen == len) {
      if (strcmp(ref_name, ml->name) == 0) {
        ml_id = i;
        break;
      }
    }

  } /* Loops on mesh locations */

  return ml_id;
}

/*----------------------------------------------------------------------------
 * \brief Define a new mesh location.
 *
 * If a list of associated elements is given (defining a subset of a main
 * location), its ownership is transferred to the mesh location.
 *
 * \param[in]  name      name of location to define
 * \param[in]  type      type of location to define
 *
 * \return   id of the newly defined mesh location
 */
/*----------------------------------------------------------------------------*/

static int
_mesh_location_define(const char               *name,
                      cs_mesh_location_type_t   type)
{
  int  i;

  cs_mesh_location_t  *ml = NULL;

  /* Check if there is already a mesh location with the same name */
  int id = _find_id_by_name(name);

  if (id == -1) /* Not found in the current list */
    id = _n_mesh_locations;

  else { /* A mesh location shares the same name */
    ml = _mesh_location + id;
    if (ml->type != type)
      bft_error(__FILE__, __LINE__, 0,
                _(" The mesh location %s is already defined as a mesh location"
                  " but with a different type.\n"
                  " Please check your settings."), name);
    return id;
  }

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
  ml->select_input = NULL;
  ml->n_sub_ids = 0;
  ml->sub_ids = NULL;
  ml->complement = false;
  ml->explicit_ids = false;

  for (i = 0; i < 3; i++)
    ml->n_elts[i] = 0;
  ml->elt_list = NULL;

  return id;
}

/*----------------------------------------------------------------------------
 * \brief Build a list of element ids from a list of existing mesh location
 *        ids
 *
 * \param[in]  name      name of location to define
 * \param[in]  type      type of location to define
 */
/*----------------------------------------------------------------------------*/

static void
_build_by_ml_ids(cs_mesh_location_t  *ml)
{
  int  i, count;

  cs_lnum_t  n_elts_max = 0;

  switch(ml->type) {
  case CS_MESH_LOCATION_CELLS:
    n_elts_max = ml->mesh->n_cells;
    break;
  case CS_MESH_LOCATION_FACES:
    n_elts_max = ml->mesh->n_b_faces + ml->mesh->n_i_faces;
    break;
  case CS_MESH_LOCATION_INTERIOR_FACES:
    n_elts_max = ml->mesh->n_i_faces;
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    n_elts_max = ml->mesh->n_b_faces;
    break;
  case CS_MESH_LOCATION_VERTICES:
    n_elts_max = ml->mesh->n_vertices;
    break;
  default:
    assert(0);
    break;
  }

  if (ml->n_sub_ids == 1 && ml->complement == false) {

    int  sub_id = ml->sub_ids[0];
    assert(sub_id < _n_mesh_locations);
    cs_mesh_location_t  *sub_ml = _mesh_location + sub_id;
    assert(sub_ml->type == ml->type);

    ml->n_elts[0] = sub_ml->n_elts[0];
    if (sub_ml->elt_list != NULL) { /* Copy */
      BFT_MALLOC(ml->elt_list, ml->n_elts[0], cs_lnum_t);
      memcpy(ml->elt_list, sub_ml->elt_list, ml->n_elts[0]*sizeof(cs_lnum_t));
    }

  }
  else {

    bool  *flag = NULL;

    /* Initialize flag */
    BFT_MALLOC(flag, n_elts_max, bool);
    for (i = 0; i < n_elts_max; i++)
      flag[i] = false;

    /* Build flag */
    for (int ii = 0; ii < ml->n_sub_ids; ii++) {

      int  sub_id = ml->sub_ids[ii];
      assert(sub_id < _n_mesh_locations);
      cs_mesh_location_t  *sub_ml = _mesh_location + sub_id;
      assert(sub_ml->type == ml->type);

      if (sub_ml->elt_list == NULL)
        for (i = 0; i < n_elts_max; i++)
          flag[i] = true;
      else
        for (i = 0; i < sub_ml->n_elts[0]; i++)
          flag[sub_ml->elt_list[i]] = true;

    } /* Loop on (sub) mesh locations */

    if (ml->complement) {
      for (i = 0; i < n_elts_max; i++) {
        if (flag[i])
          flag[i] = false;
        else
          flag[i] = true;
      }
    }

    /* Compute n_elts[0] */
    count = 0;
    for (i = 0; i < n_elts_max; i++)
      if (flag[i]) count++;
    ml->n_elts[0] = count;

    /* Build elt_list */
    if (ml->n_elts[0] != 0 && ml->n_elts[0] != n_elts_max) {
      BFT_MALLOC(ml->elt_list, ml->n_elts[0], cs_lnum_t);
      count = 0;
      for (i = 0; i < n_elts_max; i++)
        if (flag[i]) ml->elt_list[count++] = i;
    }

    BFT_FREE(flag);

  } /* If simple case (n_sub_ids = 1 and no complement) or not */

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
 * By default, 7 mesh locations are built, matching the 7 first values of
 * the cs_mesh_location_type_t enum: CS_MESH_LOCATION_NONE for global
 * values, CS_MESH_LOCATION_CELLS for the cells of the (default) global mesh,
 * CS_MESH_LOCATION_INTERIOR_FACES and CS_MESH_LOCATION_BOUNDARY_FACES for
 * its faces, and CS_MESH_LOCATION_VERTICES for its vertices.
 * CS_MESH_LOCATION_FACES and a placeholder for CS_MESH_LOCATION_EDGES are
 * also added for CDO discretizations.
 *
 * Locations should then be built once the global mesh is complete, and
 * its halo structures completed.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_initialize(void)
{
  cs_mesh_location_add(N_("global"),
                       CS_MESH_LOCATION_NONE,
                       NULL);
  cs_mesh_location_add(N_("cells"),
                       CS_MESH_LOCATION_CELLS,
                       NULL);
  cs_mesh_location_add(N_("interior_faces"),
                       CS_MESH_LOCATION_INTERIOR_FACES,
                       NULL);
  cs_mesh_location_add(N_("boundary_faces"),
                       CS_MESH_LOCATION_BOUNDARY_FACES,
                       NULL);
  cs_mesh_location_add(N_("vertices"),
                       CS_MESH_LOCATION_VERTICES,
                       NULL);
  cs_mesh_location_add(N_("faces"),
                       CS_MESH_LOCATION_FACES,
                       NULL);
  cs_mesh_location_add(N_("edges (type)"),
                       CS_MESH_LOCATION_EDGES,
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

  BFT_FREE(_explicit_ids);

  for (i = 0; i < _n_mesh_locations; i++) {
    cs_mesh_location_t  *ml = _mesh_location + i;
    BFT_FREE(ml->elt_list);
    BFT_FREE(ml->select_str);
    BFT_FREE(ml->sub_ids);
  }

  _explicit_ids_size = 0;
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
  int  ml_id;
  int id_start = 0, id_end = _n_mesh_locations;

  cs_lnum_t  explicit_ids_size = 0;

  assert(mesh != NULL);
  assert(   mesh->halo != NULL
         || (cs_glob_n_ranks == 1 && mesh->periodicity == NULL));

  if (id >= 0) {
    id_start = id;
    if (id < _n_mesh_locations)
      id_end = id + 1;
  }

  for (ml_id = id_start; ml_id < id_end; ml_id++) {

    int n_elts_max = 0;
    fvm_selector_t *selector = NULL;
    cs_mesh_location_t  *ml = _mesh_location + ml_id;

    ml->mesh = mesh;

    if (ml->elt_list != NULL)
      BFT_FREE(ml->elt_list);

    switch(ml->type) {
    case CS_MESH_LOCATION_CELLS:
      selector = mesh->select_cells;
      n_elts_max = ml->mesh->n_cells;
      break;
    case CS_MESH_LOCATION_FACES:
      n_elts_max = ml->mesh->n_b_faces + ml->mesh->n_i_faces;
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
                                         0,
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
        bft_error(__FILE__, __LINE__, 0,
                  _("A selection criteria is given but no associated selector\n"
                    "is available for mesh location %d of type %d."),
                  ml_id, (int)ml->type);
    }
    else if (ml->select_fp != NULL) {
      ml->select_fp(ml->select_input,
                    ml->mesh,
                    ml_id,
                    ml->n_elts,
                    &(ml->elt_list));
      if (ml->elt_list != NULL)
        cs_sort_lnum(ml->elt_list, ml->n_elts[0]);
    }
    else if (ml->n_sub_ids > 0 && ml->sub_ids != NULL)
      _build_by_ml_ids(ml);
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

    if (ml->explicit_ids && explicit_ids_size < ml->n_elts[0])
      explicit_ids_size = ml->n_elts[0];

  } /* Loop on mesh locations */

  if (explicit_ids_size != _explicit_ids_size) {
    if (id == 0 || explicit_ids_size > _explicit_ids_size) {
      cs_lnum_t s_id = (id == 0) ? 0 : _explicit_ids_size;
      _explicit_ids_size = explicit_ids_size;
      BFT_REALLOC(_explicit_ids, _explicit_ids_size, cs_lnum_t);
      for (cs_lnum_t i = s_id; i < _explicit_ids_size; i++)
        _explicit_ids[i] = i;
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
cs_mesh_location_add(const char                *name,
                     cs_mesh_location_type_t    type,
                     const char                *criteria)
{
  int  ml_id = _mesh_location_define(name, type);
  cs_mesh_location_t  *ml = _mesh_location + ml_id;

  if (criteria != NULL) {
    BFT_MALLOC(ml->select_str, strlen(criteria) + 1, char);
    strcpy(ml->select_str, criteria);
  }

  return ml_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location with an associated selection function.
 *
 * So as to define a subset of mesh entities of a given type, a pointer
 * to a selection function may be given.
 *
 * This requires more programming but allows finer control than selection
 * criteria, as the function has access to the complete mesh structure.
 *
 * \param[in]  name        name of location to define
 * \param[in]  type        type of location to define
 * \param[in]  func        pointer to selection function for associated
 *                         elements, or NULL
 * \param[in, out]  input  pointer to optional (untyped) value
 *                         or structure.
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_add_by_func(const char                 *name,
                             cs_mesh_location_type_t     type,
                             cs_mesh_location_select_t  *func,
                             void                       *input)
{
  int  ml_id = _mesh_location_define(name, type);
  cs_mesh_location_t  *ml = _mesh_location + ml_id;

  ml->select_fp = func;
  ml->select_input = input;

  return ml_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location.
 *
 * So as to define a subset of mesh entities of a given type, a list of ids
 * related to existing mesh locations may be given
 *
 * \param[in]  name        name of location to define
 * \param[in]  type        type of location to define
 * \param[in]  n_ml_ids    number of mesh location ids
 * \param[in]  ml_ids      list of mesh location ids
 * \param[in]  complement  take the complement of the selected entities if true
 *
 * \return  id of newly created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_add_by_union(const char               *name,
                              cs_mesh_location_type_t   type,
                              int                       n_ml_ids,
                              const int                *ml_ids,
                              bool                      complement)
{
  int  ml_id = _mesh_location_define(name, type);
  cs_mesh_location_t  *ml = _mesh_location + ml_id;

  ml->complement = complement;
  ml->n_sub_ids = n_ml_ids;
  if (ml->n_sub_ids > 0) {
    BFT_MALLOC(ml->sub_ids, ml->n_sub_ids, int);
    for (int i = 0; i < ml->n_sub_ids; i++)
      ml->sub_ids[i] = ml_ids[i];
  }

  return ml_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the related location id from the location name
 *
 * \param[in]  ref_name    name of the location to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_location_get_id_by_name(const char  *ref_name)
{
  return _find_id_by_name(ref_name);
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
cs_mesh_location_get_name(int  id)
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
cs_mesh_location_get_type(int  id)
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
 * \deprecated  Use \ref cs_mesh_location_get_elt_ids_try or
 * \ref cs_mesh_location_get_elt_ids instead.
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to elements array (0 to n-1 numbering).
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_elt_list(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->elt_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's element ids, if present.
 *
 * An array of element ids is returned if the location is a subset of a main
 * main location type, and NULL is returned when the id array would map to
 * the identity function (i.e. {0, 1, 2, ..., n_elts-1}).
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to elements array (0 to n-1 numbering).
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_elt_ids_try(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->elt_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's element ids.
 *
 * This function may only be used with a given location if
 * \ref cs_mesh_location_set_explicit_ids has been used to indicate
 * explicit ids are needed for this location type.
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to elements array (0 to n-1 numbering).
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_mesh_location_get_elt_ids(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  if (! (   ml->explicit_ids
         || (_mesh_location + ml->type)->explicit_ids))
    bft_error(__FILE__, __LINE__, 0,
              _("Explicit ids have not been built for mesh location %d\n"
                "or its base type.\n"
                "Use cs_mesh_location_set_explicit_ids."), id);

  const cs_lnum_t *retval = ml->elt_list;
  if (retval == NULL)
    retval = _explicit_ids;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's selection criteria string
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to mesh location selection criteria, or NULL
 */
/*----------------------------------------------------------------------------*/

const char *
cs_mesh_location_get_selection_string(int  id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->select_str;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a mesh location's selection function pointer
 *
 * \param[in]  id  id of mesh location
 *
 * \return  pointer to mesh location selection function pointer, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_mesh_location_select_t *
cs_mesh_location_get_selection_function(int  id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);

  return ml->select_fp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if \ref cs_mesh_location_get_elt_ids always returns
 *         explicit element ids for a given mesh location type.
 *
 * \param[in]  id  id or type of location
 *
 * \return true if explicit element ids are needed, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_mesh_location_get_explicit_ids(int id)
{
  const cs_mesh_location_t  *ml = _const_mesh_location_by_id(id);
  return ml->explicit_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set behavior of \ref cs_mesh_location_get_elt_ids for mesh
 *         mesh location type.
 *
 * \param[in]  id                id or type of location
 * \param[in]  explicit_elt_ids  indicate if explicit element ids are needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_location_set_explicit_ids(int   id,
                                  bool  explicit_elt_ids)
{
  cs_mesh_location_t  *ml = _mesh_location_by_id(id);
  ml->explicit_ids = explicit_elt_ids;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
