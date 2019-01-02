/*============================================================================
 * Handle the "physical" boundary conditions attached to a computational domain
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_boundary.c

  \brief Handle the "physical" boundary conditions attached to a computational
         domain
*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local variables
 *============================================================================*/

static const char
cs_boundary_name[CS_BOUNDARY_N_TYPES][CS_BASE_STRING_LEN] =
  { N_("wall"),
    N_("sliding wall"),
    N_("inlet"),
    N_("outlet"),
    N_("symmetry")
  };

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_boundary_t  *cs_glob_boundaries = NULL ; /* Pointer to the shared boundaries
                                             * on the computational domain */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of boundary faces attached to a wall boundary
 *         condition
 *         Function pointer to mesh location elements selection definition.
 *
 * If non-empty and not containing all elements, a list of elements
 * of the parent mesh belonging to the location should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * This list's lifecycle is then managed by the mesh location object.
 *
 * \param[in]   input        pointer to a structure cast on-the-fly
 * \param[in]   m            pointer to associated mesh structure.
 * \param[in]   location_id  id of associated location.
 * \param[out]  n_elts       number of selected elements
 * \param[out]  elt_list     list of selected elements.
 */
/*----------------------------------------------------------------------------*/

static void
_wall_boundary_selection(void              *input,
                         const cs_mesh_t   *m,
                         int                location_id,
                         cs_lnum_t         *n_elts,
                         cs_lnum_t        **elt_ids)
{
  CS_UNUSED(location_id);

  const cs_boundary_t  *bdy = (const cs_boundary_t *)input;

  bool  *is_wall = NULL;
  BFT_MALLOC(is_wall, m->n_b_faces, bool);

  if (bdy->default_type == CS_BOUNDARY_WALL) {

#   pragma omp parallel for if (m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = true;

    for (int i = 0; i < bdy->n_boundaries; i++) {
      if (bdy->types[i] != CS_BOUNDARY_WALL) {

        const int  z_id = bdy->zone_ids[i];
        const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);

        /* At this stage, zone are not defined contrary to the mesh location
         * So, we retrieve the mesh location information
         */
        const int  ml_id = z->location_id;
        const cs_lnum_t  _n_elts = cs_mesh_location_get_n_elts(ml_id)[0];
        const cs_lnum_t  *_elt_ids = cs_mesh_location_get_elt_ids(ml_id);

        if (_elt_ids == NULL)
          for (cs_lnum_t j = 0; j < _n_elts; j++) is_wall[j] = false;
        else
          for (cs_lnum_t j = 0; j < _n_elts; j++) is_wall[_elt_ids[j]] = false;

      }
    } /* Loop on boundary definitions */

  }
  else { /* Wall is not the default boundary */

#   pragma omp parallel for if (m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = false;

    for (int i = 0; i < bdy->n_boundaries; i++) {
      if (bdy->types[i] == CS_BOUNDARY_WALL) {

        const int  z_id = bdy->zone_ids[i];
        const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);

        /* At this stage, zone are not defined contrary to the mesh location
         * So, we retrieve the mesh location information
         */
        const int  ml_id = z->location_id;
        const cs_lnum_t  _n_elts = cs_mesh_location_get_n_elts(ml_id)[0];
        const cs_lnum_t  *_elt_ids = cs_mesh_location_get_elt_ids(ml_id);

        if (_elt_ids == NULL)
          for (cs_lnum_t j = 0; j < _n_elts; j++) is_wall[j] = true;
        else
          for (cs_lnum_t j = 0; j < _n_elts; j++) is_wall[_elt_ids[j]] = true;

      }
    } /* Loop on boundary definitions */

  } /* Which kind of default boundary ? */

  /* Count the number of boundary faces attached to a wall boundary */
  cs_lnum_t  n_wall_elts = 0;
  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    if (is_wall[i]) n_wall_elts++;

  cs_lnum_t *wall_elts = NULL;
  if (n_wall_elts < m->n_b_faces) {

    /* Fill list  */
    BFT_MALLOC(wall_elts, n_wall_elts, cs_lnum_t);

    cs_lnum_t shift = 0;
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      if (is_wall[i]) wall_elts[shift++] = i;

    assert(shift == n_wall_elts);

  } /* Build elt_ids */

  BFT_FREE(is_wall);

  /* Return pointers */
  *n_elts = n_wall_elts;
  *elt_ids = wall_elts;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the domain boundary condition
 *          This name is also used as a name for zone definition
 *
 * \param[in] type     type of boundary
 *
 * \return the associated boundary name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_boundary_get_name(cs_boundary_type_t  type)
{
  if (type == CS_BOUNDARY_N_TYPES)
    return "Undefined";
  else
    return cs_boundary_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to the given \ref cs_boundary_t
 *         structure
 *
 * \param[in, out]   boundaries   pointer to a structure storing boundary info
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_set_default(cs_boundary_t        *boundaries,
                        cs_boundary_type_t    type)
{
  if (boundaries == NULL)
    return;

  if (type != CS_BOUNDARY_WALL && type != CS_BOUNDARY_SYMMETRY)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid type of default boundary.\n"
                " A valid choice is either \"CS_BOUNDARY_WALL\" or"
                " \"CS_BOUNDARY_SYMMETRY\"."), __func__);

  boundaries->default_type = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a default boundary structure for the computational domain
 *
 * \param[in]        type         default type of boundary to set
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_boundary_t *
cs_boundary_create(cs_boundary_type_t    type)
{
  cs_boundary_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_boundary_t);

  b->default_type = type;
  b->n_boundaries = 0;
  b->zone_ids = NULL;
  b->types = NULL;

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 *
 * \param[in, out]   p_boundaries   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_free(cs_boundary_t    **p_boundaries)
{
  if (*p_boundaries == NULL)
    return;

  cs_boundary_t  *bdy = *p_boundaries;

  BFT_FREE(bdy->types);
  BFT_FREE(bdy->zone_ids);
  BFT_FREE(bdy);
  bdy = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in, out] bdy          pointer to a structure storing boundary info
 * \param[in]      type         type of boundary to set
 * \param[in]      zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_add(cs_boundary_t        *bdy,
                cs_boundary_type_t    type,
                const char           *zone_name)
{
  if (bdy == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Empty boundary structure", __func__);

  const cs_zone_t  *zone = cs_boundary_zone_by_name(zone_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid zone name %s.\n"
                " This zone is not already defined.\n"), __func__, zone_name);

  int  new_id = bdy->n_boundaries;

  /* Add a new boundary for this zone */
  bdy->n_boundaries += 1;

  BFT_REALLOC(bdy->zone_ids, bdy->n_boundaries, int);
  BFT_REALLOC(bdy->types, bdy->n_boundaries, cs_boundary_type_t);

  bdy->zone_ids[new_id] = zone->id;
  bdy->types[new_id] = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_BOUNDARY_WALL zone type
 *
 * \param[in, out]  boundaries    pointer to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_def_wall_zones(cs_boundary_t   *bdy)
{
  if (bdy == NULL)
    return;

  /* Add a new boundary zone (and also a new mesh location) related to all
   * wall boundary faces */
  const char  zone_name[] = CS_BOUNDARY_WALLS_NAME;

  int  flag = CS_BOUNDARY_ZONE_WALL | CS_BOUNDARY_ZONE_PRIVATE;
  int  z_id = cs_boundary_zone_define_by_func(zone_name,
                                              _wall_boundary_selection,
                                              bdy,
                                              flag);

  /* Allow overlay with other boundary zones used to set BCs on transport
     equations for instance (not really needed since zone is private) */
  cs_boundary_zone_set_overlay(z_id, true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 *
 * \param[in] bdy          pointer to a structure storing boundary info
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_log_setup(const cs_boundary_t     *bdy)
{
  if (bdy == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n  Domain boundary by default: %s\n",
                cs_boundary_get_name(bdy->default_type));

  for (int i = 0; i < bdy->n_boundaries; i++) {

    const int  z_id = bdy->zone_ids[i];
    const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);

    /* Count the number of boundary faces related to this definition */
    cs_gnum_t  n_g_elts = (cs_gnum_t)z->n_elts;
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&n_g_elts, 1);

    cs_log_printf(CS_LOG_SETUP, " %s -- %s: %u boundary faces,",
                  z->name, cs_boundary_get_name(bdy->types[i]),
                  (unsigned int)n_g_elts);

  } /* Loop on domain boundaries */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
