/*============================================================================
 * Handle the "physical" boundary conditions attached to a computational domain
 *============================================================================*/

/* VERS */

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

#include "cs_domain_boundary.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_domain_boundary.c

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
cs_domain_boundary_name[CS_DOMAIN_N_BOUNDARY_TYPES][CS_BASE_STRING_LEN] =
  { N_("wall"),
    N_("sliding wall"),
    N_("inlet"),
    N_("outlet"),
    N_("symmetry")
  };

static cs_domain_boundary_type_t cs_domain_boundary_default_type =
  CS_DOMAIN_N_BOUNDARY_TYPES;   /* No boundary set by default */

static int  cs_domain_n_boundaries = 0;
static int  *cs_domain_boundary_zone_ids = NULL;
static cs_domain_boundary_type_t  *cs_domain_boundary_zone_types = NULL;

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
  CS_UNUSED(input);

  cs_lnum_t  n_wall_elts = 0;
  cs_lnum_t *wall_elts = NULL;
  bool  *is_wall = NULL;

  BFT_MALLOC(is_wall, m->n_b_faces, bool);

  if (cs_domain_boundary_default_type == CS_DOMAIN_BOUNDARY_WALL) {

#   pragma omp parallel for if (m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = true;

    for (int i = 0; i < cs_domain_n_boundaries; i++) {
      if (cs_domain_boundary_zone_types[i] != CS_DOMAIN_BOUNDARY_WALL) {

        const int  z_id = cs_domain_boundary_zone_ids[i];
        const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  _n_faces = z->n_elts;
        const cs_lnum_t  *_face_ids = z->elt_ids;

        for (cs_lnum_t j = 0; j < _n_faces; j++)
          is_wall[_face_ids[j]] = false;

      }
    } /* Loop on boundary definitions */

  }
  else { /* Wall is not the default boundary */

#   pragma omp parallel for if (m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = false;

    for (int i = 0; i < cs_domain_n_boundaries; i++) {
      if (cs_domain_boundary_zone_types[i] == CS_DOMAIN_BOUNDARY_WALL) {

        const int  z_id = cs_domain_boundary_zone_ids[i];
        const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  _n_faces = z->n_elts;
        const cs_lnum_t  *_face_ids = z->elt_ids;

        for (cs_lnum_t j = 0; j < _n_faces; j++)
          is_wall[_face_ids[j]] = true;

      }
    } /* Loop on boundary definitions */

  } /* Which default ? */

  /* Count the number of boundary faces attached to a wall boundary */
  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    if (is_wall[i]) n_wall_elts++;

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
cs_domain_boundary_get_name(cs_domain_boundary_type_t  type)
{
  if (type == CS_DOMAIN_N_BOUNDARY_TYPES)
    return "Not set";
  else
    return cs_domain_boundary_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to this domain
 *
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_set_default(cs_domain_boundary_type_t    type)
{
  cs_domain_boundary_default_type = type;

  if (type != CS_DOMAIN_BOUNDARY_WALL && type != CS_DOMAIN_BOUNDARY_SYMMETRY)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid type of boundary by default.\n"
                " A valid choice is CS_DOMAIN_BOUNDARY_WALL or"
                " CS_DOMAIN_BOUNDARY_SYMMETRY."), __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the default boundary associated to the computational domain
 *
 * \return the type of the domain boundary defined by default
 */
/*----------------------------------------------------------------------------*/

cs_domain_boundary_type_t
cs_domain_boundary_get_default(void)
{
  return cs_domain_boundary_default_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_free(void)
{
  BFT_FREE(cs_domain_boundary_zone_types);
  BFT_FREE(cs_domain_boundary_zone_ids);
  cs_domain_n_boundaries = 0;
  cs_domain_boundary_default_type = CS_DOMAIN_N_BOUNDARY_TYPES;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in] type         type of boundary to set
 * \param[in] zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_add(cs_domain_boundary_type_t    type,
                       const char                  *zone_name)
{
  const cs_zone_t  *zone = cs_boundary_zone_by_name(zone_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid zone name %s.\n"
                " This zone is not already defined.\n"), __func__, zone_name);

  int  new_id = cs_domain_n_boundaries;

  /* Add a new boundary for this zone */
  cs_domain_n_boundaries += 1;
  BFT_REALLOC(cs_domain_boundary_zone_ids, cs_domain_n_boundaries, int);
  BFT_REALLOC(cs_domain_boundary_zone_types, cs_domain_n_boundaries,
              cs_domain_boundary_type_t);

  cs_domain_boundary_zone_ids[new_id] = zone->id;
  cs_domain_boundary_zone_types[new_id] = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_DOMAIN_BOUNDARY_WALL zone type
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_def_wall_zones(void)
{
  /* Add a new boundary zone (and also a new mesh location) related to all
     wall boundary faces */
  const char  zone_name[] = "cs_domain_boundary_walls";

  int flag = CS_BOUNDARY_ZONE_WALL | CS_BOUNDARY_ZONE_PRIVATE;

  int  z_id = cs_boundary_zone_define_by_func(zone_name,
                                              _wall_boundary_selection,
                                              NULL, /* No context */
                                              flag);

  /* Allow overlay with other boundary zones used to set BCs on transport
     equations for instance (not really needed since zone is private) */
  cs_boundary_zone_set_overlay(z_id, true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_log_setup(void)
{
  cs_domain_boundary_type_t  db_default_type = cs_domain_boundary_get_default();
  cs_log_printf(CS_LOG_SETUP, "\n  Domain boundary by default: %s\n",
                cs_domain_boundary_get_name(db_default_type));

  for (int i = 0; i < cs_domain_n_boundaries; i++) {

    const int  z_id = cs_domain_boundary_zone_ids[i];
    const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);

    /* Count the number of boundary faces related to this definition */
    cs_gnum_t  n_g_elts = (cs_gnum_t)z->n_elts;
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&n_g_elts, 1);

    cs_domain_boundary_type_t  db_type = cs_domain_boundary_zone_types[i];
    cs_log_printf(CS_LOG_SETUP, " %s -- %s: %u boundary faces,",
                  z->name, cs_domain_boundary_get_name(db_type),
                  (unsigned int)n_g_elts);

  } /* Loop on domain boundaries */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
