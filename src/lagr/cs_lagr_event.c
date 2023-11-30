/*============================================================================
 * Lagrangian particle event model
 *============================================================================*/

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

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_timer_stats.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_event.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_event.c
        Particle event management.

Particle events allow keeping track of particle events suh as boundary
interactions using a structure and API similar to that used for particles.

Particles events could actually be stored as particle data, given the
necessary additional attributes defined here (i.e. using a different
attribute map), and could be refactored in the future based on feedback.

As of now, it was chosen to use an independent structure, so as to avoid
confusion regarding the more restricted uses and (temporary) lifetime
of particle events, which should need to be maintained only as a buffer
mechanism for event-based statistics.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  N_GEOL 13
#define  CS_LAGR_MIN_COMM_BUF_SIZE  8

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/* keys to sort attributes by type. */

typedef enum {
  CS_LAGR_P_NULLV,
  CS_LAGR_P_RV,
  CS_LAGR_P_IV,
} _array_map_id_t;

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* Private data associated to each particle event */
/* -----------------------------------------------*/

/* Event data value */
/*------------------*/

union cs_lagr_value_t {
  cs_lnum_t      l; /* v_lnum_t */
  cs_gnum_t      g; /* v_gnum_t */
  cs_real_t      f; /* v_real_t */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Enumerator names */
/* Enumerator names */

const char *_event_attribute_name[] = {
  "flag",
  "cell_id",
  "face_id",
  "velocity_post",
  "<none>"};

/* Global particle event attributes map */

static cs_lagr_event_attribute_map_t  *_e_attr_map = NULL;

/* Quick mapping from particle attributes to event attributes
   to allow for quick copy from particle to event */

static int   _n_mapped_part_attr = 0;
static int  *_mapped_part_attr = NULL;

static cs_lagr_event_set_t  *_boundary_events = NULL;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * Compute new extents to ensure alignment of data
 *
 * returns:
 *   padded extents ensuring alignement
 *----------------------------------------------------------------------------*/

static size_t
_align_extents(size_t  size)
{
  size_t retval = size;

  size_t align_size = sizeof(union cs_lagr_value_t);

  size_t r = size % align_size;
  if (r > 0)
    retval += (align_size - r);

  return retval;
}

/*----------------------------------------------------------------------------*
 * Map particle event attributes for a given configuration.
 *
 * parameters:
 *   attr_keys   <-> keys to sort attributes by Fortran array and index
 *                   for each attribute: array, index in array, count
 *
 * returns:
 *   pointer to structure mapping particle event attributes
 *----------------------------------------------------------------------------*/

static cs_lagr_event_attribute_map_t *
_create_attr_map(cs_lnum_t attr_keys[CS_LAGR_N_E_ATTRIBUTES][3])
{
  cs_lagr_event_attribute_t attr;
  cs_lnum_t *order;

  cs_lagr_event_attribute_map_t  *e_am;

  BFT_MALLOC(e_am, 1, cs_lagr_event_attribute_map_t);

  e_am->lb = 0;
  e_am->extents = e_am->lb;

  for (attr = 0; attr < CS_LAGR_N_E_ATTRIBUTES; attr++) {
    e_am->size[attr] = 0;
    e_am->datatype[attr] = CS_REAL_TYPE;
    e_am->displ[attr] = -1;
    e_am->count[attr] = 1;
  }

  BFT_MALLOC(order, CS_LAGR_N_E_ATTRIBUTES, cs_lnum_t);

  cs_order_lnum_allocated_s(NULL,
                            (const cs_lnum_t *)attr_keys,
                            3,
                            order,
                            CS_LAGR_N_E_ATTRIBUTES);

  int array_prev = 0;

  /* Now loop on ordered attributes */

  for (int i = 0; i < CS_LAGR_N_E_ATTRIBUTES; i++) {

    cs_datatype_t datatype = CS_REAL_TYPE;

    attr = order[i];

    e_am->datatype[attr] = CS_DATATYPE_NULL;
    e_am->displ[attr] =-1;
    e_am->count[attr] = 0;

    if (attr_keys[attr][0] < 1) continue;

    /* Behavior depending on array */

    switch(attr_keys[attr][0]) {
    case CS_LAGR_P_RV:
      break;
    case CS_LAGR_P_IV:
      datatype = CS_LNUM_TYPE;
      break;
    default:
      continue;
    }

    /* Add padding for alignment when changing array */

    if (attr_keys[attr][0] != array_prev) {
      e_am->extents = _align_extents(e_am->extents);
      array_prev = attr_keys[attr][0];
    }

    /* Add attribute to map */

    e_am->displ[attr] = e_am->extents;
    e_am->count[attr] = attr_keys[attr][2];
    e_am->datatype[attr] = datatype;
    e_am->size[attr] =   e_am->count[attr]
                       * cs_datatype_size[e_am->datatype[attr]];

    e_am->extents += e_am->size[attr];

  }

  e_am->extents = _align_extents(e_am->extents);

  BFT_FREE(order);

  return e_am;
}

/*----------------------------------------------------------------------------*
 * Free particle event attributes for a given configuration.
 *----------------------------------------------------------------------------*/

static void
_destroy_attr_map(cs_lagr_event_attribute_map_t  **e_am)
{
  if (*e_am != NULL) {
    BFT_FREE(*e_am);
  }
}

/*----------------------------------------------------------------------------
 * Allocate a cs_lagr_particle_event_set_t structure.
 *
 * parameters:
 *   n_events_max <-- local max. number of particles
 *   e_am         <-- particle event attributes map
 *
 * returns:
 *   a new allocated cs_lagr_event_set_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_event_set_t *
_create_event_set(cs_lnum_t                             n_events_max,
                  const cs_lagr_event_attribute_map_t  *e_am)

{
  cs_lagr_event_set_t  *new_set = NULL;

  if (n_events_max == 0)
    return NULL;

  BFT_MALLOC(new_set, 1, cs_lagr_event_set_t);

  BFT_MALLOC(new_set->e_buffer, n_events_max * e_am->extents, unsigned char);

  new_set->n_events     = 0;
  new_set->n_events_max = n_events_max;

  assert(n_events_max >= 1);

  new_set->e_am = e_am;

  return new_set;
}

/*----------------------------------------------------------------------------
 * Dump an events structure
 *
 * parameter
 *   particles   <-- cs_lagr_event_set_t structure to dump
 *   particle_id <-- id of particle to dump
 *----------------------------------------------------------------------------*/

static void
_dump_event(const cs_lagr_event_set_t  *events,
            cs_lnum_t                   event_id)
{
  const cs_lagr_event_attribute_map_t *am = events->e_am;

  bft_printf("  event: %lu\n", (unsigned long)event_id);

  bft_printf("    values:\n");

  for (cs_lagr_event_attribute_t attr = 0;
       attr < CS_LAGR_N_E_ATTRIBUTES;
       attr++) {
    if (am->count[attr] > 0) {
      const char *attr_name = cs_lagr_event_get_attr_name(attr);
      switch (am->datatype[attr]) {
      case CS_LNUM_TYPE:
        {
          const cs_lnum_t *v
            = cs_lagr_events_attr_const(events, event_id, attr);
          bft_printf("      %24s: %10ld\n", attr_name, (long)v[0]);
          for (int i = 1; i < am->count[attr]; i++)
            bft_printf("      %24s: %10ld\n", " ", (long)v[i]);
        }
        break;
      case CS_REAL_TYPE:
        {
          const cs_real_t *v
            = cs_lagr_events_attr_const(events, event_id, attr);
          bft_printf("      %24s: %10.3g\n", attr_name, v[0]);
          for (int i = 1; i < am->count[attr]; i++)
            bft_printf("      %24s: %10.3g\n", " ", v[i]);
        }
        break;
      default:
        break;
      }
    }
  }
  bft_printf("\n");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle map based on defined options.
 *
 * This function should only be called after
 * \ref cs_lagr_particle_attr_initialize,
 * as it may use elements from the main particle attributes map.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_initialize(void)
{
  int  i;

  int loc_count = 0;

  cs_lnum_t attr_keys[CS_LAGR_N_E_ATTRIBUTES][3];

  /* Initialize global parameter relative to the lagrangian module */

  /* Set indexes */

  for (i = 0; i < CS_LAGR_N_E_ATTRIBUTES; i++) {
    attr_keys[i][0] = CS_LAGR_P_NULLV;
    attr_keys[i][1] = 0;
    attr_keys[i][2] = 0;
  }

  /* Copy some selected attributes from particle data */

  cs_lagr_attribute_t default_p_attrs[]
    = {CS_LAGR_STAT_WEIGHT,
       CS_LAGR_RESIDENCE_TIME,
       CS_LAGR_MASS,
       CS_LAGR_DIAMETER,
       CS_LAGR_SHAPE,
       CS_LAGR_ORIENTATION,
       CS_LAGR_QUATERNION,
       CS_LAGR_RADII,
       CS_LAGR_ANGULAR_VEL,
       CS_LAGR_EULER,
       CS_LAGR_SHAPE_PARAM,
       CS_LAGR_TAUP_AUX,
       CS_LAGR_COORDS,
       CS_LAGR_VELOCITY,
       CS_LAGR_YPLUS,
       CS_LAGR_INTERF,
       CS_LAGR_MARKO_VALUE,
       CS_LAGR_FOULING_INDEX,
       CS_LAGR_TEMPERATURE,
       CS_LAGR_FLUID_TEMPERATURE,
       CS_LAGR_CP,
       CS_LAGR_WATER_MASS,
       CS_LAGR_COAL_MASS,
       CS_LAGR_COKE_MASS,
       CS_LAGR_SHRINKING_DIAMETER,
       CS_LAGR_STAT_CLASS,
       CS_LAGR_USER};

  _n_mapped_part_attr = 0;

  int n_attrs = sizeof(default_p_attrs) / sizeof(cs_lagr_attribute_t);

  const cs_lagr_attribute_map_t *p_am = cs_lagr_particle_get_attr_map();

  for (i = 0; i < n_attrs; i++) {
    int j = default_p_attrs[i];
    int count = p_am->count[0][j];
    if (count > 0) {
      if (p_am->datatype[j] == CS_REAL_TYPE)
        attr_keys[j][0] = CS_LAGR_P_RV;
      else if (p_am->datatype[j] == CS_LNUM_TYPE)
        attr_keys[j][0] = CS_LAGR_P_IV;
      if (attr_keys[j][0] != CS_DATATYPE_NULL) {
        attr_keys[j][1] = ++loc_count;
        attr_keys[j][2] = count;
        _n_mapped_part_attr += 1;
      }
    }
  }

  BFT_REALLOC(_mapped_part_attr, _n_mapped_part_attr, int);
  _n_mapped_part_attr = 0;

  for (i = 0; i < n_attrs; i++) {
    int j = default_p_attrs[i];
    if (attr_keys[j][0] != CS_DATATYPE_NULL) {
      _mapped_part_attr[_n_mapped_part_attr] = j;
      _n_mapped_part_attr += 1;
    }
  }

  /* Now handle event-specific attributes */

  attr_keys[CS_LAGR_E_FLAG][0] = CS_LAGR_P_IV;
  attr_keys[CS_LAGR_E_FLAG][1] = ++loc_count;
  attr_keys[CS_LAGR_E_FLAG][2] = 1;

  attr_keys[CS_LAGR_E_CELL_ID][0] = CS_LAGR_P_IV;
  attr_keys[CS_LAGR_E_CELL_ID][1] = ++loc_count;
  attr_keys[CS_LAGR_E_CELL_ID][2] = 1;

  attr_keys[CS_LAGR_E_FACE_ID][0] = CS_LAGR_P_IV;
  attr_keys[CS_LAGR_E_FACE_ID][1] = ++loc_count;
  attr_keys[CS_LAGR_E_FACE_ID][2] = 1;

  attr_keys[CS_LAGR_E_VELOCITY][0] = CS_LAGR_P_RV;
  attr_keys[CS_LAGR_E_VELOCITY][1] = ++loc_count;
  attr_keys[CS_LAGR_E_VELOCITY][2] = 3;

  /* Default count of 1 */

  for (i = 0; i < CS_LAGR_N_E_ATTRIBUTES; i++) {
    if (attr_keys[i][1] > 0 && attr_keys[i][2] == 0)
      attr_keys[i][2] = 1;
    else if (attr_keys[i][1] < 1)
      attr_keys[i][0] = 0;
  }

  /* Build mappings
     (in the future, they should be created first, then marked,
     then built) */

  _e_attr_map = _create_attr_map(attr_keys);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy main particle set and map if they exist.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_finalize(void)
{
  if (_boundary_events != NULL)
    cs_lagr_event_set_destroy(&_boundary_events);

  BFT_FREE(_mapped_part_attr);
  _n_mapped_part_attr = 0;

  _destroy_attr_map(&_e_attr_map);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return const pointer to the main particle event attribute
 *         map structure.
 *
 * \return pointer to current particle event attrbute map, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lagr_event_attribute_map_t *
cs_lagr_event_get_attr_map(void)
{
  const cs_lagr_event_attribute_map_t *e_am = _e_attr_map;
  return e_am;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return name associated with a given attribute.
 *
 * \param[in]   attr   event attribute
 */
/*----------------------------------------------------------------------------*/

const char *
cs_lagr_event_get_attr_name(cs_lagr_event_attribute_t   attr)
{
  const char *retval = _event_attribute_name[  CS_LAGR_N_E_ATTRIBUTES
                                             - CS_LAGR_N_ATTRIBUTES];

  if (attr >= 0) {
    if ((int)attr < CS_LAGR_N_ATTRIBUTES)
      retval = cs_lagr_attribute_name[attr];
    else {
      if (attr < CS_LAGR_N_E_ATTRIBUTES)
        retval = _event_attribute_name[attr - CS_LAGR_N_ATTRIBUTES];
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * Create a cs_lagr_event_set_t structure.
 *
 * \return pointer to event set
 */
/*----------------------------------------------------------------------------*/

cs_lagr_event_set_t  *
cs_lagr_event_set_create(void)
{
  cs_lagr_event_set_t  *events = _create_event_set(256, _e_attr_map);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n EVENT SET AFTER CREATION\n");
  cs_lagr_event_set_dump(events);
#endif

  return events;
}

/*----------------------------------------------------------------------------*/
/*!
 * Destroy a cs_lagr_event_set_t structure.
 *
 * \param[in, out]  events  pointer to pointer to event set to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_set_destroy(cs_lagr_event_set_t  **events)
{
  if (events != NULL) {

    cs_lagr_event_set_t *_set = *events;
    BFT_FREE(_set->e_buffer);

    BFT_FREE(*events);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data extents for a given event attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * \param[in]   events     associated event set
 * \param[in]   attr       event attribute
 * \param[out]  extents    size (in bytes) of event structure, or NULL
 * \param[out]  size       size (in bytes) of attribute in event structure,
 *                         or NULL
 * \param[out]  displ      displacement (in bytes) in event structure,
 *                         or NULL
 * \param[out]  datatype   datatype of associated attribute, or NULL
 * \param[out]  count      number of type values associated with attribute,
 *                         or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_get_attr_info(const cs_lagr_event_set_t  *events,
                            cs_lagr_event_attribute_t   attr,
                            size_t                     *extents,
                            size_t                     *size,
                            ptrdiff_t                  *displ,
                            cs_datatype_t              *datatype,
                            int                        *count)
{
  if (extents)
    *extents = events->e_am->extents;
  if (size)
    *size = events->e_am->size[attr];
  if (displ)
    *displ = events->e_am->displ[attr];
  if (datatype)
    *datatype = events->e_am->datatype[attr];
  if (count)
    *count = events->e_am->count[attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if an event attribute is in a valid range.
 *
 * If this is not the case, a fatal error is provoked.

 * \param[in]   attr       event attribute
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_attr_in_range(int  attr)
{
  if (attr < 0 || attr >= CS_LAGR_N_E_ATTRIBUTES)
    bft_error(__FILE__, __LINE__,0,
              _("Out-of range attribute type: %d"),
              (int)attr);
}

/*----------------------------------------------------------------------------
 * Resize event set buffers if needed.
 *
 * \param[in, out]  event_set  pointer to event set
 * \param[in]       minimum required
 *----------------------------------------------------------------------------*/

void
cs_lagr_event_set_resize(cs_lagr_event_set_t  *event_set,
                         cs_lnum_t             min_size)
{
  if (min_size == event_set->n_events_max)
    return;

  assert(min_size >= event_set->n_events);

  event_set->n_events_max = min_size;

  BFT_REALLOC(event_set->e_buffer,
              event_set->n_events_max * event_set->e_am->extents,
              unsigned char);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_lagr_event_set_t structure
 *
 * \param[in]  events  cs_lagr_event_set_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_set_dump(const cs_lagr_event_set_t  *events)
{
  if (events != NULL) {

    bft_printf("Particle events set\n");
    bft_printf("-------------------\n");
    bft_printf("  n_events:      %10ld\n", (long)events->n_events);
    bft_printf("  n_events_max:  %10ld\n", (long)events->n_events_max);

    bft_printf_flush();

    for (cs_lnum_t i = 0; i < events->n_events; i++) {
      _dump_event(events, i);
    }

  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Resize event set buffers if needed.
 *
 * \param[in, out]  events       pointer to event set
 * \param[in, out]  particles    pointer to particle set
 * \param[in]       event_id     event id
 * \param[in]       particle_id  particle id
 *----------------------------------------------------------------------------*/

void
cs_lagr_event_init_from_particle(cs_lagr_event_set_t     *events,
                                 cs_lagr_particle_set_t  *particles,
                                 cs_lnum_t                event_id,
                                 cs_lnum_t                particle_id)
{
  memset(events->e_buffer + events->e_am->extents*event_id,
         0,
         events->e_am->extents);

  for (cs_lnum_t i = 0; i < _n_mapped_part_attr; i++) {
    int attr = _mapped_part_attr[i];

    const unsigned char *p_attr = cs_lagr_particles_attr(particles,
                                                         particle_id,
                                                         attr);

    unsigned char *e_attr = cs_lagr_events_attr(events,
                                                event_id,
                                                attr);

    size_t size = particles->p_am->size[attr];

    for (size_t j = 0; j < size; j++)
      e_attr[j] = p_attr[j];
  }

  cs_lnum_t cell_id = cs_lagr_particles_get_lnum(particles, particle_id,
                                                 CS_LAGR_CELL_ID);
  cs_lagr_events_set_lnum(events,
                          event_id,
                          CS_LAGR_E_CELL_ID,
                          cell_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * Return a cs_lagr_event_set_t structure for particle/boundary interactions.
 *
 * The event set is created if not present yet.
 *
 * This event set is automatically freed and destroyed at the end of the
 * computation.
 *
 * \return pointer to event set
 */
/*----------------------------------------------------------------------------*/

cs_lagr_event_set_t  *
cs_lagr_event_set_boundary_interaction(void)
{
  if (_boundary_events == NULL)
    _boundary_events = _create_event_set(256, _e_attr_map);

  return _boundary_events;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
