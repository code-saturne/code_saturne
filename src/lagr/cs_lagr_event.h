#ifndef __CS_LAGR_EVENT_H__
#define __CS_LAGR_EVENT_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_lagr_particle.h"

#include "assert.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * Flags specifying general event attributes
 */

/*! inflow event */
#define CS_EVENT_INFLOW                (1 << 0)

/*! outflow event */
#define CS_EVENT_OUTFLOW               (1 << 1)

/*! rebound on wall */
#define CS_EVENT_REBOUND               (1 << 2)

/*! deposistion event */
#define CS_EVENT_DEPOSITION            (1 << 3)

/*! resuspension event */
#define CS_EVENT_RESUSPENSION          (1 << 4)

/*! roll off */
#define CS_EVENT_ROLL_OFF              (1 << 5)

/*! roll on */
#define CS_EVENT_ROLL_ON               (1 << 6)

/*! fouling  */
#define CS_EVENT_FOULING               (1 << 7)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Prefedined event attributes */
/* ----------------------------- */

/* Numbering of predefined attributes starts after particle attributes,
   so that particle attribute can be added to mapped variables traced
   a particle events */

typedef enum {

  CS_LAGR_E_FLAG    = CS_LAGR_N_ATTRIBUTES,   /*!< local flag */
  CS_LAGR_E_CELL_ID,                          /*!< local cell id */
  CS_LAGR_E_FACE_ID,                          /*!< local face id
                                                (-1 if not on boundary) */

  CS_LAGR_E_VELOCITY,                         /*!< velocity following event */

  /* End of event attributes */

  CS_LAGR_N_E_ATTRIBUTES

} cs_lagr_event_attribute_t;

/*! Event attribute structure mapping */
/* ------------------------------------- */

typedef struct {

  size_t         extents;                         /* size (in bytes) of event
                                                     structure */
  size_t         lb;                              /* size (in bytes) of lower
                                                     bounds of event data
                                                     (work area before) */

  size_t         size[CS_LAGR_N_E_ATTRIBUTES];    /* size (in bytes) of
                                                     attributes in event
                                                     structure for a given
                                                     time value */
  cs_datatype_t  datatype[CS_LAGR_N_E_ATTRIBUTES]; /* datatype of associated
                                                      attributes */
  int            count[CS_LAGR_N_E_ATTRIBUTES];    /* number of values for each
                                                      attribute */
  ptrdiff_t      displ[CS_LAGR_N_E_ATTRIBUTES];    /* displacement (in bytes)
                                                      of attributes in event
                                                      data */

} cs_lagr_event_attribute_map_t;

/* Event set */
/* ------------ */

typedef struct {

  cs_lnum_t  n_events;                              /* number of events */
  cs_lnum_t  n_events_max;

  const cs_lagr_event_attribute_map_t  *e_am;       /*!< event attributes map */
  unsigned char                        *e_buffer;   /*!< Events data buffer */

} cs_lagr_event_set_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define event map based on defined options.
 *
 * This function should only be called after
 * \ref cs_lagr_particle_attr_initialize,
 * as it may use elements from the main particle attributes map.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy event set map if it exists.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return const pointer to the main event attribute map structure.
 *
 * \return pointer to current event attribute map, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lagr_event_attribute_map_t *
cs_lagr_event_get_attr_map(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return name associated with a given attribute.
 *
 * \param[in]   attr   event attribute
 */
/*----------------------------------------------------------------------------*/

const char *
cs_lagr_event_get_attr_name(cs_lagr_event_attribute_t   attr);

/*----------------------------------------------------------------------------*/
/*!
 * Create a cs_lagr_event_set_t structure.
 *
 * \return pointer to event set
 */
/*----------------------------------------------------------------------------*/

cs_lagr_event_set_t  *
cs_lagr_event_set_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * Destroy a cs_lagr_event_set_t structure.
 *
 * \param[in, out]  events  pointer to pointer to event set to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_set_destroy(cs_lagr_event_set_t  **events);

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
                            int                        *count);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if an event attribute is in a valid range.
 *
 * If this is not the case, a fatal error is provoked.

 * \param[in]   attr       event attribute
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_attr_in_range(int  attr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to a current attribute of a given event in a set.
 *
 * \param[in]  event_set  pointer to event set
 * \param[in]  event_id   event id
 * \param[in]  attr       requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static void *
cs_lagr_events_attr(cs_lagr_event_set_t  *event_set,
                    cs_lnum_t             event_id,
                    int                   attr)
{
  assert(event_set->e_am->count[attr] > 0);

  return   (unsigned char *)event_set->e_buffer
         + event_set->e_am->extents*event_id
         + event_set->e_am->displ[attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to current attribute data
 *        of a given event in a set.
 *
 * \param[in]  event_set  pointer to event set
 * \param[in]  event_id   event id
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static const void *
cs_lagr_events_attr_const(const cs_lagr_event_set_t  *event_set,
                          cs_lnum_t                   event_id,
                          int                         attr)
{
  assert(event_set->e_am->count[attr] > 0);

  return   event_set->e_buffer
         + event_set->e_am->extents*event_id
         + event_set->e_am->displ[attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_lnum_t of a given event in a set.
 *
 * \param[in]  event_set  pointer to event set
 * \param[in]  event_id   event id
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_lagr_events_get_lnum(const cs_lagr_event_set_t  *event_set,
                        cs_lnum_t                   event_id,
                        int                         attr)
{
  assert(event_set->e_am->count[attr] > 0);

  return *((const cs_lnum_t *)(  event_set->e_buffer
                               + event_set->e_am->extents*event_id
                               + event_set->e_am->displ[attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_lnum_t of a given event in a set.
 *
 * \param[in, out]   event_set  pointer to event set
 * \param[in]        event_id   event id
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_events_set_lnum(cs_lagr_event_set_t        *event_set,
                        cs_lnum_t                   event_id,
                        int                         attr,
                        cs_lnum_t                   value)
{
  assert(event_set->e_am->count[attr] > 0);

  *((cs_lnum_t *)(  event_set->e_buffer
                  + event_set->e_am->extents*event_id
                  + event_set->e_am->displ[attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_real_t of a given event in a set.
 *
 * \param[in]  event_set  pointer to event set
 * \param[in]  event_id   event id
 * \param[in]  attr       requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_lagr_events_get_real(const cs_lagr_event_set_t  *event_set,
                        cs_lnum_t                   event_id,
                        int                         attr)
{
  assert(event_set->e_am->count[attr] > 0);

  return *((const cs_real_t *)(  event_set->e_buffer
                               + event_set->e_am->extents*event_id
                               + event_set->e_am->displ[attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_real_t of a given event in a set.
 *
 * \param[in, out]   event_set  pointer to event set
 * \param[in]        event_id   event id
 * \param[in]        attr       requested attribute id
 * \param[in]        value      value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_events_set_real(cs_lagr_event_set_t  *event_set,
                        cs_lnum_t             event_id,
                        int                   attr,
                        cs_real_t             value)
{
  assert(event_set->e_am->count[attr] > 0);

  *((cs_real_t *)(  event_set->e_buffer
                  + event_set->e_am->extents*event_id
                  + event_set->e_am->displ[attr])) = value;
}

/*----------------------------------------------------------------------------
 * Resize event set buffers if needed.
 *
 * \param[in, out]  event_set  pointer to event set
 * \param[in]       mini mum required
 *----------------------------------------------------------------------------*/

void
cs_lagr_event_set_resize(cs_lagr_event_set_t  *event_set,
                         cs_lnum_t             min_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_lagr_event_set_t structure
 *
 * \param[in]  events  cs_lagr_event_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_event_set_dump(const cs_lagr_event_set_t  *events);

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
                                 cs_lnum_t                particle_id);

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
cs_lagr_event_set_boundary_interaction(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_EVENT_H__ */
