/*============================================================================
 * Extract information from lagrangian particles.
 *============================================================================*/

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_lagr_tracking.h"
#include "cs_selector.h"
#include "cs_timer.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_extract.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_extract.c

  \brief Extract information from lagrangian particles.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the local number of particles.
 *
 * \return  current number of particles.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_get_n_particles(void)
{
  cs_lnum_t retval = 0;

  const cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  if (p_set != NULL)
    retval = p_set->n_particles;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract a list of particles using an optional cell filter and
 *        statistical density filter.
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * \param[in]   n_cells          number of cells in filter
 * \param[in]   cell_list        optional list of containing cells filter
 *                               (1 to n numbering)
 * \param[in]   density          if < 1, fraction of particles to select
 * \param[out]  n_particles      number of selected particles, or NULL
 * \param[out]  particle_list    particle_list (1 to n numbering), or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_particle_list(cs_lnum_t         n_cells,
                          const cs_lnum_t   cell_list[],
                          double            density,
                          cs_lnum_t        *n_particles,
                          cs_lnum_t        *particle_list)
{
  cs_lnum_t i;

  ptrdiff_t  displ = 0;

  cs_lnum_t p_count = 0;

  bool *cell_flag = NULL;

  const cs_mesh_t *mesh = cs_glob_mesh;

  const cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  assert(p_set != NULL);

  size_t extents = p_set->p_am->extents;

  if (density < 1) {

    size_t  _extents, size;
    cs_datatype_t  datatype;
    int  count;

    cs_lagr_get_attr_info(p_set,
                          0,
                          CS_LAGR_RANDOM_VALUE,
                          &_extents, &size, &displ,
                          &datatype, &count);

    assert(   (displ > 0 && count == 1 && datatype == CS_REAL_TYPE)
           || (displ < 0 && count == 0 && datatype == CS_DATATYPE_NULL));

  }

  /* Case where we have a filter */

  if (n_cells < mesh->n_cells) {

    /* Build cel filter */

    BFT_MALLOC(cell_flag, mesh->n_cells, bool);

    for (i = 0; i < mesh->n_cells; i++)
      cell_flag[i] = false;

    if (cell_list != NULL) {
      for (i = 0; i < n_cells; i++)
        cell_flag[cell_list[i] - 1] = true;
    }
    else {
      for (i = 0; i < n_cells; i++)
        cell_flag[i] = true;
    }

  }

  /* Loop on particles */

  for (i = 0; i < p_set->n_particles; i++) {

    /* If density < 1, randomly select which particles are added;
       normally, particles maintain a random_value for this purpose,
       but we also plan for the case where this could be optional. */

    if (density < 1) {
      double r;
      if (displ < 0)
        r = (double)rand() / RAND_MAX;
      else {
        const unsigned char
          *p = (const unsigned char *)(p_set->p_buffer + i*extents) + displ;
        r = *(const cs_real_t *)p;
      }
      if (r > density)
        continue;
    }

    /* Check for filter cell */

    if (cell_flag != NULL) {
      cs_lnum_t cur_cell_num
        = cs_lagr_particles_get_lnum(p_set, i, CS_LAGR_CELL_NUM);
      cs_lnum_t  cell_id = CS_ABS(cur_cell_num) - 1;
      if (cell_flag[cell_id] == false)
        continue;
    }

    if (particle_list != NULL)
      particle_list[p_count] = i+1;

    p_count += 1;

  }

  if (cell_flag != NULL)
    BFT_FREE(cell_flag);

  *n_particles = p_count;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract values for a set of particles.
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * \param[in]   particles             associated particle set
 * \param[in]   attr                  attribute whose values are required
 * \param[in]   datatype              associated value type
 * \param[in]   stride                number of values per particle
 * \param[in]   component_id          if -1 : extract the whole attribute
 *                                    if >0 : id of the component to extract
 * \param[in]   n_particles           number of particles in filter
 * \param[in]   particle_list         particle_list (1 to n numbering), or NULL
 * \param[out]  values                particle values for given attribute
 *
 * \return 0 in case of success, 1 if attribute is not present
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_get_particle_values(const cs_lagr_particle_set_t  *particles,
                            cs_lagr_attribute_t            attr,
                            cs_datatype_t                  datatype,
                            int                            stride,
                            int                            component_id,
                            cs_lnum_t                      n_particles,
                            const cs_lnum_t                particle_list[],
                            void                          *values)
{
  size_t j;
  cs_lnum_t i;

  size_t  extents, size, _length;
  ptrdiff_t  displ;
  cs_datatype_t _datatype;
  int  _count;
  unsigned char *_values = values;

  assert(particles != NULL);

  cs_lagr_get_attr_info(particles, 0, attr,
                        &extents, &size, &displ, &_datatype, &_count);

  if (_count == 0)
    return 1;
  else {
    if (component_id == -1)
      _length = size;
    else
      _length = size/_count;
  }

  /* Check consistency */

  if (cs_lagr_check_attr_query(particles,
                               attr,
                               datatype,
                               stride,
                               component_id) != 0)
    return 1;

  /* No offset if export of the whole attribute */
  if (component_id == -1)
    component_id = 0;

  /* Case where we have no filter */

  if (particle_list == NULL) {
    for (i = 0; i < n_particles; i++) {
      unsigned char *dest = _values + i*_length;
      const unsigned char
        *src = (const unsigned char *)(particles->p_buffer + i*extents)
               + displ
               + component_id * _length;
      for (j = 0; j < _length; j++)
        dest[j] = src[j];
    }
  }

  /* Case where we have a filter list */
  else {
    for (i = 0; i < n_particles; i++) {
      cs_lnum_t p_id = particle_list[i] - 1;
      unsigned char *dest = _values + i*_length;
      const unsigned char
        *src = (const unsigned char *)(particles->p_buffer + p_id*extents)
               + displ
               + component_id * _length;
      for (j = 0; j < _length; j++)
        dest[j] = src[j];
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract trajectory values for a set of particles.
 *
 * Trajectories are defined as a mesh of segments, whose start and end
 * points are copied in an interleaved manner in the segment_values array
 * (p1_old, p1_new, p2_old, p2_new, ... pn_old, pn_new).
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * \param[in]   particles        associated particle set
 * \param[in]   attr             attribute whose values are required
 * \param[in]   datatype         associated value type
 * \param[in]   stride           number of values per particle
 * \param[in]   component_id     if -1 : extract the whole attribute
 *                               if >0 : id of the component to extract
 * \param[in]   n_particles      number of particles in filter
 * \param[in]   particle_list    particle_list (1 to n numbering), or NULL
 * \param[out]  segment_values   particle segment values
 *
 * \return 0 in case of success, 1 if attribute is not present
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_get_trajectory_values(const cs_lagr_particle_set_t  *particles,
                              cs_lagr_attribute_t            attr,
                              cs_datatype_t                  datatype,
                              int                            stride,
                              int                            component_id,
                              cs_lnum_t                      n_particles,
                              const cs_lnum_t                particle_list[],
                              void                          *segment_values)
{
  size_t j;
  cs_lnum_t i;

  size_t  extents, size, _length;
  ptrdiff_t  displ, displ_p;
  cs_datatype_t _datatype;
  int  _count;
  unsigned char *_values = segment_values;

  assert(particles != NULL);

  const unsigned char *p_buffer = particles->p_buffer;

  cs_lagr_get_attr_info(particles, 0, attr,
                        &extents, &size, &displ, &_datatype, &_count);

  if (_count == 0)
    return 1;
  else {
    if (component_id == -1)
      _length = size;
    else
      _length = size/_count;
  }

  if (particles->p_am->count[1][attr] > 0)
    cs_lagr_get_attr_info(particles, 1, attr,
                          &extents, NULL, &displ_p, NULL, NULL);

  /* Check consistency */

  if (cs_lagr_check_attr_query(particles,
                               attr,
                               datatype,
                               stride,
                               component_id) != 0)
    return 1;

  /* No offset if export of the whole attribute */
  if (component_id == -1)
    component_id = 0;

  /* Case where we have no filter */

  if (particle_list == NULL) {

    if (particles->p_am->count[1][attr] > 0) {

      for (i = 0; i < n_particles; i++) {
        unsigned char *dest = _values + i*_length*2;
        const unsigned char
          *src = (const unsigned char *)(p_buffer + i*extents)
                  + displ
                  + component_id * _length;
        const unsigned char
          *srcp = (const unsigned char *)(p_buffer + i*extents)
                   + displ_p
                   + component_id * _length;
        for (j = 0; j < _length; j++) {
          dest[j] = src[j];
          dest[j + _length] = srcp[j];
        }

      }

    }
    else { /* With no previous value available; copy current value */

      for (i = 0; i < n_particles; i++) {
        unsigned char *dest = _values + i*_length*2;
        const unsigned char
          *src = (const unsigned char *)(p_buffer + i*extents)
                  + displ
                  + component_id * _length;
        for (j = 0; j < _length; j++) {
          dest[j] = src[j];
          dest[j + _length] = src[j];
        }

      }

    }
  }

  /* Case where we have a filter list */
  else {

    if (particles->p_am->count[1][attr] > 0) {

      for (i = 0; i < n_particles; i++) {
        cs_lnum_t p_id = particle_list[i] - 1;
        unsigned char *dest = _values + i*_length*2;
        const unsigned char
          *src = (const unsigned char *)(p_buffer + p_id*extents)
                  + displ
                  + component_id * _length;
        const unsigned char
          *srcp = (const unsigned char *)(p_buffer + p_id*extents)
                   + displ_p
                   + component_id * _length;
        for (j = 0; j < _length; j++) {
          dest[j] = src[j];
          dest[j + _length] = srcp[j];
        }
      }

    }

    else { /* With no previous value available; copy current value */

      for (i = 0; i < n_particles; i++) {
        cs_lnum_t p_id = particle_list[i] - 1;
        unsigned char *dest = _values + i*_length*2;
        const unsigned char
          *src = (const unsigned char *)(p_buffer + p_id*extents)
                  + displ
                  + component_id * _length;
        for (j = 0; j < _length; j++) {
          dest[j] = src[j];
          dest[j + _length] = src[j];
        }
      }

    }

  }

  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
