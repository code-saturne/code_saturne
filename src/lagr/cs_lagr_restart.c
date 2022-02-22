/*============================================================================
 * Checkpoint/restart handling for Lagrangian module.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_lagr.h"
#include "cs_lagr_extract.h"
#include "cs_lagr_tracking.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_mesh_location.h"
#include "cs_restart.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_restart.c
        Checkpoint/restart handling for Lagrangian module.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build particle values section name
 *
 * parameters:
 *   attr     <-- attribute
 *   comp_id  <-- component id, or -1 for all
 *   sec_name --> associated section name
 *----------------------------------------------------------------------------*/

static void
_lagr_section_name(cs_lagr_attribute_t  attr,
                   int                  comp_id,
                   char                 sec_name[128])
{
  if (comp_id < 0)
    snprintf(sec_name, 127, "particle_%s::vals::0",
             cs_lagr_attribute_name[attr]);
  else
    snprintf(sec_name, 127, "particle_%s::vals::%d::0",
             cs_lagr_attribute_name[attr], comp_id);
}

/*----------------------------------------------------------------------------
 * Build legacy particle values section name or formats for names
 *
 * parameters:
 *   attr     <-- attribute
 *   sec_name --> associated section name, component 0
 *----------------------------------------------------------------------------*/

static void
_legacy_section_name(cs_lagr_attribute_t  attr,
                     char                 sec_name[128])
{
  sec_name[0] = '\0';

  switch(attr) {
  case CS_LAGR_RANDOM_VALUE:
    strcpy(sec_name, "random_value");
    break;
  default:
    break;
  }
}

/*----------------------------------------------------------------------------
 * Assign values for a set of particles.
 *
 * parameters:
 *   particles    <-> associated particle set
 *   attr         <-- attribute whose values are set
 *   datatype     <-- associated value type
 *   stride       <-- number of values per particle
 *   component_id <-- if -1 : set the whole attribute
 *                    if >0 : id of the component to set
 *   values                particle values for given attribute
 *
 * returns:
 *   0 in case of success, 1 if attribute is not present
 *----------------------------------------------------------------------------*/

static int
_set_particle_values(cs_lagr_particle_set_t  *particles,
                     cs_lagr_attribute_t      attr,
                     cs_datatype_t            datatype,
                     int                      stride,
                     int                      component_id,
                     void                    *values)
{
  size_t j;
  cs_lnum_t i;

  size_t  extents, size, _length;
  ptrdiff_t  displ;
  cs_datatype_t _datatype;
  int  _count;
  unsigned char *_values = values;

  cs_lnum_t n_particles = particles->n_particles;

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

  for (i = 0; i < n_particles; i++) {
    const unsigned char *src = _values + i*_length;
    unsigned char
      *dest = (unsigned char *)(particles->p_buffer + i*extents)
               + displ
               + component_id * _length;
    for (j = 0; j < _length; j++)
      dest[j] = src[j];
  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Initialize values for a set of particles.
 *
 * parameters:
 *   particles    <-> associated particle set
 *   attr         <-- attribute whose values are set
 *   component_id <-- if -1 : set the whole attribute
 *                    if >0 : id of the component to set
 *
 * returns:
 *   0 in case of success, 1 if attribute is not present
 *----------------------------------------------------------------------------*/

static int
_init_particle_values(cs_lagr_particle_set_t  *particles,
                      cs_lagr_attribute_t      attr,
                      int                      component_id)
{
  size_t  extents, size;
  ptrdiff_t  displ;
  cs_datatype_t datatype;
  int  stride;

  cs_lnum_t n_particles = particles->n_particles;

  assert(particles != NULL);

  cs_lagr_get_attr_info(particles, 0, attr,
                        &extents, &size, &displ, &datatype, &stride);

  if (stride == 0)
    return 1;

  /* Check component_id */
  if (component_id < -1 || component_id >= stride) {
    char attr_name[128] = "";
    attr_name[127] = '\0';
    const char *_attr_name = attr_name;
    if (attr < CS_LAGR_N_ATTRIBUTES) {
      snprintf(attr_name, 127, "CS_LAGR_%s", cs_lagr_attribute_name[attr]);
      size_t l = strlen(attr_name);
      for (size_t i = 0; i < l; i++)
        attr_name[i] = toupper(attr_name[i]);
    }
    else {
      snprintf(attr_name, 127, "%d", (int)attr);
    }
    bft_error(__FILE__, __LINE__, 0,
              _("Attribute %s has a number of components equal to %d\n"
                "but component %d is requested."),
              _attr_name,
              stride,
              component_id);
    return 1;
  }

  /* No offset if export of the whole attribute */
  if (component_id == -1)
    component_id = 0;

  switch(attr) {

  case CS_LAGR_RANDOM_VALUE:
    {
      assert(datatype == CS_REAL_TYPE);
      assert(stride == 1);
      for (cs_lnum_t i = 0; i < n_particles; i++)
        cs_lagr_particles_set_real(particles, i, attr, rand()/RAND_MAX);
    }
    break;

  case CS_LAGR_TEMPERATURE:
  case CS_LAGR_FLUID_TEMPERATURE:
    {
      cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
      double c_kelvin = 0;
      const double *t = NULL;

      if (extra->temperature != NULL)
        t = extra->temperature->val;
      if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_KELVIN)
        c_kelvin = 273.15;

      assert(datatype == CS_REAL_TYPE);

      /* initialize particle temperature with fluid temperature */
      if (t != NULL) {
        for (cs_lnum_t i = 0; i < n_particles; i++) {
          cs_lnum_t cell_id
            = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_ID);
          cs_real_t *part_val
            = cs_lagr_particles_attr_n(particles, i, 0, attr);
          for (cs_lnum_t j = 0; j < stride; j++)
            part_val[j] = t[cell_id] - c_kelvin;
        }
      }

      /* set to reference temperature */
      else {
        for (cs_lnum_t i = 0; i < n_particles; i++) {
          cs_real_t *part_val
            = cs_lagr_particles_attr_n(particles, i, 0, attr);
          for (cs_lnum_t j = 0; j < stride; j++)
            part_val[j] = cs_glob_fluid_properties->t0 - c_kelvin;
        }
      }
    }
    break;

  default:

    if (datatype == CS_LNUM_TYPE) {
      assert(stride == 1);
      for (cs_lnum_t i = 0; i < n_particles; i++) {
        cs_lnum_t *dest = cs_lagr_particles_attr_n(particles, i, 0, attr);
        for (cs_lnum_t j = 0; j < stride; j++)
          dest[j] = 0;
      }
    }
    else if (datatype == CS_GNUM_TYPE) {
      assert(stride == 1);
      for (cs_lnum_t i = 0; i < n_particles; i++) {
        cs_gnum_t *dest = cs_lagr_particles_attr_n(particles, i, 0, attr);
        for (cs_lnum_t j = 0; j < stride; j++)
          dest[j] = 0;
      }
    }
    else if (datatype == CS_REAL_TYPE) {
      for (cs_lnum_t i = 0; i < n_particles; i++) {
        cs_real_t *dest = cs_lagr_particles_attr_n(particles, i, 0, attr);
        for (cs_lnum_t j = 0; j < stride; j++)
          dest[j] = 0.;
      }
    }

  }

  return 0;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read particle data from checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 *
 * \return  number of particle arrays read
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_restart_read_particle_data(cs_restart_t  *r)
{
  int retval = 0;

  size_t  extents, size;
  ptrdiff_t  displ;
  cs_datatype_t datatype;
  int  stride, sec_code;
  cs_restart_val_type_t restart_type;
  unsigned char *vals = NULL;
  size_t max_size = 0;

  char sec_name[128], old_name[128];

  /* Initialization */

  cs_lnum_t n_particles;
  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();

  if (p_set == NULL)
    return retval;

  p_set->n_particles = 0;
  p_set->weight = 0;
  p_set->n_part_out = 0;
  p_set->weight_out = 0;
  p_set->n_part_dep = 0;
  p_set->weight_dep = 0;
  p_set->n_part_fou = 0;
  p_set->weight_fou = 0;
  p_set->n_failed_part = 0;
  p_set->weight_failed = 0;

  /* Read particles location and info */
  /*----------------------------------*/

  _lagr_section_name(CS_LAGR_COORDS, -1, sec_name);

  const int particles_location_id
    = cs_restart_read_particles_info(r, sec_name, &n_particles);

  if (particles_location_id < 0)
    return retval;

  /* Read coordinates and get mesh location */
  /*-----------------------------------------*/

  cs_lnum_t  *p_cell_id;
  cs_real_t  *p_coords;

  BFT_MALLOC(p_cell_id, n_particles, cs_lnum_t);
  BFT_MALLOC(p_coords, n_particles*3, cs_real_t);

  sec_code = cs_restart_read_particles(r,
                                       particles_location_id,
                                       p_cell_id,
                                       p_coords);

  if (sec_code == CS_RESTART_SUCCESS) {

    p_set->n_particles = n_particles;

    if (p_set->n_particles_max < p_set->n_particles)
      cs_lagr_particle_set_resize(p_set->n_particles);

    _set_particle_values(p_set, CS_LAGR_COORDS, CS_REAL_TYPE,
                         3, -1, p_coords);

    _set_particle_values(p_set, CS_LAGR_CELL_ID, CS_LNUM_TYPE,
                         1, -1, p_cell_id);

  }

  BFT_FREE(p_cell_id);
  BFT_FREE(p_coords);

  if (sec_code == CS_RESTART_SUCCESS)
    retval = 1;

  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Error reading main particle restart data."));
    return retval;
  }

  /* Loop on all other attributes, handling special cases */
  /*------------------------------------------------------*/

  for (cs_lagr_attribute_t attr = 0; attr < CS_LAGR_N_ATTRIBUTES; attr++) {

    cs_lagr_get_attr_info(p_set, 0, attr,
                          &extents, &size, &displ, &datatype, &stride);

    if (stride == 0)
      continue;

    if (datatype == CS_LNUM_TYPE)
      restart_type = CS_TYPE_int;
    else if (datatype == CS_GNUM_TYPE)
      restart_type = CS_TYPE_cs_gnum_t;
    else {
      assert(datatype == CS_REAL_TYPE);
      restart_type = CS_TYPE_cs_real_t;
    }

    switch(attr) {

    case CS_LAGR_COORDS:
    case CS_LAGR_CELL_ID:
    case CS_LAGR_RANK_ID:
      break;

    case CS_LAGR_NEIGHBOR_FACE_ID:
      {
        cs_lnum_t *face_id = (cs_lnum_t *)vals;

        /* Initialize to default */

        for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
          face_id[i] = -1;

        /* Now read global numbers (1 to n, 0 for none) */

        cs_gnum_t *gnum_read = NULL;
        BFT_MALLOC(gnum_read, p_set->n_particles, cs_gnum_t);

        snprintf(sec_name, 127, "particle_%s::vals::0", "neighbor_face_num");
        _legacy_section_name(attr, old_name);

        sec_code = cs_restart_read_section_compat(r,
                                                  sec_name,
                                                  old_name,
                                                  particles_location_id,
                                                  stride,
                                                  CS_TYPE_cs_gnum_t,
                                                  gnum_read);

        if (sec_code == CS_RESTART_SUCCESS)
          retval += 1;
        else {
          for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
            gnum_read[i] = 0;
        }

        /* Now assign values */

        const cs_lnum_t   *cell_b_faces_idx = NULL;
        const cs_lnum_t   *cell_b_faces = NULL;
        const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

        if (ma != NULL) {
          cell_b_faces_idx = ma->cell_b_faces_idx;
          cell_b_faces = ma->cell_b_faces;
        }

        if (cell_b_faces_idx != NULL) {
          const cs_gnum_t *g_b_face_num = cs_glob_mesh->global_b_face_num;
          if (g_b_face_num != NULL) {
            for (cs_lnum_t i = 0; i < p_set->n_particles; i++) {
              if (gnum_read[i] > 0) {
                cs_lnum_t cell_id = cs_lagr_particles_get_lnum(p_set, i,
                                                               CS_LAGR_CELL_ID);
                cs_lnum_t s_id = cell_b_faces_idx[cell_id];
                cs_lnum_t e_id = cell_b_faces_idx[cell_id+1];
                for (cs_lnum_t j = s_id; j < e_id; j++) {
                  if (g_b_face_num[cell_b_faces[j]] == gnum_read[i]) {
                    face_id[i] = cell_b_faces[j];
                    break;
                  }
                }
              }
            }
          }
          else {
            for (cs_lnum_t i = 0; i < p_set->n_particles; i++) {
              if (gnum_read[i] > 0) {
                cs_lnum_t cell_id = cs_lagr_particles_get_lnum(p_set, i,
                                                               CS_LAGR_CELL_ID);
                cs_lnum_t s_id = cell_b_faces_idx[cell_id];
                cs_lnum_t e_id = cell_b_faces_idx[cell_id+1];
                for (cs_lnum_t j = s_id; j < e_id; j++) {
                  if ((cs_gnum_t)(j+1) == gnum_read[i]) {
                    face_id[i] = cell_b_faces[j];
                    break;
                  }
                }
              }
            }
          }
        }
        else {
          assert(cell_b_faces != NULL);
        }

        BFT_FREE(gnum_read);

        _set_particle_values(p_set,
                             attr,
                             CS_LNUM_TYPE,
                             1,
                             -1,
                             vals);
      }
      break;

    default:

      if (size > max_size) {
        max_size = size;
        BFT_REALLOC(vals, max_size*n_particles, unsigned char);
      }

      int n_sections = stride;

      if (   attr == CS_LAGR_VELOCITY
          || attr == CS_LAGR_VELOCITY_SEEN)
        n_sections = 1;

      for (int s_id = 0; s_id < n_sections; s_id++) {

        int comp_id = (n_sections == 1) ? -1 : s_id;
        int c_stride = (n_sections == 1) ? stride : 1;

        _lagr_section_name(attr, comp_id, sec_name);
        _legacy_section_name(attr, old_name);

        /* Try to read data */

        sec_code = cs_restart_read_section_compat(r,
                                                  sec_name,
                                                  old_name,
                                                  particles_location_id,
                                                  c_stride,
                                                  restart_type,
                                                  vals);

        /* Finish setting values */

        if (sec_code == CS_RESTART_SUCCESS) {

          if (attr == CS_LAGR_STAT_WEIGHT) {
            cs_real_t *w = (cs_real_t *)vals;
            for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
              p_set->weight += w[i];
          }

          _set_particle_values(p_set,
                               attr,
                               datatype,
                               stride,
                               -1,
                               vals);

          retval += 1;

        }
        else
          _init_particle_values(p_set, attr, -1);

      }

    }

  }

  BFT_FREE(vals);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write particle data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 *
 * \return  number of particle arrays written
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_restart_write_particle_data(cs_restart_t  *r)
{
  int retval = 0;

  size_t  extents, size;
  ptrdiff_t  displ;
  cs_datatype_t datatype;
  int  stride;
  cs_restart_val_type_t restart_type;
  unsigned char *vals = NULL;
  size_t max_size = 0;

  char sec_name[128];

  /* Initialization */

  const cs_lnum_t n_particles = cs_lagr_get_n_particles();
  const cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();

  if (p_set == NULL)
    return retval;

  /* Write coordinates and get mesh location */
  /*-----------------------------------------*/

  _lagr_section_name(CS_LAGR_COORDS, -1, sec_name);

  cs_lnum_t  *p_cell_id;
  cs_real_t  *p_coords;

  BFT_MALLOC(p_cell_id, n_particles, cs_lnum_t);
  BFT_MALLOC(p_coords, n_particles*3, cs_real_t);

  cs_lagr_get_particle_values(p_set, CS_LAGR_COORDS, CS_REAL_TYPE,
                              3, -1, n_particles, NULL, p_coords);

  cs_lagr_get_particle_values(p_set, CS_LAGR_CELL_ID, CS_LNUM_TYPE,
                              1, -1, n_particles, NULL, p_cell_id);

  const int particles_location_id
    = cs_restart_write_particles(r,
                                 sec_name,
                                 false, /* number_by_coords */
                                 n_particles,
                                 p_cell_id,
                                 p_coords);

  BFT_FREE(p_cell_id);
  BFT_FREE(p_coords);

  retval = 1;

  /* Loop on all other attributes, handling special cases */
  /*------------------------------------------------------*/

  for (cs_lagr_attribute_t attr = 0; attr < CS_LAGR_N_ATTRIBUTES; attr++) {

    cs_lagr_get_attr_info(p_set, 0, attr,
                          &extents, &size, &displ, &datatype, &stride);

    if (stride == 0)
      continue;

    if (datatype == CS_LNUM_TYPE)
      restart_type = CS_TYPE_int;
    else if (datatype == CS_GNUM_TYPE)
      restart_type = CS_TYPE_cs_gnum_t;
    else {
      assert(datatype == CS_REAL_TYPE);
      restart_type = CS_TYPE_cs_real_t;
    }

    switch(attr) {

    case CS_LAGR_COORDS:
    case CS_LAGR_CELL_ID:
    case CS_LAGR_RANK_ID:
      break;

    case CS_LAGR_NEIGHBOR_FACE_ID:
      {
        cs_lagr_get_particle_values(p_set,
                                    attr,
                                    CS_LNUM_TYPE,
                                    1,
                                    -1,
                                    n_particles,
                                    NULL,
                                    vals);

        cs_gnum_t *gnum_write = NULL;
        BFT_MALLOC(gnum_write, p_set->n_particles, cs_gnum_t);

        const cs_gnum_t *g_b_face_num = cs_glob_mesh->global_b_face_num;
        if (g_b_face_num != NULL) {
          for (cs_lnum_t i = 0; i < p_set->n_particles; i++) {
            cs_lnum_t face_id = ((cs_lnum_t *)vals)[i];
            if (face_id < 0)
              gnum_write[i] = 0;
            else
              gnum_write[i] = g_b_face_num[face_id];
          }
        }
        else {
          for (cs_lnum_t i = 0; i < p_set->n_particles; i++) {
            cs_lnum_t face_id = ((cs_lnum_t *)vals)[i];
            if (face_id < 0)
              gnum_write[i] = 0;
            else
              gnum_write[i] = face_id + 1;
          }
        }

        snprintf(sec_name, 127, "particle_%s::vals::0", "neighbor_face_num");

        cs_restart_write_section(r,
                                 sec_name,
                                 particles_location_id,
                                 1,
                                 CS_TYPE_cs_gnum_t,
                                 gnum_write);

        BFT_FREE(gnum_write);

        retval += 1;
      }
      break;

    default:

      if (size > max_size) {
        max_size = size;
        BFT_REALLOC(vals, max_size*n_particles, unsigned char);
      }

      int n_sections = stride;

      if (   attr == CS_LAGR_VELOCITY
          || attr == CS_LAGR_VELOCITY_SEEN)
        n_sections = 1;

      for (int s_id = 0; s_id < n_sections; s_id++) {

        int comp_id = (n_sections == 1) ? -1 : s_id;
        int c_stride = (n_sections == 1) ? stride : 1;

        cs_lagr_get_particle_values(p_set,
                                    attr,
                                    datatype,
                                    stride,
                                    comp_id,
                                    n_particles,
                                    NULL,
                                    vals);

        _lagr_section_name(attr, comp_id, sec_name);

        cs_restart_write_section(r,
                                 sec_name,
                                 particles_location_id,
                                 c_stride,
                                 restart_type,
                                 vals);

        retval += 1;

      }

    }

  }

  BFT_FREE(vals);

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
