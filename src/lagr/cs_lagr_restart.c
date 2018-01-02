/*============================================================================
 * Checkpoint/restart handling for Lagrangian module.
 *============================================================================*/

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
#include "cs_lagr_extract.h"
#include "cs_lagr_tracking.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
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
 * Build particle values section name
 *
 * parameters:
 *   base     <-- base name
 *   comp_id  <-- component id, or -1 for all
 *   sec_name --> associated section name
 *----------------------------------------------------------------------------*/

static void
_section_name(const char  *base,
              int          comp_id,
              char         sec_name[128])
{
  size_t l;

  if (comp_id < 0)
    l= snprintf(sec_name, 127, "%s::vals::0", base);
  else
    l= snprintf(sec_name, 127, "%s::vals::%d::0", base, comp_id);

  sec_name[127] = '\0';
  for (size_t i = 0; i < l; i++)
    sec_name[i] = tolower(sec_name[i]);
}

/*----------------------------------------------------------------------------
 * Build legacy particle values section names or formats for names
 *
 * parameters:
 *   attr       <-- attribute
 *   comp_id    <-- component id, or -1 for all
 *   sec_name_0 --> associated section name, component 0
 *   sec_name_0 --> associated section name, component 1
 *   sec_name_0 --> associated section name, component 2
 *----------------------------------------------------------------------------*/

static void
_legacy_section_names(cs_lagr_attribute_t  attr,
                      int                  comp_id,
                      char                 sec_name_0[128],
                      char                 sec_name_1[128],
                      char                 sec_name_2[128])
{
  sec_name_0[0] = '\0';
  sec_name_1[0] = '\0';
  sec_name_2[0] = '\0';

  switch(attr) {
  case CS_LAGR_RANDOM_VALUE:
    strcpy(sec_name_0, "random_value");
    break;
  case CS_LAGR_STAT_WEIGHT:
    strcpy(sec_name_0, "poids_statistiques_particules");
    break;
  case CS_LAGR_RESIDENCE_TIME:
    strcpy(sec_name_0, "temps_sejour_particules");
    break;
  case CS_LAGR_MASS:
    strcpy(sec_name_0, "variable_masse_particule");
    break;
  case CS_LAGR_DIAMETER:
    strcpy(sec_name_0, "variable_diametre_particule");
    break;
  case CS_LAGR_VELOCITY:
    strcpy(sec_name_0, "variable_vitesseU_particule");
    strcpy(sec_name_1, "variable_vitesseV_particule");
    strcpy(sec_name_2, "variable_vitesseW_particule");
    break;
  case CS_LAGR_VELOCITY_SEEN:
    strcpy(sec_name_0, "variable_vitesseU_fluide_vu");
    strcpy(sec_name_1, "variable_vitesseV_fluide_vu");
    strcpy(sec_name_2, "variable_vitesseW_fluide_vu");
    break;
  case CS_LAGR_YPLUS:
    strcpy(sec_name_0, "yplus_particules");
    break;
  case CS_LAGR_INTERF:
    strcpy(sec_name_0, "dx_particules");
    break;
  case CS_LAGR_NEIGHBOR_FACE_ID:
    strcpy(sec_name_0, "dfac_particules");
    break;
  case CS_LAGR_MARKO_VALUE:
    strcpy(sec_name_0, "indicateur_de_saut");
    break;
  case CS_LAGR_DEPOSITION_FLAG:
    strcpy(sec_name_0, "part_depo");
    break;
  case CS_LAGR_N_LARGE_ASPERITIES:
    strcpy(sec_name_0, "nb_ls_aspe");
    break;
  case CS_LAGR_N_SMALL_ASPERITIES:
    strcpy(sec_name_0, "nb_sms_aspe");
    break;
  case CS_LAGR_ADHESION_FORCE:
    strcpy(sec_name_0, "force_adhesion");
    break;
  case CS_LAGR_ADHESION_TORQUE:
    strcpy(sec_name_0, "moment_adhesion");
    break;
  case CS_LAGR_DISPLACEMENT_NORM:
    strcpy(sec_name_0, "disp_norm");
    break;
  case CS_LAGR_TEMPERATURE:
    if (comp_id < 0)
      strcpy(sec_name_0, "variable_temperature_particule");
    else
      sprintf(sec_name_0, "variable_temperature_particule_couche_%04d", comp_id+1);
    break;
  case CS_LAGR_FLUID_TEMPERATURE:
    strcpy(sec_name_0, "variable_temperature_fluide_vu");
    break;
  case CS_LAGR_CP:
    strcpy(sec_name_0, "variable_chaleur_specifique_particule");
    break;
  case CS_LAGR_WATER_MASS:
    strcpy(sec_name_0, "variable_masse_humidite");
    break;
  case CS_LAGR_COAL_MASS:
    sprintf(sec_name_0, "variable_masse_charbon_reactif_%04d", comp_id+1);
    break;
  case CS_LAGR_COKE_MASS:
    sprintf(sec_name_0, "variable_masse_coke_couche_%04d", comp_id+1);
    break;
  case CS_LAGR_SHRINKING_DIAMETER:
    strcpy(sec_name_0, "diametre_coeur_retrecissant_charbon");
    break;
  case CS_LAGR_INITIAL_DIAMETER:
    strcpy(sec_name_0, "diametre_initial_charbon");
    break;
  case CS_LAGR_COAL_NUM:
    strcpy(sec_name_0, "numero_charbon");
    break;
  case CS_LAGR_COAL_DENSITY:
    sprintf(sec_name_0, "masse_volumique_coke_couche_%04d", comp_id+1);
    break;
  case CS_LAGR_EMISSIVITY:
    strcpy(sec_name_0, "emissivite_particules");
    break;
  case CS_LAGR_STAT_CLASS:
    strcpy(sec_name_0, "numero_groupe_statistiques");
    break;
  case CS_LAGR_USER:
    sprintf(sec_name_0, "variable_supplementaire_%04d", comp_id+1);
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
      /* Determine temperature or enthalpy */
      double c_kelvin = 0;
      const double *t = NULL;
      double *h = NULL;
      if (t == NULL) {
        const cs_field_t *f = cs_field_by_name_try("t_gas");
        if (f != NULL) {
          t = f->val;
          c_kelvin = 273.15;
        }
      }
      if (t == NULL) {
        const cs_field_t *f = cs_field_by_name_try("t_celcius");
        if (f != NULL)
          t = f->val;
      }
      if (t == NULL) {
        const cs_field_t *f = cs_field_by_name_try("temperature");
        if (f != NULL) {
          t = f->val;
          if (cs_glob_thermal_model->itpscl == 1)
            c_kelvin = 273.15;
        }
      }
      if (t == NULL) {
        const cs_field_t *f = cs_field_by_name_try("enthalpy");
        if (f != NULL)
          h = f->val;
      }

      assert(datatype == CS_REAL_TYPE);

      /* initialize particle temperature with fluid temperature */
      if (t != NULL) {
        for (cs_lnum_t i = 0; i < n_particles; i++) {
          cs_lnum_t cell_num
            = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM);
          cs_lnum_t cell_id = CS_ABS(cell_num) - 1;
          cs_real_t *part_val
            = cs_lagr_particles_attr_n(particles, i, 0, attr);
          for (cs_lnum_t j = 0; j < stride; j++)
            part_val[j] = t[cell_id] - c_kelvin;
        }
      }

      /* initialize particle temperature with enthalpy */
      else if (h != NULL) {
        int mode = 1;
        for (cs_lnum_t i = 0; i < n_particles; i++) {
          cs_lnum_t cell_num
            = cs_lagr_particles_get_lnum(particles, i, CS_LAGR_CELL_NUM);
          cs_lnum_t cell_id = CS_ABS(cell_num) - 1;
          cs_real_t *part_val
            = cs_lagr_particles_attr_n(particles, i, 0, attr);
          for (cs_lnum_t j = 0; j < stride; j++)
            CS_PROCF(usthht, USTHHT)(&mode, h + cell_id, part_val + j);
        }
      }

      /* set to zero */
      else {
        for (cs_lnum_t i = 0; i < n_particles; i++) {
          cs_real_t *part_val
            = cs_lagr_particles_attr_n(particles, i, 0, attr);
          for (cs_lnum_t j = 0; j < stride; j++)
            part_val[j] = 0.0;
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

  char sec_name[128], old_name_0[128], old_name_1[128], old_name_2[128];

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

  cs_lnum_t  *p_cell_num;
  cs_real_t  *p_coords;

  BFT_MALLOC(p_cell_num, n_particles, cs_lnum_t);
  BFT_MALLOC(p_coords, n_particles*3, cs_real_t);

  sec_code = cs_restart_read_particles(r,
                                       particles_location_id,
                                       p_cell_num,
                                       p_coords);

  if (sec_code == CS_RESTART_SUCCESS) {

    p_set->n_particles = n_particles;

    if (p_set->n_particles_max < p_set->n_particles)
      cs_lagr_particle_set_resize(p_set->n_particles);

    _set_particle_values(p_set, CS_LAGR_COORDS, CS_REAL_TYPE,
                         3, -1, p_coords);

    _set_particle_values(p_set, CS_LAGR_CELL_NUM, CS_LNUM_TYPE,
                         1, -1, p_cell_num);

  }

  BFT_FREE(p_cell_num);
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
      restart_type = CS_TYPE_cs_int_t;
    else if (datatype == CS_GNUM_TYPE)
      restart_type = CS_TYPE_cs_gnum_t;
    else {
      assert(datatype == CS_REAL_TYPE);
      restart_type = CS_TYPE_cs_real_t;
    }

    switch(attr) {

    case CS_LAGR_COORDS:
    case CS_LAGR_RANK_ID:
      break;

    case CS_LAGR_NEIGHBOR_FACE_ID:
      _lagr_section_name(attr, -1, sec_name);
      sec_code = cs_restart_read_ids(r,
                                     sec_name,
                                     particles_location_id,
                                     CS_MESH_LOCATION_BOUNDARY_FACES,
                                     1, /* numbering base */
                                     (cs_lnum_t *)vals);
      for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
        ((cs_lnum_t *)vals)[i] -= 1;
      if (sec_code == CS_RESTART_SUCCESS) {
        _set_particle_values(p_set,
                             attr,
                             CS_LNUM_TYPE,
                             1,
                             -1,
                             vals);
        if (attr == CS_LAGR_STAT_WEIGHT) {
          cs_real_t *w = (cs_real_t *)vals;
          for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
            p_set->weight += w[i];
        }
      }
      retval += 1;
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

        _lagr_section_name(attr, comp_id, sec_name);
        _legacy_section_names(attr, comp_id,
                              old_name_0, old_name_1, old_name_2);

        /* Special case for cell number: sign used for deposition */

        if (attr == CS_LAGR_CELL_NUM) {
          _section_name("particle_status_flag", comp_id, sec_name);
          strcpy(old_name_0, "particle_status_flag");
        }

        /* Try to read data */

        if (   attr == CS_LAGR_VELOCITY
            || attr == CS_LAGR_VELOCITY_SEEN) {
          sec_code = cs_restart_read_real_3_t_compat(r,
                                                     sec_name,
                                                     old_name_0,
                                                     old_name_1,
                                                     old_name_2,
                                                     particles_location_id,
                                                     (cs_real_3_t *)vals);
        }
        else
          sec_code = cs_restart_read_section_compat(r,
                                                    sec_name,
                                                    old_name_0,
                                                    particles_location_id,
                                                    stride,
                                                    restart_type,
                                                    vals);

        /* Finish setting values */

        if (sec_code == CS_RESTART_SUCCESS && attr == CS_LAGR_CELL_NUM) {

          cs_lnum_t  *flag = (cs_lnum_t *)vals;
          for (cs_lnum_t j = 0; j < n_particles; j++) {
            if (flag[j] == 1) {
              cs_lnum_t cell_num
                = cs_lagr_particles_get_lnum(p_set, j, CS_LAGR_CELL_NUM);
              cs_lagr_particles_set_lnum(p_set, j, CS_LAGR_CELL_NUM, -cell_num);
            }
          }

          retval += 1;

        }

        else if (sec_code == CS_RESTART_SUCCESS) {

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

  cs_lnum_t  *p_cell_num;
  cs_real_t  *p_coords;

  BFT_MALLOC(p_cell_num, n_particles, cs_lnum_t);
  BFT_MALLOC(p_coords, n_particles*3, cs_real_t);

  cs_lagr_get_particle_values(p_set, CS_LAGR_COORDS, CS_REAL_TYPE,
                              3, -1, n_particles, NULL, p_coords);

  cs_lagr_get_particle_values(p_set, CS_LAGR_CELL_NUM, CS_LNUM_TYPE,
                              1, -1, n_particles, NULL, p_cell_num);

  const int particles_location_id
    = cs_restart_write_particles(r,
                                 sec_name,
                                 false, /* number_by_coords */
                                 n_particles,
                                 p_cell_num,
                                 p_coords);

  BFT_FREE(p_cell_num);
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
      restart_type = CS_TYPE_cs_int_t;
    else if (datatype == CS_GNUM_TYPE)
      restart_type = CS_TYPE_cs_gnum_t;
    else {
      assert(datatype == CS_REAL_TYPE);
      restart_type = CS_TYPE_cs_real_t;
    }

    switch(attr) {

    case CS_LAGR_COORDS:
    case CS_LAGR_RANK_ID:
      break;

    case CS_LAGR_NEIGHBOR_FACE_ID:
      cs_lagr_get_particle_values(p_set,
                                  attr,
                                  CS_LNUM_TYPE,
                                  1,
                                  -1,
                                  n_particles,
                                  NULL,
                                  vals);
      _lagr_section_name(attr, -1, sec_name);
      for (cs_lnum_t i = 0; i < p_set->n_particles; i++)
        ((cs_lnum_t *)vals)[i] += 1;
      cs_restart_write_ids(r,
                           sec_name,
                           particles_location_id,
                           CS_MESH_LOCATION_BOUNDARY_FACES,
                           1, /* numbering base */
                           (const cs_lnum_t *)vals);
      retval += 1;
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

        /* Special case for cell number: sign used for deposition */

        if (attr == CS_LAGR_CELL_NUM) {
          cs_lnum_t  *flag = (cs_lnum_t *)vals;
          for (cs_lnum_t j = 0; j < n_particles; j++) {
            if (flag[j] < 0)
              flag[j] = 1;
            else
              flag[j] = 0;
          }
          _section_name("particle_status_flag", comp_id, sec_name);
        }

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
