#ifndef __CS_LAGR_PARTICLE_H__
#define __CS_LAGR_PARTICLE_H__

/*============================================================================
 * Lagrangian module particle model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "assert.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * Flags specifying general event attributes
 */

/*! particle will be deleted */
#define CS_LAGR_PART_TO_DELETE         (1 << 0)

/*! particle is fixed, no resuspension possible */
#define CS_LAGR_PART_FIXED             (1 << 1)

/*! particle is deposited */
#define CS_LAGR_PART_DEPOSITED         (1 << 2)

/*! particle is rolling */
#define CS_LAGR_PART_ROLLING           (1 << 3)

/*! particle motion is prescribed */
#define CS_LAGR_PART_IMPOSED_MOTION    (1 << 4)

/*! Flag sets (useful for cancelling nay flag of a set) */

/*! particle is in flow */
#define CS_LAGR_PART_DEPOSITION_FLAGS \
  (  CS_LAGR_PART_TO_DELETE | CS_LAGR_PART_FIXED | CS_LAGR_PART_DEPOSITED \
   | CS_LAGR_PART_ROLLING | CS_LAGR_PART_IMPOSED_MOTION)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Prefedined particle attributes */
/* ------------------------------- */

typedef enum {

  CS_LAGR_P_FLAG,             /*!< local metadata flag */

  CS_LAGR_CELL_ID,             /*!< local cell id (0 to n-1) */
  CS_LAGR_RANK_ID,             /*!< local parallel rank id */

  CS_LAGR_REBOUND_ID,          /*!< number of time steps since rebound, or -1 */

  CS_LAGR_RANDOM_VALUE,        /*!< associated random value (for sampling) */

  CS_LAGR_STAT_WEIGHT,         /*!< statistical weight */
  CS_LAGR_RESIDENCE_TIME,      /*!< time since particle injection */
  CS_LAGR_MASS,
  CS_LAGR_DIAMETER,
  CS_LAGR_TAUP_AUX,
  CS_LAGR_COORDS,
  CS_LAGR_VELOCITY,
  CS_LAGR_VELOCITY_SEEN,

  CS_LAGR_TR_TRUNCATE,         /*!< portion of trajectory truncated */
  CS_LAGR_TR_REPOSITION,       /*!< number of times the particle is repositioned
                                    at the cell center, or -1 when lost (to be
                                    discarded at the next displacement step) */

  /* Arrays for 2nd order scheme */

  CS_LAGR_TURB_STATE_1,        /* turbulence characteristics of first pass */
  CS_LAGR_PRED_VELOCITY,       /* 1st step prediction for particle velocity */
  CS_LAGR_PRED_VELOCITY_SEEN,  /* 1st step prediction for relative velocity */
  CS_LAGR_V_GAUSS,             /* 1st step Gaussian variable */
  CS_LAGR_BR_GAUSS,            /* 1st step Brownian motion Gaussian variable */

  /* Non-spherical particles submodel additoinal parameters */

  CS_LAGR_SHAPE,               /*!< Shape of spheroids */
  CS_LAGR_ORIENTATION,         /*!< Orientation of spheroids */
  CS_LAGR_RADII,               /*!< Radii a, b, c for ellispoids */
  CS_LAGR_ANGULAR_VEL,         /*!< Angular velocity */
  CS_LAGR_EULER,               /*!< Euler four parameters */
  CS_LAGR_SHAPE_PARAM,         /*!< Shape parameters computed from a b c */

  /* Deposition submodel additional parameters */

  CS_LAGR_YPLUS,
  CS_LAGR_INTERF,
  CS_LAGR_NEIGHBOR_FACE_ID,
  CS_LAGR_MARKO_VALUE,
  CS_LAGR_FOULING_INDEX,

  /* Resuspension model additional parameters */

  CS_LAGR_N_LARGE_ASPERITIES,
  CS_LAGR_N_SMALL_ASPERITIES,
  CS_LAGR_ADHESION_FORCE,
  CS_LAGR_ADHESION_TORQUE,
  CS_LAGR_DISPLACEMENT_NORM,

  /* Clogging model additional parameters */

  CS_LAGR_HEIGHT,
  CS_LAGR_CLUSTER_NB_PART,
  CS_LAGR_DEPO_TIME,
  CS_LAGR_CONSOL_HEIGHT,

  /* Thermal model additional parameters */

  CS_LAGR_TEMPERATURE,
  CS_LAGR_FLUID_TEMPERATURE,
  CS_LAGR_CP,

  /* Coal combustion additional parameters */

  CS_LAGR_WATER_MASS,
  CS_LAGR_COAL_MASS,
  CS_LAGR_COKE_MASS,

  CS_LAGR_SHRINKING_DIAMETER,
  CS_LAGR_INITIAL_DIAMETER,

  CS_LAGR_COAL_ID,
  CS_LAGR_COAL_DENSITY,

  /* Radiative model additional parameters */

  CS_LAGR_EMISSIVITY,

  /* Statistical class */

  CS_LAGR_STAT_CLASS,

  CS_LAGR_PARTICLE_AGGREGATE,

  /* User attributes */

  CS_LAGR_USER,

  /* End of attributes */

  CS_LAGR_N_ATTRIBUTES

} cs_lagr_attribute_t;

/*! Particle attribute structure mapping */
/* ------------------------------------- */

typedef struct {

  size_t          extents;                         /* size (in bytes) of particle
                                                      structure */
  size_t          lb;                              /* size (in bytes) of lower
                                                      bounds of particle data
                                                      (work area before) */

  int             n_time_vals;                     /* number of time values
                                                      handled */

  size_t          size[CS_LAGR_N_ATTRIBUTES];      /* size (in bytes) of
                                                      attributes in particle
                                                      structure for a given
                                                      time value */
  cs_datatype_t   datatype[CS_LAGR_N_ATTRIBUTES];  /* datatype of associated
                                                      attributes */
  int            (*count)[CS_LAGR_N_ATTRIBUTES];   /* number of values for each
                                                      attribute, per associated
                                                      time_id */
  ptrdiff_t      (*displ)[CS_LAGR_N_ATTRIBUTES];   /* displacement (in bytes) of
                                                      attributes in particle data,
                                                      per associated time_id*/

  ptrdiff_t      *source_term_displ;               /* displacement (in bytes) of
                                                      source term values
                                                      for second-order scheme,
                                                      or NULL */

} cs_lagr_attribute_map_t;

/* Particle set */
/* ------------ */

typedef struct {

  cs_lnum_t  n_particles;                     /* number of particle in domain */
  cs_lnum_t  n_part_new;
  cs_lnum_t  n_part_out;
  cs_lnum_t  n_part_merged;                   /* number of merged particles */
  cs_lnum_t  n_part_dep;
  cs_lnum_t  n_part_fou;
  cs_lnum_t  n_part_resusp;
  cs_lnum_t  n_failed_part;

  cs_real_t  weight;
  cs_real_t  weight_new;
  cs_real_t  weight_out;
  cs_real_t  weight_merged;
  cs_real_t  weight_dep;
  cs_real_t  weight_fou;
  cs_real_t  weight_resusp;
  cs_real_t  weight_failed;

  cs_lnum_t  n_particles_max;

  const cs_lagr_attribute_map_t  *p_am;       /*!< particle attributes maps
                                                   (p_am + i for time n-i) */
  unsigned char                  *p_buffer;   /*!< Particles data buffer */

} cs_lagr_particle_set_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*! Names of particle attributes */

extern const char *cs_lagr_attribute_name[];

/*! Pointer to the main particle set */

extern cs_lagr_particle_set_t  *cs_glob_lagr_particle_set;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle map based on defined options.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_attr_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return const pointer to the main particle attribute map structure.
 *
 * \return pointer to current particle attribute map, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_lagr_attribute_map_t *
cs_lagr_particle_get_attr_map(void);

/*----------------------------------------------------------------------------*/
/*!
 * Allocate main cs_lagr_particle_set_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy main particle set and map if they exist.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy attributes from one particle to another
 *
 * The random value associated with the particle is modified.
 *
 * \param  dest  id (0-based) of destination particle
 * \param  src   id (0-based) of source particle
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_part_copy(cs_lnum_t  dest,
                  cs_lnum_t  src);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data extents for a given particle attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * \param[in]   particles  associated particle set
 * \param[in]   time_id    associated time id (0: current, 1: previous)
 * \param[in]   attr       particle attribute
 * \param[out]  extents    size (in bytes) of particle structure, or NULL
 * \param[out]  size       size (in bytes) of attribute in particle structure,
 *                         or NULL
 * \param[out]  displ      displacement (in bytes) in particle structure,
 *                         or NULL
 * \param[out]  datatype   datatype of associated attribute, or NULL
 * \param[out]  count      number of type values associated with attribute,
 *                         or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_attr_info(const cs_lagr_particle_set_t  *particles,
                      int                            time_id,
                      cs_lagr_attribute_t            attr,
                      size_t                        *extents,
                      size_t                        *size,
                      ptrdiff_t                     *displ,
                      cs_datatype_t                 *datatype,
                      int                           *count);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that query for a given particle attribute is valid.
 *
 * \param[in]   particles             associated particle set
 * \param[in]   attr                  attribute whose values are required
 * \param[in]   datatype              associated value type
 * \param[in]   stride                number of values per particle
 * \param[in]   component_id          if -1 : extract the whole attribute
 *                                    if >0 : id of the component to extract
 *
 * \return 0 in case of success, 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_check_attr_query(const cs_lagr_particle_set_t  *particles,
                         cs_lagr_attribute_t            attr,
                         cs_datatype_t                  datatype,
                         int                            stride,
                         int                            component_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a particle attribute is in a valid range.
 *
 * If this is not the case, a fatal error is provoked.

 * \param[in]   attr       particle attribute
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_attr_in_range(int  attr);

/*----------------------------------------------------------------------------
 * Return pointer to the main cs_lagr_particle_set_t structure.
 *
 * returns:
 *   pointer to current particle set, or NULL
 *----------------------------------------------------------------------------*/

cs_lagr_particle_set_t  *
cs_lagr_get_particle_set(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to a current attribute of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static void *
cs_lagr_particles_attr(cs_lagr_particle_set_t  *particle_set,
                       cs_lnum_t                particle_id,
                       cs_lagr_attribute_t      attr)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  return   (unsigned char *)particle_set->p_buffer
         + particle_set->p_am->extents*particle_id
         + particle_set->p_am->displ[0][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to current attribute data
 *        of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static const void *
cs_lagr_particles_attr_const(const cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                      particle_id,
                             cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  return   particle_set->p_buffer
         + particle_set->p_am->extents*particle_id
         + particle_set->p_am->displ[0][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to attribute data of a given particle in a set
 *        at a given time.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  time_id       0 for current, 1 for previous
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to attribute data at given time
 */
/*----------------------------------------------------------------------------*/

inline static void *
cs_lagr_particles_attr_n(cs_lagr_particle_set_t  *particle_set,
                         cs_lnum_t                particle_id,
                         int                      time_id,
                         cs_lagr_attribute_t      attr)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  return   particle_set->p_buffer
         + particle_set->p_am->extents*particle_id
         + particle_set->p_am->displ[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to attribute data of a given particle in a set
 *        at a given time.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  time_id       0 for current, 1 for previous
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to attribute data at given time
 */
/*----------------------------------------------------------------------------*/

inline static const void *
cs_lagr_particles_attr_n_const(const cs_lagr_particle_set_t  *particle_set,
                               cs_lnum_t                      particle_id,
                               int                            time_id,
                               cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  return   particle_set->p_buffer
         + particle_set->p_am->extents*particle_id
         + particle_set->p_am->displ[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get flag value with a selected mask for a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  mask          requested attribute flag value
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static int
cs_lagr_particles_get_flag(const cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                      particle_id,
                           int                            mask)
{
  int flag
    = *((const cs_lnum_t *)(  particle_set->p_buffer
                            + particle_set->p_am->extents*particle_id
                            + particle_set->p_am->displ[0][CS_LAGR_P_FLAG]));

  return (flag & mask);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set flag value with a selected mask for a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  mask          attribute flag value to set
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_flag(const cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                      particle_id,
                           int                            mask)
{
  int flag
    = *((const cs_lnum_t *)(  particle_set->p_buffer
                            + particle_set->p_am->extents*particle_id
                            + particle_set->p_am->displ[0][CS_LAGR_P_FLAG]));

  flag = flag | mask;

  *((cs_lnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[0][CS_LAGR_P_FLAG])) = flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Unset flag value with a selected mask for a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  mask          attribute flag value to set
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_unset_flag(const cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                      particle_id,
                             int                            mask)
{
  int flag
    = *((const cs_lnum_t *)(  particle_set->p_buffer
                            + particle_set->p_am->extents*particle_id
                            + particle_set->p_am->displ[0][CS_LAGR_P_FLAG]));

  flag = (flag | mask) - mask;

  *((cs_lnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[0][CS_LAGR_P_FLAG])) = flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_lnum_t of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_lagr_particles_get_lnum(const cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                      particle_id,
                           cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  return *((const cs_lnum_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_lnum_t of a given particle in a set
 *        at a given time.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  time_id       0 for current, 1 for previous
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_lagr_particles_get_lnum_n(const cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                      particle_id,
                             int                            time_id,
                             cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  return *((const cs_lnum_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_lnum_t of a given particle in a set.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_lnum(cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                particle_id,
                           cs_lagr_attribute_t      attr,
                           cs_lnum_t                value)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  *((cs_lnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[0][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_lnum_t of a given particle in a set
 *        at a given time.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        time_id       0 for current, 1 for previous
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_lnum_n(cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                particle_id,
                             int                      time_id,
                             cs_lagr_attribute_t      attr,
                             cs_lnum_t                value)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  *((cs_lnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_gnum_t of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_gnum_t
cs_lagr_particles_get_gnum(const cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                      particle_id,
                           cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  return *((const cs_gnum_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_gnum_t of a given particle in a set
 *        at a given time.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  time_id       0 for current, 1 for previous
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_gnum_t
cs_lagr_particles_get_gnum_n(const cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                      particle_id,
                             int                            time_id,
                             cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  return *((const cs_gnum_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_gnum_t of a given particle in a set.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_gnum(cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                particle_id,
                           cs_lagr_attribute_t      attr,
                           cs_gnum_t                value)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  *((cs_gnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[0][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_gnum_t of a given particle in a set
 *        at a given time.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        time_id       0 for current, 1 for previous
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_gnum_n(cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                particle_id,
                             int                      time_id,
                             cs_lagr_attribute_t      attr,
                             cs_gnum_t                value)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  *((cs_gnum_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_real_t of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_lagr_particles_get_real(const cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                      particle_id,
                           cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  return *((const cs_real_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_real_t of a given particle in a set
 *        at a given time.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  time_id       0 for current, 1 for previous
 * \param[in]  attr          requested attribute id
 *
 * \return  attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_lagr_particles_get_real_n(const cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                      particle_id,
                             int                            time_id,
                             cs_lagr_attribute_t            attr)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  return *((const cs_real_t *)(  particle_set->p_buffer
                               + particle_set->p_am->extents*particle_id
                               + particle_set->p_am->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_real_t of a given particle in a set.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_real(cs_lagr_particle_set_t  *particle_set,
                           cs_lnum_t                particle_id,
                           cs_lagr_attribute_t      attr,
                           cs_real_t                value)
{
  assert(particle_set->p_am->count[0][attr] > 0);

  *((cs_real_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[0][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_real_t of a given particle in a set
 *        at a given time.
 *
 * \param[in, out]   particle_set  pointer to particle set
 * \param[in]        particle_id   particle id
 * \param[in]        time_id       0 for current, 1 for previous
 * \param[in]        attr          requested attribute id
 * \param[in]        value         value to assign
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particles_set_real_n(cs_lagr_particle_set_t  *particle_set,
                             cs_lnum_t                particle_id,
                             int                      time_id,
                             cs_lagr_attribute_t      attr,
                             cs_real_t                value)
{
  assert(particle_set->p_am->count[time_id][attr] > 0);

  *((cs_real_t *)(  particle_set->p_buffer
                  + particle_set->p_am->extents*particle_id
                  + particle_set->p_am->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to 2nd order scheme source terms for an attribute
 *        of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t *
cs_lagr_particles_source_terms(cs_lagr_particle_set_t  *particle_set,
                               cs_lnum_t                particle_id,
                               cs_lagr_attribute_t      attr)
{
  assert(particle_set->p_am->source_term_displ != NULL);
  assert(particle_set->p_am->source_term_displ[attr] >= 0);

  return (cs_real_t *)(  (unsigned char *)particle_set->p_buffer
                       + particle_set->p_am->extents*particle_id
                       + particle_set->p_am->source_term_displ[attr]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to 2nd order scheme source terms an attribute
 *        of a given particle in a set.
 *
 * \param[in]  particle_set  pointer to particle set
 * \param[in]  particle_id   particle id
 * \param[in]  attr          requested attribute id
 *
 * \return    pointer to current attribute data
 */
/*----------------------------------------------------------------------------*/

inline static const cs_real_t *
cs_lagr_particles_source_terms_const(cs_lagr_particle_set_t  *particle_set,
                                     cs_lnum_t                particle_id,
                                     cs_lagr_attribute_t      attr)
{
  assert(particle_set->p_am->source_term_displ != NULL);
  assert(particle_set->p_am->source_term_displ[attr] >= 0);

  return (const cs_real_t *)(  (unsigned char *)particle_set->p_buffer
                             + particle_set->p_am->extents*particle_id
                             + particle_set->p_am->source_term_displ[attr]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to current attribute data of a particle.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr      requested attribute id
 *
 * \return  pointer to attribute data
 */
/*----------------------------------------------------------------------------*/

inline static void *
cs_lagr_particle_attr(void                           *particle,
                      const cs_lagr_attribute_map_t  *attr_map,
                      cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[0][attr] > 0);

  return  (unsigned char *)particle + attr_map->displ[0][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to current attribute data of a particle.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr      requested attribute id
 *
 * \return  const pointer to attribute
 */
/*----------------------------------------------------------------------------*/

inline static const void *
cs_lagr_particle_attr_const(const void                     *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[0][attr] > 0);

  return  (const unsigned char *)particle + attr_map->displ[0][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to attribute data of a particle at a given time.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  time_id   0 for current, 1 for previous
 * \param[in]  attr      requested attribute id
 *
 * \return  pointer to attribute data
 */
/*----------------------------------------------------------------------------*/

inline static void *
cs_lagr_particle_attr_n(void                           *particle,
                        const cs_lagr_attribute_map_t  *attr_map,
                        int                             time_id,
                        cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[time_id][attr] > 0);

  return  (unsigned char *)particle + attr_map->displ[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get const pointer to attribute data of a particle at a given time.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  time_id   0 for current, 1 for previous
 * \param[in]  attr      requested attribute id
 *
 * \return  pointer to attribute data
 */
/*----------------------------------------------------------------------------*/

inline static const void *
cs_lagr_particle_attr_n_const(const void                     *particle,
                              const cs_lagr_attribute_map_t  *attr_map,
                              int                             time_id,
                              cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[time_id][attr] > 0);

  return    (const unsigned char *)particle
          + attr_map->displ[time_id][attr];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_lnum_t of a given particle in a set
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr     requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_lagr_particle_get_lnum(const void                     *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[0][attr] > 0);

  return  *((const cs_lnum_t *)(  (const unsigned char *)particle
                                + attr_map->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_lnum_t of a given particle
 *        at a given time.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  time_id   0 for current, 1 for previous
 * \param[in]  attr      requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_lagr_particle_get_lnum_n(const void                     *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[time_id][attr] > 0);

  return  *((const cs_lnum_t *)(  (const unsigned char *)particle
                                + attr_map->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_lnum_t of a given particle.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_lnum(void                           *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr,
                          cs_lnum_t                       value)
{
  assert(attr_map->count[0][attr] > 0);

  *((cs_lnum_t *)((unsigned char *)particle + attr_map->displ[0][attr]))
    = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_lnum_t of a given particle
 *        at a given time.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       time_id   0 for current, 1 for previous
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_lnum_n(void                           *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr,
                            cs_lnum_t                       value)
{
  assert(attr_map->count[time_id][attr] > 0);

  *((cs_lnum_t *)(  (unsigned char *)particle
                  + attr_map->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_gnum_t of a given particle in a set
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr     requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_gnum_t
cs_lagr_particle_get_gnum(const void                     *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[0][attr] > 0);

  return  *((const cs_gnum_t *)(  (const unsigned char *)particle
                                + attr_map->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_gnum_t of a given particle
 *        at a given time.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  time_id   0 for current, 1 for previous
 * \param[in]  attr      requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_gnum_t
cs_lagr_particle_get_gnum_n(const void                     *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[time_id][attr] > 0);

  return  *((const cs_gnum_t *)(  (const unsigned char *)particle
                                + attr_map->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_gnum_t of a given particle.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_gnum(void                           *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr,
                          cs_gnum_t                       value)
{
  assert(attr_map->count[0][attr] > 0);

  *((cs_gnum_t *)((unsigned char *)particle + attr_map->displ[0][attr]))
    = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_gnum_t of a given particle
 *        at a given time.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       time_id   0 for current, 1 for previous
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_gnum_n(void                           *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr,
                            cs_gnum_t                       value)
{
  assert(attr_map->count[time_id][attr] > 0);

  *((cs_gnum_t *)(  (unsigned char *)particle
                  + attr_map->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_real_t of a given particle in a set
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr     requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_lagr_particle_get_real(const void                     *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[0][attr] > 0);

  return  *((const cs_real_t *)(  (const unsigned char *)particle
                                + attr_map->displ[0][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get attribute value of type cs_real_t of a given particle
 *        at a given time.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  time_id   0 for current, 1 for previous
 * \param[in]  attr      requested attribute id
 *
 * \return attribute value
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
cs_lagr_particle_get_real_n(const void                     *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr)
{
  assert(attr_map->count[time_id][attr] > 0);

  return  *((const cs_real_t *)(  (const unsigned char *)particle
                                + attr_map->displ[time_id][attr]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_real_t of a given particle.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_real(void                           *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lagr_attribute_t             attr,
                          cs_real_t                       value)
{
  assert(attr_map->count[0][attr] > 0);

  *((cs_real_t *)((unsigned char *)particle + attr_map->displ[0][attr]))
    = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set attribute value of type cs_real_t of a given particle
 *        at a given time.
 *
 * \param[in, out]  particle  pointer to particle data
 * \param[in]       attr_map  pointer to attribute map
 * \param[in]       time_id   0 for current, 1 for previous
 * \param[in]       attr      requested attribute id
 * \param[in]       value     value to assign
 */
 /*----------------------------------------------------------------------------*/

inline static void
cs_lagr_particle_set_real_n(void                           *particle,
                            const cs_lagr_attribute_map_t  *attr_map,
                            int                             time_id,
                            cs_lagr_attribute_t             attr,
                            cs_real_t                       value)
{
  assert(attr_map->count[time_id][attr] > 0);

  *((cs_real_t *)(  (unsigned char *)particle
                  + attr_map->displ[time_id][attr])) = value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to 2nd order scheme attribute source terms of a particle.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr      requested attribute id
 *
 * \return  pointer to attribute source terms
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t *
cs_lagr_particle_source_term(void                           *particle,
                             const cs_lagr_attribute_map_t  *attr_map,
                             cs_lagr_attribute_t             attr)
{
  assert(attr_map->source_term_displ != NULL);
  assert(attr_map->source_term_displ[attr] >= 0);

  return  (cs_real_t *)(  (unsigned char *)particle
                        + attr_map->source_term_displ[attr]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to 2nd order scheme attribute source terms of a particle.
 *
 * \param[in]  particle  pointer to particle data
 * \param[in]  attr_map  pointer to attribute map
 * \param[in]  attr      requested attribute id
 *
 * \return  pointer to attribute source terms
 */
/*----------------------------------------------------------------------------*/

inline static const cs_real_t *
cs_lagr_particle_source_term_const(void                           *particle,
                                   const cs_lagr_attribute_map_t  *attr_map,
                                   cs_lagr_attribute_t             attr)
{
  assert(attr_map->source_term_displ != NULL);
  assert(attr_map->source_term_displ[attr] >= 0);

  return  (const cs_real_t *)(  (unsigned char *)particle
                              + attr_map->source_term_displ[attr]);
}

/*----------------------------------------------------------------------------
 * Resize particle set buffers if needed.
 *
 * parameters:
 *   n_particles <-- minumum number of particles required
 *
 *
 * returns:
 *   1 if resizing was required, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_lagr_particle_set_resize(cs_lnum_t  n_min_particles);

/*----------------------------------------------------------------------------
 * Set reallocation factor for particle sets.
 *
 * This factor determines the multiplier used for reallocations when
 * the particle set's buffers are too small to handle the new number of
 * particles.
 *
 * parameters:
 *  f <-- reallocation size multiplier
 *----------------------------------------------------------------------------*/

void
cs_lagr_set_reallocation_factor(double f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global maximum number of particles.
 *
 * By default, the number is limited only by local \ref cs_lnum_t and global
 * \ref cs_gnum_t data representation limits.
 *
 * \return  global maximum number of particles
 */
/*----------------------------------------------------------------------------*/

unsigned long long
cs_lagr_get_n_g_particles_max(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global maximum number of particles.
 *
 * By default, the number is limited only by local \ref cs_lnum_t and global
 * \ref cs_gnum_t data representation limits.
 *
 * \param[in]  n_g_particles_max  global maximum number of particles
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_n_g_particles_max(unsigned long long  n_g_particles_max);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy current attributes to previous attributes.
 *
 * \param[in, out]  particles     associated particle set
 * \param[in]       particle_id  id of particle
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particles_current_to_previous(cs_lagr_particle_set_t  *particles,
                                      cs_lnum_t                particle_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_lagr_particle_set_t structure
 *
 * \param[in]  particles  cs_lagr_particle_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_dump(const cs_lagr_particle_set_t  *particles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set number of user particle variables.
 *
 * \param[in]  n_user_variables  number of user variables
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_n_user_variables(int  n_user_variables);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_PARTICLE_H__ */
