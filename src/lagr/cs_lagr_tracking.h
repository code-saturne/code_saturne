#ifndef __CS_LAGR_TRACKING_H__
#define __CS_LAGR_TRACKING_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_LAGR_CELL_NUM,         /* local cell number */
  CS_LAGR_RANK_ID,          /* local parallel rank id */

  CS_LAGR_SWITCH_ORDER_1,

  CS_LAGR_RANDOM_VALUE,     /* random value associated with the particle */

  CS_LAGR_STAT_WEIGHT,
  CS_LAGR_RESIDENCE_TIME,
  CS_LAGR_MASS,
  CS_LAGR_DIAMETER,
  CS_LAGR_TAUP_AUX,
  CS_LAGR_COORDS,
  CS_LAGR_VELOCITY,
  CS_LAGR_VELOCITY_SEEN,

  /* Arrays for 2nd order scheme */

  CS_LAGR_TURB_STATE_1,        /* turbulence characteristics of first pass */
  CS_LAGR_PRED_VELOCITY,       /* 1st step prediction for particle velocity */
  CS_LAGR_PRED_VELOCITY_SEEN,  /* 1st step prediction for relative velocity */

  /* Deposition submodel additional parameters */

  CS_LAGR_YPLUS,
  CS_LAGR_INTERF,
  CS_LAGR_NEIGHBOR_FACE_ID,
  CS_LAGR_MARKO_VALUE,
  CS_LAGR_DEPOSITION_FLAG,

  /* Resuspension model additional parameters */

  CS_LAGR_N_LARGE_ASPERITIES,
  CS_LAGR_N_SMALL_ASPERITIES,
  CS_LAGR_ADHESION_FORCE,
  CS_LAGR_ADHESION_TORQUE,
  CS_LAGR_DISPLACEMENT_NORM,

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

  CS_LAGR_COAL_NUM,
  CS_LAGR_COAL_DENSITY,

  /* Radiative model additional parameters */

  CS_LAGR_EMISSIVITY,

  /* Statistical class */

  CS_LAGR_STAT_CLASS,

  /* User variables */

  CS_LAGR_USER,

  /* End of attributes */

  CS_LAGR_N_ATTRIBUTES

} cs_lagr_attribute_t;

/* Particle structure mapping */
/* -------------------------- */

typedef struct {

  size_t          extents;                         /* size (in bytes) of particle
                                                      structure */

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

/* Linked list */
/* ----------- */

typedef struct _cs_lagr_tracking_list_t cs_lagr_tracking_list_t;

/* Particle set */
/* ------------ */

typedef struct {

  int        time_id;                         /* 0 for current time,
                                                 -1 for previous */
  cs_lnum_t  n_particles;
  cs_lnum_t  n_part_out;
  cs_lnum_t  n_part_dep;
  cs_lnum_t  n_part_fou;
  cs_lnum_t  n_failed_part;

  cs_real_t  weight;
  cs_real_t  weight_out;
  cs_real_t  weight_dep;
  cs_real_t  weight_fou;
  cs_real_t  weight_failed;

  cs_lnum_t  n_particles_max;

  cs_lnum_t  first_used_id;
  cs_lnum_t  first_free_id;

  const cs_lagr_attribute_map_t  *p_am;       /* particle attributes maps
                                                 (p_am + i for time n-i) */
  unsigned char                  *p_buffer;   /* Particles data buffer */

  cs_lagr_tracking_list_t        *used_id;    /* active particles list,
                                                 or NULL for secondary sets */

} cs_lagr_particle_set_t;

/* Global parameters for Lagrangian module */
/*-----------------------------------------*/

typedef struct {

  int  physical_model;  /* FIXME: => enum: CS_LAGR_PHYS_STD,
                                           CS_LAGR_PHYS_COAL,
                                           CS_LAGR_PHYS_HEAT... */
  int  n_temperature_layers;

  int  deposition;
  int  roughness;
  int  resuspension;
  int  clogging;

  int  n_stat_classes;
  int  n_user_variables;

  int  t_order;          /* Algorithm order in time */

} cs_lagr_param_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

extern const char *cs_lagr_attribute_name[];

/* Pointer to global Lagragian module parameters */

extern const cs_lagr_param_t  *cs_glob_lagr_params;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate cs_lagr_particle_set_t structure and initialize useful buffers
 * and indexes
 *
 * parameters:
 *   nordre          <--  time algorithm order (1 or 2)
 *   iphyla          <--  kind of physics used for the lagrangian approach
 *   nvls            <--  number of user-defined variables
 *   nbclst          <--  number of stat. class to study sub-set of particles
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagbeg, LAGBEG)(const cs_int_t    *nordre,
                          const cs_int_t    *nlayer,
                          const cs_int_t    *iphyla,
                          const cs_int_t    *idepst,
                          const cs_int_t    *irough,
                          const cs_int_t    *ireent,
                          const cs_int_t    *iclogst,
                          const cs_int_t    *nvls,
                          const cs_int_t    *nbclst,
                          cs_lnum_t          icocel[],
                          cs_lnum_t          itycel[],
                          cs_int_t          *jisor,
                          cs_int_t          *jisora,
                          cs_int_t          *jirka,
                          cs_int_t          *jord1,
                          cs_int_t          *jrval,
                          cs_int_t          *jrpoi,
                          cs_int_t          *jrtsp,
                          cs_int_t          *jdp,
                          cs_int_t          *jmp,
                          cs_int_t          *jxp,
                          cs_int_t          *jyp,
                          cs_int_t          *jzp,
                          cs_int_t          *jup,
                          cs_int_t          *jvp,
                          cs_int_t          *jwp,
                          cs_int_t          *juf,
                          cs_int_t          *jvf,
                          cs_int_t          *jwf,
                          cs_int_t          *jtaux,
                          cs_int_t           jbx1[3],
                          cs_int_t           jtsup[3],
                          cs_int_t           jtsuf[3],
                          cs_int_t          *jryplu,
                          cs_int_t          *jrinpf,
                          cs_int_t          *jdfac,
                          cs_int_t          *jimark,
                          cs_int_t          *jtp,
                          cs_int_t           jhp[],
                          cs_int_t          *jtf,
                          cs_int_t          *jmwat,
                          cs_int_t           jmch[],
                          cs_int_t           jmck[],
                          cs_int_t          *jcp,
                          cs_int_t          *jrdck,
                          cs_int_t          *jrd0p,
                          cs_int_t          *jinch,
                          cs_int_t           jrhock[],
                          cs_int_t          *jreps,
                          cs_int_t          *jdepo,
                          cs_int_t          *jnbasg,
                          cs_int_t          *jnbasp,
                          cs_int_t          *jfadh,
                          cs_int_t          *jmfadh,
                          cs_int_t          *jndisp,
                          cs_int_t          *jclst,
                          cs_int_t          *jvls);

/*----------------------------------------------------------------------------
 * Get variables and parameters associated to each particles and keep it in
 * a new structure
 *
 * parameters:
 *   ...
 *----------------------------------------------------------------------------*/

void
CS_PROCF (getbdy, GETBDY)(const cs_int_t    *nflagm,
                          const cs_int_t    *nfrlag,
                          const cs_int_t    *injcon,
                          const cs_int_t     ilflag[],
                          const cs_int_t     iusncl[],
                          const cs_int_t     iusclb[],
                          const cs_real_t    deblag[],
                          const cs_int_t     ifrlag[]);

/*----------------------------------------------------------------------------
 * Displacement of particles.
 *
 * parameters:
 *   scheme_order  <-- current order of the scheme used for Lagragian
 *----------------------------------------------------------------------------*/

void
CS_PROCF (dplprt, DPLPRT)(cs_int_t        *p_scheme_order,
                          cs_real_t        boundary_stat[],
                          const cs_int_t  *iensi3,
                          const cs_int_t  *inbr,
                          const cs_int_t  *inbrbd,
                          const cs_int_t  *iflm,
                          const cs_int_t  *iflmbd,
                          const cs_int_t  *iang,
                          const cs_int_t  *iangbd,
                          const cs_int_t  *ivit,
                          const cs_int_t  *ivitbd,
                          const cs_int_t  *iencnd,
                          const cs_int_t  *iencma,
                          const cs_int_t  *iencdi,
                          const cs_int_t  *iencck,
                          const cs_int_t  *iencnbbd,
                          const cs_int_t  *iencmabd,
                          const cs_int_t  *iencdibd,
                          const cs_int_t  *iencckbd,
                          const cs_int_t  *inclg,
                          const cs_int_t  *iscovc,
                          const cs_int_t  *nusbor,
                          cs_int_t         iusb[],
                          cs_real_t        visc_length[],
                          cs_real_t        dlgeo[],
                          cs_real_t        energt[],
                          const cs_real_t  tprenc[],
                          const cs_real_t  visref[],
                          const cs_real_t  enc1[],
                          const cs_real_t  enc2[],
                          const cs_real_t  *tkelvi);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get data access information for a given particle attribute.
 *
 * For attributes not currently present, the displacement and data
 * size should be -1 and 0 respectively.
 *
 * parameters:
 *   particles <-- associated particle set
 *   time_id   <-- associated time id (0: current, 1: previous)
 *   attr      <-- particle attribute
 *   extents   --> size (in bytes) of particle structure, or NULL
 *   size      --> size (in bytes) of attribute in particle structure, or NULL
 *   displ     --> displacement (in bytes) in particle structure, or NULL
 *   datatype  --> associated datatype, or NULL
 *   count     --> associated elements count, or NULL
 *----------------------------------------------------------------------------*/

void
cs_lagr_get_attr_info(const cs_lagr_particle_set_t  *particles,
                      int                            time_id,
                      cs_lagr_attribute_t            attr,
                      size_t                        *extents,
                      size_t                        *size,
                      ptrdiff_t                     *displ,
                      cs_datatype_t                 *datatype,
                      int                           *count);

/*----------------------------------------------------------------------------
 * Return pointer to the main cs_lagr_particle_set_t structure.
 *
 * returns:
 *   pointer to current particle set, or NULL
 *----------------------------------------------------------------------------*/

cs_lagr_particle_set_t  *
cs_lagr_get_particle_set(void);

/*----------------------------------------------------------------------------
 * Delete particle set structure and other useful buffers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_destroy(void);

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
cs_lagr_resize_particle_set(cs_lnum_t  n_min_particles);

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

/*----------------------------------------------------------------------------
 * \brief Set global maximum number of particles.
 *
 * By default, the number is limited only by local cs_lnum_t and global
 * cs_gnum_t data representation limits.
 *
 * parameters:
 *   n_g_particles_max <-- global maximum number of particles
*----------------------------------------------------------------------------*/

void
cs_lagr_set_n_g_particles_max(unsigned long long  n_g_particles_max);

/*----------------------------------------------------------------------------
 * Dump a cs_lagr_particle_t structure
 *
 * parameters:
 *   particles <-- cs_lagr_particle_t structure to dump
 *----------------------------------------------------------------------------*/

void
cs_lagr_particle_set_dump(const cs_lagr_particle_set_t  *particles);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_TRACKING_H__ */
