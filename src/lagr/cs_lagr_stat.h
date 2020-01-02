#ifndef __CS_LAGR_STAT_H__
#define __CS_LAGR_STAT_H__

/*============================================================================
 * Functions and types for the Lagrangian module
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
#include "cs_base.h"
#include "cs_field.h"
#include "cs_restart.h"

#include "cs_lagr.h"
#include "cs_lagr_event.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

 /*! Particle statistics moment type */

typedef enum {

  CS_LAGR_MOMENT_MEAN,
  CS_LAGR_MOMENT_VARIANCE

} cs_lagr_stat_moment_t;

 /*! Particle statistics moment data type */

typedef enum {

  CS_LAGR_STAT_GROUP_PARTICLE,
  CS_LAGR_STAT_GROUP_TRACKING_EVENT,

  CS_LAGR_STAT_GROUP_N_GROUPS

} cs_lagr_stat_group_t;

/*! Moment restart behavior */

typedef enum {

  CS_LAGR_MOMENT_RESTART_RESET,
  CS_LAGR_MOMENT_RESTART_AUTO,
  CS_LAGR_MOMENT_RESTART_EXACT

} cs_lagr_stat_restart_t;

/*! Predefined particle and event statistics */
/* ----------------------------------------- */

typedef enum {

  /* Volume statistics */

  CS_LAGR_STAT_CUMULATIVE_WEIGHT,        /*!< cumulative particle statistical
                                           weight (active if any particle
                                           attribute statistics are active,
                                           may be activated separately) */
  CS_LAGR_STAT_VOLUME_FRACTION,          /*!< particle volume fraction */

  /* Boundary statistics */

  /*! cumulative particle event  statistical weight */
  CS_LAGR_STAT_E_CUMULATIVE_WEIGHT,

  /*! cumulative fouling event statistical weight */
  CS_LAGR_STAT_RESUSPENSION_CUMULATIVE_WEIGHT,

  /*! cumulative fouling event statistical weight */
  CS_LAGR_STAT_FOULING_CUMULATIVE_WEIGHT,

  CS_LAGR_STAT_MASS_FLUX,                /*!< particle mass flux */
  CS_LAGR_STAT_RESUSPENSION_MASS_FLUX,   /*!< particle resuspension mass flux */
  CS_LAGR_STAT_FOULING_MASS_FLUX,        /*!< particle fouling mass flux */

  CS_LAGR_STAT_IMPACT_ANGLE,             /*!< particle impact angle
                                           (in radians) */
  CS_LAGR_STAT_IMPACT_VELOCITY,          /*!< particle impact velocity */
  CS_LAGR_STAT_FOULING_DIAMETER,         /*!< fouled particle diameter */
  CS_LAGR_STAT_FOULING_COKE_FRACTION,    /*!< fouled particle coke fraction */

  /* Particle or event-based attributes */

  CS_LAGR_STAT_ATTR                      /*!< particle or event attribute; add
                                           attribute id for given attribute */

} cs_lagr_stat_type_t;

/*----------------------------------------------------------------------------
 * Function pointer for computation of particle data values for
 * Lagrangian statistics.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   input    <-- pointer to optional (untyped) value or structure.
 *   particle <-- pointer to particle data
 *   p_am     <-- pointer to particle attribute map
 *   vals     --> pointer to values
 *----------------------------------------------------------------------------*/

typedef void
(cs_lagr_moment_p_data_t) (const void                     *input,
                           const void                     *particle,
                           const cs_lagr_attribute_map_t  *p_am,
                           cs_real_t                       vals[]);

/*----------------------------------------------------------------------------
 * Function pointer for computation of event data values for
 * Lagrangian statistics.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   input     <-- pointer to optional (untyped) value or structure.
 *   events    <-- pointer to events
 *   event_id  <-- event id range (first to past-last)
 *   vals      --> pointer to values
 *----------------------------------------------------------------------------*/

typedef void
(cs_lagr_moment_e_data_t) (const void                 *input,
                           const cs_lagr_event_set_t  *events,
                           cs_lnum_t                   id_range[2],
                           cs_real_t                   vals[]);

/*----------------------------------------------------------------------------
 * Function pointer for computation of data values for particle statistics
 * based on mesh
 *
 * If the matching values are multidimensional, they must be interleaved.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   input       <-- pointer to optional value or structure, or NULL
 *   events      <-- pointer to optional events set, or NULL.
 *   location_id <-- associated mesh location id
 *   class_id    <-- associated particle class id (0 for all)
 *   vals        --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

typedef void
(cs_lagr_moment_m_data_t) (const void                 *input,
                           const cs_lagr_event_set_t  *events,
                           int                         location_id,
                           int                         class_id,
                           cs_real_t                   vals[]);

/*! Structure defining Lagrangian statistics options */

typedef struct {

  /*! during a Lagrangian calculation restart, indicates whether the particle
    statistics (volume and boundary) and two-way coupling terms are to be read
    from a restart file (=1) or reinitialized (=0).
    Useful if \ref isuila = 1 */
  int  isuist;

  /*! absolute time step number (including the restarts) after
    which the calculation of the volume statistics is activated. */
  int  idstnt;

  /*! absolute time step number (includings the restarts) after
    which the volume statistics are cumulated over time (they are then said
    to be steady).
    if the absolute time step number is lower than \ref nstist,
    or if the flow is unsteady (\ref isttio=0), the statistics are reset
    to zero at every time step (the volume statistics are then said
    to be non-steady).
    Useful if \ref isttio=1 */
  int  nstist;

  /*! threshold for statistical meaning when used by other model
    features (such as the Poisson correction) */
  cs_real_t  threshold;

} cs_lagr_stat_options_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointer to global statistic options structure */

extern cs_lagr_stat_options_t   *cs_glob_lagr_stat_options;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particle-based statistic.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_type      predefined statistics type, or -1
 * \param[in]  stat_group     statistics group (particle or event)
 * \param[in]  m_type         moment type
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  component_id   attribute component id, or < 0 for all
 * \param[in]  data_func      pointer to function to compute statistics
 *                            (if stat_type < 0)
 * \param[in]  data_input     associated input
 * \param[in]  w_data_func    pointer to function to compute weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  w_data_input   associated input for w_data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_particle_define(const char                *name,
                             int                        location_id,
                             int                        stat_type,
                             cs_lagr_stat_moment_t      m_type,
                             int                        class_id,
                             int                        dim,
                             int                        component_id,
                             cs_lagr_moment_p_data_t   *data_func,
                             void                      *data_input,
                             cs_lagr_moment_p_data_t   *w_data_func,
                             void                      *w_data_input,
                             int                        nt_start,
                             double                     t_start,
                             cs_lagr_stat_restart_t     restart_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an event-based statistic.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_type      predefined statistics type, or -1
 * \param[in]  stat_group     statistics group (event type)
 * \param[in]  m_type         moment type
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  component_id   attribute component id, or < 0 for all
 * \param[in]  data_func      pointer to function to compute statistics
 *                            (if stat_type < 0)
 * \param[in]  data_input     associated input
 * \param[in]  w_data_func    pointer to function to compute weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  w_data_input   associated input for w_data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_event_define(const char                *name,
                          int                        location_id,
                          int                        stat_type,
                          cs_lagr_stat_group_t       stat_group,
                          cs_lagr_stat_moment_t      m_type,
                          int                        class_id,
                          int                        dim,
                          int                        component_id,
                          cs_lagr_moment_e_data_t   *data_func,
                          void                      *data_input,
                          cs_lagr_moment_e_data_t   *w_data_func,
                          void                      *w_data_input,
                          int                        nt_start,
                          double                     t_start,
                          cs_lagr_stat_restart_t     restart_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a particle weight type statistic.
 *
 * Weights are automatically associated to general statitistics, but defining
 * them explicitely allows activation of standard logging and postprocessing
 * for those weights, as well as defining specific weights.
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_group     statistics group (particle or event)
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  p_data_func    pointer to function to compute particle weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  e_data_func    pointer to function to compute event weight
 *                            (if NULL, statistic weight assumed)
 * \param[in]  data_input     associated input for data_func
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_accumulator_define(const char                *name,
                                int                        location_id,
                                cs_lagr_stat_group_t       stat_group,
                                int                        class_id,
                                cs_lagr_moment_p_data_t   *p_data_func,
                                cs_lagr_moment_e_data_t   *e_data_func,
                                void                      *data_input,
                                int                        nt_start,
                                double                     t_start,
                                cs_lagr_stat_restart_t     restart_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh-based statistic based on particles or particle events.
 *
 * This type of statistic is reinitialized and evaluated during each time step,
 * but may be computed incrementally when based on particle events, so the
 * associated data function must uptate the statistics without reinitializing
 * them at each call.
 *
 * As this type of statistic does not need to keep state between time steps,
 * it is ignored by the lagragian statistics checkpoint/restart mechanism.
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_group     statistics group (particle or event)
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  data_func      pointer to function to compute statistics
 * \param[in]  data_input     associated input
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_mesh_define(const char                *name,
                         int                        location_id,
                         cs_lagr_stat_group_t       stat_group,
                         int                        class_id,
                         int                        dim,
                         cs_lagr_moment_m_data_t   *data_func,
                         void                      *data_input,
                         int                        nt_start,
                         double                     t_start);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time moment associated to particle statistics.
 *
 * This is similar to general time moments (see \ref cs_time_moment.c),
 * with restart, logging, and unsteady reinitialization behavior
 * aligned with other particle statistics.
 *
 * Time moments must be based on values available at the end of each
 * time step, so they cannot be based directly on events (though they
 * can be based on fields defined through \ref cs_lagr_stat_mesh_define,
 * as the matching event-based fields will be updated first).
 *
 * If dimension > 1, the val array is interleaved
 *
 * \param[in]  name           statistics base name
 * \param[in]  location_id    id of associated mesh location
 * \param[in]  stat_type      predefined statistics type, or -1
 * \param[in]  m_type         moment type
 * \param[in]  class_id       particle class id, or 0 for all
 * \param[in]  dim            dimension associated with element data
 * \param[in]  component_id   attribute component id, or < 0 for all
 * \param[in]  data_func      pointer to function to compute statistics
 *                            (if stat_type < 0)
 * \param[in]  data_input     associated input
 * \param[in]  nt_start       starting time step (or -1 to use t_start,
 *                            0 to use idstnt)
 * \param[in]  t_start        starting time
 * \param[in]  restart_mode   behavior in case of restart (reset,
 *                            automatic, or strict)
 *
 * \return id of new moment in case of success, -1 in case of error.
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_time_moment_define(const char                *name,
                                int                        location_id,
                                int                        stat_type,
                                cs_lagr_stat_moment_t      m_type,
                                int                        class_id,
                                int                        dim,
                                int                        component_id,
                                cs_lagr_moment_m_data_t   *data_func,
                                void                      *data_input,
                                int                        nt_start,
                                double                     t_start,
                                cs_lagr_stat_restart_t     restart_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate Lagrangian statistics for a given statistics type.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  stat_type   particle statistics type
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_activate(int  stat_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate time moment for some predefined Lagrangian statistics types.
 *
 * By default, statistics such as mass flows are based on a current time step,
 * and time moments are not computed by default. This function allows forcing
 * the associated moment level so that it is computed also.
 *
 * Note that requesting a higher order moment will automatically include lower
 * order moments, so activating the variance also activates the mean.
 *
 * \param[in]  stat_type   particle statistics type
 * \param[in]  moment      associated time moment level
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_activate_time_moment(int                    stat_type,
                                  cs_lagr_stat_moment_t  moment);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate Lagrangian statistics for a given statistics type.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  stat_type   particle statistics type
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_deactivate(int  stat_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate Lagrangian statistics for a given particle attribute.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  attr_id   particle attribute id
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_activate_attr(int  attr_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate Lagrangian statistics for a given particle attribute.
 *
 * This function is ignored if called after \ref cs_lagr_stat_initialize.
 *
 * \param[in]  attr_id   particle attribute id
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_deactivate_attr(int  attr_id);

/*---------------------------------------------------------------------------*/
/*!
 * \brief Return statistics type associated with a given particle
 *        attribute id.
 *
 * \param[in]  attr_id  particle  attribute id
 *
 * \return  associated  particle statistics type id
 */
/*---------------------------------------------------------------------------*/

int
cs_lagr_stat_type_from_attr_id(int attr_id);

/*---------------------------------------------------------------------------*/
/*!
 * \brief Return attribute id associated with a given statistics type.
 *
 * \param[in]   stat_type     particle statistics type
 *
 * \return attribute id, or -1 if not applicable
 */
/*---------------------------------------------------------------------------*/

int
cs_lagr_stat_type_to_attr_id(int  stat_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine a basic statistic type by its base name.
 *
 * \param[in]  name  particle statistics base name (without class id)
 *
 * \return  matching stat type id, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_stat_type_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map time step values array for Lagrangian statistics.
 *
 * If this function is not called, the field referenced by field pointer
 * CS_F_(dt) will be used instead.
 *
 * \param[in]   dt   pointer to time step values array
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_map_cell_dt(const cs_real_t  *dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian statistics initialization.
 *
 * Statistics activated or deactivated by previous calls to
 * \ref cs_lagr_stat_activate, \ref cs_lagr_stat_deactivate,
 * \ref cs_lagr_stat_activate_attr, and \ref cs_lagr_stat_deactivate_attr
 * will be initialized here.
 *
 * Restart info will be used after to fill in the moments structure
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if a given statistics type has active statistics.
 *
 * \param[in]  group   event group to update
 *
 * \return true if statistics are active for the given group
 */
/*----------------------------------------------------------------------------*/

bool
cs_lagr_stat_is_active( cs_lagr_stat_group_t   group);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read particle statistics restart info if needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_restart_read(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare particle statistics for a given time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_prepare(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update particle statistics for a given time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_update(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update event-based moment accumulators.
 *
 * Partial updates are allowed, so as to balance memory cost for storing
 * events and repetition of mesh-location-based weight updates.
 *
 * \param[in]  events  pointer to event set
 * \param[in]  group   event group to update
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_update_event(cs_lagr_event_set_t   *events,
                          cs_lagr_stat_group_t   group);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all moments management metadata.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log moment definition setup information
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log moment definition information for a given iteration.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_log_iteration(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Checkpoint moment data
 *
 * \param[in]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_stat_restart_write(cs_restart_t  *restart);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return field associated with a given Lagrangian statistic,
 *        given a statistics type (i.e. variable), group (particles or event),
 *        moment order, statistical class, and component id.
 *
 * \param[in]  stat_type     statistics type
 * \param[in]  stat_group    statistics group (particle or event)
 * \param[in]  m_type        moment type (mean or variance)
 * \param[in]  class_id      particle statistical class
 * \param[in]  component_id  component id, or -1 for all
 *
 * \returns pointer to the field associated to the corresponding moment
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_lagr_stat_get_moment(int                    stat_type,
                        cs_lagr_stat_group_t   stat_group,
                        cs_lagr_stat_moment_t  m_type,
                        int                    class_id,
                        int                    component_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return statistical weight
 *
 * \param[in]  class_id    particle statistical class
 *
 * \returns pointer to the field associated to the corresponding weight
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_lagr_stat_get_stat_weight(int  class_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return global volume statistics age
 *
 * \returns age of volume statistics, or -1 if statistics not active yet
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_stat_get_age(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return statistics age for a given moment
 *
 * \param[in]  f  field associated with given statistic
 *
 * \returns age of given statistic, or -1 if not active yet
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_stat_get_moment_age(cs_field_t  *f);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_STAT_H__ */
