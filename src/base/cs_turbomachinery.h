#ifndef __CS_TURBOMACHINERY_H__
#define __CS_TURBOMACHINERY_H__

/*============================================================================
 * Turbomachinery modeling features.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Turbomachinery model type */

typedef enum {

  CS_TURBOMACHINERY_NONE,          /* No turbomachinery modeling */
  CS_TURBOMACHINERY_FROZEN,        /* Frozen rotor model */
  CS_TURBOMACHINERY_TRANSIENT      /* full transient simulation */

} cs_turbomachinery_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define rotor/stator model.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_model(cs_turbomachinery_model_t  model);

/*----------------------------------------------------------------------------
 * return rotor/stator model.
 *----------------------------------------------------------------------------*/

cs_turbomachinery_model_t
cs_turbomachinery_get_model(void);

/*----------------------------------------------------------------------------
 * Define a rotor by its axis and cell selection criteria.
 *
 * parameters:
 *   cell_criteria         <-- cell selection criteria string
 *   rotation_velocity     <-- rotation velocity, in radians/second
 *   rotation_axis         <-- rotation axis vector
 *   rotation_invariant    <-- rotation invariant point
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_add_rotor(const char    *cell_criteria,
                            double         rotation_velocity,
                            const double   rotation_axis[3],
                            const double   rotation_invariant[3]);

/*----------------------------------------------------------------------------
 * Add a cs_join_t structure to the list of rotor/stator joinings.
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *
 * returns:
 *   number (1 to n) associated with new joining
 *----------------------------------------------------------------------------*/

int
cs_turbomachinery_join_add(const char  *sel_criteria,
                           float        fraction,
                           float        plane,
                           int          verbosity,
                           int          visualization);

/*----------------------------------------------------------------------------
 * Add a cs_join_t structure to the list of rotor/stator couplings.
 *
 * parameters:
 *   sel_criteria <-- boundary face selection criteria
 *   tolerance    <-- value of the search tolerance
 *   verbosity    <-- level of verbosity required
 *
 * returns:
 *   number (1 to n) associated with new coupling
 *----------------------------------------------------------------------------*/

int
cs_turbomachinery_coupling_add(const char  *sel_criteria,
                               float        tolerance,
                               int          verbosity);

/*----------------------------------------------------------------------------
 * Definitions for turbomachinery computation.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_define(void);

/*----------------------------------------------------------------------------
 * Initializations for turbomachinery computation
 *
 * Note: this function should be called after the mesh is built,
 *       but before cs_post_init_meshes() so that postprocessing meshes are
 *       updated correctly in the transient case.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_initialize(void);

/*----------------------------------------------------------------------------
 * Free turbomachinery info
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_finalize(void);

/*----------------------------------------------------------------------------
 * Update mesh for unsteady rotor/stator computation
 *
 * parameters:
 *   t_cur_mob    <-- current rotor time
 *   t_elapsed    --> elapsed computation time
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_update_mesh(double   t_cur_mob,
                              double  *t_elapsed);

/*----------------------------------------------------------------------------
 * Update mesh for unsteady rotor/stator computation in case of restart.
 *
 * Reads mesh from checkpoint when available.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_mesh(void);

/*----------------------------------------------------------------------------
 * Reinitialize interior face-based fields.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_reinit_i_face_fields(void);

/*----------------------------------------------------------------------------
 * Resize cell-based fields.
 *
 * This function only handles fields owning their values.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_resize_cell_fields(void);

/*----------------------------------------------------------------------------
 * Compute rotation matrix
 *
 * parameters:
 *   rotor_num <-- rotor number (1 to n numbering)
 *   theta     <-- rotation angle, in radians
 *   matrix    --> resulting rotation matrix
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_rotation_matrix(int        rotor_num,
                                  double     theta,
                                  cs_real_t  matrix[3][4]);

/*----------------------------------------------------------------------------
 * Return number of rotors.
 *
 * Note that the number of associated rotations is n_rotors + 1, as the
 * first rotation id is reserved for the fixed portion of the domain.
 *
 * return:
 *   number of rotors
 *----------------------------------------------------------------------------*/

int
cs_turbomachinery_n_rotors(void);

/*----------------------------------------------------------------------------
 * Return cell rotor number.
 *
 * Each cell may be associated with a given rotor, or rotation, with 0
 * indicating that that cell does not rotate.
 *
 * returns:
 *   array defining rotor number associated with each cell
 *   (0 for none, 1 to n otherwise)
 *----------------------------------------------------------------------------*/

const int *
cs_turbomachinery_get_cell_rotor_num(void);

/*----------------------------------------------------------------------------
 * Return rotation velocity
 *
 * parameters:
 *   rotor_num <-- rotor number (1 to n numbering)
 *----------------------------------------------------------------------------*/

double
cs_turbomachinery_get_rotation_velocity(int  rotor_num);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set rotation velocity
 *
 * param[in]  rotor_num  rotor number (1 to n numbering)
 * param[in]  omega      rotation velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_rotation_velocity(int     rotor_num,
                                        double  omega);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build rotation matrices for a given time interval.
 *
 * The caller is responsible for freeing the array when not needed.
 *
 * \param[in]  dt  associated time delta (0 for current, unmodified time)
 *
 * \return  array of rotation matrices.
 */
/*----------------------------------------------------------------------------*/

cs_real_34_t *
cs_turbomachinery_get_rotation_matrices(double dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set turbomachinery joining retry parameters.
 *
 * When a joing leads to a different number of boundary faces from the
 * previous position, the rotor positions may be perturbed by a small
 * quantity to try to obtain a better joining.
 *
 * param[in]  n_max_join_retries   maximum number of retries before considering
 *                                 the joining has failed
 * param[in]  dt_retry_multiplier  time step multiplier for new position retry
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_rotation_retry(int     n_max_join_retries,
                                     double  dt_retry_multiplier);

/*----------------------------------------------------------------------------
 * Rotation of vector and tensor fields.
 *
 * parameters:
 *   dt <-- cell time step values
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_rotate_fields(const cs_real_t dt[]);

/*----------------------------------------------------------------------------
 * Compute velocity relative to fixed coordinates at a given point
 *
 * Deprecated:
 * Use cs_rotation_velocity for more consistent naming of this reference
 * frame velocity.
 *
 * parameters:
 *   rotor_num <-- associated rotor number (1 to n numbering)
 *   coords    <-- point coordinates
 *   velocity  --> velocity relative to fixed coordinates
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_relative_velocity(int              rotor_num,
                                    const cs_real_t  coords[3],
                                    cs_real_t        velocity[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read turbomachinery metadata from restart file.
 *
 * The mesh is handled separately.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_read(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write turbomachinery metadata to checkpoint file.
 *
 * The mesh is handled separately.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_write(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBOMACHINERY_H__ */

