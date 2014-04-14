#ifndef __CS_TURBOMACHINERY_H__
#define __CS_TURBOMACHINERY_H__

/*============================================================================
 * Turbomachinery modeling features.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 *   cell_criteria     <-- cell selection criteria string
 *   rotation_velocity <-- rotation velocity, in radians/second
 *   rotation_axis     <-- rotation axis vector
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_add_rotor(const char    *cell_criteria,
                            double         rotation_velocity,
                            const double   rotation_axis[3]);

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
 * Initializations for turbomachinery computation
 *
 * Note: this function should be called before once the mesh is built,
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
 *   t_cur_mob <-- current rotor time
 *   t_elapsed --> elapsed computation time
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_update_mesh(double   t_cur_mob,
                              double  *t_elapsed);

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
 * Return cell rotation axis number
 *----------------------------------------------------------------------------*/

int *
cs_turbomachinery_get_rotation_axis_number(void);

/*----------------------------------------------------------------------------
 * Return rotation velocity
 *----------------------------------------------------------------------------*/

double
cs_turbomachinery_get_rotation_velocity(void);

/*----------------------------------------------------------------------------
 * Rotation of vector and tensor fields.
 *----------------------------------------------------------------------------*/

void
cs_turbomachinery_rotate_fields(const cs_real_t dt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBOMACHINERY_H__ */

