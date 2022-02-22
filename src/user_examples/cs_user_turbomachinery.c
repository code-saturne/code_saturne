/*============================================================================
 * Definition of turbomachinery related options.
 *
 * Turbomachinery-related user functions (called in this order):
 *   1) Define rotor cells and associated axis
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_turbomachinery.c
 *
 * \brief Definition of turbomachinery related options.
 *
 * See \ref turbomachinery for examples.
*/
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define rotor/stator model.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery(void)
{
  /*! [user_tbm_set_model] */

  /* Set turbomachinery model type:

     CS_TURBOMACHINERY_NONE,          No turbomachinery modeling
     CS_TURBOMACHINERY_FROZEN,        Frozen rotor model
     CS_TURBOMACHINERY_TRANSIENT      Full transient simulation
  */

  cs_turbomachinery_set_model(CS_TURBOMACHINERY_TRANSIENT);

  /*! [user_tbm_set_model] */
}

/*----------------------------------------------------------------------------
 * Define rotor axes, associated cells, and rotor/stator faces.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_rotor(void)
{
  /* Define rotor axis and cells, with rotor/stator interface face joining */
  /* --------------------------------------------------------------------- */

  /*! [user_tbm_set_rotor] */
  {
    /* Define cells belonging to rotor and associated axis */

    double rotation_velocity = 2.;
    double rotation_axis[3] = {0., 0., 1.};

    double rotation_invariant[3] = {0., 0., 0.};

    const char cell_criteria[]  = "cylinder[0.0, 0.0, 0.0,"
                                          " 0.0, 0.0, 1.0,"
                                          " 2.0]";

    cs_turbomachinery_add_rotor(cell_criteria,
                                rotation_velocity,
                                rotation_axis,
                                rotation_invariant);

  }
  /*! [user_tbm_set_rotor] */


  /*! [user_tbm_set_interface] */
  {
    /* Define joining associated with rotor/stator interface */

    const char faces_criteria[] = "rotor_interface or stator_interface";

    int    verbosity = 0;     /* per-task dump if > 1, debug level if >= 3 */
    int    visualization = 0; /* debug level if >= 3 */
    float  fraction = 0.10, plane = 25.;

    int join_num = cs_turbomachinery_join_add(faces_criteria,
                                              fraction,
                                              plane,
                                              verbosity,
                                              visualization);

    /* Note that advanced parameters may be defined
       using cs_join_set_advanced_param(),
       just as for regular joinings or periodicities. */

  }
  /*! [user_tbm_set_interface] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotation velocity of rotor.
*/
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_set_rotation_velocity(void)
{
  /*! [user_tbm_set_linear_rotation_velocity] */
  {
    /* Linearly increase the rotation velocity from 0 to 1470 rd/min in 0.2 s */
    /* ---------------------------------------------------------------------- */

    int rotor_num = 1;
    double two_pi = 2. * acos(-1.);
    double rotation_velocity = -1470. * two_pi / 60.;
    double rotor_vel = rotation_velocity * CS_MIN(cs_glob_time_step->t_cur / 0.2, 1.);

    cs_turbomachinery_set_rotation_velocity(rotor_num,
                                            rotor_vel);
  }
  /*! [user_tbm_set_linear_rotation_velocity] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
