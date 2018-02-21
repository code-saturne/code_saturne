#ifndef __CS_ROTATION_H__
#define __CS_ROTATION_H__

/*============================================================================
 * Rotation modeling features.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Rotation structure */

typedef struct {

  double  omega;
  double  angle;
  double  axis[3];
  double  invariant[3];

} cs_rotation_t;

/*============================================================================
 * Global variables
 *============================================================================*/

extern cs_rotation_t  *cs_glob_rotation;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a global rotation.
 *
 * The rotation vector's length determines the angular velocity (in rad/s).
 *
 * parameters:
 *   omega_x     <-- rotation vector x component
 *   omega_y     <-- rotation vector y component
 *   omega_z     <-- rotation vector z component
 *   invariant_x <-- invariant point x component
 *   invariant_y <-- invariant point y component
 *   invariant_z <-- invariant point z component
 *----------------------------------------------------------------------------*/

void
cs_rotation_define(double  omega_x,
                   double  omega_y,
                   double  omega_z,
                   double  invariant_x,
                   double  invariant_y,
                   double  invariant_z);

/*----------------------------------------------------------------------------
 * Compute rotation matrix
 *
 * parameters:
 *   theta           <-- rotation angle, in radians
 *   axis            <-- rotation axis direction vector
 *   invariant_point <-- invariant point coordinates
 *   matrix          -> resulting rotation matrix
 *---------------------------------------------------------------------------*/

void
cs_rotation_matrix(double        theta,
                   const double  axis[3],
                   const double  invariant_point[3],
                   double        matrix[3][4]);

/*----------------------------------------------------------------------------
 * Update coordinates based on a global rotation and time.
 *
 * parameters:
 *   n_coords <-- number of coordinates
 *   t_rot    <-- time since rotation start
 *   coords   <-> coordinates array
 *----------------------------------------------------------------------------*/

void
cs_rotation_update_coords(cs_lnum_t    n_coords,
                          double       t_rot,
                          cs_real_3_t  coords[]);

/*----------------------------------------------------------------------------
 * Compute rotation velocity relative to fixed coordinates at a given point.
 *
 * parameters:
 *   r      <-- pointer to rotation structure
 *   coords <-- coordinates at point
 *   vr     --> relative velocity
 *---------------------------------------------------------------------------*/

static inline void
cs_rotation_velocity(const cs_rotation_t  *r,
                     const cs_real_t       coords[3],
                     cs_real_t             vr[3])
{
  vr[0] = (- r->axis[2] * (coords[1] - r->invariant[1])
           + r->axis[1] * (coords[2] - r->invariant[2])) * r->omega;
  vr[1] = (  r->axis[2] * (coords[0] - r->invariant[0])
           - r->axis[0] * (coords[2] - r->invariant[2])) * r->omega;
  vr[2] = (- r->axis[1] * (coords[0] - r->invariant[0])
           + r->axis[0] * (coords[1] - r->invariant[1])) * r->omega;
}

/*----------------------------------------------------------------------------
 * Add a Coriolis term to a vector.
 *
 * parameters:
 *   r  <-- pointer to rotation structure
 *   c  <-- multiplicative coefficient
 *   v  <-- velocity
 *   vr <-> vector to which Coriolis term is added
 *---------------------------------------------------------------------------*/

static inline void
cs_rotation_add_coriolis_v(const cs_rotation_t  *r,
                           cs_real_t             c,
                           const cs_real_t       v[3],
                           cs_real_t             vr[3])
{
  double f = r->omega * c;

  vr[0] += (- r->axis[2]*v[1] + r->axis[1]*v[2]) * f;
  vr[1] += (- r->axis[0]*v[2] + r->axis[2]*v[0]) * f;
  vr[2] += (- r->axis[1]*v[0] + r->axis[0]*v[1]) * f;
}

/*----------------------------------------------------------------------------
 * Compute a vector Coriolis term.
 *
 * parameters:
 *   r  <-- pointer to rotation structure
 *   c  <-- multiplicative coefficient
 *   v  <-- velocity
 *   vr --> vector associted to Coriolis term
 *---------------------------------------------------------------------------*/

static inline void
cs_rotation_coriolis_v(const cs_rotation_t  *r,
                       cs_real_t             c,
                       const cs_real_t       v[3],
                       cs_real_t             vr[3])
{
  double f = r->omega * c;

  vr[0] = (- r->axis[2]*v[1] + r->axis[1]*v[2]) * f;
  vr[1] = (- r->axis[0]*v[2] + r->axis[2]*v[0]) * f;
  vr[2] = (- r->axis[1]*v[0] + r->axis[0]*v[1]) * f;
}

/*----------------------------------------------------------------------------
 * Add the dual tensor of a rotation vector to a tensor
 * The dual tensor is such that:
 *  tr[i][j] * v[j] = (omage ^ v)_(ij)
 *
 * parameters:
 *   r  <-- pointer to rotation structure
 *   c  <-- multiplicative coefficient
 *   tr <-> tensor to which dual tensor of rotation is added
 *---------------------------------------------------------------------------*/

static inline void
cs_rotation_add_coriolis_t(const cs_rotation_t  *r,
                           cs_real_t             c,
                           cs_real_t             tr[3][3])
{
  double f = r->omega * c;

  tr[0][1] -= r->axis[2]*f;
  tr[0][2] += r->axis[1]*f;

  tr[1][0] += r->axis[2]*f;
  tr[1][2] -= r->axis[0]*f;

  tr[2][0] -= r->axis[1]*f;
  tr[2][1] += r->axis[0]*f;
}

/*----------------------------------------------------------------------------
 * Compute the dual tensor of a rotation vector
 * The dual tensor is such that:
 *  tr[i][j] * v[j] = (omage ^ v)_(ij)
 *
 * parameters:
 *   r  <-- pointer to rotation structure
 *   c  <-- multiplicative coefficient
 *   tr --> dual tensor of rotation
 *---------------------------------------------------------------------------*/

static inline void
cs_rotation_coriolis_t(const cs_rotation_t  *r,
                       cs_real_t             c,
                       cs_real_t             tr[3][3])
{
  double f = r->omega * c;

  tr[0][0] =   0.;
  tr[0][1] = - r->axis[2]*f;
  tr[0][2] =   r->axis[1]*f;

  tr[1][0] =   r->axis[2]*f;
  tr[1][1] =   0.;
  tr[1][2] = - r->axis[0]*f;

  tr[2][0] = - r->axis[1]*f;
  tr[2][1] =   r->axis[0]*f;
  tr[2][2] =   0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Express a vector in the cyclindrical system associated to a rotation.
 *
 * \param[in]   r   pointer to rotation structure
 * \param[in]   p   cartesian coordinates of the location point
 * \param[in]   v   vector components in cartesian coordinates system
 * \param[out]  vc  vector components in cylindrical coordinates system
 */
/*----------------------------------------------------------------------------*/

void
cs_rotation_cyl_v(const cs_rotation_t  *r,
                  const cs_real_t       coords[3],
                  const cs_real_t       v[3],
                  cs_real_t             vc[3]);

/*----------------------------------------------------------------------------
 * Copy rotation structure values to an array
 *
 * This may be useful to avoid requiring specific type mappings for MPI or
 * other programming languages.
 *
 * parameters:
 *   r_num <-- rotation number (1 to n numbering, 0 for none)
 *   fra   --> flat rotation array: axis (0-2), invariant(3-5),
 *             omega (6), angle(7)
 *---------------------------------------------------------------------------*/

void
cs_rotation_to_array(int        r_num,
                     cs_real_t  fra[8]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ROTATION_H__ */
