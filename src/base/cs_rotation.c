/*============================================================================
 * Functions associated to ALE formulation
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_math.h"
#include "cs_time_step.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rotation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*! \struct cs_rotation_t

  \brief Subdomain rotation description.

  Members of this structure are publicly accessible, to allow for
  concise syntax, and for use by inline functions.

  double  omega;
  double  angle;
  double  axis[3];
  double  invariant[3];


  \var  cs_rotation_t::omega
        rotation velocity
  \var  cs_rotation_t::angle
        cumulated rotation
  \var  cs_rotation_t::axis
        rotation vector
  \var  cs_rotation_t::invariant
        coordinates of invariant point
*/

/*! \fn inline static void \
        cs_rotation_velocity(const cs_rotation_t  *r, \
                             const cs_real_t       coords[3], \
                             cs_real_t             vr[3])
 *
 * \brief Compute velocity relative to a fixed frame at a given point.
 *
 * \param[in]   r       pointer to rotation structure
 * \param[in]   coords  point coordinates
 * \param[out]  vr      resulting rotation frame velocity
 */

/*! \fn inline static void \
        cs_rotation_add_coriolis_v(const cs_rotation_t  *r, \
                                   cs_real_t             c, \
                                   const cs_real_t       v[3], \
                                   cs_real_t             vr[3])
 *
 * \brief Add a Coriolis term to a vector.
 *
 * \param[in]       r   pointer to rotation structure
 * \param[in]       c   multiplicative coefficient
 * \param[in]       v   velocity
 * \param[in, out]  vr  resulting Coriolis term
 */

/*! \fn inline static void \
        cs_rotation_coriolis_v(const cs_rotation_t  *r, \
                               cs_real_t             c, \
                               const cs_real_t       v[3], \
                               cs_real_t             vr[3])
 *
 * \brief Compute a vector Coriolis term
 *
 * \param[in]   r   pointer to rotation structure
 * \param[in]   c   multiplicative coefficient
 * \param[in]   v   velocity
 * \param[out]  vr  resulting Coriolis term
 */

/*! \fn inline static void \
        cs_rotation_add_coriolis_t(const cs_rotation_t  *r, \
                                   cs_real_t             c, \
                                   cs_real_t             tr[3][3])
 *
 * \brief Add the dual tensor of a rotation vector to a tensor.
 *
 * \param[in]       r   pointer to rotation structure
 * \param[in]       c   multiplicative coefficient
 * \param[in, out]  tr  tensor to which dual tensor of rotation is added
 */

/*! \fn inline static void \
        cs_rotation_coriolis_t(const cs_rotation_t  *r, \
                               cs_real_t             c, \
                               cs_real_t             tr[3][3])
 *
 * \brief Compute the dual tensor of a rotation vector.
 *
 * \param[in]   r   pointer to rotation structure
 * \param[in]   c   multiplicative coefficient
 * \param[out]  tr  dual tensor of rotation is added
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_rotation_t   _glob_rotation_0[2] = {{0, 0, {0, 0, 0}, {0, 0, 0}},
                                              {0, 0, {0, 0, 0}, {0, 0, 0}}};

/*============================================================================
 * Global variables
 *============================================================================*/

cs_rotation_t  *cs_glob_rotation = _glob_rotation_0;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_rotation_angular_velocity(int         r_num,
                               cs_real_t  *omega);

void
cs_f_rotation_velocity(int              r_num,
                       const cs_real_t  coords[3],
                       cs_real_t        vr[3]);

void
cs_f_rotation_add_coriolis_v(int              r_num,
                             cs_real_t        c,
                             const cs_real_t  v[3],
                             cs_real_t        vr[3]);

void
cs_f_rotation_coriolis_v(int              r_num,
                         cs_real_t        c,
                         const cs_real_t  v[3],
                         cs_real_t        vr[3]);

void
cs_f_rotation_add_coriolis_t(int              r_num,
                             const cs_real_t  c,
                             cs_real_t        tr[3][3]);

void
cs_f_rotation_coriolis_t(int              r_num,
                         const cs_real_t  c,
                         cs_real_t        tr[3][3]);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a vector.
 *
 * parameters:
 *   m[3][4] <-- matrix of the transformation in homogeneous coord.
 *               last line = [0; 0; 0; 1] (Not used here)
 *   c[3]    <-> coordinates
 *----------------------------------------------------------------------------*/

static inline void
_apply_vector_transfo(double     matrix[3][4],
                      cs_real_t  c[3])
{
  int  i, j;

  double  c_a[4] = {c[0], c[1], c[2], 1.}; /* homogeneous coords */
  double  c_b[3] = {0, 0, 0};

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      c_b[i] += matrix[i][j]*c_a[j];

  for (i = 0; i < 3; i++)
    c[i] = c_b[i];
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return angular velocity associated with a rotation.
 *
 * parameters:
 *   r_num  <-- rotation number (1 to n numbering, 0 for none)
 *   omega  --> angular velocity
 *----------------------------------------------------------------------------*/

void
cs_f_rotation_angular_velocity(int         r_num,
                               cs_real_t  *omega)
{
  *omega = (cs_glob_rotation + r_num)->omega;
}

/*----------------------------------------------------------------------------
 * Compute velocity relative to fixed coordinates at a given point.
 *
 * parameters:
 *   r_num  <-- rotation number (1 to n numbering, 0 for none)
 *   coords <-- point coordinates
 *   vr     --> velocity relative to fixed coordinates
 *----------------------------------------------------------------------------*/

void
cs_f_rotation_velocity(int              r_num,
                       const cs_real_t  coords[3],
                       cs_real_t        vr[3])
{
  cs_rotation_velocity(cs_glob_rotation + r_num, coords, vr);
}

/*----------------------------------------------------------------------------
 * Add a Coriolis term to a vector.
 *
 * parameters:
 *   r_num <-- rotation number (1 to n numbering, 0 for none)
 *   c     <-- multiplicative coefficient
 *   v     <-- velocity
 *   vr    <-> vector to which Coriolis term is added
 *---------------------------------------------------------------------------*/

void
cs_f_rotation_add_coriolis_v(int              r_num,
                             cs_real_t        c,
                             const cs_real_t  v[3],
                             cs_real_t        vr[3])
{
  cs_rotation_add_coriolis_v(cs_glob_rotation + r_num, c, v, vr);
}

/*----------------------------------------------------------------------------
 * Compute a vector Coriolis term.
 *
 * parameters:
 *   r_num <-- rotation number (1 to n numbering, 0 for none)
 *   c     <-- multiplicative coefficient
 *   v     <-- velocity
 *   vr    --> vector associted to Coriolis term
 *---------------------------------------------------------------------------*/

void
cs_f_rotation_coriolis_v(int              r_num,
                         cs_real_t        c,
                         const cs_real_t  v[3],
                         cs_real_t        vr[3])
{
  cs_rotation_coriolis_v(cs_glob_rotation + r_num, c, v, vr);
}

/*----------------------------------------------------------------------------
 * Add the dual tensor of a rotation vector to a tensor
 *
 * parameters:
 *   r_num <-- rotation number (1 to n numbering, 0 for none)
 *   c     <-- multiplicative coefficient
 *   tr    <-> tensor to which dual tensor of rotation is added
 *---------------------------------------------------------------------------*/

void
cs_f_rotation_add_coriolis_t(int              r_num,
                             const cs_real_t  c,
                             cs_real_t        tr[3][3])
{
  cs_rotation_add_coriolis_t(cs_glob_rotation + r_num, c, tr);
}

/*----------------------------------------------------------------------------
 * Compute the dual tensor of a rotation vector
 *
 * parameters:
 *   r_num <-- rotation number (1 to n numbering, 0 for none)
 *   c     <-- multiplicative coefficient
 *   tr    --> dual tensor of rotation
 *---------------------------------------------------------------------------*/

void
cs_f_rotation_coriolis_t(int              r_num,
                         const cs_real_t  c,
                         cs_real_t        tr[3][3])
{
  cs_rotation_coriolis_t(cs_glob_rotation + r_num, c, tr);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a global rotation.
 *
 * The rotation vector's length determines the angular velocity (in rad/s).
 *
 * \param[in]  omega_x      rotation vector x component
 * \param[in]  omega_y      rotation vector y component
 * \param[in]  omega_z      rotation vector z component
 * \param[in]  invariant_x  invariant point x component
 * \param[in]  invariant_y  invariant point y component
 * \param[in]  invariant_z  invariant point z component
 */
/*----------------------------------------------------------------------------*/

void
cs_rotation_define(double  omega_x,
                   double  omega_y,
                   double  omega_z,
                   double  invariant_x,
                   double  invariant_y,
                   double  invariant_z)
{
  cs_rotation_t  *r = _glob_rotation_0;

  r->axis[0] = omega_x;
  r->axis[1] = omega_y;
  r->axis[2] = omega_z;
  r->invariant[0] = invariant_x;
  r->invariant[1] = invariant_y;
  r->invariant[2] = invariant_z;

  r->omega = sqrt(cs_math_3_square_norm(r->axis));

  r->angle = 0;

  for (int i = 0; i < 3; i++) {
    r->axis[i] /= r->omega;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rotation matrix
 *
 * \param[in]   theta            rotation angle, in radians
 * \param[in]   axis             rotation axis direction vector
 * \param[in]   invariant_point  invariant point coordinates
 * \param[out]  matrix           resulting rotation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_rotation_matrix(double        theta,
                   const double  axis[3],
                   const double  invariant_point[3],
                   double        matrix[3][4])
{
  int  i, j;
  double norm;
  double direction[3];
  double rot[3][3];

  if (fabs(theta) > 0) {

    const double cost = cos(theta);
    const double sint = sin(theta);
    const double onemcost = (1.0 - cost);

    /* Compute the rotation matrix, using formula:
     *  R = (1-cos(theta))axis.transp(axis) + cos(theta)I + sin(theta)V
     *
     *           [ 0            -direction(3)  direction(2)]
     *  with V = [ direction(3)       0       -direction(1)]
     *           [-direction(2)  direction(1)       0      ]
     */

    norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

    direction[0] = axis[0] / norm;
    direction[1] = axis[1] / norm;
    direction[2] = axis[2] / norm;

    /* first row of rotation maxtrix */
    rot[0][0] = onemcost*direction[0]*direction[0] + cost;
    rot[0][1] = onemcost*direction[0]*direction[1] - sint*direction[2];
    rot[0][2] = onemcost*direction[0]*direction[2] + sint*direction[1];

    /* second row of rotation maxtrix */
    rot[1][0] = onemcost*direction[1]*direction[0] + sint*direction[2];
    rot[1][1] = onemcost*direction[1]*direction[1] + cost;
    rot[1][2] = onemcost*direction[1]*direction[2] - sint*direction[0];

    /* third row of rotation maxtrix */
    rot[2][0] = onemcost*direction[2]*direction[0] - sint*direction[1];
    rot[2][1] = onemcost*direction[2]*direction[1] + sint*direction[0];
    rot[2][2] = onemcost*direction[2]*direction[2] + cost;

    /* Now compute full rotation matrix in homogeneous coordinates,
     * accounting for invariant point of coordiantes t[], with the formula:
     *
     *     [1 0 0 t[0]] [r[0][0] r[0][1] r[0][3] 0] [1 0 0 -t[0]]
     * M = [0 1 0 t[1]].[r[1][0] r[1][1] r[1][3] 0].[0 1 0 -t[1]]
     *     [0 0 1 t[2]] [r[2][0] r[2][1] r[2][3] 0] [0 0 1 -t[2]]
     *     [0 0 0 1   ] [0       0       0       1] [0 0 0  1]
     */

    for (i = 0; i < 3; i++) {       /* rotation part of matrix */
      for (j = 0; j < 3; j++) {
        matrix[i][j] = rot[i][j];
      }
    }

    for (i = 0; i < 3; i++) {
      matrix[i][3] = invariant_point[i];
      for (j = 0; j < 3; j++)
        matrix[i][3] -= rot[i][j]*invariant_point[j];
    }

  }

  /* Zero rotation angle case */

  else {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 4; j++)
        matrix[i][j] = 0;
      matrix[i][i] = 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update coordinates based on a global rotation and time.
 *
 * \param[in]       n_coords  number of coordinates
 * \param[in]       t_rot     time since rotation start
 * \param[in, out]  coords    coordinates array
 */
/*----------------------------------------------------------------------------*/

void
cs_rotation_update_coords(cs_lnum_t    n_coords,
                          double       t_rot,
                          cs_real_3_t  coords[])
{
  assert(cs_glob_rotation == _glob_rotation_0);

  double   matrix[3][4];

  cs_rotation_matrix((cs_glob_rotation + 1)->omega * t_rot,
                     (cs_glob_rotation + 1)->axis,
                     (cs_glob_rotation + 1)->invariant,
                     matrix);

# pragma omp parallel for
  for (cs_lnum_t i = 0; i < n_coords; i++)
    _apply_vector_transfo(matrix, coords[i]);
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
                  cs_real_t             vc[3])
{
  /* Axial unit vector */

  const cs_real_t *e_ax = r->axis;

  /* Tangential unit vector */

  cs_real_t e_th[3];

  e_th[0] =  e_ax[1] * (coords[2] - r->invariant[2])
           - e_ax[2] * (coords[1] - r->invariant[1]);
  e_th[1] =  e_ax[2] * (coords[0] - r->invariant[0])
           - e_ax[0] * (coords[2] - r->invariant[2]);
  e_th[2] =  e_ax[0] * (coords[1] - r->invariant[1])
           - e_ax[1] * (coords[0] - r->invariant[0]);

  cs_real_t xnrm =  sqrt(cs_math_3_square_norm(e_th));

  e_th[0] /= xnrm;
  e_th[1] /= xnrm;
  e_th[2] /= xnrm;

  /* Radial unit vector */

  cs_real_t e_r[3];

  e_r[0] = - e_ax[1]*e_th[2] + e_ax[2]*e_th[1];
  e_r[1] = - e_ax[2]*e_th[0] + e_ax[0]*e_th[2];
  e_r[2] = - e_ax[0]*e_th[1] + e_ax[1]*e_th[0];

  /* Transformation into cylindrical coordinates */

  vc[0] = v[0]*e_r[0]  + v[1]*e_r[1]  + v[2]*e_r[2];
  vc[1] = v[0]*e_th[0] + v[1]*e_th[1] + v[2]*e_th[2];
  vc[2] = v[0]*e_ax[0] + v[1]*e_ax[1] + v[2]*e_ax[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy rotation structure values to an array
 *
 * This may be useful to avoid requiring specific type mappings for MPI or
 * other programming languages.
 *
 * \param[in]  r_num  rotation number (1 to n numbering, 0 for none)
 * \param[in]  fra    flat rotation array: axis (0-2), invariant(3-5),
 *                    omega (6), angle(7)
 */
/*----------------------------------------------------------------------------*/

void
cs_rotation_to_array(int        r_num,
                     cs_real_t  fra[8])
{
  const cs_rotation_t *r = cs_glob_rotation + r_num;

  fra[0] = r->axis[0];
  fra[1] = r->axis[1];
  fra[2] = r->axis[2];
  fra[3] = r->invariant[0];
  fra[4] = r->invariant[1];
  fra[5] = r->invariant[2];
  fra[6] = r->omega;
  fra[7] = r->angle;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
