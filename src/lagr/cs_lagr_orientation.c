/*============================================================================
 * Methods for particle deposition
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

/*============================================================================
 * Functions dealing with particle deposition
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_random.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_event.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_orientation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function prototypes (definitions follow)
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid stretching by symetric part of the mean
 *        velocity gradient
 *
 * Principle: exact exponential scheme
 *
 * \param[in]  dtp           time step
 * \param[in]  orientation   spheroid orientation (vector)
 * \param[in]  Lambda        spheroid shape parameter
 * \param[in]  beta          coefficient to be added to diagonal elements (if needed)
 * \param[in]  gradvf_sym    symetric part of the fluid velocity gradient (tensor)
 */
/*----------------------------------------------------------------------------*/

static void
_mean_stretching_phase_spheroid(const cs_real_t dtp,
                                      cs_real_t orientation[3],
                                const cs_real_t Lambda,
                                const cs_real_t beta,
                                const cs_real_t gradvf_sym[3][3])
{
  /* Get the eigenvalues and eigenvectors of gradvf_sym */
  cs_real_3_t eig_val;
  // Initialize eig_vec to identity
  cs_real_33_t eig_vec = {{1.0, 0.0, 0.0},
                          {0.0, 1.0, 0.0},
                          {0.0, 0.0, 1.0}};
  cs_real_t tol_err = 1.0e-12;
  cs_real_33_t mat_loc ;
  mat_loc[0][0] = Lambda * (gradvf_sym[0][0]) + beta;
  mat_loc[0][1] = Lambda *  gradvf_sym[0][1] ;
  mat_loc[0][2] = Lambda *  gradvf_sym[0][2] ;
  mat_loc[1][0] = Lambda *  gradvf_sym[1][0] ;
  mat_loc[1][1] = Lambda * (gradvf_sym[1][1]) + beta;
  mat_loc[1][2] = Lambda *  gradvf_sym[1][2] ;
  mat_loc[2][0] = Lambda *  gradvf_sym[2][0] ;
  mat_loc[2][1] = Lambda *  gradvf_sym[2][1] ;
  mat_loc[2][2] = Lambda * (gradvf_sym[2][2]) + beta;
  cs_math_33_eig_val_vec(mat_loc, tol_err, eig_val, eig_vec);

  /* Get the orientation in the local frame of reference */
  cs_real_3_t orient_new;
  cs_math_33t_3_product(eig_vec, orientation, orient_new);

  /* Update the orientation*/
  for (int id = 0; id < 3; id++) {
    orient_new[id] = orient_new[id] * exp( eig_val[id] * dtp);
  }

  /* Get the orient_new on the global frame of reference */
  cs_math_33_3_product(eig_vec, orient_new, orientation);

  /* Renormalize the orientation (for security) */
  cs_math_3_normalise(orientation, orientation);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid rotation increment by antisymetric mean
 *        velocity gradient
 *
 * Principle: exact matrix rotation (Rodrigues' formula)
 *
 * TODO: use directly the antisymetric mean velocity gradient matrix Omega
 *
 * \param[in]  dtp           time step
 * \param[in]  orientation   spheroid orientation (vector)
 * \param[in]  gradvf_ant    antisymetric part of the fluid velocity gradient (tensor)
 */
/*----------------------------------------------------------------------------*/
// TODO: use directly the antisymetric mean velocity gradient matrix Omega

static void
_mean_rotation_phase_spheroid(const cs_real_t  dtp,
                                    cs_real_t  orientation[3],
                              const cs_real_t  gradvf_ant[3][3])
{
  /* Calculate the vector n and its norm */
  cs_real_t n_norm
    = sqrt(  cs_math_pow2( gradvf_ant[0][1] )
           + cs_math_pow2( gradvf_ant[0][2] )
           + cs_math_pow2( gradvf_ant[1][2] ) );

  cs_real_t n[3];
  n[0] = gradvf_ant[1][2] / n_norm;
  n[1] = gradvf_ant[0][2] / n_norm;
  n[2] = gradvf_ant[0][1] / n_norm;

  /* Calculate the rotation of the spheroid */
  cs_real_t orientation_new[3];
  cs_real_t t_n = dtp * n_norm;

  // TODO: use multiplication of matrix and vector (to simplify)
  orientation_new[0]
    =    ( cos(t_n) + n[0]*n[0]*(1-cos(t_n))      ) * orientation[0]
       + ( n[0]*n[1]*(1-cos(t_n)) - n[2]*sin(t_n) ) * orientation[1]
       + ( n[0]*n[2]*(1-cos(t_n)) + n[1]*sin(t_n) ) * orientation[2];
  orientation_new[1]
    =    ( n[0]*n[1]*(1-cos(t_n)) + n[2]*sin(t_n) ) * orientation[0]
       + ( cos(t_n) + n[1]*n[1]*(1-cos(t_n))      ) * orientation[1]
       + ( n[1]*n[2]*(1-cos(t_n)) - n[0]*sin(t_n) ) * orientation[2];
  orientation_new[2]
    =    ( n[0]*n[2]*(1-cos(t_n)) - n[1]*sin(t_n) ) * orientation[0]
       + ( n[1]*n[2]*(1-cos(t_n)) + n[0]*sin(t_n) ) * orientation[1]
       + ( cos(t_n) + n[2]*n[2]*(1-cos(t_n))      ) * orientation[2];

  /* Update values */
  for (int i = 0; i < 3; i++)
    orientation[i] = orientation_new[i];

  /* Renormalize the orientation (for security) */
  cs_math_3_normalise(orientation, orientation);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid stretching by Brownian effect
 *
 * Principle: Euler scheme
 *
 * \param[in]  ip           particle number
 * \param[in]  dtp          time step
 * \param[in]  d_plus      reduced parameter for the (isotropic) fluctuation tensor
 * \param[in]  orientation  spheroid orientation (vector)
 * \param[in]  Lambda       spheroid shape parameter
 * \param[in]  brown        normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_stretching_phase_spheroid(const cs_lnum_t ip,
                              const cs_real_t dtp,
                              const cs_real_t d_plus,
                                    cs_real_t orientation[3],
                              const cs_real_t Lambda,
                              const cs_real_t brown[])
{
  /* Auxiliary parameters for calculation */
  cs_real_t aux1 = orientation[0] * orientation[1] * orientation[2];
  cs_real_t r1_w2 = cs_math_pow2(orientation[0]);
  cs_real_t r2_w2 = cs_math_pow2(orientation[1]);
  cs_real_t r3_w2 = cs_math_pow2(orientation[2]);
  cs_real_t aux4 = Lambda * d_plus;
  cs_real_t aux5 = 1. / (1. + 0.5 * cs_math_pow2(aux4) * dtp);

  /* W11 W12 W13 W21 W22 W23 W31 W32 W33 */
  cs_real_t  dw11 = sqrt(dtp) * brown[ip*9+0];
  cs_real_t  dw12 = sqrt(dtp) * brown[ip*9+1];
  cs_real_t  dw13 = sqrt(dtp) * brown[ip*9+2];
  cs_real_t  dw21 = sqrt(dtp) * brown[ip*9+3];
  cs_real_t  dw22 = sqrt(dtp) * brown[ip*9+4];
  cs_real_t  dw23 = sqrt(dtp) * brown[ip*9+5];
  cs_real_t  dw31 = sqrt(dtp) * brown[ip*9+6];
  cs_real_t  dw32 = sqrt(dtp) * brown[ip*9+7];
  cs_real_t  dw33 = sqrt(dtp) * brown[ip*9+8];

  cs_real_t orientation_new[3];

  /* solve the orientation dynamics produced by Brownian stretching */
  orientation_new[0]
    = aux5 * (orientation[0] + aux4
              * (- aux1 * (dw23+dw32)
                 + orientation[0] * ((1. - r1_w2)*dw11 - r2_w2*dw22 - r3_w2*dw33)
                 + 0.5 * (1.-2.*r1_w2)*(orientation[1]
                                        *(dw12+dw21)+orientation[2]*(dw13+dw31))));

  orientation_new[1]
    = aux5 * (orientation[1] + aux4
              * (- aux1 * (dw13+dw31)
                 + orientation[1] * ((1. - r2_w2)*dw22 - r1_w2*dw11 - r3_w2*dw33)
                 + 0.5 * (1.-2.*r2_w2)*(orientation[0]
                                        *(dw12+dw21)+orientation[2]*(dw23+dw32))));

  orientation_new[2]
    = aux5 * (orientation[2] + aux4
              * (- aux1 * (dw12+dw21)
                 + orientation[2] * ((1. - r3_w2)*dw33 - r1_w2*dw11 - r2_w2*dw22)
                 + 0.5 * (1.-2.*r3_w2)*(orientation[1]
                                        *(dw23+dw32)+orientation[0]*(dw13+dw31))));

  /* Renormalise for security */
  cs_math_3_normalise(orientation_new, orientation);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid rotation by Brownian on the sphere
 *        using quaternions
 *
 * Principle:
 *            Use of quaternions to get the rotation angle
 *
 * \param[in]  ip           particle number
 * \param[in]  dtp          time step
 * \param[in]  d_minus      reduced parameter for the (isotropic)
 *                          fluctuation tensor
 * \param[in]  orientation  spheroid orientation (vector)
 * \param[in]  brown        normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_rotation_phase_spheroid_by_quaternion(const cs_lnum_t  ip,
                                          const cs_real_t  dtp,
                                          const cs_real_t  d_minus,
                                          cs_real_t        orient_loc[3],
                                          cs_real_t        quater_loc[4],
                                          const cs_real_t  brown[])
{
  /* Define the nine Brownian increments W11 W12 W13 W21 W22 W23 W31 W32 W33 */
  cs_real_t  dw12 = sqrt(dtp) * brown[ip*9+1];
  cs_real_t  dw13 = sqrt(dtp) * brown[ip*9+2];
  cs_real_t  dw21 = sqrt(dtp) * brown[ip*9+3];
  cs_real_t  dw23 = sqrt(dtp) * brown[ip*9+5];
  cs_real_t  dw31 = sqrt(dtp) * brown[ip*9+6];
  cs_real_t  dw32 = sqrt(dtp) * brown[ip*9+7];

  /* Define the Brownian rotation dw */
  cs_real_t dw[3];
  dw[0] =  dw23-dw32;
  dw[1] =  dw31-dw13;
  dw[2] =  dw12-dw21;

  /* Increment quaternions */
  cs_real_t quater_fix[4] = {1.0, 0.0, 0.0, 0.0};
  cs_real_t aux1 = 1. / ( 1.0 + 0.09375*cs_math_pow2(d_minus) * dtp );
  cs_real_t quater_new[4]
    = {aux1 * (quater_fix[0] + 0.25 * d_minus * (- dw[0]*quater_fix[1]
                                                 - dw[1]*quater_fix[2]
                                                 - dw[2]*quater_fix[3])),
       aux1 * (quater_fix[1] + 0.25 * d_minus * (  dw[0]*quater_fix[0]
                                                 - dw[1]*quater_fix[3]
                                                 + dw[2]*quater_fix[2])),
       aux1 * (quater_fix[2] + 0.25 * d_minus * (  dw[0]*quater_fix[3]
                                                 + dw[1]*quater_fix[0]
                                                 - dw[2]*quater_fix[1])),
       aux1 * (quater_fix[3] + 0.25 * d_minus * ( -dw[0]*quater_fix[2]
                                                 + dw[1]*quater_fix[1]
                                                 + dw[2]*quater_fix[0]))};

  cs_real_t quater_pow2[4];
  quater_pow2[0] = cs_math_pow2(quater_new[0]);
  quater_pow2[1] = cs_math_pow2(quater_new[1]);
  quater_pow2[2] = cs_math_pow2(quater_new[2]);
  quater_pow2[3] = cs_math_pow2(quater_new[3]);

  cs_real_t quater_norm =   quater_pow2[0] + quater_pow2[1]
                          + quater_pow2[2] + quater_pow2[3];
  for (int dim = 0; dim < 4; dim++) {
    quater_loc[dim] = quater_new[dim] / sqrt(quater_norm);
    quater_pow2[dim] = quater_pow2[dim] / quater_norm;
  }

  /* Rotation matrix */
  cs_real_33_t rot_m = {
    {quater_pow2[0] + quater_pow2[1] - quater_pow2[2] - quater_pow2[3],    // [0][0]
     2.0 * (quater_loc[1]*quater_loc[2] - quater_loc[0]*quater_loc[3]),    // [0][1]
     2.0 * (quater_loc[1]*quater_loc[3] + quater_loc[0]*quater_loc[2])} ,  // [0][2]

    {2.0 * (quater_loc[1]*quater_loc[2] + quater_loc[0]*quater_loc[3]),    // [1][0]
     quater_pow2[0] - quater_pow2[1] + quater_pow2[2] - quater_pow2[3],    // [1][1]
     2.0 * (quater_loc[2]*quater_loc[3] - quater_loc[0]*quater_loc[1])} ,  // [1][2]

    {2.0 * (quater_loc[1]*quater_loc[3] - quater_loc[0]*quater_loc[2]),    // [2][0]
     2.0 * (quater_loc[2]*quater_loc[3] + quater_loc[0]*quater_loc[1]),    // [2][1]
     quater_pow2[0] - quater_pow2[1] - quater_pow2[2] + quater_pow2[3]} }; // [2][2]

  /* Rotate orientation */
  cs_real_3_t orient_new;
  cs_math_33_3_product(rot_m, orient_loc, orient_new);

  /* Renormalise */
  cs_math_3_normalise(orient_new, orient_loc);

}

#if 0

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid rotation by Brownian on the sphere
 *        surface increment
 *
 * Principle:
 * ALGO 1: Trajectorial stochastic Euler scheme :
 * ALGO 2: semi exact law increment
 *         BM_rotation_phase_spheroid_by_wright_fisher
 *
 * \param[in]  ip           particle number
 * \param[in]  dtp          time step
 * \param[in]  d_minus      reduced parameter for the (isotropic)
 *                          fluctuation tensor
 * \param[in]  orientation  spheroid orientation (vector)
 * \param[in]  brown        normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_rotation_phase_spheroid_by_spherical_coordinates(const cs_lnum_t ip,
                                                     const cs_real_t dtp,
                                                     const cs_real_t d_minus,
                                                     cs_real_t   orient_loc[3],
                                                     const cs_real_t brown[])
{
  /* Time factor d_minus*/
  cs_real_t half_d_minus = 0.5 * d_minus;

  /* Define the nine Brownian increments W11 W12 W13 W21 W22 W23 W31 W32 W33 */
  //cs_real_t  dw11 = sqrt(dtp) * brown[ip*9+0];
  cs_real_t  dw12 = sqrt(dtp) * brown[ip*9+1];
  cs_real_t  dw13 = sqrt(dtp) * brown[ip*9+2];
  cs_real_t  dw21 = sqrt(dtp) * brown[ip*9+3];
  //cs_real_t  dw22 = sqrt(dtp) * brown[ip*9+4];
  cs_real_t  dw23 = sqrt(dtp) * brown[ip*9+5];
  cs_real_t  dw31 = sqrt(dtp) * brown[ip*9+6];
  cs_real_t  dw32 = sqrt(dtp) * brown[ip*9+7];
  //cs_real_t  dw33 = sqrt(dtp) * brown[ip*9+8];

  /* Define the angles theta and phi (3d case) */
  cs_real_t theta;
  cs_real_t phi;

  /* Define the Brownian rotation in global frame */
  cs_real_t dw_glob_frame[3];
  dw_glob_frame[0] =  dw32-dw23;
  dw_glob_frame[1] =  dw13-dw31;
  dw_glob_frame[2] =  dw21-dw12;

  /* Detect singularity to be handled
   * (division by 0 when the fiber is aligned with z direction) */
  cs_real_t singularity_threshold = 1.0e-6 ;
  cs_real_t singularity_value = 1.0 - cs_math_pow2(orient_loc[2]);

  /* Handle the case where the orientation is nearly aligned with z
     (division by 0) */
  if (cs_math_fabs(singularity_value) < singularity_threshold) {
    // Define the reference to singularity
    cs_real_t ax_singularity[3] = {0.0, 1.0, 0.0};
    // Get vector for rotation
    cs_real_t n_rot[3];
    n_rot[0] = orient_loc[1]*ax_singularity[2] -orient_loc[2]*ax_singularity[1];
    n_rot[1] = orient_loc[2]*ax_singularity[0] -orient_loc[0]*ax_singularity[2];
    n_rot[2] = orient_loc[0]*ax_singularity[1] -orient_loc[1]*ax_singularity[0];
    cs_math_3_normalise(n_rot, n_rot);
    // Compute rotation angle
    cs_real_t rot_angle = acos(  orient_loc[0]*ax_singularity[0]
                               + orient_loc[1]*ax_singularity[1]
                               + orient_loc[2]*ax_singularity[2]);
    // Compute the rotation matrix
    cs_real_33_t rot_m = {
       {cos(rot_angle) + cs_math_pow2(n_rot[0])*(1.0 - cos(rot_angle)),
        n_rot[0]*n_rot[1]*(1.0 - cos(rot_angle)) + n_rot[2]*sin(rot_angle),
        n_rot[0]*n_rot[2]*(1.0 - cos(rot_angle)) - n_rot[1]*sin(rot_angle)},

       {n_rot[0]*n_rot[1]*(1.0 - cos(rot_angle)) - n_rot[2]*sin(rot_angle),
        cos(rot_angle) + cs_math_pow2(n_rot[1])*(1.0 - cos(rot_angle)),
        n_rot[1]*n_rot[2]*(1.0 - cos(rot_angle)) + n_rot[0]*sin(rot_angle)},

       {n_rot[0]*n_rot[2]*(1.0 - cos(rot_angle)) + n_rot[1]*sin(rot_angle),
        n_rot[1]*n_rot[2]*(1.0 - cos(rot_angle)) - n_rot[0]*sin(rot_angle),
        cos(rot_angle) + cs_math_pow2(n_rot[2])*(1.0 - cos(rot_angle))}
      };

    // Compute the Brownian rotation in local frame
    cs_real_t dw_loc_frame[3];
    cs_math_33t_3_product(rot_m, dw_glob_frame, dw_loc_frame);

    // Compute the angles theta and phi
    theta = acos(ax_singularity[2])
            + (cs_math_pow2(half_d_minus)*ax_singularity[2]*dtp
               + half_d_minus*(  ax_singularity[1]*dw_loc_frame[0]
                               - ax_singularity[0]*dw_loc_frame[1]))
            / sqrt(1.0 - cs_math_pow2(ax_singularity[2]));
    phi =   atan2(ax_singularity[1],ax_singularity[0])
          + half_d_minus * (ax_singularity[2] * (ax_singularity[1]*dw_loc_frame[1]
          + ax_singularity[0]*dw_loc_frame[0])
                            / (1.0 - cs_math_pow2(ax_singularity[2]))
          - dw_loc_frame[2]);

    /* From spherical to cartesian coordinates theta in [0,pi], phi in [0,2pi) */

    // Start by a modulo operation (to have theta and phi in [-2pi,2pi]
    theta = fmod(theta, 2.0*cs_math_pi) ;
    phi = fmod(phi, 2.0*cs_math_pi) ;

    // Handle cases where theta is not in [0,pi]
    if (theta > cs_math_pi) {
         theta = 2.0*cs_math_pi - theta ;
         phi   = fmod(phi + cs_math_pi, 2.0*cs_math_pi) ;
    }
    else if ( (-cs_math_pi < theta) & (theta < 0.0) ){
         theta = -theta ;
         phi   = fmod(phi + cs_math_pi, 2.0*cs_math_pi) ;
    }
    else if ( (-2.0*cs_math_pi < theta) & (theta < -cs_math_pi) ){
         theta = 2.0*cs_math_pi - cs_math_fabs(theta) ;
    }
    if (phi < 0.0)
      phi = 2*cs_math_pi + phi;

    /* Compute the orientation from angles (theta and phi) */
    cs_real_3_t orient_new;
    orient_new[0] = sin(theta)*cos(phi);
    orient_new[1] = sin(theta)*sin(phi);
    orient_new[2] = cos(theta);
    // Project on the global coordinate frame
    cs_math_33_3_product(rot_m, orient_new, orient_loc);

  }
  else{ // Handle the general case
    theta = acos(orient_loc[2])
        + (cs_math_pow2(half_d_minus)*orient_loc[2]*dtp
        + half_d_minus*(  orient_loc[1]*dw_glob_frame[0]
                        - orient_loc[0]*dw_glob_frame[1]))
        / sqrt(singularity_value);

    phi = atan2(orient_loc[1],orient_loc[0])
        + half_d_minus*(orient_loc[2]*(  orient_loc[1]*dw_glob_frame[1]
                                       + orient_loc[0]*dw_glob_frame[0])
                        / singularity_value
        - dw_glob_frame[2]);

    /* From spherical to cartesian coordinates theta in [0,pi], phi in [0,2pi) */

    // Start by a modulo operation (to have theta and phi in [-2pi,2pi]
    theta = fmod(theta, 2.0*cs_math_pi) ;
    phi = fmod(phi, 2.0*cs_math_pi) ;

    // Handle cases where theta is not in [0,pi]
    if (theta > cs_math_pi) {
      theta = 2*cs_math_pi - theta ;
      phi   = fmod(phi + cs_math_pi, 2.0*cs_math_pi);
    }
    else if ( (-cs_math_pi < theta) & (theta < 0.0) ){
      theta = -theta ;
      phi   = fmod(phi + cs_math_pi, 2.0*cs_math_pi);
    }
    else if ( (-2.0*cs_math_pi < theta) & (theta < -cs_math_pi) ){
      theta = 2.0*cs_math_pi - cs_math_fabs(theta) ;
    }
    if (phi < 0.0)
      phi = 2*cs_math_pi + phi;

    /* Compute the orientation from angles (theta and phi) */
    orient_loc[0] = sin(theta)*cos(phi);
    orient_loc[1] = sin(theta)*sin(phi);
    orient_loc[2] = cos(theta);
  }
}

#endif

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of the equations of motion for particle orientation:
 *
 * Details in publication by M. Bossy and L. Campana (INRIA)
 *
 * \param[in] iprev     time step indicator for fields
 *                        0: use fields at current time step
 *                        1: use fields at previous time step
 * \param[in] dt_p      lagrangian time step
 * \param[in] gradvf    fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_orientation_dyn_spheroids(int              iprev,
                                  const cs_real_t  dt_p,
                                  const cs_real_t  gradvf[][3][3])
{
  /*===================
   * 1. Initializations
   *===================*/

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_real_t *brown = NULL;
  BFT_MALLOC(brown, p_set->n_particles*9, cs_real_t);

  const cs_real_t *cvar_k = NULL, *cvar_ep = NULL, *cvar_omg = NULL;
  if (extra->cvar_k != NULL) {
    cvar_k = (extra->cvar_k->n_time_vals > 1) ?
      extra->cvar_k->vals[iprev] : extra->cvar_k->val;
  }
  if (extra->cvar_ep != NULL) {
    cvar_ep = (extra->cvar_ep->n_time_vals > 1) ?
      extra->cvar_ep->vals[iprev] : extra->cvar_ep->val;
  }
  if (extra->cvar_omg != NULL) {
    cvar_omg = (extra->cvar_omg->n_time_vals > 1) ?
      extra->cvar_omg->vals[iprev] : extra->cvar_omg->val;
  }

  if (! (   (extra->itytur >= 2 && extra->itytur <= 50)
         || extra->iturb == 60))
    bft_error
      (__FILE__, __LINE__, 0,
       _("The lagrangian turbulent dispersion model is not compatible\n"
         "with the selected turbulence model.\n"
         "\n"
         "Turbulent dispersion is taken into account with IDISTU = %d.\n"
         " Activated turbulence model is %d, when only k-eps, LES, Rij-eps,\n"
         " V2f or k-omega are compatible lagrangian turbulent dispersion."),
       (int)cs_glob_lagr_model->idistu,
       (int)extra->iturb);

  /*===============================================
   * 2. Integration of the (S)DE on the orientation
   *===============================================*/

  /* Loop on particles */
  for (cs_lnum_t p_id = 0; p_id < p_set->n_particles; p_id++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

    /* Generation of random numbers */
    /* 9 gaussian increments numered as 3x3 tensor */
    /* W11 W12 W13 W21 W22 W23 W31 W32 W33 */
    cs_random_normal(9, &(brown[9 * p_id]));

    /* Get local flow properties
       cell_id :    id of the cell
       romf    :    fluid density
       visccf  :    kinematic viscosity
       epsilon :    dissipation rate of turbulent kinetic energy
       tau_eta :    Kolmogorov timescale
    */

    cs_lnum_t cell_id
      = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

    cs_real_t romf = extra->cromf->val[cell_id];

    cs_real_t visccf = extra->viscl->val[cell_id] / romf;

    cs_real_t epsilon = 0;
    if (cvar_ep != NULL) {
      epsilon = cvar_ep[cell_id];
    }
    else if (extra->iturb == 60) {
      epsilon = extra->cmu * cvar_k[cell_id] * cvar_omg[cell_id];
    }

    cs_real_t tau_eta = sqrt(visccf / epsilon);
    /*
      Parameters for the (isotropic) fluctuation tensor (d_1,....,d_9)
      Assuming isotropy, the only non zero components are :
      d_1 = -1.0 / (sqrt(3.0 * tau_eta);
      d_2 = 0.5 * (sqrt(3.0/tau_eta) + sqrt(5.0/tau_eta));
      d_3 = 0.5 * (sqrt(3.0/tau_eta) - sqrt(5.0/tau_eta));
      Only this reduced information is needed :
      d_plus = d_2 + d_3
    */
    cs_real_t d_1 = -1.0 / sqrt(3.0 * tau_eta);

    cs_real_t d_plus = sqrt(3.0/tau_eta);
    cs_real_t d_minus = sqrt(5.0/tau_eta);

    /* Get particle orientation*/
    cs_real_t *orientation  = cs_lagr_particles_attr(p_set, p_id,
                                                     CS_LAGR_ORIENTATION);

    /* Get particle quaternion*/
    cs_real_t *quaternion  = cs_lagr_particles_attr(p_set, p_id,
                                                    CS_LAGR_QUATERNION);

    /* Get particle shape_param*/
    cs_real_t *radii = cs_lagr_particles_attr(p_set, p_id,
                                              CS_LAGR_RADII);
    cs_real_t Lambda = (cs_math_pow2(radii[2]/radii[1]) - 1.0)
                     / (cs_math_pow2(radii[2]/radii[1]) + 1.0);

    /* Extract symmetric and anti-symmetric parts
     * of the velocity gradients (matrix) */
    cs_real_t gradvf_sym[3][3];
    cs_real_t gradvf_ant[3][3];
    cs_math_33_extract_sym_ant(gradvf[cell_id], gradvf_sym, gradvf_ant);

    /* First step:
       Stretching of spheroids by the mean velocity gradient */
    // Add a value to the diagonal elements (if needed)
    cs_real_t beta = 0.5 * cs_math_pow2(Lambda) *
                    (2. * cs_math_pow2(d_1) +
                     cs_math_pow2(d_plus)  +
                     cs_math_pow2(d_1+d_plus))
                   + 0.5 * cs_math_pow2(d_minus);

    _mean_stretching_phase_spheroid(dt_p,
                                    orientation,
                                    Lambda,
                                    beta,
                                    gradvf_sym);

    /* Second step:
       Rotation of spheroids by the mean velocity gradient */
    _mean_rotation_phase_spheroid(dt_p,
                                  orientation,
                                  gradvf_ant);

    /* Third step:
       Stretching of spheroids by Brownian motion */
    _bm_stretching_phase_spheroid(p_id,
                                  dt_p,
                                  d_plus,
                                  orientation,
                                  Lambda,
                                  brown);

    /* Fourth step:
       Rotation of spheroids by Brownian motion and renormalisation */
    /*
    _bm_rotation_phase_spheroid_by_spherical_coordinates(p_id,
                                                         dt_p,
                                                         d_minus,
                                                         orientation,
                                                         brown);
    */
    _bm_rotation_phase_spheroid_by_quaternion(p_id,
                                              dt_p,
                                              d_minus,
                                              orientation,
                                              quaternion,
                                              brown);

    cs_math_3_normalise(orientation, orientation);
  }

  /* Free memory */
  BFT_FREE(brown);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of the Jeffey equations in DNS mode
 *
 * \param[in] dt_p      lagrangian time step
 * \param[in] gradvf    fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_orientation_dyn_jeffery(cs_real_t        dt_p,
                                const cs_real_t  gradvf[][3][3])
{
  /* 1. Initializations
     ================== */

  cs_lagr_particle_set_t         *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  /* 2. Integration of the (S)DE on the angular velocity
     =================================================== */

  /* Loop on particles */

  for (cs_lnum_t p_id = 0; p_id < p_set->n_particles; p_id++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

    /* Get local flow properties
       cell_id :    id of the cell
       romf    :    fluid density
       visccf  :    kinematic viscosity
       epsilon :    dissipation rate of turbulent kinetic energy
       tau_eta :    Kolmogorov timescale
    */

    cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

    /* Euler parameters */
    cs_real_t *euler = cs_lagr_particle_attr(particle, p_am,
                                             CS_LAGR_EULER);

    cs_real_33_t trans_m = {
      {2.*(euler[0]*euler[0]+euler[1]*euler[1]-0.5),   /* (0,0) */
       2.*(euler[1]*euler[2]+euler[0]*euler[3])    ,   /* (0,1) */
       2.*(euler[1]*euler[3]-euler[0]*euler[2])    },  /* (0,2) */

      {2.*(euler[1]*euler[2]-euler[0]*euler[3])    ,  /* (1,0) */
       2.*(euler[0]*euler[0]+euler[2]*euler[2]-0.5),  /* (1,1) */
       2.*(euler[2]*euler[3]+euler[0]*euler[1])    }, /* (1,2) */

      {2.*(euler[1]*euler[3]+euler[0]*euler[2])    ,  /* (2,0) */
       2.*(euler[2]*euler[3]-euler[0]*euler[1])    ,  /* (2,1) */
       2.*(euler[0]*euler[0]+euler[3]*euler[3]-0.5)}  /* (2,2) */
    };

    /* Fluid velocity gradient in the relative reference frame of the particle */
    cs_real_t grad_vf_r[3][3];

    cs_math_33_transform_a_to_r(gradvf[cell_id], trans_m, grad_vf_r);

    /* Ellipsoid radii */
    cs_real_t *radii = cs_lagr_particle_attr(particle, p_am,
                                             CS_LAGR_RADII);

    /* Corresponding shape parameters */
    cs_real_t *s_p = cs_lagr_particle_attr(particle, p_am,
                                           CS_LAGR_SHAPE_PARAM);

    /* Tau_p */
    cs_real_t taup = dt_p; //FIXME: need of small time steps !! (e-8 in dns)

    /* Inverse of a time scale for angular velocity
     * in the relative reference frame */
    cs_real_t zeta[3] =
      {   (40.0/9.0)*pow(radii[0]*radii[1]*radii[2], 2.0/3.0)
        / (radii[1]*radii[1]*s_p[1]+radii[2]*radii[2]*s_p[2])/taup,
          (40.0/9.0)*pow(radii[0]*radii[1]*radii[2],2.0/3.0)
        / (radii[0]*radii[0]*s_p[0]+radii[2]*radii[2]*s_p[2])/taup,
          (40.0/9.0)*pow(radii[0]*radii[1]*radii[2],2.0/3.0)
        / (radii[0]*radii[0]*s_p[0]+radii[1]*radii[1]*s_p[1])/taup};

    cs_real_t ezeta[3] = {exp(-dt_p * zeta[0]),
                          exp(-dt_p * zeta[1]),
                          exp(-dt_p * zeta[2])};

    /* Shape factors (diagonal matrix) */
    cs_real_t aom[3] = {
      (radii[2]*radii[2]-radii[1]*radii[1])/(radii[2]*radii[2]+radii[1]*radii[1]),
      (radii[0]*radii[0]-radii[2]*radii[2])/(radii[0]*radii[0]+radii[2]*radii[2]),
      (radii[1]*radii[1]-radii[0]*radii[0])/(radii[1]*radii[1]+radii[0]*radii[0])};

    /* 3. Integration of the (S)DE on the Euler parameters and angular velocity
       ======================================================================== */

    cs_real_t *ang_vel = cs_lagr_particle_attr(particle, p_am,
                                               CS_LAGR_ANGULAR_VEL);

    /* Integration of the Euler parameters:
     * Equation (26) of P. H. Mortensen (2008)*/
    cs_real_t d_euler[4] = {
      0.5 * (-ang_vel[0]*euler[1] - ang_vel[1]*euler[2] - ang_vel[2]*euler[3]),
      0.5 * ( ang_vel[0]*euler[0] - ang_vel[1]*euler[3] + ang_vel[2]*euler[2]),
      0.5 * ( ang_vel[0]*euler[3] + ang_vel[1]*euler[0] - ang_vel[2]*euler[1]),
      0.5 * (-ang_vel[0]*euler[2] + ang_vel[1]*euler[1] + ang_vel[2]*euler[0]),
    };

    for (int dim = 0; dim < 4; dim++)
      euler[dim] += dt_p * d_euler[dim];

    cs_real_t euler_norm = sqrt( cs_math_pow2(euler[0])
                               + cs_math_pow2(euler[1])
                               + cs_math_pow2(euler[2])
                               + cs_math_pow2(euler[3]) );

    for (int dim = 0; dim < 4; dim++)
      euler[dim] = euler[dim] / euler_norm;

    /* Time integration for the angular velocity
     * (correspond to equation (9) of C. Siewert et al 2013) */

    /* Angular velocity, exponential scheme */
    cs_real_3_t d_ang_vel = {
      aom[0]*(  ang_vel[1]*ang_vel[2]/zeta[0]
              - 0.5*(grad_vf_r[2][1] + grad_vf_r[1][2]))
      +  0.5*(grad_vf_r[2][1] - grad_vf_r[1][2]),

      aom[1]*(  ang_vel[0]*ang_vel[2]/zeta[1]
              - 0.5*(grad_vf_r[0][2] + grad_vf_r[2][0]))
      + 0.5*(grad_vf_r[0][2] - grad_vf_r[2][0]),

      aom[2]*(  ang_vel[0]*ang_vel[1]/zeta[2]
              - 0.5*(grad_vf_r[1][0] + grad_vf_r[0][1]))
      + 0.5*(grad_vf_r[1][0] - grad_vf_r[0][1])};

    for (int dim = 0; dim < 3; dim++)
      ang_vel[dim] = ezeta[dim] * ang_vel[dim]
        + (1. - ezeta[dim]) * d_ang_vel[dim];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
