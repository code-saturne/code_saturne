/*============================================================================
 * Methods for particle deposition
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

/*============================================================================
 * Functions dealing with particle deposition
 *============================================================================*/

#include "cs_defs.h"
#include "cs_math.h"

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

/* Boltzmann constant */
static const double _k_boltz = 1.38e-23;

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
  cs_real_33_t eig_vec = { 1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0 };
  cs_real_t tol_err = 1.0e-12;
  cs_real_33_t mat_loc ;
  mat_loc[0][0] = Lambda * (gradvf_sym[0][0] + beta);
  mat_loc[0][1] = Lambda *  gradvf_sym[0][1] ;
  mat_loc[0][2] = Lambda *  gradvf_sym[0][2] ;
  mat_loc[1][0] = Lambda *  gradvf_sym[1][0] ;
  mat_loc[1][1] = Lambda * (gradvf_sym[1][1] + beta);
  mat_loc[1][2] = Lambda *  gradvf_sym[1][2] ;
  mat_loc[2][0] = Lambda *  gradvf_sym[2][0] ;
  mat_loc[2][1] = Lambda *  gradvf_sym[2][1] ;
  mat_loc[2][2] = Lambda * (gradvf_sym[2][2] + beta);
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
  n[0] = - gradvf_ant[1][2] / n_norm;
  n[1] =   gradvf_ant[0][2] / n_norm;
  n[2] = - gradvf_ant[0][1] / n_norm;

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
 * \param[in]  d_prime      reduced parameter for the (isotropic) fluctuation tensor
 * \param[in]  orientation  spheroid orientation (vector)
 * \param[in]  Lambda       spheroid shape parameter
 * \param[in]  brown        normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_stretching_phase_spheroid(const cs_lnum_t ip,
                              const cs_real_t dtp,
                              const cs_real_t d_prime,
                                    cs_real_t orientation[3],
                              const cs_real_t Lambda,
                              const cs_real_t brown[])
{
  /* Auxiliary parameters for calculation */
  cs_real_t aux1 = orientation[0] * orientation[1] * orientation[2];
  cs_real_t r1_w2 = cs_math_pow2(orientation[0]);
  cs_real_t r2_w2 = cs_math_pow2(orientation[1]);
  cs_real_t r3_w2 = cs_math_pow2(orientation[2]);
  cs_real_t aux5 = Lambda * d_prime / (1. + 0.5 * cs_math_pow2(Lambda*d_prime) * dtp);

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

  cs_real_3_t orientation_new;

  /* solve the orientation dynamics produced by Brownian stretching */
  orientation_new[0] =  orientation[0] +
    aux5 * ( - aux1 * (dw23+dw32)
             + orientation[0] * ( (1. - r1_w2)*dw11 - r2_w2*dw22 - r3_w2*dw33 )
             + 0.5 * (1.-2.*r1_w2)*(orientation[1]*(dw12+dw21)+orientation[2]*(dw13+dw31)) );

  orientation_new[1] =  orientation[1] +
    aux5 * ( - aux1 * (dw13+dw31)
             + orientation[1] * ( (1. - r2_w2)*dw22 - r1_w2*dw11 - r3_w2*dw33 )
             + 0.5 * (1.-2.*r2_w2)*(orientation[0]*(dw12+dw21)+orientation[2]*(dw23+dw32)) );

  orientation_new[2] =  orientation[2] +
    aux5 * ( - aux1 * (dw12+dw21)
             + orientation[2] * ( (1. - r3_w2)*dw33 - r1_w2*dw11 - r2_w2*dw22 )
             + 0.5 * (1.-2.*r3_w2)*(orientation[1]*(dw23+dw32)+orientation[0]*(dw13+dw31)) );

  /* Renormalise for security */
  cs_math_3_normalise(orientation_new, orientation);

}



/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid rotation by Brownian on the sphere surface increment
 *
 * Principle:
 * ALGO 1  :  Trajectorial stochastic Euler scheme :
 * ALGO 2  :  semi exact law increment  BM_rotation_phase_spheroid_by_wright_fisher
 *
 *
 * \param[in]  ip                 particle number
 * \param[in]  dtp                time step
 * \param[in]  d_prime            reduced parameter for the (isotropic) fluctuation tensor
 * \param[in]  orientation        spheroid orientation (vector)
 * \param[in]  brown              normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_rotation_phase_spheroid_by_spherical_coordinates(const cs_lnum_t ip,
                                                     const cs_real_t dtp,
                                                     const cs_real_t d_prime,
                                                           cs_real_t orient_loc[3],
                                                     const cs_real_t brown[])
{

  /* Time factor d_prime*/
  cs_real_t half_d_prime = 0.5 * d_prime;

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
  cs_real_t dW_glob_frame[3];
  dW_glob_frame[0] =  dw23-dw32;
  dW_glob_frame[1] =  dw31-dw13;
  dW_glob_frame[2] =  dw12-dw21;

  /* Detect singularity to be handled
   * (division by 0 when the fiber is aligned with z direction) */
  cs_real_t singularity_threshold = 1.0e-6 ;
  cs_real_t singularity_value = 1.0 - cs_math_pow2(orient_loc[2]);

  if (singularity_value < singularity_threshold) {// Handle the case where the orientation is nearly aligned with z (division by 0)
    // Define the reference to singularity
    cs_real_t axe_singularity[3] = { 1.0, 0.0, 0.0 };
    // Get vector for rotation
    cs_real_t n_rot[3];
    n_rot[0] = orient_loc[1]*axe_singularity[2] -orient_loc[2]*axe_singularity[1];
    n_rot[1] = orient_loc[2]*axe_singularity[0] -orient_loc[0]*axe_singularity[2];
    n_rot[2] = orient_loc[0]*axe_singularity[1] -orient_loc[1]*axe_singularity[0];
    cs_math_3_normalise(n_rot, n_rot);
    // Compute rotation angle
    cs_real_t rot_angle = acos( orient_loc[0]*axe_singularity[0]
                              + orient_loc[1]*axe_singularity[1]
                              + orient_loc[2]*axe_singularity[2]);
    // Compute the rotation matrix
    cs_real_33_t rot_m = {
      cos(rot_angle) + cs_math_pow2(n_rot[0])*(1.0 - cos(rot_angle)),     // [0][0]
      n_rot[0]*n_rot[1]*(1.0 - cos(rot_angle)) + n_rot[2]*sin(rot_angle), // [0][1]
      n_rot[0]*n_rot[2]*(1.0 - cos(rot_angle)) - n_rot[1]*sin(rot_angle), // [0][2]
      n_rot[0]*n_rot[1]*(1.0 - cos(rot_angle)) - n_rot[2]*sin(rot_angle), // [1][0]
      cos(rot_angle) + cs_math_pow2(n_rot[1])*(1.0 - cos(rot_angle)),     // [1][1]
      n_rot[1]*n_rot[2]*(1.0 - cos(rot_angle)) + n_rot[0]*sin(rot_angle), // [1][2]
      n_rot[0]*n_rot[2]*(1.0 - cos(rot_angle)) + n_rot[1]*sin(rot_angle), // [2][0]
      n_rot[1]*n_rot[2]*(1.0 - cos(rot_angle)) - n_rot[0]*sin(rot_angle), // [2][1]
      cos(rot_angle) + cs_math_pow2(n_rot[2])*(1.0 - cos(rot_angle)) };   // [2][2]

    // Compute the Brownian rotation in local frame
    cs_real_t dW_loc_frame[3];
    cs_math_33_3_product(rot_m, dW_glob_frame, dW_loc_frame);

    // Compute the angles theta and phi
    theta = acos(axe_singularity[2])
            + (cs_math_pow2(half_d_prime)*axe_singularity[2]*dtp
               + half_d_prime*(axe_singularity[1]*dW_loc_frame[0]
                             - axe_singularity[0]*dW_loc_frame[1]))
            / sqrt(1.0 - cs_math_pow2(axe_singularity[2]));
    phi   = atan2(axe_singularity[1],axe_singularity[0])
            + half_d_prime * (axe_singularity[2] * (axe_singularity[1]*dW_loc_frame[1]
            + axe_singularity[0]*dW_loc_frame[0]) / (1.0 - cs_math_pow2(axe_singularity[2]))
            - dW_loc_frame[2]);

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
    cs_math_33t_3_product(rot_m, orient_new, orient_loc);

  }
  else{ // Handle the general case
    theta = acos(orient_loc[2])
        + (cs_math_pow2(half_d_prime)*orient_loc[2]*dtp
        + half_d_prime*(orient_loc[1]*dW_glob_frame[0] - orient_loc[0]*dW_glob_frame[1]))
        / sqrt(singularity_value);

    phi = atan2(orient_loc[1],orient_loc[0])
        + half_d_prime*(orient_loc[2]*(orient_loc[1]*dW_glob_frame[1] + orient_loc[0]*dW_glob_frame[0])/ singularity_value
        - dW_glob_frame[2]);

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
cs_lagr_orientation_dyn_spheroids(int iprev,
                                  const cs_real_t  dt_p,
                                  const cs_real_33_t gradvf[])
{
  /* ==============================================================================*/
  /* 1. Initialisations                                                            */
  /* ==============================================================================*/

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_real_t *brown = NULL;
  BFT_MALLOC(brown, p_set->n_particles*9, cs_real_t);

  /* =============================================================================
   * 2. Integration of the (S)DE on the orientation
   * =============================================================================*/

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


    cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

    cs_real_t romf = extra->cromf->val[cell_id];

    cs_real_t visccf = extra->viscl->val[cell_id] / romf;

    cs_real_t epsilon;
    if (extra->itytur == 2 || extra->itytur == 3 || extra->iturb == 50 ) {
      epsilon = extra->cvar_ep->vals[iprev][cell_id];
    }
    else if ( extra->iturb == 60) {
      epsilon = extra->cmu * extra->cvar_k->vals[iprev][cell_id]
                           * extra->cvar_omg->vals[iprev][cell_id];
    }
    else {
      bft_printf(_("\n WARNING: STOP AT LAGRANGIAN MODULE EXECUTION\n"));
      bft_printf
        (_("The lagrangian module is not compatible with the selected turbulence model.\n"
           "\n"
           "Turbulent dispersion is taken into account with IDISTU = %d\n"
           " Activated turbulence model is %d, when only k-eps, Rij-eps,\n"
           " V2f or k-omega are compatible with turbulent dispersion and Lagrangian module.\n"
           "\n"),
         (int)cs_glob_lagr_time_scheme->idistu,
         (int)extra->iturb);
      cs_exit(1);
    }

    cs_real_t tau_eta = sqrt(visccf / epsilon);
    /*
      Parameters for the (isotropic) fluctuation tensor (d_1,....,d_9)
      Assuming isotropy, the only non zero components are :
      d_1 = -1.0 / (sqrt(3.0 * tau_eta);
      d_2 = 0.5 * (sqrt(3.0/tau_eta) + sqrt(5.0/tau_eta));
      d_3 = 0.5 * (sqrt(3.0/tau_eta) - sqrt(5.0/tau_eta));
      Only this reduced information is needed :
      d_prime = d_2 + d_3
    */
    cs_real_t d_1 = -1.0 / sqrt(3.0 * tau_eta);
    cs_real_t d_2 = 0.5 * (sqrt(3.0/tau_eta) + sqrt(5.0/tau_eta));
    cs_real_t d_3 = 0.5 * (sqrt(3.0/tau_eta) - sqrt(5.0/tau_eta));

    cs_real_t d_prime = sqrt(3.0/tau_eta);

    /* Get particle orientation*/
    cs_real_t *orientation  = cs_lagr_particles_attr(p_set, p_id,
                                                     CS_LAGR_ORIENTATION);

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
    cs_real_t beta = 0.5 * ( cs_math_pow2(d_1+d_2+d_3) +
                             2.0*cs_math_pow2(d_1) + cs_math_pow2(d_prime) );

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
                                  d_prime,
                                  orientation,
                                  Lambda,
                                  brown);

    /* Fourth step:
       Rotation of spheroids by Brownian motion and renormalisation */
    _bm_rotation_phase_spheroid_by_spherical_coordinates(p_id,
                                                         dt_p,
                                                         d_prime,
                                                         orientation,
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
 * \param[in] iprev     time step indicator for fields
 *                        0: use fields at current time step
 *                        1: use fields at previous time step
 * \param[in] dt_p      lagrangian time step
 * \param[in] gradvf    fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_orientation_dyn_jeffery(int              iprev,
                                cs_real_t        dt_p,
                                const cs_real_t  gradvf[][3][3])
{
  /* 1. Initializations
     ================== */

  cs_lagr_particle_set_t         *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

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
      2.*(euler[0]*euler[0]+euler[1]*euler[1]-0.5), /* (0,0) */
      2.*(euler[1]*euler[2]-euler[0]*euler[3]),     /* (0,1) */
      2.*(euler[1]*euler[3]+euler[0]*euler[2]),     /* (0,2) */
      2.*(euler[1]*euler[2]+euler[0]*euler[3]),     /* (1,0) */
      2.*(euler[0]*euler[0]+euler[2]*euler[2]-0.5), /* (1,1) */
      2.*(euler[2]*euler[3]-euler[0]*euler[1]),     /* (1,2) */
      2.*(euler[1]*euler[3]-euler[0]*euler[2]),     /* (2,0) */
      2.*(euler[2]*euler[3]+euler[0]*euler[1]),     /* (2,1) */
      2.*(euler[0]*euler[0]+euler[3]*euler[3]-0.5)  /* (2,2) */
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
      aom[0]*(ang_vel[1]*ang_vel[2]/zeta[0] - 0.5*(grad_vf_r[2][1] + grad_vf_r[1][2]))
                                            + 0.5*(grad_vf_r[2][1] - grad_vf_r[1][2]),
      aom[1]*(ang_vel[0]*ang_vel[2]/zeta[1] - 0.5*(grad_vf_r[0][2] + grad_vf_r[2][0]))
                                            + 0.5*(grad_vf_r[0][2] - grad_vf_r[2][0]),
      aom[2]*(ang_vel[0]*ang_vel[1]/zeta[2] - 0.5*(grad_vf_r[1][0] + grad_vf_r[0][1]))
                                            + 0.5*(grad_vf_r[1][0] - grad_vf_r[0][1])};

    for (int dim = 0; dim < 3; dim++)
      ang_vel[dim] = ezeta[dim] * ang_vel[dim] + (1. - ezeta[dim]) * d_ang_vel[dim];

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
