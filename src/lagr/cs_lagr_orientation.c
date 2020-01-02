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
 * \param[in]  cell_id       cell where the spheroid is
 * \param[in]  orientation   spheroid orientation (vector)
 * \param[in]  gradvf        fluid velocity gradient (tensor)
 */
/*----------------------------------------------------------------------------*/

static void
_mean_stretching_phase_spheroid(cs_real_t        dtp,
                                cs_lnum_t        cell_id,
                                cs_real_t        orientation[3],
                                const cs_real_t  gradvf[][3][3])
{
  /* Calculate stretching */
  //TODO : add the analylitcal algorithm
  for (int id = 0; id < 3; id++) {
    // orientation[id] = orientation[id] * exp( gradvf[cell_id][id][id] * dtp);
    orientation[id] = orientation[id] + 1 -1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renormalization of the spheroid orientation
 *
 * Principle: ensure that the length of the spheroid (norm) stays equal to 1
 *
 * \param[in]  orientation             Rod orientation (vector)
 */
/*----------------------------------------------------------------------------*/

static void
_renormalize_spheroid(cs_real_t  orientation[3])
{
  /* Calculate the norm of the spheroid */
  cs_real_t norm_orient = cs_math_3_norm(orientation);

  /* Renormalise the spheroid */
  for (int id = 0; id < 3; id++) {
    orientation[id] = orientation[id] / norm_orient;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid rotation increment by antisymetric mean
 *        velocity gradient
 *
 * Principle: exact matrix rotation (Rodriges' formula)
 *
 * TODO: use directly the antisymetric mean velocity gradient matrix Omega
 *
 * \param[in]  dtp                time step
 * \param[in]  cell_id            cell where the spheroid is
 * \param[in]  orientation        spheroid orientation (vector)
 * \param[in]  gradvf             fluid velocity gradient (tensor)
 */
/*----------------------------------------------------------------------------*/

static void
_mean_rotation_phase_spheroid(cs_real_t        dtp,
                              cs_lnum_t        cell_id,
                              cs_real_t        orientation[3],
                              const cs_real_t  gradvf[][3][3])
{
  /* Calculate the vector n and its norm */
  cs_real_t n_norm
    = sqrt(  cs_math_pow2(gradvf[cell_id][0][1] - gradvf[cell_id][1][0])
           + cs_math_pow2(gradvf[cell_id][1][2] - gradvf[cell_id][2][1])
           + cs_math_pow2(gradvf[cell_id][0][2] - gradvf[cell_id][2][0]));

  cs_real_t n[3];
  n[0] = (gradvf[cell_id][0][1] - gradvf[cell_id][1][0]) / n_norm;
  n[1] = (gradvf[cell_id][2][0] - gradvf[cell_id][0][2]) / n_norm;
  n[2] = (gradvf[cell_id][1][2] - gradvf[cell_id][2][1]) / n_norm;

  /* Calculate the rotation of the spheroid */
  cs_real_t orientation_new[3];
  cs_real_t t_n = dtp * n_norm;

  // TODO: use multiplication of matrix and vector (to simplify)
  orientation_new[0]
    =    (cos(t_n) + n[0]*n[0]*(1-cos(t_n)))      * orientation[0]
       + (n[0]*n[1]*(1-cos(t_n)) - n[2]*sin(t_n)) * orientation[1]
       + (n[0]*n[2]*(1-cos(t_n)) + n[1]*sin(t_n)) * orientation[2];
  orientation_new[1]
    =    (n[0]*n[1]*(1-cos(t_n)) + n[2]*sin(t_n)) * orientation[0]
       + (cos(t_n) + n[1]*n[1]*(1-cos(t_n)))      * orientation[1]
       + (n[1]*n[2]*(1-cos(t_n)) + n[0]*sin(t_n)) * orientation[2];
  orientation_new[2]
    =    (n[0]*n[2]*(1-cos(t_n)) - n[1]*sin(t_n)) * orientation[0]
       + (n[1]*n[2]*(1-cos(t_n)) + n[1]*sin(t_n)) * orientation[1]
       + (cos(t_n) + n[2]*n[2]*(1-cos(t_n)))      * orientation[2];

  /* Update values */
  for (int i = 0; i < 3; i++)
    orientation = orientation_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Treatment of spheroid stretching by Brownian effect
 *
 * Principle: Euler scheme
 *
 * \param[in]  ip           particle number
 * \param[in]  dtp          time step
 * \param[in]  d_s          parameter for symmetric part
 * \param[in]  orientation  spheroid orientation (vector)
 * \param[in]  brown        normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_stretching_phase_spheroid(cs_lnum_t        ip,
                              cs_real_t        dtp,
                              cs_real_t        d_s,
                              cs_real_t        orientation[3],
                              const cs_real_t  brown[])
{
  /* Time factor b_prime*/
  cs_real_t b_prime = sqrt(d_s / 5.0);

  /* Auxiliary parameters for calculation */
  cs_real_t aux1 = b_prime * orientation[0] * orientation[1] * orientation[2];
  cs_real_t r1_w2 = cs_math_pow2(orientation[0]);
  cs_real_t r2_w2 = cs_math_pow2(orientation[1]);
  cs_real_t r3_w2 = cs_math_pow2(orientation[2]);
  cs_real_t aux5 = 1. / (1. + 0.5 * cs_math_pow2(b_prime) * dtp);

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

  orientation_new[0] = aux5 * ( orientation[0] - aux1 * (dw23+dw32)
      + orientation[0] * b_prime * ( (1. - r1_w2)*dw11 - r2_w2*dw22 - r3_w2*dw33 )
      + 0.5*b_prime*(1.-2.*r1_w2)*(orientation[1]*(dw12+dw21)+orientation[2]*(dw13+dw31)));

  orientation_new[1] = aux5 * ( orientation[1] - aux1 * (dw13+dw31)
      + orientation[1] * b_prime * ( (1. - r2_w2)*dw22 - r1_w2*dw11 - r3_w2*dw33 )
      + 0.5*b_prime*(1.-2.*r2_w2)*(orientation[0]*(dw12+dw21)+orientation[2]*(dw23+dw32)));

  orientation_new[2] = aux5 * ( orientation[2] - aux1 * (dw12+dw21)
      + orientation[2] * b_prime * ( (1. - r3_w2)*dw33 - r1_w2*dw11 - r2_w2*dw22 )
      + 0.5*b_prime*(1.-2.*r3_w2)*(orientation[1]*(dw23+dw32)+orientation[0]*(dw13+dw31)));

  /* Renormalise for security */
  _renormalize_spheroid(orientation_new);

  /* Update spheroid orientation */
  for (int i = 0; i < 3; i++)
    orientation  = orientation_new;
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
 * \param[in]  D_s                Parameter for symmetric part
 * \param[in]  orientation        spheroid orientation (vector)
 * \param[in]  brown              normal random number (9 values by particles)
 */
/*----------------------------------------------------------------------------*/

static void
_bm_rotation_phase_spheroid_by_spherical_coordinates(cs_lnum_t        ip,
                                                     cs_real_t        dtp,
                                                     cs_real_t        d_s,
                                                     cs_real_t        orient_loc[3],
                                                     const cs_real_t  brown[])
{
  /* Time factor b_prime*/
  cs_real_t half_b_prime = 0.5 * sqrt(d_s / 5.0);

  /* W11 W12 W13 W21 W22 W23 W31 W32 W33 */
  cs_real_t  dw11 = brown[ip*9+0];
  cs_real_t  dw12 = brown[ip*9+1];
  cs_real_t  dw13 = brown[ip*9+2];
  cs_real_t  dw21 = brown[ip*9+3];
  cs_real_t  dw22 = brown[ip*9+4];
  cs_real_t  dw23 = brown[ip*9+5];
  cs_real_t  dw31 = brown[ip*9+6];
  cs_real_t  dw32 = brown[ip*9+7];
  cs_real_t  dw33 = brown[ip*9+8];

  cs_real_t theta0;
  cs_real_t phi0;
  cs_random_uniform(1, &theta0);
  cs_random_uniform(1, &phi0);
  theta0   = acos(0.8*(2.0*theta0-1.0));
  phi0     = phi0*2.0*cs_math_pi;
  cs_real_t axe0[3];
  axe0[0] = sin(theta0)*cos(phi0);
  axe0[1] = sin(theta0)*sin(phi0);
  axe0[2] = cos(theta0);

  cs_real_t brown_old[3];
  brown_old[0] = sqrt(dtp) * dw23-dw32;
  brown_old[1] = sqrt(dtp) * dw31-dw13;
  brown_old[2] = sqrt(dtp) * dw12-dw21;

  cs_real_t n_rot[3];
  n_rot[0] = orient_loc[1]*axe0[2] -orient_loc[2]*axe0[1];
  n_rot[1] = orient_loc[2]*axe0[0] -orient_loc[0]*axe0[2];
  n_rot[2] = orient_loc[0]*axe0[1] -orient_loc[1]*axe0[0];
  cs_real_t norm_rot = cs_math_3_norm(n_rot);
  n_rot[0] = n_rot[0] / norm_rot;
  n_rot[1] = n_rot[1] / norm_rot;
  n_rot[2] = n_rot[2] / norm_rot;

  cs_real_t rot_angle = acos(  orient_loc[0]*axe0[0]
                             + orient_loc[1]*axe0[1]
                             + orient_loc[2]*axe0[2]);

  double rot_m[3][3];
  rot_m[0][0] = cos(rot_angle)+ cs_math_pow2(n_rot[0])*(1.0 -cos(rot_angle));
  rot_m[1][0] = n_rot[0]*n_rot[1]*(1.0 -cos(rot_angle)) +n_rot[2]*sin(rot_angle);
  rot_m[2][0] = n_rot[0]*n_rot[2]*(1.0 -cos(rot_angle)) -n_rot[1]*sin(rot_angle);
  rot_m[0][1] = n_rot[0]*n_rot[1]*(1.0 -cos(rot_angle)) -n_rot[2]*sin(rot_angle);
  rot_m[1][1] = cos(rot_angle)+ cs_math_pow2(n_rot[1])*(1.0 -cos(rot_angle));
  rot_m[2][1] = n_rot[1]*n_rot[2]*(1.0 -cos(rot_angle)) +n_rot[0]*sin(rot_angle);
  rot_m[0][2] = n_rot[0]*n_rot[2]*(1.0 -cos(rot_angle)) +n_rot[1]*sin(rot_angle);
  rot_m[1][2] = n_rot[1]*n_rot[2]*(1.0 -cos(rot_angle)) -n_rot[0]*sin(rot_angle);
  rot_m[2][2] = cos(rot_angle)+ cs_math_pow2(n_rot[2])*(1.0 -cos(rot_angle));

  cs_real_t brown_rot[3];
  brown_rot[0] =  rot_m[0][0]*brown_old[0] + rot_m[0][1]*brown_old[1]
                + rot_m[0][2]*brown_old[2];
  brown_rot[1] =  rot_m[1][0]*brown_old[0] + rot_m[1][1]*brown_old[1]
                + rot_m[1][2]*brown_old[2];
  brown_rot[2] =  rot_m[2][0]*brown_old[0] + rot_m[2][1]*brown_old[1]
                + rot_m[2][2]*brown_old[2];

  cs_real_t theta  = theta0 + (cs_math_pow2(half_b_prime)*axe0[2]*dtp
                               + half_b_prime*(axe0[1]*brown_rot[0] -axe0[0]*brown_rot[0]))
                     /(sqrt(1.0 -cs_math_pow2(axe0[2])));
  cs_real_t phi    = phi0  + half_b_prime*(axe0[2]*(axe0[1]*brown_rot[1]
                      - axe0[0]*brown_rot[0])/(1.0 -cs_math_pow2(axe0[2])) -brown_rot[2]);

  cs_real_t orient_new[3];
  orient_new[0] = sin(theta)*cos(phi);
  orient_new[1] = sin(theta)*sin(phi);
  orient_new[2] = cos(theta);

  /* Update orientation */
  orient_loc[0] =   rot_m[0][0]*orient_new[0]
                  + rot_m[1][0]*orient_new[1]
                  + rot_m[2][0]*orient_new[2];
  orient_loc[1] =   rot_m[0][1]*orient_new[0]
                  + rot_m[1][1]*orient_new[1]
                  + rot_m[2][1]*orient_new[2];
  orient_loc[2] =   rot_m[0][2]*orient_new[0]
                  + rot_m[1][2]*orient_new[1]
                  + rot_m[2][2]*orient_new[2];
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
cs_lagr_orientation_dyn_spheroids(int       iprev,
                                  cs_real_t dt_p,
                                  const cs_real_33_t gradvf[])
{
  /* ==============================================================================*/
  /* 1. Initialisations                                                            */
  /* ==============================================================================*/

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_real_t *brown = NULL;
  cs_real_t *unif = NULL;
  BFT_MALLOC(brown, p_set->n_particles*9, cs_real_t);
  BFT_MALLOC(unif, p_set->n_particles, cs_real_t);

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
    cs_random_uniform(1, &(unif[p_id]));

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

    cs_real_t tau_eta = sqrt(visccf / epsilon) + 1.0;

    /* Parameters for the correlated noise
       D_s:        symmetric part for the forcing
       D_a:        antisymmetric part for the forcing
       b1, b2, b3: diagonal terms for the tensor
    */

    cs_real_t d_s = 15.0 / tau_eta;
    cs_real_t d_a = 15.0 / tau_eta;

    cs_real_t b_1 = (-1.0 / 3.0) * sqrt(d_s / 5.0);
    cs_real_t b_2 = 0.5 * (sqrt(d_s/5.0) + sqrt(d_a/3.0));
    cs_real_t b_3 = 0.5 * (sqrt(d_s/5.0) - sqrt(d_a/3.0));

    /* Get particle orientation*/
    cs_real_t *orientation  = cs_lagr_particles_attr(p_set, p_id,
                                                     CS_LAGR_ORIENTATION);

    /* First step:
       Stretching of spheroids by the mean velocity gradient

    _mean_stretching_phase_spheroid(dt_p,
                                    cell_id,
                                    orientation,
                                    gradvf);
    */

    /* Second step:
       Rotation of spheroids by the mean velocity gradient
    _mean_rotation_phase_spheroid(dt_p,
                                  cell_id,
                                  orientation,
                                  gradvf); */

    /* Third step:
       Stretching of spheroids by Brownian motion

    _bm_stretching_phase_spheroid(p_id,
                                  dt_p,
                                  d_s,
                                  orientation,
                                  brown);
    */

    /* Fourth step:
       Rotation of spheroids by Brownian motion and renormalisation
    */

    _bm_rotation_phase_spheroid_by_spherical_coordinates(p_id,
                                                         dt_p,
                                                         d_s,
                                                         orientation,
                                                         brown);

    _renormalize_spheroid(orientation);
  }

  /* Free memory */
  BFT_FREE(brown);
  BFT_FREE(unif);
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
cs_lagr_orientation_dyn_ellipsoids(int              iprev,
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
    cs_real_t taup = dt_p;

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
