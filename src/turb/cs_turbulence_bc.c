/*============================================================================
 * Base turbulence boundary conditions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_map.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_bc.c
        Base turbulence boundary conditions.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Variable ids for main turbulence variables

   CAUTION:

   A redesign of the icodcl/rcodcl boundary condition logic
   and simpler mapping to fields should remove the need for this
   function in the future.

   So this structure allows us an efficient access to current "rcodcl"
   values for turbulent variables in the meantime, but should not
   be a starting point for more general developments.

*/

typedef struct {

  cs_field_bc_coeffs_t  *bc_k;      /* BC coefficients for k */
  cs_field_bc_coeffs_t  *bc_eps;    /* BC coefficients for epsilon */

  cs_field_bc_coeffs_t  *bc_rij;    /* BC coefficients for R_ij tensor */

  cs_field_bc_coeffs_t  *bc_phi;    /* BC coefficients for phi */
  cs_field_bc_coeffs_t  *bc_f_bar;  /* BC coefficients for f_bar */
  cs_field_bc_coeffs_t  *bc_alp_bl; /* BC coefficients for blending alpha */

  cs_field_bc_coeffs_t  *bc_omg;    /* BC coefficients for omega */
  cs_field_bc_coeffs_t  *bc_nusa;   /* BC coefficients for nu_t (SA model) */

  int                    size_ut;   /* number of turbulent flux fields */
  int                    size_alp_bl_t;  /* number of blending alpha fields */
  cs_field_t  **f_ut;               /* Fields for turbulent fluxes */
  cs_field_t  **f_alp_bl_t;         /* Fields for blending alpha */

} cs_turb_bc_ptr_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* BC coefficientss for main turbulence variables */

static cs_turb_bc_ptr_t
_turb_bc_id =
  {
   NULL,  /* bc_k */
   NULL,  /* bc_eps */
   NULL,  /* bc_rij */
   NULL,  /* bc_phi */
   NULL,  /* bc_f_bar */
   NULL,  /* bc_alp_bl */
   NULL,  /* bc_omg */
   NULL,  /* bc_nusa */

   0,    /* size of ut */
   0,    /* size of alp_bl_t */
   NULL,  /* f_id_ut */
   NULL   /* f_id_alp_bl_t */
};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_num,
                                  double      uref2,
                                  double      dh,
                                  double      rho,
                                  double      mu);

void
cs_f_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_num,
                                        double      uref2,
                                        double      t_intensity,
                                        double      dh);

void
cs_f_turbulence_bc_inlet_k_eps(cs_lnum_t   face_num,
                               double      k,
                               double      eps,
                               double      vel_dir[],
                               double      shear_dir[]);

void
cs_f_turbulence_bc_init_inlet_k_eps(cs_lnum_t   face_num,
                                    double      k,
                                    double      eps);

void
cs_f_turbulence_bc_set_uninit_inlet_k_eps(cs_lnum_t   face_num,
                                          double      k,
                                          double      eps,
                                          double      vel_dir[],
                                          double      shear_dir[]);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Calculation of \f$ u^\star \f$, \f$ k \f$ and \f$\varepsilon \f$
 * from a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 * for a circular duct flow with smooth wall
 * (for inlet boundary conditions).
 *
 * Both \f$ u^\star \f$ and\f$ (k,\varepsilon )\f$ are returned, so that
 * the user may compute other values of \f$ k \f$ and \f$ \varepsilon \f$
 * with \f$ u^\star \f$.
 *
 * We use the laws from Idel'Cik, i.e.
 * the head loss coefficient \f$ \lambda \f$ is defined by:
 * \f[ |\dfrac{\Delta P}{\Delta x}| =
 *                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
 *
 * then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
 * \f$\lambda \f$ depends on the hydraulic Reynolds number
 * \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
 *  - for \f$ Re < 2000 \f$
 *      \f[ \lambda = \dfrac{64}{Re} \f]
 *
 *  - for \f$ Re > 4000 \f$
 *      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
 *
 *  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
 *      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
 *
 *  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
 *  from the well known formulae of developped turbulence
 *
 * \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
 * \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
 *
 * parameters:
 *   uref2  <-- square of the reference flow velocity
 *   dh     <-- hydraulic diameter \f$ D_H \f$
 *   rho    <-- mass density \f$ \rho \f$
 *   mu     <-- dynamic viscosity \f$ \nu \f$
 *   ustar2 --> square of friction speed
 *   k      --> calculated turbulent intensity \f$ k \f$
 *   eps    --> calculated turbulent dissipation
 *                          \f$ \varepsilon \f$
 *----------------------------------------------------------------------------*/

static inline void
_ke_hyd_diam(double   uref2,
             double   dh,
             double   rho,
             double   mu,
             double  *ustar2,
             double  *k,
             double  *eps)
{
  double re = sqrt(uref2)*dh*rho/mu;

  if (re < 2000) {

    /* in this case we calculate directly \f$u*^2\f$ to avoid an issue with
       \f$ xlmbda= \dfrac{64}{Re} \f$ when Re->0 */

    *ustar2 = 8.0*mu*sqrt(uref2)/rho/dh;

  }
  else if (re < 4000) {

    double xlmbda = 0.021377 + 5.3115e-6*re;
    *ustar2 = uref2*xlmbda/8.0;

  }
  else {

    /* xlmbda = 1.0/(1.8*log(re) / log(10.0)-1.64)^2;
       ustar2 = uref2*xlmbda/8.0; */

    double a = 1.8*log(re) / log(10.0)-1.64;
    *ustar2 = uref2*0.125/(a*a);

  }

  *k   = *ustar2 / sqrt(cs_turb_cmu);
  *eps = pow(*ustar2, 1.5) / (cs_turb_xkappa*dh*0.1);
}

/*----------------------------------------------------------------------------
 * Calculation of \f$ k \f$ and \f$\varepsilon\f$
 * from a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 * and the reference velocity \f$ U_{ref} \f$
 * for a circular duct flow with smooth wall
 * (for inlet boundary conditions).
 *
 * \f[ k = 1.5 I {U_{ref}}^2 \f]
 * \f[ \varepsilon = 10 \dfrac{{C_\mu}^{0.75} k^{1.5}}{ \kappa D_H} \f]
 *
 * parameters:
 *   uref2       <-- square of the reference flow velocity
 *   t_intensity <-- turbulent intensity \f$ I \f$
 *   dh          <-- hydraulic diameter \f$ D_H \f$
 *   k           --> calculated turbulent intensity \f$ k \f$
 *   eps         --> calculated turbulent dissipation
 *                              \f$ \varepsilon \f$
 *----------------------------------------------------------------------------*/

static inline void
_ke_turb_intensity(double   uref2,
                   double   t_intensity,
                   double   dh,
                   double  *k,
                   double  *eps)
{
  *k   = 1.5*uref2*t_intensity*t_intensity;
  *eps = 10.0*pow(cs_turb_cmu, 0.75)*pow(*k, 1.5)/(cs_turb_xkappa*dh);
}

/*----------------------------------------------------------------------------*
 * Assign turbulent boundary condition to a given face
 *
 * parameters:
 *   face_id     <-- face id
 *   k           <-- k
 *   eps         <-- epsilon
 *   vel_dir     <-- velocity direction
 *   shear_dir   <-- shear direction, it also contains the level of
 *                   anisotropy (Rnt = a_nt k)
 *----------------------------------------------------------------------------*/

static inline void
_inlet_bc(cs_lnum_t   face_id,
          double      k,
          double      eps,
          double      vel_dir[],
          double      shear_dir[])
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_turb_model_t  *turb_model = cs_glob_turb_model;
  assert(turb_model != NULL);

  if (turb_model->itytur == 2) {

    _turb_bc_id.bc_k->rcodcl1[face_id] = k;
    _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

  }

  else if (turb_model->order == CS_TURB_SECOND_ORDER) {

    double d2s3 = 2./3.;

    for (int ii = 0; ii < 3; ii++)
      _turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id] = d2s3 * k;
    if (vel_dir != NULL) {
      cs_math_3_normalize(vel_dir, vel_dir);
      /* Note: do not normalize shear_dir because it contains the
         level of anisotropy */
      /* Rxy */
      _turb_bc_id.bc_rij->rcodcl1[3*n_b_faces + face_id]
        = k * (vel_dir[0]*shear_dir[1] + vel_dir[1]*shear_dir[0]);
      /* Ryz */
      _turb_bc_id.bc_rij->rcodcl1[4*n_b_faces + face_id]
        = k * (vel_dir[1]*shear_dir[2] + vel_dir[2]*shear_dir[1]);
      /* Rxz */
      _turb_bc_id.bc_rij->rcodcl1[5*n_b_faces + face_id]
        =  k * (vel_dir[0]*shear_dir[2] + vel_dir[2]*shear_dir[0]);
    }
    else {
      for (int ii = 3; ii < 6; ii++)
        _turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id] = 0.;
    }
    _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
      _turb_bc_id.bc_alp_bl->rcodcl1[face_id] = 1.;

    /* Initialization of the turbulent fluxes to 0 if DFM or
     * EB-DFM are used for scalars (turbulence_flux_model = 30 or 31)
     * Alpha_theta for EB-DFM / EB-AFM / EB-GGDH is
     * initialize to 1 */

    if (_turb_bc_id.size_ut > 0) {
      for (int id = 0; id < _turb_bc_id.size_ut; id++) {
        cs_field_t *fld = _turb_bc_id.f_ut[id];
        for (int ii = 0; ii < fld->dim; ii++)
          fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id] = 0.;
      }
    }

    if (_turb_bc_id.size_alp_bl_t > 0) {
      for (int id = 0; id < _turb_bc_id.size_alp_bl_t; id++) {
        cs_field_t *fld = _turb_bc_id.f_alp_bl_t[id];
        for (int ii = 0; ii < fld->dim; ii++)
          fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id] = 1.;
      }
    }

  }
  else if (turb_model->itytur == 5) {

    _turb_bc_id.bc_k->rcodcl1[face_id] = k;
    _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

    _turb_bc_id.bc_phi->rcodcl1[face_id] = 2./3.;
    if (turb_model->iturb == CS_TURB_V2F_PHI) {
      _turb_bc_id.bc_f_bar->rcodcl1[face_id] = 0.;
    }
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      _turb_bc_id.bc_alp_bl->rcodcl1[face_id] = 0.;
    }

  }
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    _turb_bc_id.bc_k->rcodcl1[face_id] = k;
    _turb_bc_id.bc_omg->rcodcl1[face_id] = eps/cs_turb_cmu/k;

  }
  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {

    _turb_bc_id.bc_nusa->rcodcl1[face_id] = cs_turb_cmu*k*k/eps;

  }
}

/*----------------------------------------------------------------------------*
 * Assign turbulent boundary condition to a given face only if unitialized.
 *
 * parameters:
 *   face_id <-- face id
 *   k       <-- k
 *   eps     <-- epsilon
 *   vel_dir     <-- velocity direction
 *   shear_dir   <-- shear direction, it also contains the level of
 *                   anisotropy (Rnt = a_nt k)
 *   rcodcl  <-> boundary condition values
 *----------------------------------------------------------------------------*/

static inline void
_set_uninit_inlet_bc(cs_lnum_t   face_id,
                     double      k,
                     double      eps,
                     double      vel_dir[],
                     double      shear_dir[])
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();

  if (turb_model->itytur == 2) {

    if (_turb_bc_id.bc_k->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_k->rcodcl1[face_id] = k;

    if (_turb_bc_id.bc_eps->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

  }

  else if (turb_model->order == CS_TURB_SECOND_ORDER) {

    double d2s3 = 2./3.;
    for (int ii = 0; ii < 3; ii++)
      if (_turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id]
          > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id] = d2s3 * k;
    if (vel_dir != NULL) {
      cs_math_3_normalize(vel_dir, vel_dir);
      if (_turb_bc_id.bc_rij->rcodcl1[3*n_b_faces + face_id]
          > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_rij->rcodcl1[3*n_b_faces + face_id]
          = k * (vel_dir[0]*shear_dir[1] + vel_dir[1]*shear_dir[0]);
      if (_turb_bc_id.bc_rij->rcodcl1[4*n_b_faces + face_id]
          > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_rij->rcodcl1[4*n_b_faces + face_id]
          = k * (vel_dir[1]*shear_dir[2] + vel_dir[2]*shear_dir[1]);
      if (_turb_bc_id.bc_rij->rcodcl1[5*n_b_faces + face_id]
          > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_rij->rcodcl1[5*n_b_faces + face_id]
          = k * (vel_dir[0]*shear_dir[2] + vel_dir[2]*shear_dir[0]);
    }
    else {
      for (int ii = 3; ii < 6; ii++)
        if (_turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id]
            > 0.5*cs_math_infinite_r)
          _turb_bc_id.bc_rij->rcodcl1[ii*n_b_faces + face_id] = 0.;
    }
    if (_turb_bc_id.bc_eps->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
      if (_turb_bc_id.bc_alp_bl->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_alp_bl->rcodcl1[face_id] = 1.;

    /* Initialization of the turbulent fluxes to 0 if DFM or
     * EB-DFM are used for scalars (turbulence_flux_model = 30 or 31)
     * Alpha_theta for EB-DFM / EB-AFM / EB-GGDH is
     * initialize to 1 */

    if (_turb_bc_id.size_ut > 0) {
      for (int id = 0; id < _turb_bc_id.size_ut; id++) {
        cs_field_t *fld = _turb_bc_id.f_ut[id];
        for (int ii = 0; ii < fld->dim; ii++)
          if (fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id]
              > 0.5*cs_math_infinite_r)
            fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id] = 0.;
      }
    }

    if (_turb_bc_id.size_alp_bl_t > 0) {
      for (int id = 0; id < _turb_bc_id.size_alp_bl_t; id++) {
        cs_field_t *fld = _turb_bc_id.f_alp_bl_t[id];
        for (int ii = 0; ii < fld->dim; ii++)
          if (  fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id]
              > 0.5*cs_math_infinite_r)
            fld->bc_coeffs->rcodcl1[ii*n_b_faces + face_id] = 1.;
      }
    }

  }
  else if (turb_model->itytur == 5) {

    if (_turb_bc_id.bc_k->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_k->rcodcl1[face_id] = k;
    if (_turb_bc_id.bc_eps->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_eps->rcodcl1[face_id] = eps;

    if (_turb_bc_id.bc_phi->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_phi->rcodcl1[face_id] = 2./3.;
    if (turb_model->iturb == CS_TURB_V2F_PHI) {
      if (_turb_bc_id.bc_f_bar->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_f_bar->rcodcl1[face_id] = 0.;
    }
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      if (_turb_bc_id.bc_alp_bl->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
        _turb_bc_id.bc_alp_bl->rcodcl1[face_id] = 0.;
    }

  }
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    if (_turb_bc_id.bc_k->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_k->rcodcl1[face_id] = k;
    if (_turb_bc_id.bc_omg->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_omg->rcodcl1[face_id] = eps/cs_turb_cmu/k;

  }
  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {
    if (_turb_bc_id.bc_nusa->rcodcl1[face_id] > 0.5*cs_math_infinite_r)
      _turb_bc_id.bc_nusa->rcodcl1[face_id] = cs_turb_cmu*k*k/eps;

  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Equivalent of cs_turbulence_bc_inlet_ke_hyd_diam for Fortran calls
 * (using 1-based face number instead of id).
 *
 * parameters:
 *   face_num <-- face number
 *   uref2    <-- square of the reference flow velocity
 *   dh       <-- hydraulic diameter \f$ D_H \f$
 *   rho      <-- mass density \f$ \rho \f$
 *   mu       <-- dynamic viscosity \f$ \nu \f$
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_num,
                                  double      uref2,
                                  double      dh,
                                  double      rho,
                                  double      mu)
{
  double ustar2, k, eps;

  _ke_hyd_diam(uref2, dh, rho, mu, &ustar2, &k, &eps);

  _inlet_bc(face_num - 1, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------
 * Equivalent of cs_turbulence_bc_inlet_ke_hyd_diam for Fortran calls
 * (using 1-based face number instead of id).
 *
 * parameters:
 *   face_num    <-- face number
 *   uref2       <-- square of the reference flow velocity
 *   t_intensity <-- turbulence intensity
 *   dh          <-- hydraulic diameter \f$ D_H \f$
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_num,
                                        double      uref2,
                                        double      t_intensity,
                                        double      dh)
{
  double k, eps;

  _ke_turb_intensity(uref2, t_intensity, dh, &k, &eps);

  _inlet_bc(face_num - 1, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------
 * Equivalent of cs_turbulence_bc_inlet_ke for Fortran calls
 * (using 1-based face number instead of id).
 *
 * parameters:
 *   face_num    <-- face number
 *   k           <-- turbulent kinetic energy
 *   eps         <-- turbulent dissipation
 *   vel_dir     <-- velocity direction
 *   shear_dir   <-- shear direction
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_k_eps(cs_lnum_t   face_num,
                               double      k,
                               double      eps,
                               double      vel_dir[],
                               double      shear_dir[])
{
  _inlet_bc(face_num - 1, k, eps, vel_dir, shear_dir);
}

/*----------------------------------------------------------------------------
 * Equivalent of cs_turbulence_bc_set_uninit_inlet_ke for Fortran calls
 * (using 1-based face number instead of id).
 *
 * parameters:
 *   face_num    <-- face number
 *   k           <-- turbulent kinetic energy
 *   eps         <-- turbulent dissipation
 *   vel_dir     <-- velocity direction
 *   shear_dir   <-- shear direction
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_set_uninit_inlet_k_eps(cs_lnum_t   face_num,
                                          double      k,
                                          double      eps,
                                          double      vel_dir[],
                                          double      shear_dir[])
{
  _set_uninit_inlet_bc(face_num - 1, k, eps, vel_dir, shear_dir);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize turbulence model boundary condition pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_init_pointers(void)
{
  /* CAUTION: see note about the cs_turb_bc_id structure above. */

  const int k_turbt          = cs_field_key_id("turbulent_flux_model");
  const int k_f_turbt        = cs_field_key_id("turbulent_flux_id");
  const int k_f_turbt_alp_bl = cs_field_key_id("alpha_turbulent_flux_id");
  const int k_sca            = cs_field_key_id("scalar_id");

  if (CS_F_(k) != NULL)
    _turb_bc_id.bc_k = CS_F_(k)->bc_coeffs;
  if (CS_F_(eps) != NULL)
    _turb_bc_id.bc_eps = CS_F_(eps)->bc_coeffs;

  if (CS_F_(rij) != NULL)
    _turb_bc_id.bc_rij = CS_F_(rij)->bc_coeffs;

  if (CS_F_(phi) != NULL)
    _turb_bc_id.bc_phi = CS_F_(phi)->bc_coeffs;
  if (CS_F_(f_bar) != NULL)
    _turb_bc_id.bc_f_bar = CS_F_(f_bar)->bc_coeffs;
  if (CS_F_(alp_bl) != NULL)
    _turb_bc_id.bc_alp_bl = CS_F_(alp_bl)->bc_coeffs;

  if (CS_F_(omg) != NULL)
    _turb_bc_id.bc_omg = CS_F_(omg)->bc_coeffs;
  if (CS_F_(nusa) != NULL)
    _turb_bc_id.bc_nusa = CS_F_(nusa)->bc_coeffs;

  int n_fields = cs_field_n_fields();
  int n_sca_ut = 0;
  int n_sca_alp_bl = 0;

  /* For scalar turbulent fluxes, loop over all scalars to determine:
   *  - number of scalars  with (EB)DFM (turbulence_flux_model = 30 or 31)
   *  - number of scalars using an elliptic blending model
   *    (turbulence_flux_model = 11 or 21 or 31) */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      int s_num = cs_field_get_key_int(f, k_sca);
      if (s_num > 0) {
        int f_turbt = cs_field_get_key_int(f, k_turbt) ;
        if (f_turbt / 10 == 3)
          n_sca_ut ++;
        if (f_turbt == 11 || f_turbt == 21 || f_turbt == 31)
          n_sca_alp_bl ++;
      }
    }
  }

  _turb_bc_id.size_ut = n_sca_ut;
  _turb_bc_id.size_alp_bl_t = n_sca_alp_bl;

  if (_turb_bc_id.size_ut > 0)
    BFT_REALLOC(_turb_bc_id.f_ut, n_sca_ut, cs_field_t *);
  if (_turb_bc_id.size_alp_bl_t > 0)
    BFT_REALLOC(_turb_bc_id.f_alp_bl_t, n_sca_alp_bl, cs_field_t *);

  n_sca_ut = 0;
  n_sca_alp_bl = 0;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      int s_num = cs_field_get_key_int(f, k_sca);
      if (s_num > 0) {
        int f_turbt = cs_field_get_key_int(f, k_turbt) ;
        if (f_turbt / 10 == 3) {
          int fid_turbt = cs_field_get_key_int(f, k_f_turbt);
          _turb_bc_id.f_ut[n_sca_ut] = cs_field_by_id(fid_turbt);
          n_sca_ut ++;
        }
        if (f_turbt == 11 || f_turbt == 21 || f_turbt == 31) {
          int fid_turbt = cs_field_get_key_int(f, k_f_turbt_alp_bl);
          _turb_bc_id.f_alp_bl_t[n_sca_alp_bl] = cs_field_by_id(fid_turbt);
          n_sca_alp_bl ++;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory allocations for turbulence boundary condition pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_free_pointers(void)
{
  if (_turb_bc_id.size_ut > 0)
    BFT_FREE(_turb_bc_id.f_ut);
  if (_turb_bc_id.size_alp_bl_t > 0)
    BFT_FREE(_turb_bc_id.f_alp_bl_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of \f$ u^\star \f$, \f$ k \f$ and \f$\varepsilon \f$
 *        from a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall
 *        (use for inlet boundary conditions).
 *
 * Both \f$ u^\star \f$ and\f$ (k,\varepsilon )\f$ are returned, so that
 * the user may compute other values of \f$ k \f$ and \f$ \varepsilon \f$
 * with \f$ u^\star \f$.
 *
 * We use the laws from Idel'Cik, i.e.
 * the head loss coefficient \f$ \lambda \f$ is defined by:
 * \f[ |\dfrac{\Delta P}{\Delta x}| =
 *                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
 *
 * then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
 * \f$\lambda \f$ depends on the hydraulic Reynolds number
 * \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
 *  - for \f$ Re < 2000 \f$
 *      \f[ \lambda = \dfrac{64}{Re} \f]
 *
 *  - for \f$ Re > 4000 \f$
 *      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
 *
 *  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
 *      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
 *
 *  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
 *  from the well known formulae of developped turbulence
 *
 * \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
 * \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
 *
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 * \param[in]     rho           mass density \f$ \rho \f$
 * \param[in]     mu            dynamic viscosity \f$ \nu \f$
 * \param[out]    ustar2        square of friction speed
 * \param[out]    k             calculated turbulent intensity \f$ k \f$
 * \param[out]    eps           calculated turbulent dissipation
 *                              \f$ \varepsilon \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_ke_hyd_diam(double   uref2,
                             double   dh,
                             double   rho,
                             double   mu,
                             double  *ustar2,
                             double  *k,
                             double  *eps)
{
  _ke_hyd_diam(uref2, dh, rho, mu, ustar2, k, eps);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of \f$ k \f$ and \f$\varepsilon\f$
 *        from a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 *        and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall
 *        (for inlet boundary conditions).
 *
 * \f[ k = 1.5 I {U_{ref}}^2 \f]
 * \f[ \varepsilon = 10 \dfrac{{C_\mu}^{0.75} k^{1.5}}{ \kappa D_H} \f]
 *
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     t_intensity   turbulent intensity \f$ I \f$
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 * \param[out]    k             calculated turbulent intensity \f$ k \f$
 * \param[out]    eps           calculated turbulent dissipation
 *                               \f$ \varepsilon \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_ke_turb_intensity(double   uref2,
                                   double   t_intensity,
                                   double   dh,
                                   double  *k,
                                   double  *eps)
{
  _ke_turb_intensity(uref2, t_intensity, dh, k, eps);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall.
 *
 * We use the laws from Idel'Cik, i.e.
 * the head loss coefficient \f$ \lambda \f$ is defined by:
 * \f[ |\dfrac{\Delta P}{\Delta x}| =
 *                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
 *
 * then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
 * \f$\lambda \f$ depends on the hydraulic Reynolds number
 * \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
 *  - for \f$ Re < 2000 \f$
 *      \f[ \lambda = \dfrac{64}{Re} \f]
 *
 *  - for \f$ Re > 4000 \f$
 *      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
 *
 *  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
 *      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
 *
 *  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
 *  from the well known formulae of developped turbulence
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     uref2      square of the reference flow velocity
 * \param[in]     dh         hydraulic diameter \f$ D_H \f$
 * \param[in]     rho        mass density \f$ \rho \f$
 * \param[in]     mu         dynamic viscosity \f$ \nu \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_id,
                                double      uref2,
                                double      dh,
                                double      rho,
                                double      mu)
{
  double ustar2, k, eps;

  _ke_hyd_diam(uref2, dh, rho, mu, &ustar2, &k, &eps);

  _inlet_bc(face_id, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 *        and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall.
 *
 * \param[in]     face_id       boundary face id
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     t_intensity   turbulent intensity \f$ I \f$
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_id,
                                      double      uref2,
                                      double      t_intensity,
                                      double      dh)
{
  double k, eps;

  _ke_turb_intensity(uref2, t_intensity, dh, &k, &eps);

  _inlet_bc(face_id, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on given k and epsilon values.
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     k          turbulent kinetic energy
 * \param[in]     eps        turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_k_eps(cs_lnum_t   face_id,
                             double      k,
                             double      eps)
{
  _inlet_bc(face_id, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign homogeneous Neumann turbulent boundary conditions to
 *        a given face.
 *
 * This is useful for outgoing flow.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_set_hmg_neumann(cs_lnum_t   face_id)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_turb_model_t  *turb_model = cs_glob_turb_model;
  assert(turb_model != NULL);

  if (turb_model->itytur == 2) {

    _turb_bc_id.bc_k->icodcl[face_id] = 3;
    _turb_bc_id.bc_k->rcodcl3[face_id] = 0;

    _turb_bc_id.bc_eps->icodcl[face_id] = 3;
    _turb_bc_id.bc_eps->rcodcl3[face_id] = 0;

  }

  else if (turb_model->order == CS_TURB_SECOND_ORDER) {

    _turb_bc_id.bc_rij->icodcl[face_id] = 3;

    for (int ii = 3; ii < 6; ii++)
      _turb_bc_id.bc_rij->rcodcl3[ii*n_b_faces + face_id] = 0.;

    _turb_bc_id.bc_eps->icodcl[face_id] = 3;
    _turb_bc_id.bc_eps->rcodcl3[face_id] = 0;

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      _turb_bc_id.bc_alp_bl->icodcl[face_id] = 0.;
      _turb_bc_id.bc_alp_bl->rcodcl3[face_id] = 0.;
    }

    /* Turbulent fluxes if DFM or
     * EB-DFM are used for scalars (turbulence_flux_model = 30 or 31)
     * Alpha_theta for EB-DFM / EB-AFM / EB-GGDH */

    if (_turb_bc_id.size_ut > 0) {
      for (int id = 0; id < _turb_bc_id.size_ut; id++) {
        cs_field_t *fld = _turb_bc_id.f_ut[id];
        fld->bc_coeffs->icodcl[face_id] = 3.;
        for (int ii = 0; ii < fld->dim; ii++)
          fld->bc_coeffs->rcodcl3[ii*n_b_faces + face_id] = 0.;
      }
    }

    if (_turb_bc_id.size_alp_bl_t > 0) {
      for (int id = 0; id < _turb_bc_id.size_alp_bl_t; id++) {
        cs_field_t *fld = _turb_bc_id.f_alp_bl_t[id];
        fld->bc_coeffs->icodcl[face_id] = 3.;
        for (int ii = 0; ii < fld->dim; ii++)
          fld->bc_coeffs->rcodcl3[ii*n_b_faces + face_id] = 0.;
      }
    }

  }
  else if (turb_model->itytur == 5) {

    _turb_bc_id.bc_k->icodcl[face_id] = 3;
    _turb_bc_id.bc_k->rcodcl3[face_id] = 0;

    _turb_bc_id.bc_eps->icodcl[face_id] = 3;
    _turb_bc_id.bc_eps->rcodcl3[face_id] = 0;

    _turb_bc_id.bc_phi->rcodcl3[face_id] = 0.;
    if (turb_model->iturb == CS_TURB_V2F_PHI) {
      _turb_bc_id.bc_f_bar->icodcl[face_id] = 3.;
      _turb_bc_id.bc_f_bar->rcodcl3[face_id] = 0.;
    }
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      _turb_bc_id.bc_alp_bl->icodcl[face_id] = 3.;
      _turb_bc_id.bc_alp_bl->rcodcl3[face_id] = 0.;
    }

  }
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    _turb_bc_id.bc_k->icodcl[face_id] = 3.;
    _turb_bc_id.bc_k->rcodcl3[face_id] = 0.;

    _turb_bc_id.bc_omg->icodcl[face_id] = 3.;
    _turb_bc_id.bc_omg->rcodcl3[face_id] = 0.;

  }
  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {

    _turb_bc_id.bc_nusa->icodcl[face_id] = 3.;
    _turb_bc_id.bc_nusa->rcodcl3[face_id] = 0.;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall, only if not already set.
 *
 * Apart from assigning values where not already initialized, this function
 * is similar to \ref cs_turbulence_bc_inlet_hyd_diam.
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     vel_dir    velocity direction (not normalized, or NULL)
 * \param[in]     shear_dir  shear direction (or NULL), it also contains the
 *                           level of anisotropy (Rnt = a_nt k)
 * \param[in]     uref2      square of the reference flow velocity
 * \param[in]     dh         hydraulic diameter \f$ D_H \f$
 * \param[in]     rho        mass density \f$ \rho \f$
 * \param[in]     mu         dynamic viscosity \f$ \nu \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_set_uninit_inlet_hyd_diam(cs_lnum_t   face_id,
                                           double      vel_dir[],
                                           double      shear_dir[],
                                           double      uref2,
                                           double      dh,
                                           double      rho,
                                           double      mu)
{
  double ustar2, k, eps;

  _ke_hyd_diam(uref2, dh, rho, mu, &ustar2, &k, &eps);

  _set_uninit_inlet_bc(face_id, k, eps, vel_dir, shear_dir);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 *        and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall, only if not already set.
 *
 * Apart from assigning values where not already initialized, this function
 * is similar to \ref cs_turbulence_bc_inlet_turb_intensity.
 *
 * \param[in]     face_id       boundary face id
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     t_intensity   turbulent intensity \f$ I \f$
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_set_uninit_inlet_turb_intensity(cs_lnum_t   face_id,
                                                 double      uref2,
                                                 double      t_intensity,
                                                 double      dh)
{
  double k, eps;

  _ke_turb_intensity(uref2, t_intensity, dh, &k, &eps);

  _set_uninit_inlet_bc(face_id, k, eps, NULL, NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on given k and epsilon values only if not already initialized.
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     k          turbulent kinetic energy
 * \param[in]     eps        turbulent dissipation
 * \param[in]     vel_dir    velocity direction
 * \param[in]     shear_dir  shear direction, it also contains the level of
 *                           anisotropy (Rnt = a_nt k)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_set_uninit_inlet_k_eps(cs_lnum_t   face_id,
                                        double      k,
                                        double      eps,
                                        double      vel_dir[],
                                        double      shear_dir[])
{
  _set_uninit_inlet_bc(face_id, k, eps, vel_dir, shear_dir);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute matrix \f$\tens{\alpha}\f$ used in the computation of the
 * Reynolds stress tensor boundary conditions.
 *
 * We note \f$\tens{R}_g\f$ the Reynolds Stress tensor in the global reference
 * frame (mesh reference frame) and \f$\tens{R}_l\f$ the Reynolds stress
 * tensor in the local reference frame (reference frame associated to the
 * boundary face).
 *
 * \f$\tens{P}_{lg}\f$ is the change of basis orthogonal matrix from local
 * to global reference frame.

 * \f$\tens{\alpha}\f$ is a 6 by 6 matrix defined such that:
 * \f[
 * \vect{R}_{g,\fib} = \tens{\alpha} \vect{R}_{g,\centip} + \vect{R}_{g}^*
 * \f]
 * where symetric tensors \f$\tens{R}_g\f$ have been unfolded as follows:
 * \f[
 * \vect{R}_g = \transpose{\left(R_{g,11},R_{g,22},R_{g,33},
 *                              R_{g,12},R_{g,13},R_{g,23}\right)}
 * \f].
 *
 * \f$ \tens{R}_{g,\fib} \f$ should be computed
 * as a function of \f$\tens{R}_{g,\centip}\f$ as follows:
 * \f[
 * \tens{R}_{g,\fib}=\tens{P}_{lg}\tens{R}_{l,\fib}\transpose{\tens{P}_{lg}}
 * \f]
 *
 * with
 * \f[
 * \tens{R}_{l,\fib} =
 * \begin{bmatrix}
 * R_{l,11,\centip}   &   0        & c R_{l,13,\centip}\\
 *   0          & R_{l,22,\centip} & 0                 \\
 * c R_{l,13,\centip} & 0                & R_{l,33,\centip}
 * \end{bmatrix} +
 * \underbrace{\begin{bmatrix}
 *                 0  &   (1-c) u^* u_k        & 0                 \\
 *   (1-c) u^* u_k          & 0                & 0                 \\
 * 0                  & 0                & 0
 * \end{bmatrix}}_{\vect{R}_l^*}
 * \f]
 *
 * and
 * \f$\tens{R}_{l,\centip}=\transpose{\tens{P}_{lg}}\tens{R}_{g,\centip}
 *                       \tens{P}_{lg}\f$.
 *
 * Constant c is chosen depending on the type of the boundary face:
 * \f$c = 0\f$ at a wall face, \f$c = 1\f$ at a symmetry face.
 *
 *
 *
 * \param[in]      is_sym  Constant c in description above
 *                         (1 at a symmetry face, 0 at a wall face)
 * \param[in]      p_lg    change of basis matrix (local to global)
 * \param[out]     alpha   transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_rij_transform(int        is_sym,
                               cs_real_t  p_lg[3][3],
                               cs_real_t  alpha[6][6])
{
  cs_real_t p_lg2[3][3];
  for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
      p_lg2[ii][jj] = cs_math_pow2(p_lg[ii][jj]);

  /* alpha(i,j)  for i in [1,3] and j in [1,3]: 9 terms */
  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      alpha[jj][ii] =  p_lg2[0][ii] * p_lg2[0][jj]
                     + p_lg2[1][ii] * p_lg2[1][jj]
                     + p_lg2[2][ii] * p_lg2[2][jj]
                     + 2. * is_sym * p_lg[0][ii] * p_lg[2][ii]
                                   * p_lg[0][jj] * p_lg[2][jj];
    }
  }

  /* alpha(i,j)  for i in [1,3] and j in [4,6]: 9 terms */

  /* Correcponding line to j */
  const int _jj_to_kk[3] = {0, 1, 0};
  /* Corresponding column to j*/
  const int _jj_to_pp[3] = {1, 2, 2};

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {

      int kk = _jj_to_kk[jj];
      int pp = _jj_to_pp[jj];

      alpha[jj + 3][ii] =
        2. * (  p_lg2[0][ii] * p_lg[0][kk] * p_lg[0][pp]
              + p_lg2[1][ii] * p_lg[1][kk] * p_lg[1][pp]
              + p_lg2[2][ii] * p_lg[2][kk] * p_lg[2][pp]
              + is_sym * p_lg[2][ii] * p_lg[0][ii]
                       * (  p_lg[0][kk]*p_lg[2][pp]
                          + p_lg[2][kk]*p_lg[0][pp]));
    }
  }

  /* alpha(i,j)  for i in [4,6] and j in [1,3]: 9 terms */

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {

      int kk = _jj_to_kk[ii];
      int pp = _jj_to_pp[ii];

      /* Note: could be simplified because it is 0.5*alpha[jj+3][ii] */
      alpha[jj][ii + 3] =
          p_lg[0][kk] * p_lg[0][pp] * p_lg2[0][jj]
        + p_lg[1][kk] * p_lg[1][pp] * p_lg2[1][jj]
        + p_lg[2][kk] * p_lg[2][pp] * p_lg2[2][jj]
        + is_sym * p_lg[2][jj] * p_lg[0][jj]
                 * (  p_lg[0][kk]*p_lg[2][pp]
                    + p_lg[2][kk]*p_lg[0][pp]);
    }
  }

  /* alpha(i,j)  for i in [4,6] and j in [4,6]: 9 terms */
  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {

      int kk = _jj_to_kk[ii];
      int pp = _jj_to_pp[ii];

      int jj1 = _jj_to_kk[jj];
      int jj2 = _jj_to_pp[jj];

      alpha[jj + 3][ii + 3] =
        2. * (  p_lg[0][kk] * p_lg[0][pp] * p_lg[0][jj1] * p_lg[0][jj2]
              + p_lg[1][kk] * p_lg[1][pp] * p_lg[1][jj1] * p_lg[1][jj2]
              + p_lg[2][kk] * p_lg[2][pp] * p_lg[2][jj1] * p_lg[2][jj2])
              + is_sym * (  p_lg[0][kk]*p_lg[2][pp]
                          + p_lg[2][kk]*p_lg[0][pp])
                       * (  p_lg[2][jj1]*p_lg[0][jj2]
                          + p_lg[0][jj1]*p_lg[2][jj2]);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
