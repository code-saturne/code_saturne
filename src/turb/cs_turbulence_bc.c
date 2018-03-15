/*============================================================================
 * Base turbulence boundary conditions.
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

  int  k;      /* variable id for k */
  int  eps;    /* variable id for epsilon */

  int  r11;    /* variable id for R_xx */
  int  r22;    /* variable id for R_yy */
  int  r33;    /* variable id for R_zz */
  int  r12;    /* variable id for R_xy */
  int  r23;    /* variable id for R_yz */
  int  r13;    /* variable id for R_xz */
  int  rij;    /* variable id for R_ij tensor (irijco=1) */

  int  phi;    /* variable id for phi */
  int  f_bar;  /* variable id for f_bar */
  int  alpha;  /* variable id for alpha */

  int  omg;    /* variable id for omega */
  int  nusa;   /* variable id for nu_t (SA model) */

} cs_turb_bc_id_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Variable ids for main turbulence variables */

static cs_turb_bc_id_t
_turb_bc_id =
{
  -1, /* k */
  -1, /* eps */

  -1, /* r11 */
  -1, /* r22 */
  -1, /* r33 */
  -1, /* r12 */
  -1, /* r23 */
  -1, /* r13 */
  -1, /* rij */

  -1, /* phi */
  -1, /* f_bar */
  -1, /* alpha */

  -1, /* omg */
  -1 /* nusa */
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
                                  double      mu,
                                  double     *rcodcl);

void
cs_f_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_num,
                                        double      uref2,
                                        double      t_intensity,
                                        double      dh,
                                        double     *rcodcl);

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
 *   face_id <-- face id
 *   k       <-- k
 *   eps     <-- epsilon
 *   rcodcl  <-> boundary condition values
 *----------------------------------------------------------------------------*/

static inline void
_inlet_bc(cs_lnum_t   face_id,
          double      k,
          double      eps,
          double     *rcodcl)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  if (cs_glob_turb_model->itytur == 2) {

    if (rcodcl[_turb_bc_id.k  *n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.k  *n_b_faces + face_id] = k;

    if (rcodcl[_turb_bc_id.eps*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.eps*n_b_faces + face_id] = eps;

  }

  else if (cs_glob_turb_model->itytur == 3) {

    double d2s3 = 2./3.;
    if (_turb_bc_id.rij == -1) {
      if (rcodcl[_turb_bc_id.r11*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r11*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[_turb_bc_id.r22*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r22*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[_turb_bc_id.r33*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r33*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[_turb_bc_id.r12*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r12*n_b_faces + face_id] = 0.;
      if (rcodcl[_turb_bc_id.r13*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r13*n_b_faces + face_id] = 0.;
      if (rcodcl[_turb_bc_id.r23*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.r23*n_b_faces + face_id] = 0.;
      if (rcodcl[_turb_bc_id.eps*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.eps*n_b_faces + face_id] = eps;
    }
    else {
      if (rcodcl[_turb_bc_id.rij*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.rij*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[(_turb_bc_id.rij + 1)*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[(_turb_bc_id.rij + 1)*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[(_turb_bc_id.rij + 2)*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[(_turb_bc_id.rij + 2)*n_b_faces + face_id] = d2s3 * k;
      if (rcodcl[(_turb_bc_id.rij + 3)*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[(_turb_bc_id.rij + 3)*n_b_faces + face_id] = 0.;
      if (rcodcl[(_turb_bc_id.rij + 4)*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[(_turb_bc_id.rij + 4)*n_b_faces + face_id] = 0.;
      if (rcodcl[(_turb_bc_id.rij + 5)*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[(_turb_bc_id.rij + 5)*n_b_faces + face_id] = 0.;
      if (rcodcl[_turb_bc_id.eps*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.eps*n_b_faces + face_id] = eps;
    }

    if (cs_glob_turb_model->iturb == 32)
      if (rcodcl[_turb_bc_id.alpha*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.alpha*n_b_faces + face_id] = 1.;

  }
  else if (cs_glob_turb_model->itytur == 5) {

    if (rcodcl[_turb_bc_id.k  *n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.k  *n_b_faces + face_id] = k;
    if (rcodcl[_turb_bc_id.eps*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.eps*n_b_faces + face_id] = eps;

    if (rcodcl[_turb_bc_id.phi*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.phi*n_b_faces + face_id] = 2./3.;
    if (cs_glob_turb_model->iturb == 50)
      if (rcodcl[_turb_bc_id.f_bar*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
        rcodcl[_turb_bc_id.f_bar*n_b_faces + face_id] = 0.;
    else if (cs_glob_turb_model->iturb == 51)
        if (rcodcl[_turb_bc_id.alpha*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
          rcodcl[_turb_bc_id.alpha*n_b_faces + face_id] = 0.;

  }
  else if (cs_glob_turb_model->itytur == 6) {

    if (rcodcl[_turb_bc_id.k  *n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.k  *n_b_faces + face_id] = k;
    if (rcodcl[_turb_bc_id.omg*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.omg*n_b_faces + face_id] = eps/cs_turb_cmu/k;

  }
  else if (cs_glob_turb_model->itytur == 7) {
    if (rcodcl[_turb_bc_id.nusa*n_b_faces + face_id] > 0.5*cs_math_infinite_r)
      rcodcl[_turb_bc_id.nusa*n_b_faces + face_id] = cs_turb_cmu*k*k/eps;

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
 *   rcodcl   <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_num,
                                  double      uref2,
                                  double      dh,
                                  double      rho,
                                  double      mu,
                                  double     *rcodcl)
{
  double ustar2, k, eps;

  _ke_hyd_diam(uref2, dh, rho, mu, &ustar2, &k, &eps);

  _inlet_bc(face_num - 1, k, eps, rcodcl);
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
 *   rcodcl      <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_num,
                                        double      uref2,
                                        double      t_intensity,
                                        double      dh,
                                        double     *rcodcl)
{
  double k, eps;

  _ke_turb_intensity(uref2, t_intensity, dh, &k, &eps);

  _inlet_bc(face_num - 1, k, eps, rcodcl);
}

/*----------------------------------------------------------------------------
 * Equivalent of cs_turbulence_bc_inlet_ke for Fortran calls
 * (using 1-based face number instead of id).
 *
 * parameters:
 *   face_num    <-- face number
 *   uref2       <-- square of the reference flow velocity
 *   t_intensity <-- turbulence intensity
 *   dh          <-- hydraulic diameter \f$ D_H \f$
 *   rcodcl      <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_f_turbulence_bc_inlet_k_eps(cs_lnum_t   face_num,
                               double      k,
                               double      eps,
                               double     *rcodcl)
{
  _inlet_bc(face_num - 1, k, eps, rcodcl);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize turbulence model boundary condition ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_model_init_bc_ids(void)
{
  /* CAUTION: see note about the cs_turb_bc_id structure above. */

  const int var_key_id = cs_field_key_id("variable_id");

  if (CS_F_(k) != NULL)
    _turb_bc_id.k = cs_field_get_key_int(CS_F_(k), var_key_id) -1;
  if (CS_F_(eps) != NULL)
    _turb_bc_id.eps = cs_field_get_key_int(CS_F_(eps), var_key_id) -1;

  if (CS_F_(r11) != NULL)
    _turb_bc_id.r11 = cs_field_get_key_int(CS_F_(r11), var_key_id) -1;
  if (CS_F_(r22) != NULL)
    _turb_bc_id.r22 = cs_field_get_key_int(CS_F_(r22), var_key_id) -1;
  if (CS_F_(r33) != NULL)
    _turb_bc_id.r33 = cs_field_get_key_int(CS_F_(r33), var_key_id) -1;
  if (CS_F_(r12) != NULL)
    _turb_bc_id.r12 = cs_field_get_key_int(CS_F_(r12), var_key_id) -1;
  if (CS_F_(r23) != NULL)
    _turb_bc_id.r23 = cs_field_get_key_int(CS_F_(r23), var_key_id) -1;
  if (CS_F_(r13) != NULL)
    _turb_bc_id.r13 = cs_field_get_key_int(CS_F_(r13), var_key_id) -1;
  if (CS_F_(rij) != NULL)
    _turb_bc_id.rij = cs_field_get_key_int(CS_F_(rij), var_key_id) -1;

  if (CS_F_(phi) != NULL)
    _turb_bc_id.phi = cs_field_get_key_int(CS_F_(phi), var_key_id) -1;
  if (CS_F_(f_bar) != NULL)
    _turb_bc_id.f_bar = cs_field_get_key_int(CS_F_(f_bar), var_key_id) -1;
  if (CS_F_(alpha) != NULL)
    _turb_bc_id.alpha = cs_field_get_key_int(CS_F_(alpha), var_key_id) -1;

  if (CS_F_(omg) != NULL)
    _turb_bc_id.omg = cs_field_get_key_int(CS_F_(omg), var_key_id) -1;
  if (CS_F_(nusa) != NULL)
    _turb_bc_id.nusa = cs_field_get_key_int(CS_F_(nusa), var_key_id) -1;
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
 * \param[out]    rcodcl     boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_id,
                                double      uref2,
                                double      dh,
                                double      rho,
                                double      mu,
                                double     *rcodcl)
{
  double ustar2, k, eps;

  _ke_hyd_diam(uref2, dh, rho, mu, &ustar2, &k, &eps);

  _inlet_bc(face_id, k, eps, rcodcl);
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
 * \param[out]    rcodcl        boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_id,
                                      double      uref2,
                                      double      t_intensity,
                                      double      dh,
                                      double     *rcodcl)
{
  double k, eps;

  _ke_turb_intensity(uref2, t_intensity, dh, &k, &eps);

  _inlet_bc(face_id, k, eps, rcodcl);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
