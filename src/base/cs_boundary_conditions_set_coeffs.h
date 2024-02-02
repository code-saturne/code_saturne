#ifndef __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__
#define __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__

/*============================================================================
 * Boundary condition management.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "bft_error.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_boundary_conditions_set_coeffs.c
 *
 * \brief Translation of the boundary conditions given by the user in a form
 *        that fits to the solver.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Translation of the boundary conditions given by the user
 * in a form that fits to the solver.
 *
 * The values at a boundary face \f$ \fib \f$ stored in the face center
 * \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
 * are written as:
 * \f[
 * P_{\face} = A_P^g + B_P^g P_{\centi}
 * \f]
 * and
 * \f[
 * Q_{\face} = A_P^f + B_P^f P_{\centi}
 * \f]
 * where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
 * neighboring cell.
 *
 * \warning
 * - If we consider an increment of a variable, the boundary conditions
 *   read:
 *   \f[
 *   \delta P_{\face} = B_P^g \delta P_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \delta Q_{\face} = -B_P^f \delta P_{\centi}
 *   \f]
 *
 * - For a vector field such as the velocity \f$ \vect{u} \f$ the boundary
 *   conditions may read:
 *   \f[
 *   \vect{u}_{\face} = \vect{A}_u^g + \tens{B}_u^g \vect{u}_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \vect{Q}_{\face} = \vect{A}_u^f + \tens{B}_u^f \vect{u}_{\centi}
 *   \f]
 *   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
 *   which coupled velocity components next to a boundary.
 *
 * Please refer to the
 * <a href="../../theory.pdf#boundary"><b>boundary conditions</b></a> section
 * of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#condli"><b>condli</b></a> section.
 *
 * \param[in]     nvar          total number of variables
 * \param[in]     iterns        iteration number on Navier-Stokes equations
 * \param[in]     isvhb         indicator to save exchange coeffient
 *                               at the walls
 * \param[in]     itrale        ALE iteration number
 * \param[in]     italim        for ALE
 * \param[in]     itrfin        for ALE
 * \param[in]     ineefl        for ALE
 * \param[in]     itrfup        for ALE
 * \param[in,out] isostd        indicator for standard outlet
 *                              and reference face index
 * \param[in]     dt            time step (per cell)
 * \param[out]    visvdr        dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]    hbord         exchange coefficient at boundary
 * \param[out]    theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 * \param[in]     nftcdt        Global indicator of condensation source terms
 *                              (ie. sum on the processors of nfbpcd) cells
 *                              associated to the face with condensation
 *                              phenomenon
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs(int        nvar,
                                  int        iterns,
                                  int        isvhb,
                                  int        itrale,
                                  int        italim,
                                  int        itrfin,
                                  int        ineefl,
                                  int        itrfup,
                                  int        isostd[],
                                  cs_real_t  visvdr[],
                                  cs_real_t  hbord[],
                                  cs_real_t  theipb[],
                                  int        nftcdt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialization of boundary condition arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Neumann BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   qimp   flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_scalar(cs_real_t  *a,
                                          cs_real_t  *af,
                                          cs_real_t  *b,
                                          cs_real_t  *bf,
                                          cs_real_t   qimp,
                                          cs_real_t   hint)
{
  /* Gradient BCs */
  *a = -qimp/cs_math_fmax(hint, 1.e-300);
  *b = 1.;

  /* Flux BCs */
  *af = qimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set homogeneous Neumann BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_scalar_hmg(cs_real_t  *a,
                                              cs_real_t  *af,
                                              cs_real_t  *b,
                                              cs_real_t  *bf)
{
  /* Gradient BCs */
  *a = 0.;
  *b = 1.;

  /* Flux BCs */
  *af = 0.;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Neumann BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   qimpv  flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector(cs_real_t        a[3],
                                          cs_real_t        af[3],
                                          cs_real_t        b[3][3],
                                          cs_real_t        bf[3][3],
                                          const cs_real_t  qimpv[3],
                                          cs_real_t        hint)
{
  /* Gradient BCs */

  for (size_t i = 0; i < 3; i++) {
    a[i] = -qimpv[i] / fmax(hint, 1.e-300);
  }

  b[0][0] = 1., b[0][1] = 0., b[0][2] = 0.;
  b[1][0] = 0., b[1][1] = 1., b[1][2] = 0.;
  b[2][0] = 0., b[2][1] = 0., b[2][2] = 1.;

  /* Flux BCs */

  for (size_t i = 0; i < 3; i++) {
    af[i] = qimpv[i];

    for (size_t j = 0; j < 3; j++)
      bf[i][j] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set neumann BC for an anisotropic vector for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   qimpv  flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector_aniso(cs_real_t        a[3],
                                                cs_real_t        af[3],
                                                cs_real_t        b[3][3],
                                                cs_real_t        bf[3][3],
                                                const cs_real_t  qimpv[3],
                                                const cs_real_t  hint[6])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};
  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  cs_real_t invdet = 1./(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  /* Gradient BCs */
  cs_math_sym_33_3_product(invh, qimpv, a);
  for (int isou = 0; isou < 3; isou++)
    a[isou] = -a[isou];

  b[0][0] = 1.0, b[0][1] = 0.0, b[0][2] = 0.0;
  b[1][0] = 0.0, b[1][1] = 1.0, b[1][2] = 0.0;
  b[2][0] = 0.0, b[2][1] = 0.0, b[2][2] = 1.0;

  for (int isou = 0; isou < 3; isou++) {

    /* Flux BCs */
    af[isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief  Set Neumann boundary conditions for a tensor for a given face.
 *
 * \param[out]  a        explicit BC coefficient for gradients
 * \param[out]  af       explicit BC coefficient for diffusive flux
 * \param[out]  b        implicit BC coefficient for gradients
 * \param[out]  bf       implicit BC coefficient for diffusive flux
 * \param[in]   qimpts   flux value to impose
 * \param[in]   hint     internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_tensor(cs_real_t        a[6],
                                          cs_real_t        af[6],
                                          cs_real_t        b[6][6],
                                          cs_real_t        bf[6][6],
                                          const cs_real_t  qimpts[6],
                                          cs_real_t        hint)
{
  for (int isou = 0; isou < 6; isou++) {

    /* Gradient BC */
    a[isou] = -qimpts[isou]/cs_math_fmax(hint, 1.e-300);
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        b[isou][jsou] = 1.0;
      else
        b[isou][jsou] = 0.0;
    }

    /* Flux BCs */
    af[isou] = qimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   dirichlet value to impose
 * \param[in]   hint   internal exchange coefficient
 * \param[in]   hext   external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_scalar(cs_real_t  *a,
                                            cs_real_t  *af,
                                            cs_real_t  *b,
                                            cs_real_t  *bf,
                                            cs_real_t   pimp,
                                            cs_real_t   hint,
                                            cs_real_t   hext)
{
  //if (fabs(hext) > cs_math_infinite_r*0.5) {
  if (hext < 0.) {

    /* Gradient BCs */
    *a = pimp;
    *b = 0.;

    /* Flux BCs */
    *af = -hint*pimp;
    *bf =  hint;

  }
  else {

    /* Gradient BCs */
    *a = hext*pimp/(hint + hext);
    *b = hint     /(hint + hext);

    /* Flux BCs */
    cs_real_t heq = hint*hext/(hint + hext);
    *af = -heq*pimp;
    *bf =  heq;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a vector for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimpv  dirichlet value to impose
 * \param[in]   hint   internal exchange coefficient
 * \param[in]   hextv  external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector(cs_real_t        a[3],
                                            cs_real_t        af[3],
                                            cs_real_t        b[3][3],
                                            cs_real_t        bf[3][3],
                                            const cs_real_t  pimpv[3],
                                            cs_real_t        hint,
                                            const cs_real_t  hextv[3])
{
  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[isou][jsou] = 0.;

      /* Flux BCs */
      af[isou] = -hint*pimpv[isou];

      bf[0][0] = hint, bf[0][1] = 0.,   bf[0][2] = 0.;
      bf[1][0] = 0.,   bf[1][1] = hint, bf[1][2] = 0.;
      bf[2][0] = 0.,   bf[2][1] = 0.,   bf[2][2] = hint;

    }
    else {

      const cs_real_t val = hint/(hint + hextv[isou]);
      const cs_real_t heq = hextv[isou]*val;

      /* Gradient BCs */
      a[isou] = hextv[isou]*pimpv[isou]/(hint + hextv[isou]);

      b[0][0] = val, b[0][1] = 0.,  b[0][2] = 0.;
      b[1][0] = 0.,  b[1][1] = val, b[1][2] = 0.;
      b[2][0] = 0.,  b[2][1] = 0.,  b[2][2] = val;

      /* Flux BCs */
      af[isou] = -heq*pimpv[isou];

      bf[0][0] = heq, bf[0][1] = 0.,  bf[0][2] = 0.;
      bf[1][0] = 0.,  bf[1][1] = heq, bf[1][2] = 0.;
      bf[2][0] = 0.,  bf[2][1] = 0.,  bf[2][2] = heq;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Set convective oulet BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   flux value to impose
 * \param[in]   cfl    local Courant number used to convect
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t  *a,
                                                    cs_real_t  *af,
                                                    cs_real_t  *b,
                                                    cs_real_t  *bf,
                                                    cs_real_t   pimp,
                                                    cs_real_t   cfl,
                                                    cs_real_t   hint);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a vector for a given face with left anisotropic
 *        diffusion.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimpv  dirichlet value to impose
 * \param[in]   hintt  internal exchange coefficient
 * \param[in]   hextv  external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector_aniso(cs_real_t        a[3],
                                                  cs_real_t        af[3],
                                                  cs_real_t        b[3][3],
                                                  cs_real_t        bf[3][3],
                                                  const cs_real_t  pimpv[3],
                                                  const cs_real_t  hintt[6],
                                                  const cs_real_t  hextv[3])
{
  /* Gradient BCs */
  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {
      a[isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[isou][jsou] = 0.;
    }
    else {
      /* FIXME: at least log error message */
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: hextv not set for component %d."),
                __func__, isou);
    }
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, pimpv, af);
  for (int isou = 0; isou < 3; isou++)
    af[isou] = -af[isou];

  bf[0][0] = hintt[0];
  bf[1][1] = hintt[1];
  bf[2][2] = hintt[2];
  bf[0][1] = hintt[3];
  bf[1][0] = hintt[3];
  bf[1][2] = hintt[4];
  bf[2][1] = hintt[4];
  bf[0][2] = hintt[5];
  bf[2][0] = hintt[5];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a tensor for a given face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpts  Dirichlet value to impose
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   hextts  external exchange coefficient (10^30 by default)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_tensor(cs_real_t        a[6],
                                            cs_real_t        af[6],
                                            cs_real_t        b[6][6],
                                            cs_real_t        bf[6][6],
                                            const cs_real_t  pimpts[6],
                                            cs_real_t        hint,
                                            const cs_real_t  hextts[6])
{
  for (int isou = 0; isou < 6; isou++) {

    if (fabs(hextts[isou]) > 0.5*cs_math_infinite_r) {
      /* Gradient BCs */
      a[isou] = pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++)
        b[isou][jsou] = 0.;

      /* Flux BCs */
      af[isou] = -hint * pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          bf[isou][jsou] = hint;
        else
          bf[isou][jsou] = 0.;
      }
    }
    else {

      const cs_real_t heq = hint * hextts[isou] / (hint + hextts[isou]);

      /* Gradient BCs */
      a[isou] = hextts[isou] * pimpts[isou] / (hint + hextts[isou]);
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          b[isou][jsou] = hint / (hint + hextts[isou]);
        else
          b[isou][jsou] = 0.;
      }

      /* Flux BCs */
      af[isou] = -heq * pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          bf[isou][jsou] = heq;
        else
          bf[isou][jsou] = 0.;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for a symmetric vector for a given face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose on the normal component
 * \param[in]   qimpv   flux value to impose on the tangential components
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   normal  normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_sym_vector(cs_real_t        a[3],
                                                  cs_real_t        af[3],
                                                  cs_real_t        b[3][3],
                                                  cs_real_t        bf[3][3],
                                                  const cs_real_t  pimpv[3],
                                                  const cs_real_t  qimpv[3],
                                                  cs_real_t        hint,
                                                  const cs_real_t  normal[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    a[isou] = - qimpv[isou]/cs_math_fmax(hint, 1.e-300);
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++) {

      a[isou] = a[isou] + normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      if (jsou == isou)
        b[isou][jsou] = 1.0 - normal[isou] * normal[jsou];
      else
        b[isou][jsou] = - normal[isou] * normal[jsou];
    }

    /* Flux BCs */
    af[isou] = qimpv[isou];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++){

      af[isou] = af[isou] - normal[isou]*normal[jsou]
                  * (hint * pimpv[jsou] + qimpv[jsou]);

      bf[isou][jsou] = hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose on the normal component
 * \param[in]   qimpv   flux value to impose on the tangential components
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   normal  normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for a vector for a given face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose on the tangential components
 * \param[in]   qimpv   flux value to impose on the normal component
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   normal  normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_dirichlet_vector
  (cs_real_t         a[3],
   cs_real_t         af[3],
   cs_real_t         b[3][3],
   cs_real_t         bf[3][3],
   const cs_real_t   pimpv[3],
   const cs_real_t   qimpv[3],
   cs_real_t         hint,
   const cs_real_t   normal[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BC*/
    /* "[1 -n(x)n] Pimp" is divided into two */
    a[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      a[isou] = a[isou] - normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      b[isou][jsou] = normal[isou] * normal[jsou];
    }

    /* Flux BC */
    /* "[1 -n(x)n] Pimp" is divided into two */
    af[isou] = -hint*pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      af[isou] = af[isou] + normal[isou]*normal[jsou]
        * (qimpv[jsou] + pimpv[jsou] * hint);

      if (jsou == isou)
        bf[isou][jsou] = hint * (1.0 - normal[isou] * normal[jsou]);
      else
        bf[isou][jsou] = - hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose on the tangential components
 * \param[in]   qimpv   flux value to impose on the normal component
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   normal  normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for a vector for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose
 * \param[in]   cflv    local Courant number used to convect
 * \param[in]   hint    internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector(cs_real_t        a[3],
                                                    cs_real_t        af[3],
                                                    cs_real_t        b[3][3],
                                                    cs_real_t        bf[3][3],
                                                    const cs_real_t  pimpv[3],
                                                    const cs_real_t  cflv[3],
                                                    cs_real_t        hint)
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou ++) {
      if (jsou == isou)
        b[isou][jsou] = cflv[isou] / (1.0 + cflv[isou]);
      else
        b[isou][jsou] = 0.0;
    }
    a[isou] = pimpv[isou] * (1.0 - b[isou][isou]);

    /* Flux BCs */
    af[isou] = -hint * a[isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        bf[isou][jsou] = hint * (1.0 - b[isou][jsou]);
    else
      bf[isou][jsou] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective outlet BC for a tensor for a given face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpts  Dirichlet value to impose
 * \param[in]   cflts   local Courant number used to convect
 * \param[in]   hint    internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_tensor(cs_real_t        a[6],
                                                    cs_real_t        af[6],
                                                    cs_real_t        b[6][6],
                                                    cs_real_t        bf[6][6],
                                                    const cs_real_t  pimpts[6],
                                                    const cs_real_t  cflts[6],
                                                    cs_real_t        hint)
{
  for (int isou = 0; isou < 6; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        b[isou][jsou] = cflts[isou] / (1.0 + cflts[isou]);
      else
        b[isou][jsou] = 0.0;
    }
    a[isou] = (1.0 - b[isou][isou]) * pimpts[isou];

    /* Flux BCs */
    af[isou] = -hint*a[isou];
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        bf[isou][jsou] = hint * (1.0 - b[isou][jsou]);
      else
        bf[isou][jsou] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for an anisotropic vector for a given face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose
 * \param[in]   cflv    local Courant number used to convect
 * \param[in]   hintt   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector_aniso
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  cflv[3],
   const cs_real_t  hintt[6])
{
  for(int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        b[isou][jsou] = cflv[isou]/(1.0+cflv[isou]);
      else
        b[isou][jsou] = 0.0;
    }
    a[isou] = (1.0-b[isou][isou])*pimpv[isou];

  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, a, af);
  for (int isou = 0; isou < 3; isou++)
    af[isou] = -af[isou];

  bf[0][0] = hintt[0]*(1.0 - b[0][0]);
  bf[1][1] = hintt[1]*(1.0 - b[1][1]);
  bf[2][2] = hintt[2]*(1.0 - b[2][2]);
  bf[0][1] = hintt[3]*(1.0 - b[0][0]);
  bf[1][0] = hintt[3]*(1.0 - b[0][0]);
  bf[1][2] = hintt[4]*(1.0 - b[1][1]);
  bf[2][1] = hintt[4]*(1.0 - b[1][1]);
  bf[0][2] = hintt[5]*(1.0 - b[2][2]);
  bf[2][0] = hintt[5]*(1.0 - b[2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set BC for an affine scalar function for a given face.
 *
 * \param[out]    a       explicit BC coefficient for gradients
 * \param[out]    af      explicit BC coefficient for diffusive flux
 * \param[out]    b       implicit BC coefficient for gradients
 * \param[out]    bf      implicit BC coefficient for diffusive flux
 * \param[in]     pinf    affine part
 * \param[in]     ratio   linear part
 * \param[in]     hint    internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_scalar(cs_real_t  *a,
                                                  cs_real_t  *af,
                                                  cs_real_t  *b,
                                                  cs_real_t  *bf,
                                                  cs_real_t   pinf,
                                                  cs_real_t   ratio,
                                                  cs_real_t   hint)
{
  /* Gradient BCs */
  *b = ratio;
  *a = pinf;

  /* Flux BCs */
  *af = -hint * *a;
  *bf =  hint * (1. - *b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, zero flux for diffusion.
 *
 * \param[out]    a       explicit BC coefficient for gradients
 * \param[out]    af      explicit BC coefficient for diffusive flux
 * \param[out]    b       implicit BC coefficient for gradients
 * \param[out]    bf      implicit BC coefficient for diffusive flux
 * \param[in]     dimp    flux value to impose
 * \param[in]     hint    internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_conv_h_neumann_diff_scalar(cs_real_t  *a,
                                                              cs_real_t  *af,
                                                              cs_real_t  *b,
                                                              cs_real_t  *bf,
                                                              cs_real_t   dimp,
                                                              cs_real_t   hint)

{
  /* Gradient BCs */
  cs_boundary_conditions_set_neumann_scalar(&*a,
                                            &*af,
                                            &*b,
                                            &*bf,
                                            dimp,
                                            hint);

  /* Flux BCs */
  *af = 0.;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, imposed flux for
 *         diffusion.
 *
 * \param[out]    a         explicit BC coefficient for gradients
 * \param[out]    af        explicit BC coefficient for diffusive flux
 * \param[out]    b         implicit BC coefficient for gradients
 * \param[out]    bf        implicit BC coefficient for diffusive flux
 * \param[in]     pinf      affine part
 * \param[in]     ratio     linear part
 * \param[in]     dimp      flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_conv_neumann_diff_scalar
  (cs_real_t  *a,
   cs_real_t  *af,
   cs_real_t  *b,
   cs_real_t  *bf,
   cs_real_t   pinf,
   cs_real_t   ratio,
   cs_real_t   dimp)
{
  /* Gradient BCs */
  *b = ratio;
  *a = pinf;

  /* Flux BCs */
  *af = dimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set total flux as a Robin condition.
 *
 * \param[out]    a         explicit BC coefficient for gradients
 * \param[out]    af        explicit BC coefficient for diffusive flux
 * \param[out]    b         implicit BC coefficient for gradients
 * \param[out]    bf        implicit BC coefficient for diffusive flux
 * \param[in]     hext      convective flux to be imposed
 * \param[in]     dimp      flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_total_flux(cs_real_t *a,
                                      cs_real_t *af,
                                      cs_real_t *b,
                                      cs_real_t *bf,
                                      cs_real_t hext,
                                      cs_real_t dimp)
{
  /* Gradients BCs */
  *a = 0.;
  *b = 1.;

  /* Flux BCs */
  *af = dimp;
  *bf = hext;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a scalar.
 *
 * \param[out]    a        explicit BC coefficient for gradients
 * \param[out]    af       explicit BC coefficient for diffusive flux
 * \param[out]    b        implicit BC coefficient for gradients
 * \param[out]    bf       implicit BC coefficient for diffusive flux
 * \param[in]     pimp     Dirichlet value to impose
 * \param[in]     dimp     flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar(cs_real_t *a,
                                                              cs_real_t *af,
                                                              cs_real_t *b,
                                                              cs_real_t *bf,
                                                              cs_real_t pimp,
                                                              cs_real_t dimp)
{
  /* Gradients BC */
  *a = pimp;
  *b = 0.;

  /* Flux BCs */
  *af = dimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a vector.
 *
 * \param[out]    a         explicit BC coefficient for gradients
 * \param[out]    af        explicit BC coefficient for diffusive flux
 * \param[out]    b         implicit BC coefficient for gradients
 * \param[out]    bf        implicit BC coefficient for diffusive flux
 * \param[in]     pimpv     Dirichlet value to impose
 * \param[in]     qimpv     flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    a[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      b[isou][jsou] = 0.0;

    /* Flux BCs */
    af[isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a tensor
 *
 * \param[out]    a         explicit BC coefficient for gradients
 * \param[out]    af        explicit BC coefficient for diffusive flux
 * \param[out]    b         implicit BC coefficient for gradients
 * \param[out]    bf        implicit BC coefficient for diffusive flux
 * \param[in]     pimpts    Dirichlet value to impose
 * \param[in]     qimpts    flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_tensor
  (cs_real_t        a[6],
   cs_real_t        af[6],
   cs_real_t        b[6][6],
   cs_real_t        bf[6][6],
   const cs_real_t  pimpts[6],
   const cs_real_t  qimpts[6])
{
  for (int isou = 0; isou < 6; isou++) {

    /* BS test sur hextv ? if (abs(hextv[isou])>rinfin*0.5d0) then */

    /* Gradient BCs */
    a[isou] = pimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      b[isou][jsou] = 0.0;

    /* Flux BCs */
    af[isou] = qimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__ */
