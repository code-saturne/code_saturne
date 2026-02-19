#ifndef __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__
#define __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__

/*============================================================================
 * Translation of the boundary conditions given by the user in a form
 * that fits the solver.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include "base/cs_profiling.h"
#include "bft/bft_error.h"
#include "base/cs_field.h"
#include "base/cs_math.h"
#include "base/cs_profiling.h"
#include "cdo/cs_equation_param.h"

#ifdef __cplusplus
#include "base/cs_dispatch.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

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
 * \param[in]       nvar        total number of variables
 * \param[in]       iterns      iteration number on Navier-Stokes equations
 * \param[in]       isvhb       id of field whose exchange coeffient should be
 *                              saved at the walls, or -1.
 * \param[in]       italim      for ALE
 * \param[in]       itrfin      Last velocity-pressure sub-iteration indicator
 * \param[in]       ineefl      for ALE
 * \param[in]       itrfup      Update after velocity-pressure sub-iterations
 * \param[in, out]  isostd      indicator for standard outlet
 *                              and reference face index
 * \param[in]       dt          time step (per cell)
 * \param[in, out]  visvdr      dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]      hbord       exchange coefficient at boundary
 * \param[out]      theipb      value of thermal scalar at \f$ \centip \f$
 *                              f boundary cells
 * \param[in]       nftcdt      Global indicator of condensation source terms
 *                              (ie. sum on the processors of nfbpcd) cells
 *                              associated to the face with condensation
 *                              phenomenon
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs(int         nvar,
                                  int         iterns,
                                  int         isvhb,
                                  int         italim,
                                  int         itrfin,
                                  int         ineefl,
                                  int         itrfup,
                                  int         isostd[],
                                  cs_real_t  *visvdr,
                                  cs_real_t   hbord[],
                                  cs_real_t   theipb[],
                                  int         nftcdt);

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
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   qimpv       flux value to impose
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_vector(cs_real_t             a[3],
                                          cs_real_t             af[3],
                                          cs_real_t             b[3][3],
                                          cs_real_t             bf[3][3],
                                          const cs_real_t       qimpv[3],
                                          cs_real_t             hint)
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
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   qimpv       flux value to impose
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_vector_aniso(cs_real_t             a[3],
                                                cs_real_t             af[3],
                                                cs_real_t             b[3][3],
                                                cs_real_t             bf[3][3],
                                                const cs_real_t       qimpv[3],
                                                const cs_real_t       hint[6])
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
  for (cs_lnum_t i = 0; i < 3; i++)
    a[i] = -a[i];

  b[0][0] = 1.0, b[0][1] = 0.0, b[0][2] = 0.0;
  b[1][0] = 0.0, b[1][1] = 1.0, b[1][2] = 0.0;
  b[2][0] = 0.0, b[2][1] = 0.0, b[2][2] = 1.0;

  for (cs_lnum_t i = 0; i < 3; i++) {
    /* Flux BCs */
    af[i] = qimpv[i];
    for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
      bf[i][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief  Set Neumann boundary conditions for a tensor for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  fb          implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   qimpts      flux value to impose
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_tensor(cs_real_t        a[6],
                                          cs_real_t        af[6],
                                          cs_real_t        b[6][6],
                                          cs_real_t        bf[6][6],
                                          const cs_real_t  qimpts[6],
                                          cs_real_t        hint)
{
  for (int i = 0; i < 6; i++) {

    /* Gradient BC */
    a[i] = -qimpts[i]/cs::max(hint, 1.e-300);
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == i)
        b[i][jsou] = 1.0;
      else
        b[i][jsou] = 0.0;
    }

    /* Flux BCs */
    af[i] = qimpts[i];
    for (int jsou = 0; jsou < 6; jsou++)
      bf[i][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a vector for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  fb          implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   hextv       external exchange coefficient
 *                          (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_vector(cs_real_t             a[3],
                                            cs_real_t             af[3],
                                            cs_real_t             b[3][3],
                                            cs_real_t             bf[3][3],
                                            const cs_real_t       pimpv[3],
                                            cs_real_t             hint,
                                            const cs_real_t       hextv[3])
{
  for (int i = 0; i < 3; i++) {
    if (fabs(hextv[i]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[i] = pimpv[i];
      for (int jsou = 0; jsou < 3; jsou++)
        b[i][jsou] = 0.;

      /* Flux BCs */
      af[i] = -hint*pimpv[i];

      bf[0][0] = hint, bf[0][1] = 0.,   bf[0][2] = 0.;
      bf[1][0] = 0.,   bf[1][1] = hint, bf[1][2] = 0.;
      bf[2][0] = 0.,   bf[2][1] = 0.,   bf[2][2] = hint;

    }
    else {

      const cs_real_t val = hint/(hint + hextv[i]);
      const cs_real_t heq = hextv[i]*val;

      /* Gradient BCs */
      a[i] = hextv[i]*pimpv[i]/(hint + hextv[i]);

      b[0][0] = val, b[0][1] = 0.,  b[0][2] = 0.;
      b[1][0] = 0.,  b[1][1] = val, b[1][2] = 0.;
      b[2][0] = 0.,  b[2][1] = 0.,  b[2][2] = val;

      /* Flux BCs */
      af[i] = -heq*pimpv[i];

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

CS_F_HOST_DEVICE void
cs_boundary_conditions_set_convective_outlet_scalar
  (cs_real_t  &a,
   cs_real_t  &af,
   cs_real_t  &b,
   cs_real_t  &bf,
   cs_real_t  pimp,
   cs_real_t  cfl,
   cs_real_t  hint);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a vector for a given face with left anisotropic
 *        diffusion.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   hintt       internal exchange coefficient
 * \param[in]   hextv       external exchange coefficient
 *                          (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_vector_aniso
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        hintt[6],
   const cs_real_t        hextv[3])
{
  /* Gradient BCs */
  for (int i = 0; i < 3; i++) {
    if (fabs(hextv[i]) > 0.5*cs_math_infinite_r) {
      a[i] = pimpv[i];
      for (int jsou = 0; jsou < 3; jsou++)
        b[i][jsou] = 0.;
    }
    else {
      /* FIXME: at least log error message */
#if defined(CS_DEVICE_COMPILE)
      assert(0);
#else
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: hextv not set for component %d."),
                __func__, i);
#endif
    }
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, pimpv, af);
  for (int i = 0; i < 3; i++)
    af[i] = -af[i];

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
 * \param[out]  coefa       explicit BC coefficient for gradients
 * \param[out]  cofaf       explicit BC coefficient for diffusive flux
 * \param[out]  coefb       implicit BC coefficient for gradients
 * \param[out]  cofbf       implicit BC coefficient for diffusive flux
 * \param[in]   pimpts      Dirichlet value to impose
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   hextts      external exchange coefficient (10^30 by default)
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_tensor(cs_real_t        a[6],
                                            cs_real_t        af[6],
                                            cs_real_t        b[6][6],
                                            cs_real_t        bf[6][6],
                                            const cs_real_t  pimpts[6],
                                            cs_real_t        hint,
                                            const cs_real_t  hextts[6])
{
  for (int i = 0; i < 6; i++) {

    if (fabs(hextts[i]) > 0.5*cs_math_infinite_r) {
      /* Gradient BCs */
      a[i] = pimpts[i];
      for (int jsou = 0; jsou < 6; jsou++)
        b[i][jsou] = 0.;

      /* Flux BCs */
      af[i] = -hint * pimpts[i];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == i)
          bf[i][jsou] = hint;
        else
          bf[i][jsou] = 0.;
      }
    }

    else {

      const cs_real_t heq = hint * hextts[i] / (hint + hextts[i]);

      /* Gradient BCs */
      a[i] = hextts[i] * pimpts[i] / (hint + hextts[i]);
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == i)
          b[i][jsou] = hint / (hint + hextts[i]);
        else
          b[i][jsou] = 0.;
      }

      /* Flux BCs */
      af[i] = -heq * pimpts[i];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == i)
          bf[i][jsou] = heq;
        else
          bf[i][jsou] = 0.;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for a symmetric vector for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpv       Dirichlet value to impose on the normal component
 * \param[in]   qimpv       flux value to impose on the tangential components
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   normal      unit normal
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_generalized_sym_vector
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   cs_real_t              hint,
   const cs_nreal_t       normal[3])
{
  for (int i = 0; i < 3; i++) {

    /* Gradient BCs */
    a[i] = - qimpv[i]/cs::max(hint, 1.e-300);
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int j = 0; j < 3; j++) {

      a[i] = a[i] + normal[i]*normal[j]
        * (pimpv[j] + qimpv[j] / cs::max(hint, 1.e-300));

      if (j == i)
        b[i][j] = 1.0 - normal[i] * normal[j];
      else
        b[i][j] = - normal[i] * normal[j];
    }

    /* Flux BCs */
    af[i] = qimpv[i];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int j = 0; j < 3; j++){

      af[i] = af[i] - normal[i]*normal[j]
                  * (hint * pimpv[j] + qimpv[j]);

      bf[i][j] = hint * normal[i] * normal[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[out]  a            explicit BC coefficient for gradients
 * \param[out]  af           explicit BC coefficient for diffusive flux
 * \param[out]  b            implicit BC coefficient for gradients
 * \param[out]  bf           implicit BC coefficient for diffusive flux
 * \param[in]   pimpv        Dirichlet value to impose on the normal component
 * \param[in]   qimpv        flux value to impose on the tangential components
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   const cs_real_t        hint[6],
   const cs_nreal_t       normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for a vector for a given face.
 *
 * \param[out]  a            explicit BC coefficient for gradients
 * \param[out]  af           explicit BC coefficient for diffusive flux
 * \param[out]  b            implicit BC coefficient for gradients
 * \param[out]  bf           implicit BC coefficient for diffusive flux
 * \param[in]   pimpv        Dirichlet value to impose on the
 *                           tangential components
 * \param[in]   qimpv        flux value to impose on the normal component
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_generalized_dirichlet_vector
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   cs_real_t              hint,
   const cs_nreal_t       normal[3])
{
  for (int i = 0; i < 3; i++) {

    /* Gradient BC*/
    /* "[1 -n(x)n] Pimp" is divided into two */
    a[i] = pimpv[i];
    for (int j = 0; j < 3; j++) {

      a[i] = a[i] - normal[i]*normal[j]
        * (pimpv[j] + qimpv[j] / cs::max(hint, 1.e-300));

      b[i][j] = normal[i] * normal[j];
    }

    /* Flux BC */
    /* "[1 -n(x)n] Pimp" is divided into two */
    af[i] = -hint*pimpv[i];
    for (int j = 0; j < 3; j++) {

      af[i] = af[i] + normal[i]*normal[j]
        * (qimpv[j] + pimpv[j] * hint);

      if (j == i)
        bf[i][j] = hint * (1.0 - normal[i] * normal[j]);
      else
        bf[i][j] = - hint * normal[i] * normal[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for a symmetric tensor for a given face.
 *
 * That is Dirichlet on shear, Neumann on "normal components"
 *
 * \param[out]  a            explicit BC coefficient for gradients
 * \param[out]  af           explicit BC coefficient for diffusive flux
 * \param[out]  b            implicit BC coefficient for gradients
 * \param[out]  bf           implicit BC coefficient for diffusive flux
 * \param[in]   pimpv        Dirichlet value to impose on the
 *                           tangential components
 * \param[in]   qimpv        flux value to impose on the normal component
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_generalized_dirichlet_sym_tensor
  (cs_real_t              a[6],
   cs_real_t              af[6],
   cs_real_t              b[6][6],
   cs_real_t              bf[6][6],
   const cs_real_t        pimpv[6],
   const cs_real_t        qimpv[6],
   cs_real_t              hint,
   const cs_nreal_t       normal[3])
{
  const int iv2t[6] = {0, 1, 2, 0, 1, 0};
  const int jv2t[6] = {0, 1, 2, 1, 2, 2};
  const cs_real_t delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

  cs_real_t d[6][6];
  cs_real_t n[6][6];

  for (int ij = 0; ij < 6; ij++) {

    int i = iv2t[ij];
    int j = jv2t[ij];
    /* Gradient BC*/
    /* "[1 -n(x)n] Pimp" is divided into two */
    //a[ij] = pimpv[ij];
    a[ij] = 0;
    for (int kl = 0; kl < 6; kl++) {

      int k = iv2t[kl];
      int l = jv2t[kl];

      d[ij][kl] = (delta[i][k] - normal[i]*normal[k]) * normal[l]*normal[j]
                 + normal[i]*normal[k] * (delta[j][l] - normal[l]*normal[j]);

      n[ij][kl] = normal[i]*normal[k]*normal[l]*normal[j]
        +   (delta[i][k] - normal[i]*normal[k])
          * (delta[j][l] - normal[j]*normal[l]);

      a[ij] += d[ij][kl]*pimpv[kl]
              - n[ij][kl]*qimpv[kl]/(cs::max(hint, 1.e-300));

      b[ij][kl] = n[ij][kl];
    }

    /* Flux BC */
    /* "[1 -n(x)n] Pimp" is divided into two */
    //af[ij] = -hint*pimpv[ij];
    af[ij] = 0;
    for (int kl = 0; kl < 6; kl++) {

      int k = iv2t[kl];
      int l = jv2t[kl];

      d[ij][kl] = (delta[i][k] - normal[i]*normal[k]) * normal[l]*normal[j]
                 + normal[i]*normal[k] * (delta[j][l] - normal[l]*normal[j]);

      n[ij][kl] = normal[i]*normal[k]*normal[l]*normal[j]
        +   (delta[i][k] - normal[i]*normal[k])
          * (delta[j][l] - normal[l]*normal[j]);

      af[ij] += n[ij][kl]*qimpv[kl]
               - hint*d[ij][kl]*pimpv[kl];

      bf[ij][kl] = hint*d[ij][kl];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[out]  a            explicit BC coefficient for gradients
 * \param[out]  af           explicit BC coefficient for diffusive flux
 * \param[out]  b            implicit BC coefficient for gradients
 * \param[out]  bf           implicit BC coefficient for diffusive flux
 * \param[in]   pimpv        Dirichlet value to impose on the
 *                           tangential components
 * \param[in]   qimpv        flux value to impose on the normal component
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   const cs_real_t        hint[6],
   const cs_nreal_t       normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for a vector for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   cflv        local Courant number used to convect
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_convective_outlet_vector
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        cflv[3],
   cs_real_t              hint)
{
  for (int i = 0; i < 3; i++) {

    /* Gradient BCs */
    for (int j = 0; j < 3; j ++) {
      if (j == i)
        b[i][j] = cflv[i] / (1.0 + cflv[i]);
      else
        b[i][j] = 0.0;
    }
    a[i] = pimpv[i] * (1.0 - b[i][i]);

    /* Flux BCs */
    af[i] = -hint * a[i];
    for (int j = 0; j < 3; j++) {
      if (j == i)
        bf[i][j] = hint * (1.0 - b[i][j]);
    else
      bf[i][j] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective outlet BC for a tensor for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpts      Dirichlet value to impose
 * \param[in]   cflts       local Courant number used to convect
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_convective_outlet_tensor(cs_real_t        a[6],
                                                    cs_real_t        af[6],
                                                    cs_real_t        b[6][6],
                                                    cs_real_t        bf[6][6],
                                                    const cs_real_t  pimpts[6],
                                                    const cs_real_t  cflts[6],
                                                    cs_real_t        hint)
{
  for (int ij = 0; ij < 6; ij++) {

    /* Gradient BCs */
    for (int kl = 0; kl < 6; kl++) {
      if (kl == ij)
        b[ij][kl] = cflts[ij] / (1.0 + cflts[ij]);
      else
        b[ij][kl] = 0.0;
    }
    a[ij] = (1.0 - b[ij][ij]) * pimpts[ij];

    /* Flux BCs */
    af[ij] = -hint*a[ij];
    for (int kl = 0; kl < 6; kl++) {
      if (kl == ij)
        bf[ij][kl] = hint * (1.0 - b[ij][kl]);
      else
        bf[ij][kl] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for an anisotropic vector for a given face.
 *
 * \param[out]  a           explicit BC coefficient for gradients
 * \param[out]  af          explicit BC coefficient for diffusive flux
 * \param[out]  b           implicit BC coefficient for gradients
 * \param[out]  bf          implicit BC coefficient for diffusive flux
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   cflv        local Courant number used to convect
 * \param[in]   hintt       internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_convective_outlet_vector_aniso
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        cflv[3],
   const cs_real_t        hintt[6])
{
  for(int i = 0; i < 3; i++) {

    /* Gradient BCs */
    for (int j = 0; j < 3; j++) {
      if (j == i)
        b[i][j] = cflv[i]/(1.0+cflv[i]);
      else
        b[i][j] = 0.0;
    }
    a[i] = (1.0-b[i][i])*pimpv[i];

  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, a, af);
  for (int i = 0; i < 3; i++)
    af[i] = -af[i];

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
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a vector.
 *
 * \param[out]    a           explicit BC coefficient for gradients
 * \param[out]    af          explicit BC coefficient for diffusive flux
 * \param[out]    b           implicit BC coefficient for gradients
 * \param[out]    bf          implicit BC coefficient for diffusive flux
 * \param[in]     pimpv       Dirichlet value to impose
 * \param[in]     qimpv       flux value to impose
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
  (cs_real_t              a[3],
   cs_real_t              af[3],
   cs_real_t              b[3][3],
   cs_real_t              bf[3][3],
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3])
{
  for (int i = 0; i < 3; i++) {

    /* Gradient BCs */
    a[i] = pimpv[i];
    for (int j = 0; j < 3; j++)
      b[i][j] = 0.0;

    /* Flux BCs */
    af[i] = qimpv[i];
    for (int j = 0; j < 3; j++)
      bf[i][j] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a tensor
 *
 * \param[out]    coefa       explicit BC coefficient for gradients
 * \param[out]    cofaf       explicit BC coefficient for diffusive flux
 * \param[out]    coefb       implicit BC coefficient for gradients
 * \param[out]    cofbf       implicit BC coefficient for diffusive flux
 * \param[in]     pimpts      Dirichlet value to impose
 * \param[in]     qimpts      flux value to impose
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_tensor
  (cs_real_t        a[6],
   cs_real_t        af[6],
   cs_real_t        b[6][6],
   cs_real_t        bf[6][6],
   const cs_real_t  pimpts[6],
   const cs_real_t  qimpts[6])
{
  for (int ij = 0; ij < 6; ij++) {

    /* BS test on hextv ? if (abs(hextv[ij]) > cs_math_infinite_r * 0.5) */

    /* Gradient BCs */
    a[ij] = pimpts[ij];
    for (int kl = 0; kl < 6; kl++)
      b[ij][kl] = 0.0;

    /* Flux BCs */
    af[ij] = qimpts[ij];
    for (int kl = 0; kl < 6; kl++)
      bf[ij][kl] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#ifdef __cplusplus

/*============================================================================
 * Public C++ function definitions
 *============================================================================*/

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

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_scalar(cs_real_t  &a,
                                          cs_real_t  &af,
                                          cs_real_t  &b,
                                          cs_real_t  &bf,
                                          cs_real_t   qimp,
                                          cs_real_t   hint)
{
  /* Gradient BCs */
  a = -qimp/cs::max(hint, 1.e-300);
  b = 1.;

  /* Flux BCs */
  af = qimp;
  bf = 0.;
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

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_scalar_hmg(cs_real_t  &a,
                                              cs_real_t  &af,
                                              cs_real_t  &b,
                                              cs_real_t  &bf)
{
  /* Gradient BCs */
  a = 0.;
  b = 1.;

  /* Flux BCs */
  af = 0.;
  bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   Dirichlet value to impose
 * \param[in]   hint   internal exchange coefficient
 * \param[in]   hext   external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_scalar(cs_real_t  &a,
                                            cs_real_t  &af,
                                            cs_real_t  &b,
                                            cs_real_t  &bf,
                                            cs_real_t   pimp,
                                            cs_real_t   hint,
                                            cs_real_t   hext)
{
  if (hext < 0.) {

    /* Gradient BCs */
    a = pimp;
    b = 0.;

    /* Flux BCs */
    af = -hint*pimp;
    bf =  hint;

  }
  else {

    /* Gradient BCs */
    a = hext*pimp/(hint + hext);
    b = hint     /(hint + hext);

    /* Flux BCs */
    cs_real_t heq = hint*hext/(hint + hext);
    af = -heq*pimp;
    bf =  heq;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, zero flux for diffusion.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   dimp   flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_neumann_conv_h_neumann_diff_scalar
  (cs_real_t  &a,
   cs_real_t  &af,
   cs_real_t  &b,
   cs_real_t  &bf,
   cs_real_t   dimp,
   cs_real_t   hint)

{
  /* Gradient BCs */
  cs_boundary_conditions_set_neumann_scalar(a,
                                            af,
                                            b,
                                            bf,
                                            dimp,
                                            hint);

  /* Flux BCs */
  af = 0.;
  bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set BC for an affine scalar function for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pinf   affine part
 * \param[in]   ratio  linear part
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_affine_function_scalar
  (cs_real_t  &a,
   cs_real_t  &af,
   cs_real_t  &b,
   cs_real_t  &bf,
   cs_real_t   pinf,
   cs_real_t   ratio,
   cs_real_t   hint)
{
  /* Gradient BCs */
  b = ratio;
  a = pinf;

  /* Flux BCs */
  af = -hint * a;
  bf =  hint * (1. - b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, imposed flux for
 *         diffusion.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pinf   affine part
 * \param[in]   ratio  linear part
 * \param[in]   dimp   flux value to impose
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_affine_function_conv_neumann_diff_scalar
  (cs_real_t  &a,
   cs_real_t  &af,
   cs_real_t  &b,
   cs_real_t  &bf,
   cs_real_t   pinf,
   cs_real_t   ratio,
   cs_real_t   dimp)
{
  /* Gradient BCs */
  b = ratio;
  a = pinf;

  /* Flux BCs */
  af = dimp;
  bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set total flux as a Robin condition.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   hext   convective flux to be imposed
 * \param[in]   dimp   flux value to impose
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_total_flux(cs_real_t  &a,
                                      cs_real_t  &af,
                                      cs_real_t  &b,
                                      cs_real_t  &bf,
                                      cs_real_t   hext,
                                      cs_real_t   dimp)
{
  /* Gradients BCs */
  a = 0.;
  b = 1.;

  /* Flux BCs */
  af = dimp;
  bf = hext;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a scalar.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   Dirichlet value to impose
 * \param[in]   dimp   flux value to impose
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
  (cs_real_t  &a,
   cs_real_t  &af,
   cs_real_t  &b,
   cs_real_t  &bf,
   cs_real_t   pimp,
   cs_real_t   dimp)
{
  /* Gradients BC */
  a = pimp;
  b = 0.;

  /* Flux BCs */
  af = dimp;
  bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Update face value for gradient and diffusion when solving
 *         in increment.
 *
 * \param[in]      ctx                   reference to dispatch context
 * \param[in]      f                     pointer to field
 * \param[in]      bc_coeffs             boundary condition structure
 * \param[in]      inc                   0 if an increment, 1 otherwise
 * \param[in]      eqp                   equation parameters
 * \param[in]      need_compute_bc_grad  val_f must be computed
 * \param[in]      need_compute_bc_flux  flux must be computed
 * \param[in]      hyd_p_flag            hydrostatic pressure indicator
 * \param[in]      f_ext                 exterior force generating pressure
 * \param[in]      visel                 viscosity by cell, or nullptr
 * \param[in]      viscel                symmetric cell tensor
                                         \f$ \tens{\mu}_\celli \f$,
                                         or nullptr
 * \param[in]      weighb                boundary face weight for cells i in
 *                                       case of tensor diffusion, or nullptr
 * \param[in]      var                   variable values at cell centers
 * \param[in,out]  var_f                 face values for the gradient computation
 * \param[in,out]  flux                  face values for the diffusion computation
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_update_bc_coeff_face_values
  (cs_dispatch_context        &ctx,
   const cs_field_t           *f,
   const cs_field_bc_coeffs_t *bc_coeffs,
   const int                   inc,
   const cs_equation_param_t  *eqp,
   const bool                  need_compute_bc_grad,
   const bool                  need_compute_bc_flux,
   int                         hyd_p_flag,
   cs_real_t                   f_ext[][3],
   cs_real_t                   visel[],
   cs_real_t                   viscel[][6],
   const cs_real_t             weighb[],
   const cs_real_t             pvar[],
   cs_real_t                   val_f[],
   cs_real_t                   flux[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Update boundary coefficient face values for gradient and diffusion
 *         when solving for a given field.
 *
 * \param[in]       ctx                   reference to dispatch context
 * \param[in, out]  f                     pointer to field
 * \param[in]       eqp                   equation parameters
 * \param[in]       need_compute_bc_grad  val_f must be computed
 * \param[in]       need_compute_bc_flux  flux must be computed
 * \param[in]       hyd_p_flag            flag for hydrostatic pressure
 * \param[in]       f_ext                 exterior force generating pressure
 * \param[in]       viscel                symmetric cell tensor
                                          \f$ \tens{\mu}_\celli \f$,
                                          or nullptr
 * \param[in]       weighb                boundary face weight for cells i in
 *                                        case of tensor diffusion, or nullptr
 * \param[in]       pvar                  variable values at cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_update_bc_coeff_face_values
  (cs_dispatch_context        &ctx,
   cs_field_t                 *f,
   const cs_equation_param_t  *eqp,
   const bool                  need_compute_bc_grad,
   const bool                  need_compute_bc_flux,
   int                         hyd_p_flag,
   cs_real_t                   f_ext[][3],
   cs_real_t                   viscel[][6],
   const cs_real_t             weighb[],
   const cs_real_t             pvar[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Update face value for gradient and diffusion when solving
 *         in increment.
 *
 * \param[in]      ctx          reference to dispatch context
 * \param[in]      f            pointer to field
 * \param[in]      bc_coeffs    boundary condition structure for the variable
 * \param[in]      inc          0 if an increment, 1 otherwise
 * \param[in]      halo_type    halo type (extended or not)
 * \param[in]      var          variable values at cell centers
 * \param[in,out]  var_ip       boundary variable values at I' position
 * \param[in,out]  var_f        face values for the gradient computation
 * \param[in,out]  var_f_d      face values for the diffusion computation
 * \param[in,out]  var_f_d_lim  face values for the diffusion computation
 *                              (with limiter)
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_boundary_conditions_update_bc_coeff_face_values_strided
  (cs_dispatch_context        &ctx,
   cs_field_t                 *f,
   const cs_field_bc_coeffs_t *bc_coeffs,
   const int                   inc,
   const cs_equation_param_t  *eqp,
   const cs_real_t             pvar[][stride],
   cs_real_t                   val_ip[][stride],
   cs_real_t                   val_f[][stride],
   cs_real_t                   flux[][stride],
   cs_real_t                   flux_lim[][stride]);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Update boundary coefficient face values for gradient and diffusion
 *         when solving for a given field.
 *
 * \param[in]       ctx    reference to dispatch context
 * \param[in, out]  f      pointer to field
 * \param[in]       pvar   variable values at cell centers
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_boundary_conditions_update_bc_coeff_face_values_strided
  (cs_dispatch_context  &ctx,
   cs_field_t           *f,
   const cs_real_t       pvar[][stride]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate BC coefficients face values if needed.
 *
 * \param[in, out]  bc_coeffs  pointer to boundary conditions coefficients.
 * \param[in]       n_b_faces  number of boundary faces
 * \param[in]       dim        associated field dimension
 * \param[in]       amode      allocation mode
 * \param[in]       limiter    is a limiter active ?
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_ensure_bc_coeff_face_values_allocated
  (cs_field_bc_coeffs_t  *bc_coeffs,
   cs_lnum_t              n_b_faces,
   cs_lnum_t              dim,
   cs_alloc_mode_t        amode,
   bool                   limiter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize boundary condition coefficients solve arrays.
 *
 * \param[in, out]  c          reference to structure to initialize.
 * \param[in]       n_b_faces  number of boundary faces
 * \param[in]       stride     variable dimension
 * \param[in]       amode      allocation mode
 * \param[in]       limiter    is a limiter active ?
 */
/*----------------------------------------------------------------------------*/

void
cs_init_bc_coeffs_solve(cs_bc_coeffs_solve_t  &c,
                        cs_lnum_t              n_b_faces,
                        cs_lnum_t              stride,
                        cs_alloc_mode_t        amode,
                        bool                   limiter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free boundary condition coefficients solve arrays.
 *
 * \param[in, out]  c          reference to structure to initialize.
 */
/*----------------------------------------------------------------------------*/

void
cs_clear_bc_coeffs_solve(cs_bc_coeffs_solve_t  &c);

#endif /* cplusplus */

/*----------------------------------------------------------------------------*/

#endif /* __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__ */
