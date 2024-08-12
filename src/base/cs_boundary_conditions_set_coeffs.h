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
#include "cs_field.h"
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
 * \param[in]     isvhb         id of field whose exchange coeffient should be
 *                               saved at the walls, or -1.
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
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary conditions structure
 * \param[in]   qimp          flux value to impose
 * \param[in]   hint          internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_scalar(cs_lnum_t             f_id,
                                          cs_field_bc_coeffs_t *bc_coeffs,
                                          cs_real_t             qimp,
                                          cs_real_t             hint)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradient BCs */
  a[f_id] = -qimp/cs_math_fmax(hint, 1.e-300);
  b[f_id] = 1.;

  /* Flux BCs */
  af[f_id] = qimp;
  bf[f_id] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set homogeneous Neumann BC for a scalar for a given face.
 *
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary conditions structure
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_scalar_hmg(cs_lnum_t             f_id,
                                              cs_field_bc_coeffs_t *bc_coeffs)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradient BCs */
  a[f_id] = 0.;
  b[f_id] = 1.;

  /* Flux BCs */
  af[f_id] = 0.;
  bf[f_id] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Neumann BC for a scalar for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   qimpv       flux value to impose
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector(cs_lnum_t             f_id,
                                          cs_field_bc_coeffs_t *bc_coeffs,
                                          const cs_real_t       qimpv[3],
                                          cs_real_t             hint)
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  /* Gradient BCs */

  for (size_t i = 0; i < 3; i++) {
    a[f_id][i] = -qimpv[i] / fmax(hint, 1.e-300);
  }

  b[f_id][0][0] = 1., b[f_id][0][1] = 0., b[f_id][0][2] = 0.;
  b[f_id][1][0] = 0., b[f_id][1][1] = 1., b[f_id][1][2] = 0.;
  b[f_id][2][0] = 0., b[f_id][2][1] = 0., b[f_id][2][2] = 1.;

  /* Flux BCs */

  for (size_t i = 0; i < 3; i++) {
    af[f_id][i] = qimpv[i];

    for (size_t j = 0; j < 3; j++)
      bf[f_id][i][j] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set neumann BC for an anisotropic vector for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   qimpv       flux value to impose
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector_aniso(cs_lnum_t             f_id,
                                                cs_field_bc_coeffs_t *bc_coeffs,
                                                const cs_real_t       qimpv[3],
                                                const cs_real_t       hint[6])
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

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
  cs_math_sym_33_3_product(invh, qimpv, a[f_id]);
  for (int isou = 0; isou < 3; isou++)
    a[f_id][isou] = -a[f_id][isou];

  b[f_id][0][0] = 1.0, b[f_id][0][1] = 0.0, b[f_id][0][2] = 0.0;
  b[f_id][1][0] = 0.0, b[f_id][1][1] = 1.0, b[f_id][1][2] = 0.0;
  b[f_id][2][0] = 0.0, b[f_id][2][1] = 0.0, b[f_id][2][2] = 1.0;

  for (int isou = 0; isou < 3; isou++) {

    /* Flux BCs */
    af[f_id][isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      bf[f_id][isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief  Set Neumann boundary conditions for a tensor for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   qimpts      flux value to impose
 * \param[in]   hint        internal exchange coefficient
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
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary conditions structure
 * \param[in]   pimp          Dirichlet value to impose
 * \param[in]   hint          internal exchange coefficient
 * \param[in]   hext          external exchange coefficient
 *                            (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_scalar(cs_lnum_t             f_id,
                                            cs_field_bc_coeffs_t *bc_coeffs,
                                            cs_real_t             pimp,
                                            cs_real_t             hint,
                                            cs_real_t             hext)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  //if (fabs(hext) > cs_math_infinite_r*0.5) {
  if (hext < 0.) {

    /* Gradient BCs */
    a[f_id] = pimp;
    b[f_id] = 0.;

    /* Flux BCs */
    af[f_id] = -hint*pimp;
    bf[f_id] =  hint;

  }
  else {

    /* Gradient BCs */
    a[f_id] = hext*pimp/(hint + hext);
    b[f_id] = hint     /(hint + hext);

    /* Flux BCs */
    cs_real_t heq = hint*hext/(hint + hext);
    af[f_id] = -heq*pimp;
    bf[f_id] =  heq;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a vector for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   hextv       external exchange coefficient
 *                          (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector(cs_lnum_t             f_id,
                                            cs_field_bc_coeffs_t *bc_coeffs,
                                            const cs_real_t       pimpv[3],
                                            cs_real_t             hint,
                                            const cs_real_t       hextv[3])
{

  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[f_id][isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[f_id][isou][jsou] = 0.;

      /* Flux BCs */
      af[f_id][isou] = -hint*pimpv[isou];

      bf[f_id][0][0] = hint, bf[f_id][0][1] = 0.,   bf[f_id][0][2] = 0.;
      bf[f_id][1][0] = 0.,   bf[f_id][1][1] = hint, bf[f_id][1][2] = 0.;
      bf[f_id][2][0] = 0.,   bf[f_id][2][1] = 0.,   bf[f_id][2][2] = hint;

    }
    else {

      const cs_real_t val = hint/(hint + hextv[isou]);
      const cs_real_t heq = hextv[isou]*val;

      /* Gradient BCs */
      a[f_id][isou] = hextv[isou]*pimpv[isou]/(hint + hextv[isou]);

      b[f_id][0][0] = val, b[f_id][0][1] = 0.,  b[f_id][0][2] = 0.;
      b[f_id][1][0] = 0.,  b[f_id][1][1] = val, b[f_id][1][2] = 0.;
      b[f_id][2][0] = 0.,  b[f_id][2][1] = 0.,  b[f_id][2][2] = val;

      /* Flux BCs */
      af[f_id][isou] = -heq*pimpv[isou];

      bf[f_id][0][0] = heq, bf[f_id][0][1] = 0.,  bf[f_id][0][2] = 0.;
      bf[f_id][1][0] = 0.,  bf[f_id][1][1] = heq, bf[f_id][1][2] = 0.;
      bf[f_id][2][0] = 0.,  bf[f_id][2][1] = 0.,  bf[f_id][2][2] = heq;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Set convective oulet BC for a scalar for a given face.
 *
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary conditions structure
 * \param[in]   pimp          flux value to impose
 * \param[in]   cfl           local Courant number used to convect
 * \param[in]   hint          internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_convective_outlet_scalar
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   cs_real_t              pimp,
   cs_real_t              cfl,
   cs_real_t              hint);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a vector for a given face with left anisotropic
 *        diffusion.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   hintt       internal exchange coefficient
 * \param[in]   hextv       external exchange coefficient
 *                          (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector_aniso
(cs_lnum_t              f_id,
 cs_field_bc_coeffs_t  *bc_coeffs,
 const cs_real_t        pimpv[3],
 const cs_real_t        hintt[6],
 const cs_real_t        hextv[3])
{

  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  /* Gradient BCs */
  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {
      a[f_id][isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[f_id][isou][jsou] = 0.;
    }
    else {
      /* FIXME: at least log error message */
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: hextv not set for component %d."),
                __func__, isou);
    }
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, pimpv, af[f_id]);
  for (int isou = 0; isou < 3; isou++)
    af[f_id][isou] = -af[f_id][isou];

  bf[f_id][0][0] = hintt[0];
  bf[f_id][1][1] = hintt[1];
  bf[f_id][2][2] = hintt[2];
  bf[f_id][0][1] = hintt[3];
  bf[f_id][1][0] = hintt[3];
  bf[f_id][1][2] = hintt[4];
  bf[f_id][2][1] = hintt[4];
  bf[f_id][0][2] = hintt[5];
  bf[f_id][2][0] = hintt[5];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a tensor for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpts      Dirichlet value to impose
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   hextts      external exchange coefficient (10^30 by default)
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
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpv       Dirichlet value to impose on the normal component
 * \param[in]   qimpv       flux value to impose on the tangential components
 * \param[in]   hint        internal exchange coefficient
 * \param[in]   normal      unit normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_sym_vector
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   cs_real_t              hint,
   const cs_real_t        normal[3])
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    a[f_id][isou] = - qimpv[isou]/cs_math_fmax(hint, 1.e-300);
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++) {

      a[f_id][isou] = a[f_id][isou] + normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      if (jsou == isou)
        b[f_id][isou][jsou] = 1.0 - normal[isou] * normal[jsou];
      else
        b[f_id][isou][jsou] = - normal[isou] * normal[jsou];
    }

    /* Flux BCs */
    af[f_id][isou] = qimpv[isou];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++){

      af[f_id][isou] = af[f_id][isou] - normal[isou]*normal[jsou]
                  * (hint * pimpv[jsou] + qimpv[jsou]);

      bf[f_id][isou][jsou] = hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[in]   f_id         face id
 * \param[out]  bc_coeffs    boundary conditions structure
 * \param[in]   pimpv        Dirichlet value to impose on the normal component
 * \param[in]   qimpv        flux value to impose on the tangential components
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   const cs_real_t        hint[6],
   const cs_real_t        normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for a vector for a given face.
 *
 * \param[in]   f_id         face id
 * \param[out]  bc_coeffs    boundary conditions structure
 * \param[in]   pimpv        Dirichlet value to impose on the
 *                           tangential components
 * \param[in]   qimpv        flux value to impose on the normal component
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_dirichlet_vector
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   cs_real_t              hint,
   const cs_real_t        normal[3])
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BC*/
    /* "[1 -n(x)n] Pimp" is divided into two */
    a[f_id][isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      a[f_id][isou] = a[f_id][isou] - normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      b[f_id][isou][jsou] = normal[isou] * normal[jsou];
    }

    /* Flux BC */
    /* "[1 -n(x)n] Pimp" is divided into two */
    af[f_id][isou] = -hint*pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      af[f_id][isou] = af[f_id][isou] + normal[isou]*normal[jsou]
        * (qimpv[jsou] + pimpv[jsou] * hint);

      if (jsou == isou)
        bf[f_id][isou][jsou] = hint * (1.0 - normal[isou] * normal[jsou]);
      else
        bf[f_id][isou][jsou] = - hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[in]   f_id         face id
 * \param[out]  bc_coeffs    boundary conditions structure
 * \param[in]   pimpv        Dirichlet value to impose on the
 *                           tangential components
 * \param[in]   qimpv        flux value to impose on the normal component
 * \param[in]   hint         internal exchange coefficient
 * \param[in]   normal       unit normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3],
   const cs_real_t        hint[6],
   const cs_real_t        normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for a vector for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   cflv        local Courant number used to convect
 * \param[in]   hint        internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        cflv[3],
   cs_real_t              hint)
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou ++) {
      if (jsou == isou)
        b[f_id][isou][jsou] = cflv[isou] / (1.0 + cflv[isou]);
      else
        b[f_id][isou][jsou] = 0.0;
    }
    a[f_id][isou] = pimpv[isou] * (1.0 - b[f_id][isou][isou]);

    /* Flux BCs */
    af[f_id][isou] = -hint * a[f_id][isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        bf[f_id][isou][jsou] = hint * (1.0 - b[f_id][isou][jsou]);
    else
      bf[f_id][isou][jsou] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective outlet BC for a tensor for a given face.
 *
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpts      Dirichlet value to impose
 * \param[in]   cflts       local Courant number used to convect
 * \param[in]   hint        internal exchange coefficient
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
 * \param[in]   f_id        face id
 * \param[out]  bc_coeffs   BC structure
 * \param[in]   pimpv       Dirichlet value to impose
 * \param[in]   cflv        local Courant number used to convect
 * \param[in]   hintt       internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector_aniso
(cs_lnum_t              f_id,
 cs_field_bc_coeffs_t  *bc_coeffs,
 const cs_real_t        pimpv[3],
 const cs_real_t        cflv[3],
 const cs_real_t        hintt[6])
{

  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for(int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        b[f_id][isou][jsou] = cflv[isou]/(1.0+cflv[isou]);
      else
        b[f_id][isou][jsou] = 0.0;
    }
    a[f_id][isou] = (1.0-b[f_id][isou][isou])*pimpv[isou];

  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, a[f_id], af[f_id]);
  for (int isou = 0; isou < 3; isou++)
    af[f_id][isou] = -af[f_id][isou];

  bf[f_id][0][0] = hintt[0]*(1.0 - b[f_id][0][0]);
  bf[f_id][1][1] = hintt[1]*(1.0 - b[f_id][1][1]);
  bf[f_id][2][2] = hintt[2]*(1.0 - b[f_id][2][2]);
  bf[f_id][0][1] = hintt[3]*(1.0 - b[f_id][0][0]);
  bf[f_id][1][0] = hintt[3]*(1.0 - b[f_id][0][0]);
  bf[f_id][1][2] = hintt[4]*(1.0 - b[f_id][1][1]);
  bf[f_id][2][1] = hintt[4]*(1.0 - b[f_id][1][1]);
  bf[f_id][0][2] = hintt[5]*(1.0 - b[f_id][2][2]);
  bf[f_id][2][0] = hintt[5]*(1.0 - b[f_id][2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set BC for an affine scalar function for a given face.
 *
 * \param[in]     f_id          face id
 * \param[out]    bc_coeffs     boundary condition structure
 * \param[in]     pinf          affine part
 * \param[in]     ratio         linear part
 * \param[in]     hint          internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_scalar
  (cs_lnum_t             f_id,
   cs_field_bc_coeffs_t *bc_coeffs,
   cs_real_t             pinf,
   cs_real_t             ratio,
   cs_real_t             hint)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradient BCs */
  b[f_id] = ratio;
  a[f_id] = pinf;

  /* Flux BCs */
  af[f_id] = -hint * a[f_id];
  bf[f_id] =  hint * (1. - b[f_id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, zero flux for diffusion.
 *
 * \param[in]     f_id          face id
 * \param[out]    bc_coeffs     boundary condition structure
 * \param[in]     dimp          flux value to impose
 * \param[in]     hint          internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_conv_h_neumann_diff_scalar
  (cs_lnum_t             f_id,
   cs_field_bc_coeffs_t *bc_coeffs,
   cs_real_t             dimp,
   cs_real_t             hint)

{
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradient BCs */
  cs_boundary_conditions_set_neumann_scalar(f_id,
                                            bc_coeffs,
                                            dimp,
                                            hint);

  /* Flux BCs */
  af[f_id] = 0.;
  bf[f_id] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, imposed flux for
 *         diffusion.
 *
 * \param[in]     f_id          face id
 * \param[out]    bc_coeffs     boundary condition structure
 * \param[in]     pinf          affine part
 * \param[in]     ratio         linear part
 * \param[in]     dimp          flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_conv_neumann_diff_scalar
  (cs_lnum_t             f_id,
   cs_field_bc_coeffs_t *bc_coeffs,
   cs_real_t             pinf,
   cs_real_t             ratio,
   cs_real_t             dimp)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradient BCs */
  b[f_id] = ratio;
  a[f_id] = pinf;

  /* Flux BCs */
  af[f_id] = dimp;
  bf[f_id] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set total flux as a Robin condition.
 *
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary condition structure
 * \param[in]   hext          convective flux to be imposed
 * \param[in]   dimp          flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_total_flux(cs_lnum_t             f_id,
                                      cs_field_bc_coeffs_t *bc_coeffs,
                                      cs_real_t             hext,
                                      cs_real_t             dimp)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradients BCs */
  a[f_id] = 0.;
  b[f_id] = 1.;

  /* Flux BCs */
  af[f_id] = dimp;
  bf[f_id] = hext;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a scalar.
 *
 * \param[in]   f_id          face id
 * \param[out]  bc_coeffs     boundary condition structure
 * \param[in]   pimp          Dirichlet value to impose
 * \param[in]   dimp          flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   cs_real_t              pimp,
   cs_real_t              dimp)
{
  cs_real_t *a = bc_coeffs->a;
  cs_real_t *b = bc_coeffs->b;
  cs_real_t *af = bc_coeffs->af;
  cs_real_t *bf = bc_coeffs->bf;

  /* Gradients BC */
  a[f_id] = pimp;
  b[f_id] = 0.;

  /* Flux BCs */
  af[f_id] = dimp;
  bf[f_id] = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a vector.
 *
 * \param[in]     f_id        face id
 * \param[out]    bc_coeffs   BC structure
 * \param[in]     pimpv       Dirichlet value to impose
 * \param[in]     qimpv       flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
  (cs_lnum_t              f_id,
   cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_real_t        pimpv[3],
   const cs_real_t        qimpv[3])
{
  cs_real_3_t  *a = (cs_real_3_t *)bc_coeffs->a;
  cs_real_33_t *b = (cs_real_33_t *)bc_coeffs->b;
  cs_real_3_t  *af = (cs_real_3_t *)bc_coeffs->af;
  cs_real_33_t *bf = (cs_real_33_t *)bc_coeffs->bf;

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    a[f_id][isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      b[f_id][isou][jsou] = 0.0;

    /* Flux BCs */
    af[f_id][isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      bf[f_id][isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a tensor
 *
 * \param[in]     f_id        face id
 * \param[out]    bc_coeffs   BC structure
 * \param[in]     pimpts      Dirichlet value to impose
 * \param[in]     qimpts      flux value to impose
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

    /* BS test on hextv ? if (abs(hextv[isou]) > cs_math_infinite_r * 0.5) */

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
