/*============================================================================
 * Godunov-type projection of the convective subsystem of second-moment
 * turbulence closures.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string>
#include <vector>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_halo.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_reducers.h"
#include "base/cs_time_step.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_turbulence_rij_projection.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Tensor to vector (t2v) and vector to tensor (v2t) mask arrays */

#undef _T2V
#undef _IV2T
#undef _JV2T
#define _T2V {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};
#define _IV2T {0, 1, 2, 0, 1, 0};
#define _JV2T {0, 1, 2, 1, 2, 2};

/* Largest number of sub-steps of the projection in one time step */

#define _N_SUB_STEPS_MAX 1000

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* One state of a face split into its normal and orthogonal parts, in the
   global components (no tangential basis) */

template <int D>
struct _proj_side_t {

  cs_real_t  un;         /* normal velocity U.n */
  cs_real_t  a;          /* sqrt(2 C_nn) */
  cs_real_t  up[D];      /* U_perp = U - u_n n */
  cs_real_t  c[D];       /* C_perp,n / sqrt(C_nn), zero if C_nn = 0 */
  cs_real_t  pi[D][D];   /* Schur complement C_perp,perp - c (x) c */

};

/* Solution of the Riemann problem of a face, frozen advection speed and
   jump relation J_T */

template <int D>
struct _proj_rp_t {

  _proj_side_t<D>  l;    /* left state */
  _proj_side_t<D>  r;    /* right state */

  cs_real_t  v;          /* frozen advection speed */
  cs_real_t  a_s;        /* intermediate sqrt(2 C_nn) */
  cs_real_t  u_s;        /* intermediate normal velocity */
  cs_real_t  up2[D];     /* U_perp in regions II and III */
  cs_real_t  c2[D];      /* C_perp,n / sqrt(C_nn) in regions II and III */
  cs_real_t  xi[7];      /* wave positions, in increasing order */

};

/* Integrals over a strip of the averaged quantities */

template <int D>
struct _proj_mean_t {

  cs_real_t  u[D];       /* U */
  cs_real_t  s;          /* sqrt(C_nn) */
  cs_real_t  c[D];       /* C_perp,n / sqrt(C_nn) */
  cs_real_t  cpp[D][D];  /* C_perp,perp */

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Time step at which the projection was last computed */

static int _nt_last = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Split a state into its normal and orthogonal parts.
 *
 * \param[in]   n   unit normal
 * \param[in]   u   mean state
 * \param[in]   c   generalised covariance
 * \param[out]  k   split state
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE void
_proj_split(const cs_real_t   n[D],
            const cs_real_t   u[D],
            const cs_real_t   c[D][D],
            _proj_side_t<D>  &k)
{
  cs_real_t cn[D];
  cs_real_t cnn = 0.;

  k.un = 0.;
  for (int i = 0; i < D; i++) {
    k.un += u[i]*n[i];
    cn[i] = 0.;
    for (int j = 0; j < D; j++)
      cn[i] += c[i][j]*n[j];
  }
  for (int i = 0; i < D; i++)
    cnn += n[i]*cn[i];

  const cs_real_t sq = sqrt(cs::max(cnn, 0.));
  const cs_real_t inv_sq = (sq > 0.) ? 1./sq : 0.;

  k.a = sqrt(2.)*sq;

  for (int i = 0; i < D; i++) {
    k.up[i] = u[i] - k.un*n[i];
    k.c[i] = (cn[i] - cnn*n[i])*inv_sq;
  }

  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++)
      k.pi[i][j] =   c[i][j] - n[i]*cn[j] - cn[i]*n[j] + cnn*n[i]*n[j]
                   - k.c[i]*k.c[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integral, over a piece of a strip, of the states connected
 *        to one side of the face by its genuinely nonlinear wave alone.
 *
 * These states depend on a = sqrt(2 C_nn) only. Over the piece, U, sqrt(C_nn)
 * and C_perp,n / sqrt(C_nn) are affine in a, and C_perp,perp in a^2.
 *
 * \param[in]       n    unit normal
 * \param[in]       k    split state of the side
 * \param[in]       sgn  -1 for the left side, +1 for the right side
 * \param[in]       am   mean of a over the piece
 * \param[in]       a2m  mean of a^2 over the piece
 * \param[in]       l    length of the piece
 * \param[in, out]  m    integrals
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE void
_proj_add_wave_state(const cs_real_t         n[D],
                     const _proj_side_t<D>  &k,
                     cs_real_t               sgn,
                     cs_real_t               am,
                     cs_real_t               a2m,
                     cs_real_t               l,
                     _proj_mean_t<D>        &m)
{
  const cs_real_t r = (k.a > 0.) ? am/k.a : 0.;
  const cs_real_t r2 = (k.a > 0.) ? a2m/(k.a*k.a) : 0.;
  const cs_real_t un = k.un + sgn*(am - k.a);
  const cs_real_t cu = sgn*sqrt(2.)*(r - 1.);

  for (int i = 0; i < D; i++) {
    m.u[i] += l*(un*n[i] + k.up[i] + cu*k.c[i]);
    m.c[i] += l*r*k.c[i];
  }
  m.s += l*am/sqrt(2.);

  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++)
      m.cpp[i][j] += l*(k.pi[i][j] + r2*k.c[i]*k.c[j]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the Riemann problem of a face.
 *
 * \param[in]   n    unit normal, from the left to the right state
 * \param[in]   v    frozen advection speed
 * \param[in]   ul   left mean state
 * \param[in]   cl   left generalised covariance
 * \param[in]   ur   right mean state
 * \param[in]   cr   right generalised covariance
 * \param[out]  rp   solution
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE void
_proj_solve(const cs_real_t  n[D],
            cs_real_t        v,
            const cs_real_t  ul[D],
            const cs_real_t  cl[D][D],
            const cs_real_t  ur[D],
            const cs_real_t  cr[D][D],
            _proj_rp_t<D>   &rp)
{
  _proj_split<D>(n, ul, cl, rp.l);
  _proj_split<D>(n, ur, cr, rp.r);

  const cs_real_t al = rp.l.a;
  const cs_real_t ar = rp.r.a;

  /* Intermediate normal state; a_s = 0 when the two sides separate */

  const cs_real_t as = cs::max(0.5*(al + ar + rp.l.un - rp.r.un), 0.);
  rp.v = v;
  rp.a_s = as;
  rp.u_s = 0.5*(rp.l.un + rp.r.un) - 0.5*(ar - al);

  /* Regions I and IV, then II and III between the shear waves */

  const cs_real_t r_1 = (al > 0.) ? as/al : 0.;
  const cs_real_t r_4 = (ar > 0.) ? as/ar : 0.;

  for (int i = 0; i < D; i++) {
    const cs_real_t up_1 = rp.l.up[i] - sqrt(2.)*(r_1 - 1.)*rp.l.c[i];
    const cs_real_t up_4 = rp.r.up[i] + sqrt(2.)*(r_4 - 1.)*rp.r.c[i];
    const cs_real_t c_1 = r_1*rp.l.c[i];
    const cs_real_t c_4 = r_4*rp.r.c[i];
    rp.up2[i] = 0.5*(up_1 + up_4) + 0.5*(c_1 - c_4);
    rp.c2[i] = 0.5*(c_1 + c_4) + 0.5*(up_1 - up_4);
  }

  /* Wave positions: genuinely nonlinear waves (head and tail of a
     rarefaction, or shock), shear waves (put on the shock when slower)
     and contact */

  cs_real_t h_1, t_1, t_4, h_4;
  if (as > al) {
    h_1 = v - 0.5*(as + al);
    t_1 = h_1;
  }
  else {
    h_1 = v - al;
    t_1 = v - as;
  }
  if (as > ar) {
    h_4 = v + 0.5*(as + ar);
    t_4 = h_4;
  }
  else {
    h_4 = v + ar;
    t_4 = v + as;
  }

  const cs_real_t s = as/sqrt(2.);

  rp.xi[0] = h_1;
  rp.xi[1] = t_1;
  rp.xi[2] = cs::max(v - s, t_1);
  rp.xi[3] = v;
  rp.xi[4] = cs::min(v + s, t_4);
  rp.xi[5] = t_4;
  rp.xi[6] = h_4;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Largest wave speed of a face.
 *
 * \param[in]  n    unit normal, from the left to the right state
 * \param[in]  v    frozen advection speed
 * \param[in]  ul   left mean state
 * \param[in]  cl   left generalised covariance
 * \param[in]  ur   right mean state
 * \param[in]  cr   right generalised covariance
 *
 * \return  largest absolute position of the heads or shocks of the
 *          genuinely nonlinear waves
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE cs_real_t
_proj_wave_speed(const cs_real_t  n[D],
                 cs_real_t        v,
                 const cs_real_t  ul[D],
                 const cs_real_t  cl[D][D],
                 const cs_real_t  ur[D],
                 const cs_real_t  cr[D][D])
{
  cs_real_t unl = 0., unr = 0., cnnl = 0., cnnr = 0.;

  for (int i = 0; i < D; i++) {
    unl += ul[i]*n[i];
    unr += ur[i]*n[i];
    for (int j = 0; j < D; j++) {
      cnnl += n[i]*cl[i][j]*n[j];
      cnnr += n[i]*cr[i][j]*n[j];
    }
  }

  const cs_real_t al = sqrt(2.*cs::max(cnnl, 0.));
  const cs_real_t ar = sqrt(2.*cs::max(cnnr, 0.));
  const cs_real_t as = cs::max(0.5*(al + ar + unl - unr), 0.);

  const cs_real_t h_1 = (as > al) ? v - 0.5*(as + al) : v - al;
  const cs_real_t h_4 = (as > ar) ? v + 0.5*(as + ar) : v + ar;

  return cs::max(cs::abs(h_1), cs::abs(h_4));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mean of the solution of a Riemann problem over a strip, and
 *        rebuilt state.
 *
 * The quantities averaged are U, sqrt(C_nn), C_perp,n / sqrt(C_nn) and
 * C_perp,perp; the means are exact (piecewise affine or quadratic in xi).
 *
 * \param[in]   n    unit normal, from the left to the right state
 * \param[in]   rp   solution of the Riemann problem
 * \param[in]   xa   lower bound of the strip in xi
 * \param[in]   xb   upper bound of the strip in xi (xb > xa)
 * \param[out]  u    rebuilt mean state
 * \param[out]  c    rebuilt generalised covariance
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE void
_proj_strip_mean(const cs_real_t       n[D],
                 const _proj_rp_t<D>  &rp,
                 cs_real_t             xa,
                 cs_real_t             xb,
                 cs_real_t             u[D],
                 cs_real_t             c[D][D])
{
  _proj_mean_t<D> m;

  for (int i = 0; i < D; i++) {
    m.u[i] = 0.;
    m.c[i] = 0.;
    for (int j = 0; j < D; j++)
      m.cpp[i][j] = 0.;
  }
  m.s = 0.;

  /* The wave positions, clipped to the strip, cut it into pieces:
     L, left fan, I, II, III, IV, right fan, R */

  cs_real_t p[9];
  p[0] = xa;
  for (int k = 0; k < 7; k++)
    p[k+1] = cs::min(cs::max(rp.xi[k], xa), xb);
  p[8] = xb;

  const cs_real_t v = rp.v;
  const cs_real_t as = rp.a_s;

  for (int k = 0; k < 8; k++) {

    const cs_real_t l = p[k+1] - p[k];
    if (!(l > 0.))
      continue;

    switch (k) {

    case 0:
      _proj_add_wave_state<D>(n, rp.l, -1., rp.l.a, rp.l.a*rp.l.a, l, m);
      break;

    case 1:
      {
        const cs_real_t a0 = v - p[k], a1 = v - p[k+1];
        _proj_add_wave_state<D>(n, rp.l, -1.,
                                0.5*(a0 + a1),
                                (a0*a0 + a0*a1 + a1*a1)/3.,
                                l, m);
      }
      break;

    case 2:
      _proj_add_wave_state<D>(n, rp.l, -1., as, as*as, l, m);
      break;

    case 3:
      [[fallthrough]];
    case 4:
      {
        const cs_real_t (*pi)[D] = (k == 3) ? rp.l.pi : rp.r.pi;
        for (int i = 0; i < D; i++) {
          m.u[i] += l*(rp.u_s*n[i] + rp.up2[i]);
          m.c[i] += l*rp.c2[i];
          for (int j = 0; j < D; j++)
            m.cpp[i][j] += l*(pi[i][j] + rp.c2[i]*rp.c2[j]);
        }
        m.s += l*as/sqrt(2.);
      }
      break;

    case 5:
      _proj_add_wave_state<D>(n, rp.r, 1., as, as*as, l, m);
      break;

    case 6:
      {
        const cs_real_t a0 = p[k] - v, a1 = p[k+1] - v;
        _proj_add_wave_state<D>(n, rp.r, 1.,
                                0.5*(a0 + a1),
                                (a0*a0 + a0*a1 + a1*a1)/3.,
                                l, m);
      }
      break;

    case 7:
      _proj_add_wave_state<D>(n, rp.r, 1., rp.r.a, rp.r.a*rp.r.a, l, m);
      break;

    default:
      assert(0);
      break;
    }

  }

  /* Rebuilt state: C_nn = (mean of sqrt(C_nn))^2,
     C_perp,n = (mean of sqrt(C_nn)) (mean of C_perp,n / sqrt(C_nn)) */

  const cs_real_t inv_w = 1./(xb - xa);
  const cs_real_t s = m.s*inv_w;

  cs_real_t cp[D];
  for (int i = 0; i < D; i++) {
    u[i] = m.u[i]*inv_w;
    cp[i] = s*m.c[i]*inv_w;
  }

  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++)
      c[i][j] =   s*s*n[i]*n[j] + n[i]*cp[j] + cp[i]*n[j]
                + m.cpp[i][j]*inv_w;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mean over a strip of the passive scalar eps k^(-C_eps1), which
 *        only jumps at the contact.
 *
 * \param[in]  v    frozen advection speed (position of the contact)
 * \param[in]  xa   lower bound of the strip in xi
 * \param[in]  xb   upper bound of the strip in xi (xb > xa)
 * \param[in]  wl   left value
 * \param[in]  wr   right value
 *
 * \return  mean value over the strip
 */
/*----------------------------------------------------------------------------*/

static CS_F_HOST_DEVICE cs_real_t
_proj_strip_mean_contact(cs_real_t  v,
                         cs_real_t  xa,
                         cs_real_t  xb,
                         cs_real_t  wl,
                         cs_real_t  wr)
{
  const cs_real_t f = (cs::min(cs::max(v, xa), xb) - xa)/(xb - xa);

  return f*wl + (1. - f)*wr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Full symmetric tensor from its 6 components (11, 22, 33, 12, 23, 13).
 *
 * \param[in]   r6   symmetric tensor, compact storage
 * \param[out]  r    full tensor
 */
/*----------------------------------------------------------------------------*/

static CS_F_HOST_DEVICE void
_proj_sym_to_full(const cs_real_t  r6[6],
                  cs_real_t        r[3][3])
{
  r[0][0] = r6[0]; r[1][1] = r6[1]; r[2][2] = r6[2];
  r[0][1] = r6[3]; r[1][0] = r6[3];
  r[1][2] = r6[4]; r[2][1] = r6[4];
  r[0][2] = r6[5]; r[2][0] = r6[5];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exterior state of a boundary face.
 *
 * Symmetry: mirror state; wall: cell state; other types: state given by
 * the boundary condition coefficients (imposed value for a Dirichlet
 * condition, cell value for a homogeneous Neumann condition).
 *
 * \param[in]   bc_type    boundary condition type
 * \param[in]   n          unit outward normal
 * \param[in]   u          cell velocity
 * \param[in]   r6         cell Reynolds stresses
 * \param[in]   coefa_u    velocity boundary coefficient a
 * \param[in]   coefb_u    velocity boundary coefficient b
 * \param[in]   coefa_r    Reynolds stress boundary coefficient a
 * \param[in]   coefb_r    Reynolds stress boundary coefficient b
 * \param[out]  u_ext      exterior velocity
 * \param[out]  r_ext      exterior Reynolds stresses (full tensor)
 */
/*----------------------------------------------------------------------------*/

static CS_F_HOST_DEVICE void
_proj_boundary_state(int              bc_type,
                     const cs_real_t  n[3],
                     const cs_real_t  u[3],
                     const cs_real_t  r6[6],
                     const cs_real_t  coefa_u[3],
                     const cs_real_t  coefb_u[3][3],
                     const cs_real_t  coefa_r[6],
                     const cs_real_t  coefb_r[6][6],
                     cs_real_t        u_ext[3],
                     cs_real_t        r_ext[3][3])
{
  if (bc_type == CS_SYMMETRY) {
    cs_real_t r[3][3], rn[3];
    _proj_sym_to_full(r6, r);
    const cs_real_t un = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
    for (int i = 0; i < 3; i++) {
      u_ext[i] = u[i] - 2.*un*n[i];
      rn[i] = r[i][0]*n[0] + r[i][1]*n[1] + r[i][2]*n[2];
    }
    const cs_real_t rnn = rn[0]*n[0] + rn[1]*n[1] + rn[2]*n[2];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        r_ext[i][j] =   r[i][j] - 2.*n[i]*rn[j] - 2.*rn[i]*n[j]
                      + 4.*rnn*n[i]*n[j];
    }
  }

  else if (bc_type == CS_SMOOTHWALL || bc_type == CS_ROUGHWALL) {
    for (int i = 0; i < 3; i++)
      u_ext[i] = u[i];
    _proj_sym_to_full(r6, r_ext);
  }

  else {
    cs_real_t r6_ext[6];
    for (int i = 0; i < 3; i++) {
      u_ext[i] = coefa_u[i];
      for (int j = 0; j < 3; j++)
        u_ext[i] += coefb_u[i][j]*u[j];
    }
    for (int ij = 0; ij < 6; ij++) {
      r6_ext[ij] = coefa_r[ij];
      for (int kl = 0; kl < 6; kl++)
        r6_ext[ij] += coefb_r[ij][kl]*r6[kl];
    }
    _proj_sym_to_full(r6_ext, r_ext);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exterior state of the scalar part (mean value, variance and
 *        turbulent flux) at a boundary face.
 *
 * Symmetry: scalars unchanged, turbulent flux mirrored like a vector; wall:
 * cell state; other types: state given by the boundary condition
 * coefficients.
 *
 * \param[in]   bc_type    boundary condition type
 * \param[in]   n          unit outward normal
 * \param[in]   t          cell mean value of the scalar
 * \param[in]   tv         cell variance
 * \param[in]   q          cell turbulent flux
 * \param[in]   coefa_t    scalar boundary coefficient a
 * \param[in]   coefb_t    scalar boundary coefficient b
 * \param[in]   coefa_tv   variance boundary coefficient a
 * \param[in]   coefb_tv   variance boundary coefficient b
 * \param[in]   coefa_q    turbulent flux boundary coefficient a
 * \param[in]   coefb_q    turbulent flux boundary coefficient b
 * \param[out]  t_ext      exterior mean value
 * \param[out]  tv_ext     exterior variance
 * \param[out]  q_ext      exterior turbulent flux
 */
/*----------------------------------------------------------------------------*/

static CS_F_HOST_DEVICE void
_proj_boundary_state_scalar(int              bc_type,
                            const cs_real_t  n[3],
                            cs_real_t        t,
                            cs_real_t        tv,
                            const cs_real_t  q[3],
                            cs_real_t        coefa_t,
                            cs_real_t        coefb_t,
                            cs_real_t        coefa_tv,
                            cs_real_t        coefb_tv,
                            const cs_real_t  coefa_q[3],
                            const cs_real_t  coefb_q[3][3],
                            cs_real_t       *t_ext,
                            cs_real_t       *tv_ext,
                            cs_real_t        q_ext[3])
{
  if (bc_type == CS_SYMMETRY) {
    const cs_real_t qn = q[0]*n[0] + q[1]*n[1] + q[2]*n[2];
    *t_ext = t;
    *tv_ext = tv;
    for (int i = 0; i < 3; i++)
      q_ext[i] = q[i] - 2.*qn*n[i];
  }

  else if (bc_type == CS_SMOOTHWALL || bc_type == CS_ROUGHWALL) {
    *t_ext = t;
    *tv_ext = tv;
    for (int i = 0; i < 3; i++)
      q_ext[i] = q[i];
  }

  else {
    *t_ext = coefa_t + coefb_t*t;
    *tv_ext = coefa_tv + coefb_tv*tv;
    for (int i = 0; i < 3; i++) {
      q_ext[i] = coefa_q[i];
      for (int j = 0; j < 3; j++)
        q_ext[i] += coefb_q[i][j]*q[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mean state U and generalised covariance C of a cell.
 *
 * U = (u, theta) and C = [[R, u'theta'], [u'theta', theta'^2]] when D = 4,
 * U = u and C = R when D = 3.
 *
 * \param[in]   u    velocity
 * \param[in]   r6   Reynolds stresses (11, 22, 33, 12, 23, 13)
 * \param[in]   t    scalar mean value (D = 4 only)
 * \param[in]   tv   scalar variance (D = 4 only)
 * \param[in]   q    turbulent flux (D = 4 only)
 * \param[out]  uu   mean state
 * \param[out]  cc   generalised covariance
 */
/*----------------------------------------------------------------------------*/

template <int D>
static CS_F_HOST_DEVICE void
_proj_cell_state(const cs_real_t  u[3],
                 const cs_real_t  r6[6],
                 cs_real_t        t,
                 cs_real_t        tv,
                 const cs_real_t  q[3],
                 cs_real_t        uu[D],
                 cs_real_t        cc[D][D])
{
  cs_real_t r[3][3];
  _proj_sym_to_full(r6, r);

  for (int i = 0; i < 3; i++) {
    uu[i] = u[i];
    for (int j = 0; j < 3; j++)
      cc[i][j] = r[i][j];
  }

  if constexpr (D == 4) {
    uu[3] = t;
    for (int i = 0; i < 3; i++) {
      cc[i][3] = q[i];
      cc[3][i] = q[i];
    }
    cc[3][3] = tv;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scalar part (mean value, turbulent flux, variance) of a rebuilt
 *        strip state, with weight s_f.
 *
 * \param[in]   u_s    rebuilt mean state (D = 4)
 * \param[in]   c_s    rebuilt generalised covariance (D = 4)
 * \param[in]   s_f    weight
 * \param[out]  a_s    weighted values (5)
 */
/*----------------------------------------------------------------------------*/

static CS_F_HOST_DEVICE void
_proj_scalar_strip_values(const cs_real_t  u_s[4],
                          const cs_real_t  c_s[4][4],
                          cs_real_t        s_f,
                          cs_real_t        a_s[5])
{
  a_s[0] = s_f*u_s[3];
  for (int i = 0; i < 3; i++)
    a_s[1 + i] = s_f*0.5*(c_s[i][3] + c_s[3][i]);
  a_s[4] = s_f*c_s[3][3];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Projected fields of a scalar with a DFM turbulent flux.
 *
 * \param[in]   f      scalar field
 * \param[out]  p_th   projected mean value
 * \param[out]  p_tv   projected variance
 * \param[out]  p_ut   projected turbulent flux
 */
/*----------------------------------------------------------------------------*/

static void
_proj_scalar_fields(const cs_field_t   *f,
                    cs_field_t        **p_th,
                    cs_field_t        **p_tv,
                    cs_field_t        **p_ut)
{
  const std::string name = std::string("algo:proj_") + f->name;
  *p_th = cs_field_by_name(name.c_str());
  *p_tv = cs_field_by_name((name + "_variance").c_str());
  *p_ut = cs_field_by_name((name + "_turbulent_flux").c_str());
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief One sub-step of the projection of a scalar with a DFM turbulent
 *        flux.
 *
 * The Riemann problem decouples scalar by scalar: the waves (positions
 * and strips) are those of (u, R), and the scalar components of the
 * solution only involve the scalar itself. The scalar is projected with
 * the Riemann problem of (u, theta, R, u'theta', theta'^2), from the
 * velocity and Reynolds stresses at the beginning of the sub-step.
 *
 * \param[in]       f        scalar field
 * \param[in]       vel      velocity at the beginning of the sub-step
 * \param[in]       rij      Reynolds stresses at the beginning of the
 *                           sub-step
 * \param[in]       i_speed  largest wave speed of the interior faces
 * \param[in]       b_speed  largest wave speed of the boundary faces
 * \param[in]       s_sum    sum of |f| speed_f over the faces of the cells
 * \param[in]       dt_s     sub-step
 * \param[in, out]  acc      work array (5 values per cell)
 */
/*----------------------------------------------------------------------------*/

static void
_proj_scalar_sub_step(const cs_field_t   *f,
                      const cs_real_3_t   vel[],
                      const cs_real_6_t   rij[],
                      const cs_real_t     i_speed[],
                      const cs_real_t     b_speed[],
                      const cs_real_t     s_sum[],
                      cs_real_t           dt_s,
                      cs_real_t           acc[])
{
  constexpr int stride = 5;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_real_t *i_face_surf = mq->i_face_surf;
  const cs_real_t *b_face_surf = mq->b_face_surf;
  const cs_real_t *weight = mq->weight;
  const cs_real_t *cell_vol = mq->cell_vol;

  const int *bc_type = cs_glob_bc_type;

  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rij = CS_F_(rij);

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *imasfl
    = cs_field_by_id(f_vel->get_key_int(kimasf))->val;
  const cs_real_t *bmasfl
    = cs_field_by_id(f_vel->get_key_int(kbmasf))->val;

  const cs_real_3_t *coefa_u = (const cs_real_3_t *)f_vel->bc_coeffs->a;
  const cs_real_33_t *coefb_u = (const cs_real_33_t *)f_vel->bc_coeffs->b;
  const cs_real_6_t *coefa_r = (const cs_real_6_t *)f_rij->bc_coeffs->a;
  const cs_real_66_t *coefb_r = (const cs_real_66_t *)f_rij->bc_coeffs->b;

  /* Scalar, variance and turbulent flux: projected fields (work arrays)
     and boundary coefficients of the transported fields */

  const cs_field_t *f_tv = cs_field_get_variance(f);
  const cs_field_t *f_ut = cs_field_by_composite_name(f->name,
                                                      "turbulent_flux");

  cs_field_t *p_th, *p_tv, *p_ut;
  _proj_scalar_fields(f, &p_th, &p_tv, &p_ut);

  cs_real_t *th = p_th->val;
  cs_real_t *tv = p_tv->val;
  cs_real_3_t *ut = (cs_real_3_t *)p_ut->val;

  const cs_real_t *coefa_t = f->bc_coeffs->a;
  const cs_real_t *coefb_t = f->bc_coeffs->b;
  const cs_real_t *coefa_tv = f_tv->bc_coeffs->a;
  const cs_real_t *coefb_tv = f_tv->bc_coeffs->b;
  const cs_real_3_t *coefa_q = (const cs_real_3_t *)f_ut->bc_coeffs->a;
  const cs_real_33_t *coefb_q = (const cs_real_33_t *)f_ut->bc_coeffs->b;

  cs_dispatch_context ctx;
  const cs_dispatch_sum_type_t i_sum_type
    = ctx.get_parallel_for_i_faces_sum_type(m);
  const cs_dispatch_sum_type_t b_sum_type
    = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (int i = 0; i < stride; i++)
      acc[stride*c_id + i] = 0.;
  });
  ctx.wait();

  /* Strip means, accumulated with the weights |f| speed_f */

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    const cs_real_t sp = i_speed[f_id];
    if (!(sp > 0.))
      return;

    const cs_lnum_t c_l = i_face_cells[f_id][0];
    const cs_lnum_t c_r = i_face_cells[f_id][1];

    const cs_real_t n[4] = {i_face_u_normal[f_id][0],
                            i_face_u_normal[f_id][1],
                            i_face_u_normal[f_id][2],
                            0.};
    const cs_real_t rho_f
      = weight[f_id]*crom[c_l] + (1. - weight[f_id])*crom[c_r];
    const cs_real_t v = imasfl[f_id]/(rho_f*i_face_surf[f_id]);
    const cs_real_t s_f = sp*i_face_surf[f_id];

    cs_real_t ul[4], cl[4][4], ur[4], cr[4][4];
    _proj_cell_state<4>(vel[c_l], rij[c_l], th[c_l], tv[c_l], ut[c_l],
                        ul, cl);
    _proj_cell_state<4>(vel[c_r], rij[c_r], th[c_r], tv[c_r], ut[c_r],
                        ur, cr);

    _proj_rp_t<4> rp;
    _proj_solve<4>(n, v, ul, cl, ur, cr, rp);

    for (int side = 0; side < 2; side++) {
      const cs_lnum_t c_id = i_face_cells[f_id][side];
      if (c_id >= n_cells)
        continue;

      /* Strip width in xi: thickness |K| speed_f / sum / dt_s */
      const cs_real_t w = cell_vol[c_id]*sp/(s_sum[c_id]*dt_s);
      const cs_real_t xa = (side == 0) ? -w : 0.;
      const cs_real_t xb = (side == 0) ? 0. : w;

      cs_real_t u_s[4], c_s[4][4], a_s[stride];
      _proj_strip_mean<4>(n, rp, xa, xb, u_s, c_s);
      _proj_scalar_strip_values(u_s, c_s, s_f, a_s);

      cs_dispatch_sum<stride>(&acc[stride*c_id], a_s, i_sum_type);
    }
  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    const cs_real_t sp = b_speed[f_id];
    if (!(sp > 0.))
      return;

    const cs_lnum_t c_id = b_face_cells[f_id];

    const cs_real_t n[4] = {b_face_u_normal[f_id][0],
                            b_face_u_normal[f_id][1],
                            b_face_u_normal[f_id][2],
                            0.};
    const cs_real_t v = bmasfl[f_id]/(brom[f_id]*b_face_surf[f_id]);
    const cs_real_t s_f = sp*b_face_surf[f_id];
    const cs_real_t w = cell_vol[c_id]*sp/(s_sum[c_id]*dt_s);

    cs_real_t u_ext[3], r_ext[3][3];
    _proj_boundary_state(bc_type[f_id], n, vel[c_id], rij[c_id],
                         coefa_u[f_id], coefb_u[f_id],
                         coefa_r[f_id], coefb_r[f_id],
                         u_ext, r_ext);

    cs_real_t t_ext, tv_ext, q_ext[3];
    _proj_boundary_state_scalar(bc_type[f_id], n, th[c_id], tv[c_id],
                                ut[c_id], coefa_t[f_id], coefb_t[f_id],
                                coefa_tv[f_id], coefb_tv[f_id],
                                coefa_q[f_id], coefb_q[f_id],
                                &t_ext, &tv_ext, q_ext);

    const cs_real_t r6_ext[6] = {r_ext[0][0], r_ext[1][1], r_ext[2][2],
                                 r_ext[0][1], r_ext[1][2], r_ext[0][2]};

    cs_real_t ul[4], cl[4][4], ur[4], cr[4][4];
    _proj_cell_state<4>(vel[c_id], rij[c_id], th[c_id], tv[c_id], ut[c_id],
                        ul, cl);
    _proj_cell_state<4>(u_ext, r6_ext, t_ext, tv_ext, q_ext, ur, cr);

    _proj_rp_t<4> rp;
    _proj_solve<4>(n, v, ul, cl, ur, cr, rp);

    cs_real_t u_s[4], c_s[4][4], a_s[stride];
    _proj_strip_mean<4>(n, rp, -w, 0., u_s, c_s);
    _proj_scalar_strip_values(u_s, c_s, s_f, a_s);

    cs_dispatch_sum<stride>(&acc[stride*c_id], a_s, b_sum_type);
  });

  ctx.wait();

  /* New cell states: convex combination of the strip states */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    if (!(s_sum[c_id] > 0.))
      return;

    const cs_real_t inv_s = 1./s_sum[c_id];
    const cs_real_t *a_c = acc + stride*c_id;

    th[c_id] = a_c[0]*inv_s;
    for (int i = 0; i < 3; i++)
      ut[c_id][i] = a_c[1 + i]*inv_s;
    tv[c_id] = a_c[4]*inv_s;
  });
  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Projection over one time step.
 *
 * The projected fields ("algo:proj_*") hold the states at the beginning of
 * the time step on input, and the projected states on output (local
 * cells). The velocity, Reynolds stresses and epsilon are projected with
 * the Riemann problem of (u, R); then, in each sub-step, every scalar with
 * a DFM turbulent flux is projected in turn (see _proj_scalar_sub_step).
 *
 * \param[in]  dfm_scalars  scalar fields with a DFM turbulent flux
 *
 * \return  number of sub-steps
 */
/*----------------------------------------------------------------------------*/

static int
_projection(const std::vector<const cs_field_t *>  &dfm_scalars)
{
  /* Accumulated values per cell: U (3), R (6), eps k^(-C_eps1) */
  constexpr int stride = 10;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_halo_t *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_nreal_3_t *i_face_u_normal = mq->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_real_t *i_face_surf = mq->i_face_surf;
  const cs_real_t *b_face_surf = mq->b_face_surf;
  const cs_real_t *weight = mq->weight;
  const cs_real_t *cell_vol = mq->cell_vol;

  const int *bc_type = cs_glob_bc_type;

  const cs_real_t ce1 = cs_turb_ce1;

  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *imasfl
    = cs_field_by_id(f_vel->get_key_int(kimasf))->val;
  const cs_real_t *bmasfl
    = cs_field_by_id(f_vel->get_key_int(kbmasf))->val;

  const cs_real_3_t *coefa_u = (const cs_real_3_t *)f_vel->bc_coeffs->a;
  const cs_real_33_t *coefb_u = (const cs_real_33_t *)f_vel->bc_coeffs->b;
  const cs_real_6_t *coefa_r = (const cs_real_6_t *)f_rij->bc_coeffs->a;
  const cs_real_66_t *coefb_r = (const cs_real_66_t *)f_rij->bc_coeffs->b;
  const cs_real_t *coefa_e = f_eps->bc_coeffs->a;
  const cs_real_t *coefb_e = f_eps->bc_coeffs->b;

  /* Projected fields, used as work arrays */

  cs_real_3_t *vel
    = (cs_real_3_t *)cs_field_by_name("algo:proj_velocity")->val;
  cs_real_6_t *rij = (cs_real_6_t *)cs_field_by_name("algo:proj_rij")->val;
  cs_real_t *eps = cs_field_by_name("algo:proj_epsilon")->val;

  /* Uniform time step (local time steps are refused at setup) */

  const cs_real_t dt = CS_F_(dt)->val[0];

  cs_dispatch_context ctx;
  const cs_dispatch_sum_type_t i_sum_type
    = ctx.get_parallel_for_i_faces_sum_type(m);
  const cs_dispatch_sum_type_t b_sum_type
    = ctx.get_parallel_for_b_faces_sum_type(m);

  cs_real_t *w_c = nullptr, *i_speed = nullptr, *b_speed = nullptr;
  cs_real_t *s_sum = nullptr, *acc = nullptr, *acc_s = nullptr;

  CS_MALLOC_HD(w_c, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(i_speed, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_speed, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(s_sum, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(acc, stride*n_cells_ext, cs_real_t, cs_alloc_mode);
  if (! dfm_scalars.empty())
    CS_MALLOC_HD(acc_s, 5*n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t t_done = 0.;
  int n_sub_steps = 0;

  while (t_done < dt && n_sub_steps < _N_SUB_STEPS_MAX) {

    if (halo != nullptr) {
      cs_halo_sync_r(halo, ctx.use_gpu(), vel);
      cs_halo_sync_r(halo, ctx.use_gpu(), rij);
      cs_halo_sync(halo, ctx.use_gpu(), eps);
      for (const cs_field_t *f : dfm_scalars) {
        cs_field_t *p_th, *p_tv, *p_ut;
        _proj_scalar_fields(f, &p_th, &p_tv, &p_ut);
        cs_halo_sync(halo, ctx.use_gpu(), p_th->val);
        cs_halo_sync(halo, ctx.use_gpu(), p_tv->val);
        cs_halo_sync_r(halo, ctx.use_gpu(), (cs_real_3_t *)p_ut->val);
      }
    }

    /* eps k^(-C_eps1), passive scalar of the subsystem */

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t k = 0.5*(rij[c_id][0] + rij[c_id][1] + rij[c_id][2]);
      w_c[c_id] = (k > 0.) ? eps[c_id]*pow(k, -ce1) : 0.;
      s_sum[c_id] = 0.;
      for (int i = 0; i < stride; i++)
        acc[stride*c_id + i] = 0.;
    });
    ctx.wait();

    /* Largest wave speed of every face (it depends on u_n and R_nn only),
       and their sums over the faces of every cell, weighted by the face
       areas */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      const cs_lnum_t c_l = i_face_cells[f_id][0];
      const cs_lnum_t c_r = i_face_cells[f_id][1];

      const cs_real_t n[3] = {i_face_u_normal[f_id][0],
                              i_face_u_normal[f_id][1],
                              i_face_u_normal[f_id][2]};
      const cs_real_t rho_f
        = weight[f_id]*crom[c_l] + (1. - weight[f_id])*crom[c_r];
      const cs_real_t v = imasfl[f_id]/(rho_f*i_face_surf[f_id]);

      cs_real_t rl[3][3], rr[3][3];
      _proj_sym_to_full(rij[c_l], rl);
      _proj_sym_to_full(rij[c_r], rr);

      const cs_real_t sp
        = _proj_wave_speed<3>(n, v, vel[c_l], rl, vel[c_r], rr);
      i_speed[f_id] = sp;

      const cs_real_t s_f = sp*i_face_surf[f_id];
      if (c_l < n_cells)
        cs_dispatch_sum(&s_sum[c_l], s_f, i_sum_type);
      if (c_r < n_cells)
        cs_dispatch_sum(&s_sum[c_r], s_f, i_sum_type);
    });

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      const cs_lnum_t c_id = b_face_cells[f_id];

      const cs_real_t n[3] = {b_face_u_normal[f_id][0],
                              b_face_u_normal[f_id][1],
                              b_face_u_normal[f_id][2]};
      const cs_real_t v = bmasfl[f_id]/(brom[f_id]*b_face_surf[f_id]);

      cs_real_t rl[3][3], u_ext[3], r_ext[3][3];
      _proj_sym_to_full(rij[c_id], rl);
      _proj_boundary_state(bc_type[f_id], n, vel[c_id], rij[c_id],
                           coefa_u[f_id], coefb_u[f_id],
                           coefa_r[f_id], coefb_r[f_id],
                           u_ext, r_ext);

      const cs_real_t sp
        = _proj_wave_speed<3>(n, v, vel[c_id], rl, u_ext, r_ext);
      b_speed[f_id] = sp;

      cs_dispatch_sum(&s_sum[c_id], sp*b_face_surf[f_id], b_sum_type);
    });

    ctx.wait();

    /* Sub-step: dt_s sum_f |f| speed_f <= |K| in every cell, which keeps
       the waves inside the strips */

    cs_real_t r_max = 0.;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      r_max = cs::max(r_max, s_sum[c_id]/cell_vol[c_id]);
    cs_parall_max(1, CS_REAL_TYPE, &r_max);

    const cs_real_t t_rem = dt - t_done;
    const bool last_sub_step = !(r_max*t_rem > 1.);
    const cs_real_t dt_s = (last_sub_step) ? t_rem : t_rem/ceil(r_max*t_rem);

    /* Velocity, Reynolds stresses and epsilon: strip means, accumulated
       with the weights |f| speed_f */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      const cs_real_t sp = i_speed[f_id];
      if (!(sp > 0.))
        return;

      const cs_lnum_t c_l = i_face_cells[f_id][0];
      const cs_lnum_t c_r = i_face_cells[f_id][1];

      const cs_real_t n[3] = {i_face_u_normal[f_id][0],
                              i_face_u_normal[f_id][1],
                              i_face_u_normal[f_id][2]};
      const cs_real_t rho_f
        = weight[f_id]*crom[c_l] + (1. - weight[f_id])*crom[c_r];
      const cs_real_t v = imasfl[f_id]/(rho_f*i_face_surf[f_id]);
      const cs_real_t s_f = sp*i_face_surf[f_id];

      cs_real_t rl[3][3], rr[3][3];
      _proj_sym_to_full(rij[c_l], rl);
      _proj_sym_to_full(rij[c_r], rr);

      _proj_rp_t<3> rp;
      _proj_solve<3>(n, v, vel[c_l], rl, vel[c_r], rr, rp);

      for (int side = 0; side < 2; side++) {
        const cs_lnum_t c_id = i_face_cells[f_id][side];
        if (c_id >= n_cells)
          continue;

        /* Strip width in xi: thickness |K| speed_f / sum / dt_s */
        const cs_real_t w = cell_vol[c_id]*sp/(s_sum[c_id]*dt_s);
        const cs_real_t xa = (side == 0) ? -w : 0.;
        const cs_real_t xb = (side == 0) ? 0. : w;

        cs_real_t u_s[3], c_s[3][3], a_s[stride];
        _proj_strip_mean<3>(n, rp, xa, xb, u_s, c_s);

        for (int i = 0; i < 3; i++)
          a_s[i] = s_f*u_s[i];
        a_s[3] = s_f*c_s[0][0];
        a_s[4] = s_f*c_s[1][1];
        a_s[5] = s_f*c_s[2][2];
        a_s[6] = s_f*0.5*(c_s[0][1] + c_s[1][0]);
        a_s[7] = s_f*0.5*(c_s[1][2] + c_s[2][1]);
        a_s[8] = s_f*0.5*(c_s[0][2] + c_s[2][0]);
        a_s[9] = s_f*_proj_strip_mean_contact(v, xa, xb,
                                              w_c[c_l], w_c[c_r]);

        cs_dispatch_sum<stride>(&acc[stride*c_id], a_s, i_sum_type);
      }
    });

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      const cs_real_t sp = b_speed[f_id];
      if (!(sp > 0.))
        return;

      const cs_lnum_t c_id = b_face_cells[f_id];

      const cs_real_t n[3] = {b_face_u_normal[f_id][0],
                              b_face_u_normal[f_id][1],
                              b_face_u_normal[f_id][2]};
      const cs_real_t v = bmasfl[f_id]/(brom[f_id]*b_face_surf[f_id]);
      const cs_real_t s_f = sp*b_face_surf[f_id];
      const cs_real_t w = cell_vol[c_id]*sp/(s_sum[c_id]*dt_s);

      cs_real_t rl[3][3], u_ext[3], r_ext[3][3];
      _proj_sym_to_full(rij[c_id], rl);
      _proj_boundary_state(bc_type[f_id], n, vel[c_id], rij[c_id],
                           coefa_u[f_id], coefb_u[f_id],
                           coefa_r[f_id], coefb_r[f_id],
                           u_ext, r_ext);

      /* eps k^(-C_eps1) of the exterior state */
      cs_real_t w_ext = w_c[c_id];
      if (   bc_type[f_id] != CS_SYMMETRY
          && bc_type[f_id] != CS_SMOOTHWALL
          && bc_type[f_id] != CS_ROUGHWALL) {
        const cs_real_t k_ext
          = 0.5*(r_ext[0][0] + r_ext[1][1] + r_ext[2][2]);
        const cs_real_t e_ext = coefa_e[f_id] + coefb_e[f_id]*eps[c_id];
        w_ext = (k_ext > 0.) ? e_ext*pow(k_ext, -ce1) : 0.;
      }

      _proj_rp_t<3> rp;
      _proj_solve<3>(n, v, vel[c_id], rl, u_ext, r_ext, rp);

      cs_real_t u_s[3], c_s[3][3], a_s[stride];
      _proj_strip_mean<3>(n, rp, -w, 0., u_s, c_s);

      for (int i = 0; i < 3; i++)
        a_s[i] = s_f*u_s[i];
      a_s[3] = s_f*c_s[0][0];
      a_s[4] = s_f*c_s[1][1];
      a_s[5] = s_f*c_s[2][2];
      a_s[6] = s_f*0.5*(c_s[0][1] + c_s[1][0]);
      a_s[7] = s_f*0.5*(c_s[1][2] + c_s[2][1]);
      a_s[8] = s_f*0.5*(c_s[0][2] + c_s[2][0]);
      a_s[9] = s_f*_proj_strip_mean_contact(v, -w, 0., w_c[c_id], w_ext);

      cs_dispatch_sum<stride>(&acc[stride*c_id], a_s, b_sum_type);
    });

    ctx.wait();

    /* Scalars with a DFM turbulent flux, one after the other, from the
       velocity and Reynolds stresses at the beginning of the sub-step */

    for (const cs_field_t *f : dfm_scalars)
      _proj_scalar_sub_step(f, vel, rij, i_speed, b_speed, s_sum, dt_s,
                            acc_s);

    /* New cell states: convex combination of the strip states; then
       eps = (mean of eps k^(-C_eps1)) k^C_eps1 with the new cell k */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      if (!(s_sum[c_id] > 0.))
        return;

      const cs_real_t inv_s = 1./s_sum[c_id];
      const cs_real_t *a_c = acc + stride*c_id;

      for (int i = 0; i < 3; i++)
        vel[c_id][i] = a_c[i]*inv_s;
      for (int ij = 0; ij < 6; ij++)
        rij[c_id][ij] = a_c[3 + ij]*inv_s;

      const cs_real_t kk = 0.5*(rij[c_id][0] + rij[c_id][1] + rij[c_id][2]);
      const cs_real_t w = a_c[9]*inv_s;
      if (kk > 0. && w > 0.)
        eps[c_id] = w*pow(kk, ce1);
    });
    ctx.wait();

    t_done = (last_sub_step) ? dt : t_done + dt_s;
    n_sub_steps += 1;

  }

  if (t_done < dt)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: the projection needs more than %d sub-steps\n"
                "to cover the time step."),
              __func__, _N_SUB_STEPS_MAX);

  CS_FREE(w_c);
  CS_FREE(i_speed);
  CS_FREE(b_speed);
  CS_FREE(s_sum);
  CS_FREE(acc);
  CS_FREE(acc_s);

  return n_sub_steps;
}

/*! (end of DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the mean velocity, the Reynolds stresses and the
 *        dissipation rate (and the scalars with a DFM turbulent flux, their
 *        variances and their turbulent fluxes) over one time step of the
 *        convective subsystem.
 *
 * The projected states are stored in the cell fields
 * "algo:proj_velocity", "algo:proj_rij", "algo:proj_epsilon" and, for each
 * DFM scalar, "algo:proj_<scalar>", "algo:proj_<scalar>_variance" and
 * "algo:proj_<scalar>_turbulent_flux" (created in cs_setup.cpp when
 * rij_discretization_scheme == CS_RIJ_SCHEME_GODUNOV_PROJECTION). They are
 * computed once per time step, from the fields at the beginning of the time
 * step; further calls in the same time step return immediately.
 *
 * The convective subsystem gathers the convection by the mass flux of the
 * mean velocity, of the scalars, of the second moments and of epsilon, the
 * divergence of the second moments in the mean equations, the production
 * of the second moments by the mean gradients and the production of
 * epsilon. The Riemann problem of each face is solved exactly with the
 * jump relation J_T; its solution is averaged over one strip in each
 * adjacent cell, whose thickness is proportional to the largest wave speed
 * of the face, and the new cell state is the convex combination of the
 * strip states. The projection is sub-cycled when the time step exceeds
 * the Courant condition of the strips. The covariances between different
 * scalars are not transported; they do not enter the projection of the
 * other components.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_projection(void)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  if (ts->nt_cur == _nt_last)
    return;
  _nt_last = ts->nt_cur;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;

  cs_dispatch_context ctx;

  /* Scalars with a DFM turbulent flux */

  std::vector<const cs_field_t *> dfm_scalars;
  {
    const int k_tfm = cs_field_key_id("turbulent_flux_model");
    const int n_fields = cs_field_n_fields();
    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->get_key_int(k_tfm) / 10 == 3)
        dfm_scalars.push_back(f);
    }
  }

  /* Projected fields initialised with the states at the beginning of the
     time step */

  cs_real_3_t *vel
    = (cs_real_3_t *)cs_field_by_name("algo:proj_velocity")->val;
  cs_real_6_t *rij = (cs_real_6_t *)cs_field_by_name("algo:proj_rij")->val;
  cs_real_t *eps = cs_field_by_name("algo:proj_epsilon")->val;

  {
    const cs_real_3_t *vel_pre = (const cs_real_3_t *)CS_F_(vel)->val_pre;
    const cs_real_6_t *rij_pre = (const cs_real_6_t *)CS_F_(rij)->val_pre;
    const cs_real_t *eps_pre = CS_F_(eps)->val_pre;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (int i = 0; i < 3; i++)
        vel[c_id][i] = vel_pre[c_id][i];
      for (int ij = 0; ij < 6; ij++)
        rij[c_id][ij] = rij_pre[c_id][ij];
      eps[c_id] = eps_pre[c_id];
    });
    ctx.wait();
  }

  for (const cs_field_t *f : dfm_scalars) {
    cs_field_t *p_th, *p_tv, *p_ut;
    _proj_scalar_fields(f, &p_th, &p_tv, &p_ut);

    cs_real_t *th = p_th->val, *tv = p_tv->val;
    cs_real_3_t *ut = (cs_real_3_t *)p_ut->val;
    const cs_real_t *th_pre = f->val_pre;
    const cs_real_t *tv_pre = cs_field_get_variance(f)->val_pre;
    const cs_real_3_t *ut_pre = (const cs_real_3_t *)
      cs_field_by_composite_name(f->name, "turbulent_flux")->val_pre;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      th[c_id] = th_pre[c_id];
      tv[c_id] = tv_pre[c_id];
      for (int i = 0; i < 3; i++)
        ut[c_id][i] = ut_pre[c_id][i];
    });
    ctx.wait();
  }

  const int n_sub_steps = _projection(dfm_scalars);

  /* Realisability: smallest eigenvalue of the projected Reynolds
     stresses, relative to k, and sign of the projected variances */

  if (cs_glob_turb_rans_model->verbosity >= 1) {

    cs_array<bool> is_neg(n_cells, cs_alloc_mode);

    struct cs_data_double_n<1> lambda_min;
    struct cs_reduce_min_nr<1> reducer_lambda;

    // FIXME: cs_math_sym_33_eigen needs to run on GPU first.
    ctx.set_use_gpu(false);

    ctx.parallel_for_reduce(n_cells, lambda_min, reducer_lambda, CS_LAMBDA
                            (cs_lnum_t c_id, cs_data_double_n<1> &res) {
      cs_real_t ev[3];
      cs_math_sym_33_eigen(rij[c_id], ev);
      const cs_real_t k = 0.5*(rij[c_id][0] + rij[c_id][1] + rij[c_id][2]);
      const cs_real_t l_min = cs::min(ev[0], cs::min(ev[1], ev[2]));
      const cs_real_t l_rel = (k > 0.) ? l_min/k : l_min;
      res.r[0] = l_rel;
      is_neg[c_id] = (l_rel < -1.e-12);
    });
    ctx.wait();

    cs_parall_min(1, CS_DOUBLE, &lambda_min.r[0]);

    for (const cs_field_t *f : dfm_scalars) {
      cs_field_t *p_th, *p_tv, *p_ut;
      _proj_scalar_fields(f, &p_th, &p_tv, &p_ut);
      const cs_real_t * val = f->val;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        if (val[c_id] < 0.)
          is_neg[c_id] = true;
      });
    }

    cs_gnum_t n_neg = 0;
    ctx.parallel_for_reduce_sum(n_cells, n_neg, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum)
    {
      if (is_neg[c_id])
        sum += 1;
    });
    ctx.wait();
    cs_parall_counter(&n_neg, 1);

    bft_printf(_(" ** Godunov-type projection of the convective subsystem\n"
                 "    ----------------------------------------------------\n"
                 "    sub-steps: %d, smallest eigenvalue of R / k: %12.5e,"
                 " non-realisable cells: %llu\n"),
               n_sub_steps, lambda_min.r[0], (unsigned long long)n_neg);
  }
}

/*----------------------------------------------------------------------------*/
