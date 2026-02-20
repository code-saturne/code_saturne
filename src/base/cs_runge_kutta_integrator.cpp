/*============================================================================
 * Explicit Runge-Kutta integrator utilities.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_dispatch.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_runge_kutta_integrator.h"
#include "base/cs_runge_kutta_integrator_priv.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Global variables
 *============================================================================*/

cs_runge_kutta_integrator_t **_rk_lst = nullptr;
int  _n_rk_integrators = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fill integrator coeffs for a given scheme.
 *
 * Reference: Accuracy analysis of explicit Runge–Kutta methods
 * applied to the incompressible Navier–Stokes equations.
 * J.Comput.Phys. 231 (2012) 3041–3063
 */
/*----------------------------------------------------------------------------*/

static void
_runge_kutta_integrator_init_scheme(cs_runge_kutta_integrator_t   *rk,
                                    cs_runge_kutta_scheme_t        scheme)
{
  double *a = rk->rk_coeff.a;
  double *c = rk->rk_coeff.c;
  memset(a, 0., RK_HIGHEST_ORDER*RK_HIGHEST_ORDER*sizeof(double));
  memset(c, 0., RK_HIGHEST_ORDER*sizeof(double));

  switch (scheme) {

  case CS_RK_NONE:
    return;

  case CS_RK1:
    rk->n_stages = 1;
    a[0] = 1.;
    c[0] = 1.;
    break;

  case CS_RK2:
    rk->n_stages = 2;

    a[0] = 1.;
    a[4] = 0.5, a[5] = 0.5;

    c[0] = 1., c[1] = 1.;
    break;

  case CS_RK3:
    rk->n_stages = 3;

    a[0] = 1.0/3.0;
    a[4] = -1.0,    a[5] = 2.0;
    a[8] = 0.0,     a[9] = 3./4.,  a[10] = 1./4.;

    c[0] = 1.0/3.0, c[1] = 1.0, c[2] = 1.0;

    break;

  case CS_RK3_WRAY:
    rk->n_stages = 3;

    a[0] = 8.0/15.0;
    a[4] = 1.0/4.0, a[5] = 5.0/12.0;
    a[8] = 1.0/4.0, a[9] = 0.,  a[10] = 3./4.;

    c[0] = 8.0/15.0, c[1] = 2./3., c[2] = 1.0;

    break;
  case CS_RK3_SSP:

    rk->n_stages = 3;

    a[0] = 1.0;
    a[4] = 1.0/4.0, a[5] = 1.0/4.0;
    a[8] = 1.0/6.0, a[9] = 1.0/6.0, a[10] = 2./3.;

    c[0] = 1.0, c[1] = 1./2., c[2] = 1.0;

    break;
  case CS_RK4:
    rk->n_stages = 4;

    a[0] = 1.0;
    a[4] = 3.0/8.0,  a[5] = 1.0/8.0;
    a[8] = -1.0/8.0, a[9] = -3.0/8.0,   a[10] = 3.0/2.0;
    a[12] = 1.0/6.0, a[13] = -1.0/18.0, a[14] = 2.0/3.0, a[15] = 2.0/9.0;
    c[0] = 1.0, c[1] = 1.0/2.0, c[2] = 1.0, c[3] = 1.0;

    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Type of Runge-Kutta not available.\n"
              "%s: Stop building Runge-Kutta integrator.\n",
              __func__, __func__);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory
 * \param[in]  rk      double pointer to a runge-kutta integrator
 */
/*----------------------------------------------------------------------------*/

static void
_runge_kutta_integrator_free(cs_runge_kutta_integrator_t  **rk)
{
  cs_runge_kutta_integrator_t *_rk = *rk;

  if (_rk == nullptr)
    return;

  CS_FREE_HD(_rk->rhs_stages);
  CS_FREE_HD(_rk->rk_coeff.a);
  CS_FREE_HD(_rk->rk_coeff.c);

  CS_FREE(_rk->name);

  CS_FREE_HD(_rk->scaled_dt);
  CS_FREE_HD(_rk->mass);

  CS_FREE_HD(_rk->u_old);
  CS_FREE_HD(_rk->u_new);
  CS_FREE(_rk);

  *rk = nullptr;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a RK integrator.
 *
 * \param[in] scheme  RK scheme type (RK_NONE, RK1, RK2, RK3, RK4)
 * \param[in] name    associated equation or field's name
 * \param[in] dt      time step
 * \param[in] dim     variable dimentsion
 * \param[in] n_elts  number of computational elements
 *
 * return the RK integrator's id in the RK list
 */
/*----------------------------------------------------------------------------*/

int
cs_runge_kutta_integrator_create(cs_runge_kutta_scheme_t      scheme,
                                 const char                  *name,
                                 const cs_real_t             *dt,
                                 int                          dim,
                                 cs_lnum_t                    n_elts)
{
  cs_runge_kutta_integrator_t *rk = nullptr;
  cs_lnum_t stride = dim;

  CS_MALLOC(rk, 1, cs_runge_kutta_integrator_t);

  rk->scheme = scheme;

  CS_MALLOC(rk->name, strlen(name) + 1, char);
  strcpy(rk->name, name);

  rk->dt = dt;

  rk->n_elts = n_elts;
  rk->n_stages = 1;
  rk->i_stage = 0;

  CS_MALLOC_HD(rk->rk_coeff.a, RK_HIGHEST_ORDER*RK_HIGHEST_ORDER,
               double, cs_alloc_mode);
  CS_MALLOC_HD(rk->rk_coeff.c, RK_HIGHEST_ORDER,
               double, cs_alloc_mode);

  _runge_kutta_integrator_init_scheme(rk, scheme);

  CS_MALLOC_HD(rk->rhs_stages, rk->n_stages*stride*n_elts,
               cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->u_old, stride*n_elts, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->scaled_dt, n_elts, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(rk->mass, n_elts, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(rk->u_new, stride*n_elts, cs_real_t, cs_alloc_mode);

  int _rk_integrator_id = _n_rk_integrators;

  CS_REALLOC(_rk_lst, _n_rk_integrators + 1,
             cs_runge_kutta_integrator_t *);

  _rk_lst[_rk_integrator_id] = rk;

  _n_rk_integrators++;

  return _rk_integrator_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create RK integrator structures
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_initialize()
{
  const int n_fields = cs_field_n_fields();
  const cs_real_t *dt = CS_F_(dt)->val;

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      if (eqp->rk_def.scheme <= CS_RK_NONE)
        continue;

      if (eqp->dim == 1 || eqp->dim == 3 || eqp->dim == 6) {
        eqp->rk_def.rk_id
          = cs_runge_kutta_integrator_create
              (eqp->rk_def.scheme, f->name, dt, eqp->dim, n_cells_ext);
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Equation's dimension for Runge-Kutta not available.\n"
                  "%s: Stop building Runge-Kutta integrator.\n",
                  __func__, __func__);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Return a Runge-Kutta integrator by id.
 */
/*----------------------------------------------------------------------------*/

cs_runge_kutta_integrator_t *
cs_runge_kutta_integrator_by_id(int  rk_id)
{
  if (   rk_id < 0
      || rk_id >= _n_rk_integrators)
    return nullptr;

  return _rk_lst[rk_id];
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Clean all Runge-Kutta integrators.
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_integrators_destroy()
{
  for (int i = 0; i < _n_rk_integrators; i++)
    _runge_kutta_integrator_free(&_rk_lst[i]);

  CS_FREE(_rk_lst);
  _n_rk_integrators = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief prepare and complete rhs per stage for Runge-Kutta integrator.
 *        Align with legacy equations' building sequence.
 *        1) Collect initial rhs done in
 *        cs_solve_navier_stokes/cs_solve_equation_scalar
 *        2) Complete rhs with adding explicit part of the
 *        convection/diffusion balance
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized
 *   before calling cs_balance_vector!
 * - mind the sign minus
 *
 * \param[in]      ctx           Reference to dispatch context
 * \param[in,out]  rk            pointer to a Runge-Kutta integrator
 * \param[in]      idtvar        indicator of the temporal scheme
 * \param[in]      f_id          field id (or -1)
 * \param[in]      imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]      eqp           pointer to a cs_equation_param_t structure which
 *                               contains variable calculation options
 * \param[in]      bc_coeffs     boundary condition structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                                at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                                at boundary faces for the r.h.s.
 * \param[in]      viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]      weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]      weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                                - 0 upwind scheme at all boundary faces
 *                                - 1 imposed flux at some boundary faces
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                                - 0 upwind scheme
 *                                - 1 imposed flux
 * \param[in]      pvar          solved variable (current time step)
 * \param[in]      xcpp          array of specific heat (Cp)
 */
/*----------------------------------------------------------------------------*/

void
cs_runge_kutta_stage_complete_scalar_rhs
  (cs_dispatch_context         &ctx,
   cs_runge_kutta_integrator_t *rk,
   int                          idtvar,
   int                          f_id,
   int                          imucpp,
   cs_equation_param_t         *eqp,
   cs_field_bc_coeffs_t        *bc_coeffs,
   const cs_real_t              i_massflux[],
   const cs_real_t              b_massflux[],
   const cs_real_t              i_visc[],
   const cs_real_t              b_visc[],
   cs_real_t                    viscel[][6],
   const cs_real_t              weighf[][2],
   const cs_real_t              weighb[],
   int                          icvflb,
   const int                    icvfli[],
   cs_real_t                    pvar[],
   const cs_real_t              xcpp[])
{
  assert(rk != nullptr);
  const int n_elts = rk->n_elts;
  // get the current stage index
  const int i_stg = rk->i_stage;

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t *rhs = rk->rhs_stages + i_stg*n_elts;

  cs_alloc_mode_t amode = ctx.alloc_mode(false);

  /* We compute the total explicit balance. */

  int  inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
  int  imasac = 0;

  eqp->theta = 1;

  cs_field_t *f = nullptr;

  if (f_id > -1)
    f = cs_field_by_id(f_id);

  cs_boundary_conditions_update_bc_coeff_face_values
      (ctx,
       f, bc_coeffs, inc,
       eqp,
       true, true,
       0, nullptr, // hyd_p_flag, f_ext
       nullptr, viscel, weighb,
       pvar);

  cs_balance_scalar(idtvar,
                    f_id,
                    imucpp,
                    imasac,
                    inc,
                    eqp,
                    nullptr,
                    pvar,
                    bc_coeffs,
                    i_massflux,
                    b_massflux,
                    i_visc,
                    b_visc,
                    viscel,
                    xcpp,
                    weighf,
                    weighb,
                    icvflb,
                    icvfli,
                    rhs,
                    nullptr,
                    nullptr);

  eqp->theta = 0;

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
