/*============================================================================
 * Building of the right hand side for a transport of a field.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "alge/cs_blas.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "alge/cs_convection_diffusion.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "alge/cs_gradient.h"
#include "base/cs_ext_neighborhood.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_prototypes.h"
#include "base/cs_timer.h"
#include "base/cs_timer_stats.h"
#include "base/cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_balance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_balance.cpp
 *
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion terms of a transport equation of a field.
 */
/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Timer statistics */

static int _balance_stat_id = -1;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize balance timers.
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_initialize(void)
{
  int stats_root = cs_timer_stats_id_by_name("operations");

  if (stats_root > -1) {
    _balance_stat_id = cs_timer_stats_create("operations",
                                             "balances",
                                             "balances");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion terms of a transport equation of
 * a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized
 *   before calling cs_balance_scalar!
 * - mind the minus sign
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: SOLU
 * - ischcp = 1: centered
 * - ischcp = 2: SOLU with upwind gradient reconstruction
 * - ischcp = 3: blending SOLU centered
 * - ischcp = 4: NVD-TVD
 * - imucpp = 0: do not multiply the convective part by \f$ C_p \f$
 * - imucpp = 1: multiply the convective part by \f$ C_p \f$
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     eqp           pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved variable (current time step)
 *                              may be nullptr if pvara != nullptr
 * \param[in]     pvara         solved variable (previous time step)
 *                              may be nullptr if pvar != nullptr
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_scalar(int                         idtvar,
                  int                         f_id,
                  int                         imucpp,
                  int                         imasac,
                  int                         inc,
                  cs_equation_param_t        *eqp,
                  cs_real_t                   pvar[],
                  const cs_real_t             pvara[],
                  const cs_field_bc_coeffs_t *bc_coeffs,
                  const cs_real_t             i_massflux[],
                  const cs_real_t             b_massflux[],
                  const cs_real_t             i_visc[],
                  const cs_real_t             b_visc[],
                  cs_real_6_t                 viscel[],
                  const cs_real_t             xcpp[],
                  const cs_real_2_t           weighf[],
                  const cs_real_t             weighb[],
                  int                         icvflb,
                  const int                   icvfli[],
                  cs_real_t                   smbrp[])
{
  cs_balance_scalar(idtvar,
                    f_id,
                    imucpp, imasac,
                    inc,
                    eqp,
                    pvar, pvara,
                    bc_coeffs,
                    i_massflux, b_massflux,
                    i_visc, b_visc,
                    viscel,
                    xcpp,
                    weighf, weighb,
                    icvflb, icvfli,
                    smbrp,
                    nullptr, nullptr);
}

END_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion terms of a transport equation of
 * a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized
 *   before calling cs_balance_scalar!
 * - mind the minus sign
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: SOLU
 * - ischcp = 1: centered
 * - ischcp = 2: SOLU with upwind gradient reconstruction
 * - ischcp = 3: blending SOLU centered
 * - ischcp = 4: NVD-TVD
 * - imucpp = 0: do not multiply the convective part by \f$ C_p \f$
 * - imucpp = 1: multiply the convective part by \f$ C_p \f$
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     eqp           pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved variable (current time step)
 *                              may be nullptr if pvara != nullptr
 * \param[in]     pvara         solved variable (previous time step)
 *                              may be nullptr if pvar != nullptr
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 * \param[in,out] i_flux        interior flux (or nullptr)
 * \param[in,out] b_flux        boundary flux (or nullptr)
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_scalar(int                         idtvar,
                  int                         f_id,
                  int                         imucpp,
                  int                         imasac,
                  int                         inc,
                  cs_equation_param_t        *eqp,
                  cs_real_t                   pvar[],
                  const cs_real_t             pvara[],
                  const cs_field_bc_coeffs_t *bc_coeffs,
                  const cs_real_t             i_massflux[],
                  const cs_real_t             b_massflux[],
                  const cs_real_t             i_visc[],
                  const cs_real_t             b_visc[],
                  cs_real_6_t                 viscel[],
                  const cs_real_t             xcpp[],
                  const cs_real_2_t           weighf[],
                  const cs_real_t             weighb[],
                  int                         icvflb,
                  const int                   icvfli[],
                  cs_real_t                   smbrp[],
                  cs_real_2_t                 i_flux[],
                  cs_real_t                   b_flux[])
{
  cs_timer_t t0 = cs_timer_time();

  /* Local variables */
  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int idftnp = eqp->idften;
  cs_equation_param_t eqp_loc = *eqp;

  if (f_id < 0) {
    eqp_loc.iwgrec = 0;  /* requires field id */
    eqp_loc.icoupl = -1; /* requires field id */
  }
  else {
    cs_field_t *f = cs_field_by_id(f_id);
    int k_id = cs_field_key_id("var_cal_opt");
    cs_field_get_key_struct(f, k_id, &eqp_loc);
    eqp_loc.theta = eqp->theta;
  }

  /* Scalar diffusivity */
  if (idftnp & CS_ISOTROPIC_DIFFUSION) {
    if (imucpp == 0) {
      cs_convection_diffusion_scalar(idtvar,
                                     f_id,
                                     eqp_loc,
                                     icvflb,
                                     inc,
                                     imasac,
                                     pvar,
                                     pvara,
                                     icvfli,
                                     bc_coeffs,
                                     i_massflux,
                                     b_massflux,
                                     i_visc,
                                     b_visc,
                                     smbrp,
                                     i_flux,
                                     b_flux);
    }
    else {
      /* The convective part is multiplied by Cp for the temperature */
      cs_convection_diffusion_thermal(idtvar,
                                      f_id,
                                      eqp_loc,
                                      inc,
                                      imasac,
                                      pvar,
                                      pvara,
                                      bc_coeffs,
                                      i_massflux,
                                      b_massflux,
                                      i_visc,
                                      b_visc,
                                      xcpp,
                                      smbrp);
    }
  }
  /* Symmetric tensor diffusivity */
  else if (idftnp & CS_ANISOTROPIC_DIFFUSION) {
    eqp_loc.idiff = 0;
    /* Convective part */
    if (imucpp == 0 && iconvp == 1) {
      cs_convection_diffusion_scalar(idtvar,
                                     f_id,
                                     eqp_loc,
                                     icvflb,
                                     inc,
                                     imasac,
                                     pvar,
                                     pvara,
                                     icvfli,
                                     bc_coeffs,
                                     i_massflux,
                                     b_massflux,
                                     i_visc,
                                     b_visc,
                                     smbrp,
                                     i_flux,
                                     b_flux);
    }
    /* The convective part is multiplied by Cp for the temperature */
    else if (imucpp == 1 && iconvp == 1) {
      cs_convection_diffusion_thermal(idtvar,
                                      f_id,
                                      eqp_loc,
                                      inc,
                                      imasac,
                                      pvar,
                                      pvara,
                                      bc_coeffs,
                                      i_massflux,
                                      b_massflux,
                                      i_visc,
                                      b_visc,
                                      xcpp,
                                      smbrp);
    }

    /* Diffusive part */
    if (idiffp == 1) {
      cs_anisotropic_diffusion_scalar(idtvar,
                                      f_id,
                                      eqp_loc,
                                      inc,
                                      pvar,
                                      pvara,
                                      bc_coeffs,
                                      i_visc,
                                      b_visc,
                                      viscel,
                                      weighf,
                                      weighb,
                                      smbrp);
    }
  }

  cs_timer_t t1 = cs_timer_time();

  if (_balance_stat_id > -1)
    cs_timer_stats_add_diff(_balance_stat_id, &t0, &t1);
}

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion
 * terms of a transport equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized
 *   before calling cs_balance_vector!
 * - mind the sign minus
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: SOLU
 * - ischcp = 1: centered
 * - ischcp = 2: SOLU with upwind gradient reconstruction
 * - ischcp = 3: blending SOLU centered
 * - ischcp = 4: NVD-TVD
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     eqp           pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_v  boundary conditions structure when solving
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     secvif        secondary viscosity at interior faces
 * \param[in]     secvib        secondary viscosity at boundary faces
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbr          right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_vector(int                         idtvar,
                  int                         f_id,
                  int                         imasac,
                  int                         inc,
                  int                         ivisep,
                  cs_equation_param_t        *eqp,
                  cs_real_3_t                 pvar[],
                  const cs_real_3_t           pvara[],
                  const cs_field_bc_coeffs_t *bc_coeffs_v,
                  const cs_bc_coeffs_solve_t *bc_coeffs_solve_v,
                  const cs_real_t             i_massflux[],
                  const cs_real_t             b_massflux[],
                  const cs_real_t             i_visc[],
                  const cs_real_t             b_visc[],
                  const cs_real_t             secvif[],
                  const cs_real_t             secvib[],
                  cs_real_6_t                 viscel[],
                  const cs_real_2_t           weighf[],
                  const cs_real_t             weighb[],
                  int                         icvflb,
                  const int                   icvfli[],
                  cs_real_3_t                 i_pvar[],
                  cs_real_3_t                 b_pvar[],
                  cs_real_3_t                 smbr[])
{
  cs_timer_t t0 = cs_timer_time();

  /* Local variables */
  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int idftnp = eqp->idften;

  assert(eqp != nullptr);
  cs_equation_param_t eqp_loc = *eqp;

  if (f_id < 0) {
    eqp_loc.iwgrec = 0;  /* requires field id */
    eqp_loc.icoupl = -1; /* requires field id */
  }
  else {
    cs_field_t *f = cs_field_by_id(f_id);
    int k_id = cs_field_key_id("var_cal_opt");
    cs_field_get_key_struct(f, k_id, &eqp_loc);
    eqp_loc.theta = eqp->theta;
  }

  /* Scalar diffusivity */
  if (idftnp & CS_ISOTROPIC_DIFFUSION) {
    cs_convection_diffusion_vector(idtvar,
                                   f_id,
                                   eqp_loc,
                                   icvflb,
                                   inc,
                                   ivisep,
                                   imasac,
                                   pvar,
                                   pvara,
                                   icvfli,
                                   bc_coeffs_v,
                                   bc_coeffs_solve_v,
                                   i_massflux,
                                   b_massflux,
                                   i_visc,
                                   b_visc,
                                   secvif,
                                   secvib,
                                   i_pvar,
                                   b_pvar,
                                   smbr);
  }
  /* Symmetric tensor diffusivity */
  else if (idftnp & CS_ANISOTROPIC_DIFFUSION) {
    /* Neither diffusive part nor secondary viscosity or transpose of gradient */
    eqp_loc.idiff = 0;
    /* Convective part */
    if (iconvp == 1) {
      cs_convection_diffusion_vector(idtvar,
                                     f_id,
                                     eqp_loc,
                                     icvflb,
                                     inc,
                                     ivisep,
                                     imasac,
                                     pvar,
                                     pvara,
                                     icvfli,
                                     bc_coeffs_v,
                                     bc_coeffs_solve_v,
                                     i_massflux,
                                     b_massflux,
                                     i_visc,
                                     b_visc,
                                     secvif,
                                     secvib,
                                     i_pvar,
                                     b_pvar,
                                     smbr);

    }

    /* Diffusive part (with a 3x3 symmetric diffusivity) */

    /* either Daly-Harlow type i.e. Nabla(v).K */
    if (idiffp == 1 && idftnp & CS_ANISOTROPIC_RIGHT_DIFFUSION) {
      /* ! Neither diffusive part neither secondary viscosity
         nor transpose of gradient */
      cs_anisotropic_right_diffusion_vector(idtvar,
                                            f_id,
                                            eqp_loc,
                                            inc,
                                            pvar,
                                            pvara,
                                            bc_coeffs_v,
                                            bc_coeffs_solve_v,
                                            i_visc,
                                            b_visc,
                                            viscel,
                                            weighf,
                                            weighb,
                                            smbr);
    }

    /* or K.Nabla(v) ) */
    else if (idiffp == 1 && idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
      cs_anisotropic_left_diffusion_vector(idtvar,
                                           f_id,
                                           eqp_loc,
                                           inc,
                                           ivisep,
                                           pvar,
                                           pvara,
                                           bc_coeffs_v,
                                           bc_coeffs_solve_v,
                                           (const cs_real_33_t *)i_visc,
                                           b_visc,
                                           secvif,
                                           smbr);
    }
  }

  cs_timer_t t1 = cs_timer_time();

  if (_balance_stat_id > -1)
    cs_timer_stats_add_diff(_balance_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper to the function which adds the explicit part of the
 * convection/diffusion
 * terms of a transport equation of a tensor field \f$ \tens{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
 *      - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \tens{Rhs} \f$ has already been initialized before calling bilscts!
 * - mind the sign minus
 *
 * Options for the convective scheme:
 * - blencp = 0: upwind scheme for the advection
 * - blencp = 1: no upwind scheme except in the slope test
 * - ischcp = 0: SOLU
 * - ischcp = 1: centered
 * - ischcp = 2: SOLU with upwind gradient reconstruction
 * - ischcp = 3: blending SOLU centered
 * - ischcp = 4: NVD-TVD
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     inc           indicator
 * \param[in]     eqp           pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_ts   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_tensor(int                         idtvar,
                  int                         f_id,
                  int                         imasac,
                  int                         inc,
                  cs_equation_param_t        *eqp,
                  cs_real_6_t                 pvar[],
                  const cs_real_6_t           pvara[],
                  const cs_field_bc_coeffs_t *bc_coeffs_ts,
                  const cs_bc_coeffs_solve_t *bc_coeffs_solve_ts,
                  const cs_real_t             i_massflux[],
                  const cs_real_t             b_massflux[],
                  const cs_real_t             i_visc[],
                  const cs_real_t             b_visc[],
                  cs_real_6_t                 viscel[],
                  const cs_real_2_t           weighf[],
                  const cs_real_t             weighb[],
                  int                         icvflb,
                  const int                   icvfli[],
                  cs_real_6_t                 smbrp[])
{
  CS_UNUSED(icvfli);

  cs_timer_t t0 = cs_timer_time();

  /* Local variables */
  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int idftnp = eqp->idften;

  assert(eqp != nullptr);
  cs_equation_param_t eqp_loc = *eqp;

  if (f_id < 0) {
    eqp_loc.iwgrec = 0;  /* requires field id */
    eqp_loc.icoupl = -1; /* requires field id */
  }
  else {
    cs_field_t *f = cs_field_by_id(f_id);
    int k_id = cs_field_key_id("var_cal_opt");
    cs_field_get_key_struct(f, k_id, &eqp_loc);
    eqp_loc.theta = eqp->theta;
  }

  /* Scalar diffusivity */
  if (idftnp & CS_ISOTROPIC_DIFFUSION) {
    cs_convection_diffusion_tensor(idtvar,
                                   f_id,
                                   eqp_loc,
                                   icvflb,
                                   inc,
                                   imasac,
                                   pvar,
                                   pvara,
                                   bc_coeffs_ts,
                                   bc_coeffs_solve_ts,
                                   i_massflux,
                                   b_massflux,
                                   i_visc,
                                   b_visc,
                                   smbrp);
  }
  /* Symmetric tensor diffusivity */
  else if (idftnp & CS_ANISOTROPIC_RIGHT_DIFFUSION) {
    /* No diffusive part */
    eqp_loc.idiff = 0;
    /* Convective part */
    if (iconvp == 1) {
      cs_convection_diffusion_tensor(idtvar,
                                     f_id,
                                     eqp_loc,
                                     icvflb,
                                     inc,
                                     imasac,
                                     pvar,
                                     pvara,
                                     bc_coeffs_ts,
                                     bc_coeffs_solve_ts,
                                     i_massflux,
                                     b_massflux,
                                     i_visc,
                                     b_visc,
                                     smbrp);
    }
    /* Diffusive part (with a 6x6 symmetric diffusivity) */
    if (idiffp == 1) {
      cs_anisotropic_diffusion_tensor(idtvar,
                                      f_id,
                                      eqp_loc,
                                      inc,
                                      pvar,
                                      pvara,
                                      bc_coeffs_ts,
                                      bc_coeffs_solve_ts,
                                      i_visc,
                                      b_visc,
                                      viscel,
                                      weighf,
                                      weighb,
                                      smbrp);
    }
  }

  cs_timer_t t1 = cs_timer_time();

  if (_balance_stat_id > -1)
    cs_timer_stats_add_diff(_balance_stat_id, &t0, &t1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
