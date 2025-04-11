/*============================================================================
 * Convection diffusion equation solver.
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_operator.h"
#include "base/cs_fp_exception.h"
#include "base/cs_halo.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_reducers.h"
#include "base/cs_timer.h"
#include "base/cs_parall.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "alge/cs_balance.h"
#include "alge/cs_blas.h"
#include "alge/cs_convection_diffusion.h"
#include "alge/cs_gradient.h"
#include "alge/cs_gradient_boundary.h"
#include "alge/cs_multigrid.h"
#include "base/cs_reducers.h"
#include "alge/cs_matrix_building.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_equation_iterative_solve.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

#ifdef __cplusplus

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

static void
_init_bc_coeffs_solve(cs_bc_coeffs_solve_t  &c,
                      cs_lnum_t              n_b_faces,
                      cs_lnum_t              stride,
                      cs_alloc_mode_t        amode,
                      bool                   limiter)
{
  c.val_ip = nullptr;
  c.val_f = nullptr;
  c.val_f_lim = nullptr;
  c.val_f_d = nullptr;
  c.val_f_d_lim = nullptr;

  CS_MALLOC_HD(c.val_ip, stride*n_b_faces, cs_real_t, amode);
  CS_MALLOC_HD(c.val_f, stride*n_b_faces, cs_real_t, amode);
  CS_MALLOC_HD(c.val_f_d, stride*n_b_faces, cs_real_t, amode);

  if (limiter == false) {
    c.val_f_lim = c.val_f;
    c.val_f_d_lim = c.val_f_d;
  }
  else {
    CS_MALLOC_HD(c.val_f_lim, stride*n_b_faces, cs_real_t, amode);
    CS_MALLOC_HD(c.val_f_d_lim, stride*n_b_faces, cs_real_t, amode);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free boundary condition coefficients solve arrays.
 *
 * \param[in, out]  c          reference to structure to initialize.
 */
/*----------------------------------------------------------------------------*/

static void
_clear_bc_coeffs_solve(cs_bc_coeffs_solve_t  &c)
{
  if (c.val_f_lim != c.val_f) {
    CS_FREE_HD(c.val_f_lim);
    CS_FREE_HD(c.val_f_d_lim);
  }
  else {
    c.val_f_lim = nullptr;
    c.val_f_d_lim = nullptr;
  }

  CS_FREE_HD(c.val_ip);
  CS_FREE_HD(c.val_f);
  CS_FREE_HD(c.val_f_d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function solves an advection diffusion equation with source terms
 * for one time step for the vector variable \f$ \vect{a} \f$.
 *
 * The equation reads:
 *
 * \f[
 * \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
 * + \divv \left( \vect{a}^{n+1} \otimes \rho \vect {u}
 *              - \mu \gradt \vect{a}^{n+1}\right)
 * = \vect{Rhs}
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}
 * + \divv \left( \delta \vect{a} \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a} \right)
 * = \vect{Rhs}^1
 * \f]
 *
 * where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
 * \f$ \vect{Rhs}^1 = \vect{Rhs}
 * - \divv \left( \vect{a}^n \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^n \right)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}^k
 * + \divv \left( \delta \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a}^k \right)
 * = \vect{Rhs}^k
 * \f]
 *
 * where \f$ \vect{Rhs}^k = \vect{Rhs}
 * - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
 * - \divv \left( \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
 *
 * \param[in]      idtvar        indicator of the temporal scheme
 * \param[in]      iterns        external sub-iteration number
 * \param[in]      f_id          field id (or -1)
 * \param[in]      name          associated name if f_id < 0, or nullptr
 * \param[in]      ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]      iescap        compute the predictor indicator if >= 1
 * \param[in]      eqp           pointer to a cs_equation_param_t structure which
 *                               contains variable calculation options
 * \param[in]      pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]      pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               \c pvara (usually \c pvar= \c pvara)
 * \param[in]      bc_coeffs_v   boundary condition structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]      b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]      i_secvis      secondary viscosity at interior faces
 * \param[in]      b_secvis      secondary viscosity at boundary faces
 * \param[in]      viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]      weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]      weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in, out] fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in, out] smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in, out] pvar          current variable
 * \param[out]     eswork        prediction-stage error estimator
 *                               (if iescap >= 0)
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_equation_iterative_solve_strided(int                   idtvar,
                                  int                   iterns,
                                  int                   f_id,
                                  const char           *name,
                                  int                   ivisep,
                                  int                   iescap,
                                  cs_equation_param_t  *eqp,
                                  const cs_real_t       pvara[][stride],
                                  const cs_real_t       pvark[][stride],
                                  cs_field_bc_coeffs_t *bc_coeffs,
                                  const cs_real_t       i_massflux[],
                                  const cs_real_t       b_massflux[],
                                  const cs_real_t       i_viscm[],
                                  const cs_real_t       b_viscm[],
                                  const cs_real_t       i_visc[],
                                  const cs_real_t       b_visc[],
                                  const cs_real_t       i_secvis[],
                                  const cs_real_t       b_secvis[],
                                  cs_real_t             viscel[][6],
                                  const cs_real_t       weighf[][2],
                                  const cs_real_t       weighb[],
                                  int                   icvflb,
                                  const int             icvfli[],
                                  cs_real_t             fimp[][stride][stride],
                                  cs_real_t             smbrp[][stride],
                                  cs_real_t             pvar[][stride],
                                  cs_real_t             eswork[][stride])
{
  /* Local variables */
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_halo_t *halo = cs_glob_mesh->halo;

  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int iwarnp = eqp->verbosity;
  int iswdyp = eqp->iswdyn;
  int idftnp = eqp->idften;
  int ndircp = eqp->ndircl;
  double epsrsp = eqp->epsrsm;
  double epsilp = eqp->epsilo;
  double thetap = eqp->theta;

  const cs_real_t  *cell_vol = mq->cell_vol;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  using var_t = cs_real_t[stride];
  using m_t = cs_real_t[stride][stride];

  int inc, niterf;
  int lvar, imasac;
  double residu, ressol;
  double nadxk, paxm1rk, paxm1ax;

  cs_real_t alph = 0., beta = 0.;

  cs_solving_info_t sinfo;

  cs_field_t *f = nullptr;

  /*============================================================================
   * Initialization
   *==========================================================================*/

  /* Name */
  const char *var_name = cs_sles_name(f_id, name);

  if (iwarnp >= 1)
    bft_printf("Equation iterative solve of: %s\n", var_name);

  /* Matrix block size */
  cs_lnum_t eb_size = 1; /* CS_ISOTROPIC_DIFFUSION
                            or CS_ANISOTROPIC_RIGHT_DIFFUSION */
  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) eb_size = stride;

  if (cs_glob_porous_model == 3 && stride == 3) { //FIXME iphydr + other option?
    if (iconvp > 0)
      eb_size = 3;
  }

  /* Parallel or device dispatch */

  cs_dispatch_context ctx, ctx_c;
  if (idtvar < 0)
    ctx.set_use_gpu(false);  /* steady case not ported to GPU */

#if defined(HAVE_CUDA)
  ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  /* Allocate temporary arrays */
  m_t *dam;
  var_t *smbini, *dpvar;
  CS_MALLOC_HD(dam, n_cells_ext, m_t, amode);
  CS_MALLOC_HD(smbini, n_cells_ext, var_t, amode);
  CS_MALLOC_HD(dpvar, n_cells_ext, var_t, amode);

  var_t *adxk = nullptr, *adxkm1 = nullptr, *dpvarm1 = nullptr, *rhs0 = nullptr;
  if (iswdyp >= 1) {
   CS_MALLOC_HD(adxk, n_cells_ext, var_t, amode);
   CS_MALLOC_HD(adxkm1, n_cells_ext, var_t, amode);
   CS_MALLOC_HD(dpvarm1, n_cells_ext, var_t, amode);
   CS_MALLOC_HD(rhs0, n_cells_ext, var_t, amode);
  }

  var_t *i_pvar = nullptr;
  var_t *b_pvar = nullptr;
  cs_field_t *i_vf = nullptr;
  cs_field_t *b_vf = nullptr;

  /* Storing face values for kinetic energy balance and initialize them */
  if (CS_F_(vel) != nullptr && CS_F_(vel)->id == f_id) {

    i_vf = cs_field_by_name_try("inner_face_velocity");
    if (i_vf != nullptr) {
      ctx.parallel_for(3*n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        i_vf->val[face_id] = 0.;
      });
      i_pvar = (var_t *)i_vf->val;
    }

    b_vf = cs_field_by_name_try("boundary_face_velocity");
    if (b_vf != nullptr) {
      ctx.parallel_for(3*n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        b_vf->val[face_id] = 0.;
      });
      b_pvar = (var_t *)b_vf->val;
    }

    ctx.wait();
  }

  /* solving info */
  int df_limiter_id = -1;
  int key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);

    df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
  }

  /* Symmetric matrix, except if advection */
  int isym = 1;
  if (iconvp > 0) isym = 2;

  bool symmetric = (isym == 1) ? true : false;

  /* Determine if we are in a case with special requirements */

  bool conv_diff_mg = false;
  if (iconvp > 0) {
    cs_sles_t *sc = cs_sles_find_or_add(f_id, name);
    const char *sles_type = cs_sles_get_type(sc);
    if (strcmp(sles_type, "cs_multigrid_t") == 0)
      conv_diff_mg = true;
  }

  /*==========================================================================
   * Building of the "simplified" matrix
   *==========================================================================*/

  cs_matrix_t *a = cs_sles_default_get_matrix
                     (f_id, var_name, stride, eb_size, symmetric);

  int tensorial_diffusion = 1;

  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION)
    tensorial_diffusion = 2;

  /* For steady computations, the diagonal is relaxed */
  cs_real_t relaxp = (idtvar < 0) ? eqp->relaxv : 1.;

  cs_matrix_compute_coeffs(a,
                           f,
                           iconvp,
                           idiffp,
                           tensorial_diffusion,
                           ndircp,
                           eb_size,
                           thetap,
                           relaxp,
                           bc_coeffs,
                           fimp,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm);

  /*===========================================================================
   * Iterative process to handle non orthogonalities (starting from the
   * second iteration).
   *===========================================================================*/

  /* Allocate non reconstructed face value only if presence of limiter */

  const int ircflp = eqp->ircflu;
  const int ircflb = (ircflp > 0) ? eqp->b_diff_flux_rc : 0;

  cs_bc_coeffs_solve_t bc_coeffs_solve;
  _init_bc_coeffs_solve(bc_coeffs_solve,
                        n_b_faces,
                        stride,
                        amode,
                        (df_limiter_id > -1 || ircflb != 1));

  var_t *val_ip = (var_t *)bc_coeffs_solve.val_ip;
  var_t *val_f = (var_t *)bc_coeffs_solve.val_f;
  var_t *val_f_lim = (var_t *)bc_coeffs_solve.val_f_lim;
  var_t *val_f_d =  (var_t *)bc_coeffs_solve.val_f_d;
  var_t *val_f_d_lim =  (var_t *)bc_coeffs_solve.val_f_d_lim;

  /* Application of the theta-scheme */

  /* We compute the total explicit balance. */
  cs_real_t thetex = 1. - thetap;

  /* If thetex=0, no need to add anything */
  if (cs::abs(thetex) > cs_math_epzero) {
    inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
    imasac = 0;

    eqp->theta = thetex;

    cs_boundary_conditions_update_bc_coeff_face_values<stride>
      (ctx, f, bc_coeffs, inc, eqp, pvara,
       val_ip, val_f, val_f_lim, val_f_d, val_f_d_lim);

    if (stride == 3)
      cs_balance_vector(idtvar,
                        f_id,
                        imasac,
                        1,       /* inc */
                        ivisep,
                        eqp,
                        nullptr, /* pvar == pvara */
                        (const cs_real_3_t *)pvara,
                        bc_coeffs,
                        &bc_coeffs_solve,
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        i_secvis,
                        b_secvis,
                        viscel,
                        weighf,
                        weighb,
                        icvflb,
                        icvfli,
                        (cs_real_3_t *)i_pvar,
                        (cs_real_3_t *)b_pvar,
                        (cs_real_3_t *)smbrp);
    else if (stride == 6)
      cs_balance_tensor(idtvar,
                        f_id,
                        imasac,
                        1,       /* inc */
                        eqp,
                        nullptr, /* pvar == pvara */
                        (const cs_real_6_t *)pvara,
                        bc_coeffs,
                        &bc_coeffs_solve,
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        viscel,
                        weighf,
                        weighb,
                        icvflb,
                        icvfli,
                        (cs_real_6_t *)smbrp);

    /* Save (1-theta)* face_value at previous time step if needed */
    if (i_vf != nullptr && b_vf != nullptr) {
      cs_field_current_to_previous(i_vf);
      cs_field_current_to_previous(b_vf);
    }

    eqp->theta = thetap;
  }

  /* Before looping, the RHS without reconstruction is stored in smbini */

  cs_lnum_t has_dc = mq->has_disable_flag;
  ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      smbini[c_id][i] = smbrp[c_id][i];
      smbrp[c_id][i] = 0.;
    }
  });

  /* pvar is initialized on n_cells_ext to avoid a synchronization */
  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++)
      pvar[c_id][i] = pvark[c_id][i];
  });

  /* Synchronize before next major operation */
  ctx.wait();
  ctx_c.wait();

  /* In the following, cs_balance is called with inc=1,
   * except for Weight Matrix (nswrsp=-1) */
  inc = 1;

  if (eqp->nswrsm == -1) {
    eqp->nswrsm = 1;
    inc = 0;
  }

  cs_boundary_conditions_update_bc_coeff_face_values<stride>
    (ctx, f, bc_coeffs, inc, eqp, pvar,
     val_ip, val_f, val_f_lim, val_f_d, val_f_d_lim);

  /*  Incrementation and rebuild of right hand side */

  /*  We enter with an explicit SMB based on PVARA.
   *  If we initialize with PVAR with something other than PVARA
   *  we need to correct SMBR (this happens when iterating over navsto) */

  /* The added convective scalar mass flux is:
   *      (thetap*Y_\face-imasac*Y_\celli)*mf.
   * When building the implicit part of the rhs, one
   * has to impose 1 on mass accumulation. */
  imasac = 1;

  if (stride == 3)
    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,
                      ivisep,
                      eqp,
                      (cs_real_3_t *)pvar,
                      (const cs_real_3_t *)pvara,
                      bc_coeffs,
                      &bc_coeffs_solve,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      i_secvis,
                      b_secvis,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      nullptr,
                      nullptr,
                      (cs_real_3_t *)smbrp);
  else if (stride == 6)
    cs_balance_tensor(idtvar,
                      f_id,
                      imasac,
                      inc,
                      eqp,
                      (cs_real_6_t *)pvar,
                      (const cs_real_6_t *)pvara,
                      bc_coeffs,
                      &bc_coeffs_solve,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      (cs_real_6_t *)smbrp);

  if (CS_F_(vel) != nullptr && CS_F_(vel)->id == f_id) {
    cs_field_t *f_ex = cs_field_by_name_try("velocity_explicit_balance");

    if (f_ex != nullptr) {
      cs_real_3_t *cpro_cv_df_v = (cs_real_3_t *)f_ex->val;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++)
          cpro_cv_df_v[c_id][i] = smbrp[c_id][i];
      });
    }
  }

  /* Dynamic relaxation */
  if (iswdyp >= 1) {
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      cs_real_t c_sum = 0;
      for (cs_lnum_t i = 0; i < stride; i++) {
        rhs0[c_id][i] = smbrp[c_id][i];

        cs_real_t diff = 0.;
        for (cs_lnum_t j = 0; j < stride; j++)
          diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

        smbrp[c_id][i] += smbini[c_id][i] - diff;
        c_sum += cs_math_pow2(smbrp[c_id][i]);

        adxkm1[c_id][i] = 0.;
        adxk[c_id][i] = 0.;
        dpvar[c_id][i] = 0.;
      }
      sum += c_sum;
    });

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      cs_real_t c_sum = 0;
      for (cs_lnum_t i = 0; i < stride; i++) {

        cs_real_t diff = 0.;
        for (cs_lnum_t j = 0; j < stride; j++)
          diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

        smbrp[c_id][i] += smbini[c_id][i] - diff;
        c_sum += cs_math_pow2(smbrp[c_id][i]);
      }
      sum += c_sum;
    });
  }

  ctx.wait();
  cs_parall_sum(1, CS_DOUBLE, &residu);

  /* --- Right hand side residual */
  residu = sqrt(residu);

  /* Normalization residual
     (L2-norm of B.C. + source terms + non-orthogonality terms)

     Caution: when calling a matrix-vector product, here for a variable
     which is not "by increments" and is assumed initialized, including
     for ghost values. */

  /* Allocate a temporary array */

  // Number of local ghost cells may be different from that of mesh
  // in case of internal coupling.
  cs_lnum_t n_cols = cs_matrix_get_n_columns(a);

  var_t *w1, *w2;
  CS_MALLOC_HD(w1, n_cols, var_t, amode);
  CS_MALLOC_HD(w2, n_cols, var_t, amode);

  /* Compute the L2 norm of the variable */

  struct cs_double_n<stride> rd;
  struct cs_reduce_sum_n<stride> reducer;

  ctx.parallel_for_reduce(n_cells, rd, reducer, [=] CS_F_HOST_DEVICE
                          (cs_lnum_t c_id, cs_double_n<stride> &res) {
    cs_real_t c_vol = cell_vol[c_id];
    for (size_t i = 0; i < stride; i++)
      res.r[i] = pvar[c_id][i] * c_vol;
  });
  ctx.wait();

  cs_parall_sum_strided<stride>(rd.r);

  cs_real_t p_mean[stride];
  for (size_t i = 0; i < stride; i++) {
    p_mean[i] = rd.r[i] / mq->tot_vol;
    if (iwarnp >= 2)
      bft_printf("Spatial average of X_%d^n = %f\n", (int)i, p_mean[i]);
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < stride; i++) {
      w2[c_id][i] = (pvar[c_id][i] - p_mean[i]);
    }
  });

  cs_matrix_vector_multiply(a,
                            (cs_real_t *)w2,
                            (cs_real_t *)w1);

  // ctx.wait(); // matrix vector multiply uses the same stream as the ctx

  if (iwarnp >= 2) {
    struct cs_data_2r rd2;
    struct cs_reduce_sum2r reducer2;

    ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                            (cs_lnum_t c_id, cs_data_2r &sum) {
      sum.r[0] = 0.;
      sum.r[1] = 0.;
      for (cs_lnum_t i = 0; i < stride; i++) {
        sum.r[0] += w1[c_id][i] * w1[c_id][i];
        sum.r[1] += smbrp[c_id][i] * smbrp[c_id][i];
      }
    });
    ctx.wait();
    cs_parall_sum(2, CS_DOUBLE, rd2.r);

    bft_printf("L2 norm ||AX^n|| = %f\n", sqrt(rd2.r[0]));
    bft_printf("L2 norm ||B^n|| = %f\n",  sqrt(rd2.r[1]));
  }

  double rnorm2;

  if (has_dc == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for_reduce_sum(n_cells, rnorm2, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      cs_real_t c_sum = 0;
      /* Remove contributions from penalized cells */
      if (c_disable_flag[c_id] != 0) {
        for (cs_lnum_t i = 0; i < stride; i++)
          w1[c_id][i] = 0.;
      }
      else {
        for (cs_lnum_t i = 0; i < stride; i++) {
          w1[c_id][i] += smbrp[c_id][i];
          c_sum += cs_math_pow2(w1[c_id][i]);
        }
      }
      sum += c_sum;
    });
  }
  else {
    ctx.parallel_for_reduce_sum(n_cells, rnorm2, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      cs_real_t c_sum = 0;
      for (cs_lnum_t i = 0; i < stride; i++) {
        w1[c_id][i] += smbrp[c_id][i];
        c_sum += cs_math_pow2(w1[c_id][i]);
      }
      sum += c_sum;
    });
  }

  ctx.wait();
  cs_parall_sum(1, CS_DOUBLE, &rnorm2);
  cs_real_t rnorm = sqrt(rnorm2);

  CS_FREE_HD(w1);
  CS_FREE_HD(w2);

  sinfo.rhs_norm = rnorm;

  /* Warning: for Weight Matrix, one and only one sweep is done. */
  int nswmod = cs::max(eqp->nswrsm, 1);

  cs_sles_t *sc = cs_sles_find_or_add(f_id, var_name);

  int isweep = 1;

  /* Reconstruction loop (beginning)
   *-------------------------------- */
  if ((iterns <= 1 && stride <= 3) || stride == 6)
    sinfo.n_it = 0;

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {
    /* --- Solving on the increment dpvar */

    /*  Dynamic relaxation of the system */
    if (iswdyp >= 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          dpvarm1[c_id][i] = dpvar[c_id][i];
        }
      });
    }

    /*  Solver residual */
    ressol = residu;

    if (conv_diff_mg) {
      cs_multigrid_t *mg
        = static_cast<cs_multigrid_t *>(cs_sles_get_context(sc));
      cs_multigrid_setup_conv_diff(mg, var_name, a, true,
                                   cs_sles_get_verbosity(sc));
    }

    cs_sles_solve_ccc_fv(sc,
                         a,
                         epsilp,
                         rnorm,
                         &niterf,
                         &ressol,
                         (cs_real_t *)smbrp,
                         (cs_real_t *)dpvar);

    /* Dynamic relaxation of the system */
    if (iswdyp >= 1) {

      /* Computation of the variable relaxation coefficient */
      lvar = -1;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          adxkm1[c_id][i] = adxk[c_id][i];
          adxk[c_id][i] = - rhs0[c_id][i];
        }
      });

      ctx.wait();

      cs_halo_sync_r(halo, ctx.use_gpu(), dpvar);

      /* update with dpvar */
      cs_boundary_conditions_update_bc_coeff_face_values<stride>
        (ctx, nullptr, bc_coeffs, inc, eqp, dpvar,
         val_ip, val_f, val_f_lim, val_f_d, val_f_d_lim);

      if (stride == 3)
        cs_balance_vector(idtvar,
                          lvar,
                          imasac,
                          inc,
                          ivisep,
                          eqp,
                          (cs_real_3_t *)dpvar,
                          nullptr, /* dpvar */
                          bc_coeffs,
                          &bc_coeffs_solve,
                          i_massflux,
                          b_massflux,
                          i_visc,
                          b_visc,
                          i_secvis,
                          b_secvis,
                          viscel,
                          weighf,
                          weighb,
                          icvflb,
                          icvfli,
                          nullptr,
                          nullptr,
                          (cs_real_3_t *)adxk);
      else if (stride == 6)
        cs_balance_tensor(idtvar,
                          lvar,
                          imasac,
                          inc,
                          eqp,
                          (cs_real_6_t *)dpvar,
                          nullptr, /* dpvar */
                          bc_coeffs,
                          &bc_coeffs_solve,
                          i_massflux,
                          b_massflux,
                          i_visc,
                          b_visc,
                          viscel,
                          weighf,
                          weighb,
                          icvflb,
                          icvfli,
                          (cs_real_6_t *)adxk);

      /* ||E.dx^(k-1)-E.0||^2 */
      cs_real_t nadxkm1 = nadxk;

      struct cs_data_2r rd2;
      struct cs_reduce_sum2r reducer2;

      ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                              (cs_lnum_t c_id, cs_data_2r &sum) {
        sum.r[0] = 0.;
        sum.r[1] = 0.;
        for (cs_lnum_t i = 0; i < stride; i++) {
          sum.r[0] += adxk[c_id][i] * adxk[c_id][i];
          sum.r[1] += smbrp[c_id][i] * adxk[c_id][i];
        }
      });

      ctx.wait();
      cs_parall_sum(2, CS_DOUBLE, rd2.r);

      nadxk = rd2.r[0];              /* ||E.dx^k-E.0||^2 */
      cs_real_t paxkrk = rd2.r[1];   /* < E.dx^k-E.0; r^k > */

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {

        ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id, cs_data_2r &sum) {
          sum.r[0] = 0.;
          sum.r[1] = 0.;
          for (cs_lnum_t i = 0; i < stride; i++) {
            sum.r[0] += smbrp[c_id][i] * adxkm1[c_id][i];
            sum.r[1] += adxk[c_id][i] * adxkm1[c_id][i];
          }
        });

        ctx.wait();
        cs_parall_sum(2, CS_DOUBLE, rd2.r);

        paxm1rk = rd2.r[0];   // < E.dx^(k-1)-E.0; r^k >
        paxm1ax = rd2.r[1];   // < E.dx^(k-1)-E.0; E.dx^k -E.0 >

        if (   (nadxkm1 > 1.e-30*rnorm2)
            && (nadxk*nadxkm1 - cs_math_pow2(paxm1ax)) > 1.e-30*rnorm2)
          beta =   (paxkrk*paxm1ax - nadxk*paxm1rk)
                 / (nadxk*nadxkm1 - cs_math_pow2(paxm1ax));
        else
          beta = 0.;

      }
      else {
        beta = 0.;
        paxm1ax = 1.;
        paxm1rk = 0.;
        paxm1ax = 0.;
      }

      /* The first sweep is not relaxed */
      if (isweep == 1) {
        alph = 1.;
        beta = 0.;
      }
      else if (isweep == 2) {
        beta = 0.;
        alph = -paxkrk/cs::max(nadxk, 1.e-30*rnorm2);
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/cs::max(nadxk, 1.e-30*rnorm2);
      }

      /* Writing */
      if (iwarnp >= 3)
        bft_printf("%s Sweep: %d Dynamic relaxation: alpha = %12.5e, "
                   "beta = %12.5e,\n< dI^k :  R^k >   = %12.5e, "
                   "||dI^k  ||^2 = %12.5e,\n< dI^k-1 : R^k >  = %12.5e, "
                   "||dI^k-1||^2 = %12.5e,\n< dI^k-1 : dI^k > = %12.5e\n",
                   var_name, isweep, alph, beta, paxkrk, nadxk, paxm1rk,
                   nadxkm1, paxm1ax);
    }

    /* Update the solution with the increment, update the face value,
       update the right hand side and compute the new residual */

    if (iswdyp <= 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++)
          pvar[c_id][i] += dpvar[c_id][i];

          /* smbini does not contain unsteady terms and mass source terms
           * of the RHS updated at each sweep */

        for (cs_lnum_t i = 0; i < stride; i++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; j++)
            diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

          smbrp[c_id][i] = smbini[c_id][i] - diff;
        }

      });
    }
    else if (iswdyp == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++)
          pvar[c_id][i] += alph*dpvar[c_id][i];

          /* smbini does not contain unsteady terms and mass source terms
           * of the RHS updated at each sweep */

        for (cs_lnum_t i = 0; i < stride; i++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; i++)
            diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

          smbrp[c_id][i] = smbini[c_id][i] - diff;
        }
      });
    }
    else if (iswdyp >= 2) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t i = 0; i < stride; i++) {
          pvar[c_id][i] +=   alph*dpvar[c_id][i]
                              + beta*dpvarm1[c_id][i];
        }

        for (cs_lnum_t i = 0; i < stride; i++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; j++)
            diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

          smbrp[c_id][i] = smbini[c_id][i] - diff;
        }
      });
    }

    /* --> Handle parallelism and periodicity */

    ctx.wait();

    cs_halo_sync_r(halo, ctx.use_gpu(), pvar);

    /* Increment face value with theta * face_value at current time step
     * if needed
     * Reinit the previous value before */
    if (i_vf != nullptr && b_vf != nullptr) {
      cs_array_real_copy(3 * n_i_faces, i_vf->val_pre, i_vf->val);
      cs_array_real_copy(3 * n_b_faces, b_vf->val_pre, b_vf->val);
    }

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the implicit part of the rhs, one
     * has to impose 1 on mass accumulation. */
    imasac = 1;

    /* Update face value for gradient and convection-diffusion */
    cs_boundary_conditions_update_bc_coeff_face_values<stride>
      (ctx, f, bc_coeffs, inc, eqp, pvar,
       val_ip, val_f, val_f_lim, val_f_d, val_f_d_lim);

    if (stride == 3)
      cs_balance_vector(idtvar,
                        f_id,
                        imasac,
                        inc,
                        ivisep,
                        eqp,
                        (cs_real_3_t *)pvar,
                        (const cs_real_3_t *)pvara,
                        bc_coeffs,
                        &bc_coeffs_solve,
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        i_secvis,
                        b_secvis,
                        viscel,
                        weighf,
                        weighb,
                        icvflb,
                        icvfli,
                        (cs_real_3_t *)i_pvar,
                        (cs_real_3_t *)b_pvar,
                        (cs_real_3_t *)smbrp);
    else if (stride == 6)
      cs_balance_tensor(idtvar,
                        f_id,
                        imasac,
                        inc,
                        eqp,
                        (cs_real_6_t *)pvar,
                        (const cs_real_6_t *)pvara,
                        bc_coeffs,
                        &bc_coeffs_solve,
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        viscel,
                        weighf,
                        weighb,
                        icvflb,
                        icvfli,
                        (cs_real_6_t *)smbrp);

    /* --- Convergence test */
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      cs_real_t c_sum = 0;
      for (cs_lnum_t i = 0; i < stride; i++) {
        c_sum += cs_math_pow2(smbrp[c_id][i]);
      }
      sum += c_sum;
    });
    ctx.wait();
    cs_parall_sum(1, CS_DOUBLE, &residu);
    residu = sqrt(residu);

    /* Writing */
    sinfo.n_it = sinfo.n_it + niterf;

    /* Writing */
    if (iwarnp >= 2) {
      bft_printf("%s: CV_DIF_TS, IT: %d, Res: %12.5e, Norm: %12.5e\n",
                 var_name, isweep, residu, rnorm);
      bft_printf("%s: Current reconstruction sweep: %d, "
                 "Iterations for solver: %d\n", var_name, isweep, niterf);
    }

    isweep++;
  }

  /* --- Reconstruction loop (end) */

  /* Writing: convergence */
  if (cs::abs(rnorm)/sqrt(cs_real_t(stride)) > cs_math_epzero)
    sinfo.res_norm = residu/rnorm;
  else
    sinfo.res_norm = 0.;

  if (iwarnp >= 1) {
    if (residu <= epsrsp*rnorm)
      bft_printf("%s : CV_DIF_TS, IT : %d, Res : %12.5e, Norm : %12.5e\n",
                 var_name, isweep-1, residu, rnorm);
    /* Writing: non-convergence */
    else if (isweep > nswmod)
      bft_printf("@\n@ @@ WARNING: %s CONVECTION-DIFFUSION-SOURCE-TERMS\n@"
                 "=========\n@  Maximum number of iterations %d reached\n",
                 var_name,nswmod);
  }

  /* Save convergence info for fields */
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_set_key_struct(f, key_sinfo_id, &sinfo);
  }

  /*==========================================================================
   * After having computed the new value, an estimator is computed for the
   * prediction step of the velocity.
   *==========================================================================*/

  if (iescap > 0 && stride == 3) {
    /* Computation of the estimator of the current component */

    /* smbini does not contain unsteady terms and mass source terms
       of the RHS updated at each sweep */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++) {
        cs_real_t diff = 0.;
        for (cs_lnum_t j = 0; j < stride; j++)
          diff += fimp[c_id][i][j]*(pvar[c_id][j] - pvara[c_id][j]);

        smbrp[c_id][i] = smbini[c_id][i] - diff;
      }
    });

    ctx.wait();

    /* need to recompute face value if below increment is zero
       else the face value is given from the last isweep iteration */
    if (inc == 0) {
      cs_boundary_conditions_update_bc_coeff_face_values<stride>
        (ctx, f, bc_coeffs,
         1, // inc
         eqp, pvar,
         val_ip, val_f, val_f_lim,
         val_f_d, val_f_d_lim);
    }
    inc = 1;

    /* Without relaxation even for a stationnary computation */

    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,
                      ivisep,
                      eqp,
                      (cs_real_3_t *)pvar,
                      (const cs_real_3_t *)pvara,
                      bc_coeffs,
                      &bc_coeffs_solve,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      i_secvis,
                      b_secvis,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      nullptr,
                      nullptr,
                      (cs_real_3_t *)smbrp);

    /* Contribution of the current component to the L2 norm stored in eswork */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        eswork[c_id][i] = cs_math_pow2(smbrp[c_id][i] / cell_vol[c_id]);
    });
    ctx.wait();

  }

  /*==========================================================================
   * Store face value for gradient and diffusion
   *==========================================================================*/

  cs_boundary_conditions_ensure_bc_coeff_face_values_allocated
    (bc_coeffs,
     n_b_faces,
     stride,
     cs_alloc_mode,
     (df_limiter_id > -1 || ircflb != 1));

  var_t *val_f_updated = (var_t *)bc_coeffs->val_f;
  var_t *val_f_lim_updated = (var_t *)bc_coeffs->val_f_lim;
  var_t *val_f_d_updated = (var_t *)bc_coeffs->val_f_d;
  var_t *val_f_d_lim_updated = (var_t *)bc_coeffs->val_f_d_lim;

  cs_boundary_conditions_update_bc_coeff_face_values<stride>
    (ctx, f, bc_coeffs, 1, eqp, pvar, val_ip,
     val_f_updated, val_f_lim_updated,
     val_f_d_updated, val_f_d_lim_updated);

  /*==========================================================================
   * Free solver setup
   *==========================================================================*/

  cs_sles_free(sc);
  cs_matrix_release(&a);

  if (stride == 3) {
    /* Save diagonal in case we want to use it */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        for (cs_lnum_t j = 0; j < stride; j++)
          fimp[c_id][i][j] = dam[c_id][i][j];
    });

    ctx.wait();
  }

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(dam);
  CS_FREE_HD(smbini);
  CS_FREE_HD(dpvar);
  if (iswdyp >= 1) {
    CS_FREE_HD(adxk);
    CS_FREE_HD(adxkm1);
    CS_FREE_HD(dpvarm1);
    CS_FREE_HD(rhs0);
  }

  _clear_bc_coeffs_solve(bc_coeffs_solve);
}

#endif /* cplusplus */

/*============================================================================
 * Public function definitions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_equation_iterative_solve.cpp
 *
 * \brief This file gathers functions that solve advection diffusion equations
 * with source terms for one time step for a scalar, vector or tensor variable.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve an advection diffusion equation with source
 * terms for one time step for the variable \f$ a \f$.
 *
 * The equation reads:
 *
 * \f[
 * f_s^{imp}(a^{n+1}-a^n)
 * + \divs \left( a^{n+1} \rho \vect{u} - \mu \grad a^{n+1} \right)
 * = Rhs
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * f_s^{imp} \delta a
 * + \divs \left( \delta a \rho \vect{u} - \mu \grad \delta a \right)
 * = Rhs^1
 * \f]
 *
 * where \f$ \delta a = a^{n+1} - a^n\f$ and
 * \f$ Rhs^1 = Rhs - \divs( a^n \rho \vect{u} - \mu \grad a^n)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * f_s^{imp} \delta a^k
 * + \divs \left(\delta a^k \rho \vect{u}-\mu\grad\delta a^k \right)
 * = Rhs^k
 * \f]
 *
 * where \f$Rhs^k=Rhs-f_s^{imp}(a^k-a^n)
 * - \divs \left( a^k\rho\vect{u}-\mu\grad a^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ f_s^{imp} \f$ here!
 *
 * Please refer to the
 * <a href="../../theory.pdf#codits"><b>codits</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     iterns        external sub-iteration number
 * \param[in]     f_id          field id (or -1)
 * \param[in]     name          associated name if f_id < 0, or nullptr
 * \param[in]     iescap        compute the predictor indicator if 1
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     normp         Reference norm to solve the system (optional)
 *                              if negative: recomputed here
 * \param[in]     eqp           pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ a^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ a^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
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
 * \param[in]     rovsdt        \f$ f_s^{imp} \f$
 * \param[in]     smbrp         Right hand side \f$ Rhs^k \f$
 * \param[in,out] pvar          current variable
 * \param[out]    dpvar         last variable increment
 * \param[in]     xcpp          array of specific heat (Cp)
 * \param[out]    eswork        prediction-stage error estimator
 *                              (if iescap > 0)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_scalar(int                   idtvar,
                                   int                   iterns,
                                   int                   f_id,
                                   const char           *name,
                                   int                   iescap,
                                   int                   imucpp,
                                   cs_real_t             normp,
                                   cs_equation_param_t   *eqp,
                                   const cs_real_t       pvara[],
                                   const cs_real_t       pvark[],
                                   const cs_field_bc_coeffs_t *bc_coeffs,
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   const cs_real_t       i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   cs_real_6_t           viscel[],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   const cs_real_t       rovsdt[],
                                   cs_real_t             smbrp[],
                                   cs_real_t             pvar[],
                                   cs_real_t             dpvar[],
                                   const cs_real_t       xcpp[],
                                   cs_real_t             eswork[])
{
  /* Local variables */
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int iwarnp = eqp->verbosity;
  int iswdyp = eqp->iswdyn;
  int ndircp = eqp->ndircl;
  cs_real_t epsrsp = eqp->epsrsm;
  cs_real_t epsilp = eqp->epsilo;
  cs_real_t thetap = eqp->theta;

  bool _need_solve = cs_time_control_is_active(eqp->time_control,
                                               cs_glob_time_step);
  if (!_need_solve)
    return;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  int *c_disable_flag = mq->c_disable_flag;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx, ctx_c;
#if defined(HAVE_CUDA)
  ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  int inc, niterf;
  int lvar, imasac;
  cs_real_t residu, rnorm, ressol;
  cs_real_t thetex, nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk;

  cs_real_t rnorm2 = 0, alph = 0, beta = 0;

  cs_solving_info_t sinfo;

  int coupling_id = -1;
  cs_field_t *f = nullptr;
  cs_real_t *smbini = nullptr;
  cs_real_t *w1 = nullptr;

  bool conv_diff_mg = false;

  /*============================================================================
   * 0.  Initialization
   *==========================================================================*/

  /* Name */
  const char *var_name = cs_sles_name(f_id, name);

  if (iwarnp >= 1)
    bft_printf("Equation iterative solve of: %s\n", var_name);

  /* solving info */
  int key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);
    coupling_id = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
  }

  /* Determine if we are in a case with special requirements */

  if (coupling_id < 0 && iconvp > 0) {
    cs_sles_t *sc = cs_sles_find_or_add(f_id, name);
    const char *sles_type = cs_sles_get_type(sc);
    if (strcmp(sles_type, "cs_multigrid_t") == 0)
      conv_diff_mg = true;
  }

  /* Allocate temporary arrays */

  cs_real_t *adxk = nullptr, *adxkm1 = nullptr;
  cs_real_t *dpvarm1 = nullptr, *rhs0 = nullptr;

  if (iswdyp >= 1) {
    CS_MALLOC_HD(adxk, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(adxkm1, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(dpvarm1, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(rhs0, n_cells_ext, cs_real_t, cs_alloc_mode);
  }

  /* Symmetric matrix, except if advection */
  int isym = 1;
  if (iconvp > 0) isym = 2;

  bool symmetric = (isym == 1) ? true : false;

  /* Periodicity has to be taken into account */

  /* Initialization for test before matrix vector product
     for computation of initial residual */

  /*============================================================================
   * 1.  Building of the "simplified" matrix
   *==========================================================================*/

  cs_matrix_t *a = cs_sles_default_get_matrix(f_id, var_name, 1, 1, symmetric);

  /* For steady computations, the diagonal is relaxed */
  cs_real_t relaxp = (idtvar < 0) ? eqp->relaxv : 1.;

  cs_matrix_compute_coeffs(a,
                           f,
                           iconvp,
                           idiffp,
                           ndircp,
                           thetap,
                           relaxp,
                           imucpp,
                           bc_coeffs,
                           rovsdt,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm,
                           xcpp);

  /*==========================================================================
   * 2. Iterative process to handle non orthogonalities (starting from the
   *    second iteration).
   *==========================================================================*/

  /* Prepare the computation of fluxes at faces if needed */

  cs_real_t *i_flux = nullptr, *b_flux = nullptr;
  cs_real_t *b_flux_k = nullptr, *b_flux_km1 = nullptr;
  cs_real_2_t *i_flux_0 = nullptr, *i_flux_k = nullptr, *i_flux_km1 = nullptr;

  if (f_id > -1) {
    f = cs_field_by_id(f_id);

    const int kiflux = cs_field_key_id("inner_flux_id");
    const int kbflux = cs_field_key_id("boundary_flux_id");

    int i_flux_id = cs_field_get_key_int(f, kiflux);
    int b_flux_id = cs_field_get_key_int(f, kbflux);

    if (i_flux_id > -1 && b_flux_id > -1) {
      assert(idtvar != -1);
      /* flux is non-conservative with the steady algorithm but flux field is of
         dimension 1. Forbidden at parameters check */

      i_flux = cs_field_by_id(i_flux_id)->val;
      b_flux = cs_field_by_id(b_flux_id)->val;

      CS_MALLOC_HD(i_flux_0, n_i_faces, cs_real_2_t, cs_alloc_mode);
      CS_MALLOC_HD(i_flux_k, n_i_faces, cs_real_2_t, cs_alloc_mode);
      CS_MALLOC_HD(i_flux_km1, n_i_faces, cs_real_2_t, cs_alloc_mode);
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        i_flux_0[face_id][0] = 0.;
        i_flux_0[face_id][1] = 0.;
        i_flux_k[face_id][0] = 0.;
        i_flux_k[face_id][1] = 0.;
      });
      CS_MALLOC_HD(b_flux_k, n_b_faces, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(b_flux_km1, n_b_faces, cs_real_t, cs_alloc_mode);
      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        b_flux_k[face_id] = 0.;
      });

      ctx_c.wait();
    }
  }

  /* Application of the theta-scheme */

  /* On calcule le bilan explicite total */
  thetex = 1. - thetap;

  /* Compute the min/ max limiter */
  inc = 1;

  if (f_id > -1)
    cs_beta_limiter_building(f_id, inc, rovsdt);

  /* If thetex = 0, no need to do more */
  if (fabs(thetex) > cs_math_epzero) {
    inc = 1;

    /* The added convective scalar mass flux is:
       (thetex*Y_\face-imasac*Y_\celli)*mf.
       When building the explicit part of the rhs, one
       has to impose 0 on mass accumulation. */
    imasac = 0;

    eqp->theta = thetex;

    /* Compute - Con-Diff((1-theta) Y^n )
     * where Y^n is pvara
     * Note: this part does not need any update, and will be stored in
     * smbini */
    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      eqp,
                      nullptr, /* pvar == pvara */
                      pvara,
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
                      smbrp,
                      i_flux_0,
                      b_flux);

    eqp->theta = thetap;

  }

  /* Before looping, the RHS without reconstruction is stored in smbini */

  CS_MALLOC_HD(smbini, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_lnum_t has_dc = mq->has_disable_flag;

  ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    smbini[cell_id] = smbrp[cell_id];
    smbrp[cell_id] = 0.;
  });

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    pvar[cell_id] = pvark[cell_id];
  });

  /* In the following, cs_balance_scalar is called with inc=1,
     except for Weight Matrix (nswrsp=-1) */
  inc = 1;

  if (eqp->nswrsm == -1) {
    eqp->nswrsm = 1;
    inc = 0;
  }

  /* Incrementation and rebuild of right hand side */

  /* The added convective scalar mass flux is:
     (thetap*Y_\face-imasac*Y_\celli)*mf.
     When building the implicit part of the rhs, one
     has to impose 1 on mass accumulation. */
  imasac = 1;

  ctx_c.wait(); /* We now need pvar, computed by ctx_c */
  ctx.wait();   /* We now need smbrp, computed by ctx */

  /* Compute - Con-Diff(theta Y^k )
   * where Y^k is pvar (possible over iteration for
   * velocity pressure coupling) */
  cs_balance_scalar(idtvar,
                    f_id,
                    imucpp,
                    imasac,
                    inc,
                    eqp,
                    pvar,
                    pvara,
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
                    smbrp,
                    i_flux_k,
                    b_flux_k);

  if (iswdyp >= 1) {
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t cell_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      rhs0[cell_id] = smbrp[cell_id];
      smbrp[cell_id]  += smbini[cell_id]
                       - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      sum += cs_math_pow2(smbrp[cell_id]);
      adxkm1[cell_id] = 0.;
      adxk[cell_id] = 0.;
      dpvar[cell_id] = 0.;
    });

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t cell_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      smbrp[cell_id]  += smbini[cell_id]
                       - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      sum += cs_math_pow2(smbrp[cell_id]);
    });
  }

  ctx.wait();
  cs_parall_sum(1, CS_DOUBLE, &residu);

  /* --- Right hand side residual */
  residu = sqrt(residu);

  if (normp > 0.)
    rnorm = normp;

  else {
    /* Normalization residual
       (L2-norm of B.C. + source terms + non-orthogonality terms)

       Caution: when calling a matrix-vector product, here for a variable
       which is not "by increments" and is assumed initialized, including
       for ghost values. */

    /* Allocate a temporary array */

    // Number of local ghost cells may be different from that of mesh
    // in case of internal coupling.
    cs_lnum_t n_cols = cs_matrix_get_n_columns(a);

    CS_MALLOC_HD(w1, n_cols, cs_real_t, cs_alloc_mode);

    cs_real_t *w2;
    CS_MALLOC_HD(w2, n_cols, cs_real_t, cs_alloc_mode);

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X^n = %f\n", p_mean);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      w2[cell_id] = (pvar[cell_id]-p_mean);
    });

    ctx.wait();

    cs_matrix_vector_multiply(a, w2, w1);

    CS_FREE_HD(w2);

    if (iwarnp >= 2) {
      struct cs_data_2r rd2;
      struct cs_reduce_sum2r reducer2;
      ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                              (cs_lnum_t c_id, cs_data_2r &sum) {
        sum.r[0] = cs_math_pow2(w1[c_id]);
        sum.r[1] = cs_math_pow2(smbrp[c_id]);
      });
      ctx.wait();
      cs_parall_sum(2, CS_DOUBLE, rd2.r);
      bft_printf("L2 norm ||AX^n|| = %f\n", sqrt(rd2.r[0]));
      bft_printf("L2 norm ||B^n|| = %f\n",  sqrt(rd2.r[1]));
    }

    ctx.parallel_for_reduce_sum(n_cells, rnorm2, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t cell_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      w1[cell_id] += smbrp[cell_id];
      /* Remove contributions from penalized cells */
      if (has_dc * c_disable_flag[has_dc * cell_id] != 0)
        w1[cell_id] = 0.;

      sum += cs_math_pow2(w1[cell_id]);
    });

    ctx.wait();
    cs_parall_sum(1, CS_DOUBLE, &rnorm2);

    rnorm = sqrt(rnorm2);
  }

  sinfo.rhs_norm = rnorm;

  /* Free memory */
  CS_FREE_HD(w1);

  /* Warning: for weight matrix, one and only one sweep is done. */
  int nswmod = cs::max(eqp->nswrsm, 1);

  /* Reconstruction loop (beginning) */
  if (iterns <= 1)
    sinfo.n_it = 0;

  cs_sles_t *sc = cs_sles_find_or_add(f_id, var_name);

  int isweep = 1;

  /* Main loop of reconstruction
   * --------------------------- */

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {

    /* Solving on the increment: dpvar */

    if (iswdyp >= 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        dpvarm1[cell_id] = dpvar[cell_id];
      });
    }

    /* Store last fluxes if needed
     * (switch pointers) */
    cs_real_2_t *_temp_i = i_flux_km1;
    i_flux_km1 = i_flux_k;
    i_flux_k = _temp_i;
    cs_real_t *_temp_b = b_flux_km1;
    b_flux_km1 = b_flux_k;
    b_flux_k = _temp_b;
    if (i_flux_k != nullptr) {
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        i_flux_k[face_id][0] = 0.;
        i_flux_k[face_id][1] = 0.;
      });
      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        b_flux_k[face_id] = 0.;
      });
      ctx.wait();
    }

    /* Solver residual */
    ressol = residu;

    if (conv_diff_mg) {
      cs_multigrid_t *mg
        = static_cast<cs_multigrid_t *>(cs_sles_get_context(sc));
      cs_multigrid_setup_conv_diff(mg, var_name, a, true,
                                   cs_sles_get_verbosity(sc));
    }

    cs_sles_solve_ccc_fv(sc,
                         a,
                         epsilp,
                         rnorm,
                         &niterf,
                         &ressol,
                         smbrp,
                         dpvar);

    /* Dynamic relaxation of the system */
    if (iswdyp >= 1) {

      /* Computation of the variable relaxation coefficient */
      lvar = -1;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        adxkm1[cell_id] = adxk[cell_id];
        adxk[cell_id] = - rhs0[cell_id];
      });

      ctx.wait(); /* We now need adxk, computed by ctx */

      cs_balance_scalar(idtvar,
                        lvar,
                        imucpp,
                        imasac,
                        inc,
                        eqp,
                        dpvar,
                        nullptr, /* dpvara == dpvar */
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
                        adxk,
                        nullptr,
                        nullptr);

      /* ||E.dx^(k-1)-E.0||^2 */
      nadxkm1 = nadxk;

      struct cs_data_2r rd2;
      struct cs_reduce_sum2r reducer2;
      ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                              (cs_lnum_t c_id, cs_data_2r &sum) {
        sum.r[0] = adxk[c_id] * adxk[c_id];
        sum.r[1] = smbrp[c_id] * adxk[c_id];
      });
      ctx.wait();
      cs_parall_sum(2, CS_DOUBLE, rd2.r);

      nadxk = rd2.r[0];    /* ||E.dx^k-E.0||^2 */
      paxkrk = rd2.r[1];   /* < E.dx^k-E.0; r^k > */

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {

        ctx.parallel_for_reduce(n_cells, rd2, reducer2, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id, cs_data_2r &sum) {
          sum.r[0] = smbrp[c_id] * adxkm1[c_id];
          sum.r[1] = adxk[c_id] * adxkm1[c_id];
        });

        ctx.wait();
        cs_parall_sum(2, CS_DOUBLE, rd2.r);

        paxm1rk = rd2.r[0];   // < E.dx^(k-1)-E.0; r^k >
        paxm1ax = rd2.r[1];   // < E.dx^(k-1)-E.0; E.dx^k -E.0 >

        if (   nadxkm1 > 1e-30*rnorm2
            && (nadxk*nadxkm1-pow(paxm1ax,2)) > 1e-30*rnorm2)
          beta =   (paxkrk*paxm1ax - nadxk*paxm1rk)
                 / (nadxk*nadxkm1-cs_math_pow2(paxm1ax));
        else
          beta = 0.;

      }
      else {
        beta = 0.;
        paxm1ax = 1.;
        paxm1rk = 0.;
        paxm1ax = 0.;
      }

      /* The first sweep is not relaxed */
      if (isweep == 1) {
        alph = 1.;
        beta = 0.;
      }
      else if (isweep == 2) {
        beta = 0.;
        cs_fp_exception_disable_trap();
        alph = -paxkrk/fmax(nadxk, 1e-30*rnorm2);
        cs_fp_exception_restore_trap();
        if (isnan(alph))
          alph = 0;
      }
      else {
        cs_fp_exception_disable_trap();
        alph = -(paxkrk + beta*paxm1ax)/fmax(nadxk, 1e-30*rnorm2);
        cs_fp_exception_restore_trap();
        if (isnan(alph))
          alph = 0;
      }

      /* Writing */
      if (iwarnp >= 3)
        bft_printf("%s Sweep: %d Dynamic relaxation: alpha = %12.5e, "
                   "beta = %12.5e,\n< dI^k :  R^k >   = %12.5e, "
                   "||dI^k  ||^2 = %12.5e,\n< dI^k-1 : R^k >  = %12.5e, "
                   "||dI^k-1||^2 = %12.5e,\n< dI^k-1 : dI^k > = %12.5e\n",
                   var_name, isweep, alph, beta, paxkrk, nadxk, paxm1rk,
                   nadxkm1, paxm1ax);
    }

    /* Update the solution with the increment, update the right hand side
       and compute the new residual */

    if (iswdyp <= 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        pvar[cell_id] += dpvar[cell_id];
        /* smbini does not contain unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbrp[cell_id] = smbini[cell_id]
                       - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      });
    }
    else if (iswdyp == 1) {
      if (alph < 0.) break;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        pvar[cell_id] += alph*dpvar[cell_id];
        /* smbini does not contain unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbrp[cell_id] = smbini[cell_id]
                       - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      });
    }
    else if (iswdyp >= 2) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        pvar[cell_id] += alph*dpvar[cell_id] + beta*dpvarm1[cell_id];
        /* smbini does not contain unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbrp[cell_id] = smbini[cell_id]
                       - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      });
    }

    ctx.wait();

    /*  ---> Handle parallelism and periodicity */
    cs_halo_sync(m->halo, ctx.use_gpu(), pvar);

    /* Compute the beta (min/max) limiter */
    if (f_id > -1)
      cs_beta_limiter_building(f_id, inc, rovsdt);

    /* The added convective scalar mass flux is:
       (theta*Y_\face-imasac*Y_\celli)*mf.
       When building the implicit part of the rhs, one
       has to impose 1 on mass accumulation. */
    imasac = 1;

    ctx.wait();

    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      eqp,
                      pvar,
                      pvara,
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
                      smbrp,
                      i_flux_k,
                      b_flux_k);

    /* --- Convergence test */
    ctx.parallel_for_reduce_sum(n_cells, residu, [=] CS_F_HOST_DEVICE
                                (cs_lnum_t c_id,
                                 CS_DISPATCH_REDUCER_TYPE(double) &sum) {
      sum += cs_math_pow2(smbrp[c_id]);
    });
    ctx.wait();
    cs_parall_sum(1, CS_DOUBLE, &residu);
    residu = sqrt(residu);

    /* Writing */
    sinfo.n_it = sinfo.n_it + niterf;

    /* Writing */
    if (iwarnp >= 2) {
      bft_printf("%s: CV_DIF_TS, IT: %d, Res: %12.5e, Norm: %12.5e\n",
                 var_name, isweep, residu, rnorm);
      bft_printf("%s: Current reconstruction sweep: %d, "
                 "Iterations for solver: %d\n", var_name, isweep, niterf);
    }

    isweep++;

  }
  /* --- Reconstruction loop (end) */

  /* Writing: convergence */
  if (fabs(rnorm) > cs_math_epzero)
    sinfo.res_norm = residu/rnorm;
  else
    sinfo.res_norm = 0.;

  /* Save convergence info */
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_set_key_struct(f, key_sinfo_id, &sinfo);
  }

  if (iwarnp >= 1) {
    if (residu <= epsrsp*rnorm)
      bft_printf("%s: CV_DIF_TS, IT : %d, Res : %12.5e, Norm : %12.5e\n",
                 var_name, isweep-1, residu, rnorm);
    /* Writing: non-convergence */
    else if (isweep > nswmod)
      bft_printf("@\n@ @@ WARNING: %s CONVECTION-DIFFUSION-SOURCE-TERMS\n@"
                 "=========\n@  Maximum number of iterations %d reached\n",
                 var_name,nswmod);
  }

  /* Finalize the face fluxes if needed */
  /* last increment in upwind to fulfill exactly the considered
     balance equation */

  if (i_flux != nullptr) {
    inc  = 1;
    imasac = 0; /* mass accumulation not taken into account */

    inc = 0;
    cs_equation_param_t eqp_loc;
    int k_id = cs_field_key_id("var_cal_opt");

    /* copy of the equation param structure in eqp_loc */
    cs_field_get_key_struct(f, k_id, &eqp_loc);

    eqp_loc.nswrgr = 0;
    eqp_loc.blencv = 0.;

    cs_face_convection_scalar(idtvar,
                              f_id,
                              eqp_loc,
                              icvflb,
                              inc,
                              imasac,
                              dpvar,
                              pvara,
                              icvfli,
                              bc_coeffs,
                              i_massflux,
                              b_massflux,
                              i_flux_0,
                              b_flux);

    /* FIXME diffusion part */

    /* Finalize the convective flux calculation ,
     * Note that the two sides are equal because imasac=0 */
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_flux[face_id] += i_flux_km1[face_id][0] + i_flux_0[face_id][0];
    });
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      b_flux[face_id] += b_flux_km1[face_id];
    });

    ctx.wait();

    CS_FREE_HD(i_flux_0);
    CS_FREE_HD(i_flux_k);
    CS_FREE_HD(i_flux_km1);
    CS_FREE_HD(b_flux_k);
    CS_FREE_HD(b_flux_km1);
  }

  /*==========================================================================
   * 3. After having computed the new value, an estimator is computed for the
   * prediction step of the velocity.
   *==========================================================================*/

  if (iescap > 0) {
    /*  Computation of the estimator of the current component */

    /* smbini does not contain unsteady terms and mass source terms
       of the RHS updated at each sweep */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      smbrp[cell_id] = smbini[cell_id]
                     - rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
    });

    inc = 1;

    ctx.wait();

    /* Without relaxation even for a stationary computation */

    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      eqp,
                      pvar,
                      pvara,
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
                      smbrp,
                      nullptr,
                      nullptr);

    /* Contribution of the current component to the L2 norm stored in eswork */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      eswork[cell_id] = cs_math_pow2(smbrp[cell_id] / cell_vol[cell_id]);
    });
    ctx.wait();

  }

  /*==========================================================================
   * 4. Free solver setup
   *==========================================================================*/

  cs_sles_free(sc);
  cs_sles_default_release_matrix(&a);

  ctx.wait();

  /*  Free memory */
  CS_FREE_HD(smbini);

  if (iswdyp >= 1) {
    CS_FREE_HD(adxk);
    CS_FREE_HD(adxkm1);
    CS_FREE_HD(dpvarm1);
    CS_FREE_HD(rhs0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function solves an advection diffusion equation with source terms
 * for one time step for the vector variable \f$ \vect{a} \f$.
 *
 * The equation reads:
 *
 * \f[
 * \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
 * + \divv \left( \vect{a}^{n+1} \otimes \rho \vect {u}
 *              - \mu \gradt \vect{a}^{n+1}\right)
 * = \vect{Rhs}
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}
 * + \divv \left( \delta \vect{a} \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a} \right)
 * = \vect{Rhs}^1
 * \f]
 *
 * where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
 * \f$ \vect{Rhs}^1 = \vect{Rhs}
 * - \divv \left( \vect{a}^n \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^n \right)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \vect{a}^k
 * + \divv \left( \delta \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \delta \vect{a}^k \right)
 * = \vect{Rhs}^k
 * \f]
 *
 * where \f$ \vect{Rhs}^k = \vect{Rhs}
 * - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
 * - \divv \left( \vect{a}^k \otimes \rho \vect{u}
 *              - \mu \gradt \vect{a}^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
 *
 * \param[in]      idtvar        indicator of the temporal scheme
 * \param[in]      iterns        external sub-iteration number
 * \param[in]      f_id          field id (or -1)
 * \param[in]      name          associated name if f_id < 0, or nullptr
 * \param[in]      ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]      iescap        compute the predictor indicator if >= 1
 * \param[in]      eqp   pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]      pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]      pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               \c pvara (usually \c pvar= \c pvara)
 * \param[in, out] bc_coeffs_v   boundary condition structure for the variable
 * \param[in]      i_massflux    mass flux at interior faces
 * \param[in]      b_massflux    mass flux at boundary faces
 * \param[in]      i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]      b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
 * \param[in]      i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]      b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the r.h.s.
 * \param[in]      i_secvis      secondary viscosity at interior faces
 * \param[in]      b_secvis      secondary viscosity at boundary faces
 * \param[in]      viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]      weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]      weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in]      icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]      icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in, out] fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in, out] smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in, out] pvar          current variable
 * \param[out]     eswork        prediction-stage error estimator
 *                               (if iescap >= 0)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_vector(int                   idtvar,
                                   int                   iterns,
                                   int                   f_id,
                                   const char           *name,
                                   int                   ivisep,
                                   int                   iescap,
                                   cs_equation_param_t  *eqp,
                                   const cs_real_t       pvara[][3],
                                   const cs_real_t       pvark[][3],
                                   cs_field_bc_coeffs_t *bc_coeffs_v,
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   const cs_real_t       i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   const cs_real_t       i_secvis[],
                                   const cs_real_t       b_secvis[],
                                   cs_real_t             viscel[][6],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   cs_real_t             fimp[][3][3],
                                   cs_real_t             smbrp[][3],
                                   cs_real_t             pvar[][3],
                                   cs_real_t             eswork[][3])
{
  _equation_iterative_solve_strided<3>(idtvar,
                                       iterns,
                                       f_id,
                                       name,
                                       ivisep,
                                       iescap,
                                       eqp,
                                       pvara,
                                       pvark,
                                       bc_coeffs_v,
                                       i_massflux,
                                       b_massflux,
                                       i_viscm,
                                       b_viscm,
                                       i_visc,
                                       b_visc,
                                       i_secvis,
                                       b_secvis,
                                       viscel,
                                       weighf,
                                       weighb,
                                       icvflb,
                                       icvfli,
                                       fimp,
                                       smbrp,
                                       pvar,
                                       eswork);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function solves an advection diffusion equation with source
 * terms for one time step for the symmetric tensor variable
 * \f$ \tens{\variat} \f$.
 *
 * The equation reads:
 *
 * \f[
 * \tens{f_s}^{imp}(\tens{\variat}^{n+1}-\tens{\variat}^n)
 * + \divt \left( \tens{\variat}^{n+1} \otimes \rho \vect {u}
 *              - \mu \gradtt \tens{\variat}^{n+1}\right)
 * = \tens{Rhs}
 * \f]
 *
 * This equation is rewritten as:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \tens{\variat}
 * + \divt \left( \delta \tens{\variat} \otimes \rho \vect{u}
 *              - \mu \gradtt \delta \tens{\variat} \right)
 * = \tens{Rhs}^1
 * \f]
 *
 * where \f$ \delta \tens{\variat} = \tens{\variat}^{n+1} - \tens{\variat}^n\f$
 * and \f$ \tens{Rhs}^1 = \tens{Rhs}
 * - \divt \left( \tens{\variat}^n \otimes \rho \vect{u}
 *              - \mu \gradtt \tens{\variat}^n \right)\f$
 *
 *
 * It is in fact solved with the following iterative process:
 *
 * \f[
 * \tens{f_s}^{imp} \delta \tens{\variat}^k
 * + \divt \left( \delta \tens{\variat}^k \otimes \rho \vect{u}
 *              - \mu \gradtt \delta \tens{\variat}^k \right)
 * = \tens{Rhs}^k
 * \f]
 *
 * where \f$ \tens{Rhs}^k = \tens{Rhs}
 * - \tens{f_s}^{imp} \left(\tens{\variat}^k-\tens{\variat}^n \right)
 * - \divt \left( \tens{\variat}^k \otimes \rho \vect{u}
 *              - \mu \gradtt \tens{\variat}^k \right)\f$
 *
 * Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     name          associated name if f_id < 0, or nullptr
 * \param[in]     eqp   pointer to a cs_equation_param_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in,out] bc_coeffs_ts  boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the matrix
 * \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at boundary faces for the matrix
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
 * \param[in,out] fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in,out] smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in,out] pvar          current variable
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_tensor(int                         idtvar,
                                   int                         f_id,
                                   const char                 *name,
                                   cs_equation_param_t        *eqp,
                                   const cs_real_t             pvara[][6],
                                   const cs_real_t             pvark[][6],
                                   cs_field_bc_coeffs_t       *bc_coeffs_ts,
                                   const cs_real_t             i_massflux[],
                                   const cs_real_t             b_massflux[],
                                   const cs_real_t             i_viscm[],
                                   const cs_real_t             b_viscm[],
                                   const cs_real_t             i_visc[],
                                   const cs_real_t             b_visc[],
                                   cs_real_t                   viscel[][6],
                                   const cs_real_2_t           weighf[],
                                   const cs_real_t             weighb[],
                                   int                         icvflb,
                                   const int                   icvfli[],
                                   cs_real_t                   fimp[][6][6],
                                   cs_real_t                   smbrp[][6],
                                   cs_real_t                   pvar[][6])
{
  _equation_iterative_solve_strided<6>(idtvar,
                                       -1, // iterns
                                       f_id,
                                       name,
                                       -1,  // ivisep
                                       -1, // iescap
                                       eqp,
                                       pvara,
                                       pvark,
                                       bc_coeffs_ts,
                                       i_massflux,
                                       b_massflux,
                                       i_viscm,
                                       b_viscm,
                                       i_visc,
                                       b_visc,
                                       nullptr, // i_secvis
                                       nullptr, // b_secvis
                                       viscel,
                                       weighf,
                                       weighb,
                                       icvflb,
                                       icvfli,
                                       fimp,
                                       smbrp,
                                       pvar,
                                       nullptr); // eswork
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
