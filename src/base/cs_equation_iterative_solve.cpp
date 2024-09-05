/*============================================================================
 * Convection diffusion equation solver.
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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_convection_diffusion.h"
#include "cs_dispatch.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_gradient.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_parall.h"
#include "cs_matrix_building.h"
#include "cs_matrix_default.h"
#include "cs_sles.h"
#include "cs_sles_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_iterative_solve.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

#ifdef __cplusplus

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
                                  const cs_field_bc_coeffs_t *bc_coeffs,
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

  int iconvp = eqp->iconv;
  int idiffp = eqp->idiff;
  int iwarnp = eqp->verbosity;
  int iswdyp = eqp->iswdyn;
  int idftnp = eqp->idften;
  int ndircp = eqp->ndircl;
  double epsrsp = eqp->epsrsm;
  double epsilp = eqp->epsilo;
  double relaxp = eqp->relaxv;
  double thetap = eqp->theta;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  using var_t = cs_real_t[stride];
  using m_t = cs_real_t[stride][stride];

  int inc, niterf;
  int lvar, imasac;
  cs_real_t residu, ressol;
  cs_real_t nadxk, paxm1rk, paxm1ax;

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
  cs_lnum_t db_size = stride;
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
  if (CS_F_(vel) != NULL && CS_F_(vel)->id == f_id) {

    i_vf = cs_field_by_name_try("inner_face_velocity");
    if (i_vf != nullptr) {
      cs_arrays_set_value<cs_real_t, 1>(3*n_i_faces, 0., i_vf->val);
      i_pvar = (var_t *)i_vf->val;
    }

    b_vf = cs_field_by_name_try("boundary_face_velocity");
    if (b_vf != nullptr) {
      cs_arrays_set_value<cs_real_t, 1>(3*n_b_faces, 0., b_vf->val);
      b_pvar = (var_t *)b_vf->val;
    }

  }

  /* solving info */
  int key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);
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

  /*  be careful here, xam is interleaved*/

  cs_lnum_t eb_stride = eb_size*eb_size;
  cs_real_t *xam;
  CS_MALLOC_HD(xam, eb_stride*isym*n_faces, cs_real_t, amode);

  /*==========================================================================
   * Building of the "simplified" matrix
   *==========================================================================*/

  int tensorial_diffusion = 1;

  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION)
    tensorial_diffusion = 2;

  if (stride == 3)
    cs_matrix_wrapper_vector(iconvp,
                             idiffp,
                             tensorial_diffusion,
                             ndircp,
                             isym,
                             eb_size,
                             thetap,
                             bc_coeffs,
                             (const cs_real_33_t *)fimp,
                             i_massflux,
                             b_massflux,
                             i_viscm,
                             b_viscm,
                             (cs_real_33_t *)dam,
                             xam);
  else if (stride == 6)
    cs_matrix_wrapper_tensor(iconvp,
                             idiffp,
                             tensorial_diffusion,
                             ndircp,
                             isym,
                             thetap,
                             bc_coeffs,
                             (const cs_real_66_t *)fimp,
                             i_massflux,
                             b_massflux,
                             i_viscm,
                             b_viscm,
                             (cs_real_66_t *)dam,
                             xam);

  /* Precaution if diagonal is 0, which may happen is all surrounding cells
   * are disabled
   * If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++)
        if (cs_math_fabs(dam[c_id][i][i]) < DBL_MIN)
          dam[c_id][i][i] += 1.;
    });
  }

  /*  For steady computations, the diagonal is relaxed */
  if (idtvar < 0) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < stride; i++) {
        for (cs_lnum_t j = 0; j < stride; j++)
          dam[c_id][i][j] /= relaxp;
      }
    });
  }

  /*===========================================================================
   * Iterative process to handle non orthogonlaities (starting from the
   * second iteration).
   *===========================================================================*/

  /* Application of the theta-scheme */

  /* We compute the total explicit balance. */
  cs_real_t thetex = 1. - thetap;

  /* If THETEX=0, no need to add anything */
  if (cs_math_fabs(thetex) > cs_math_epzero) {
    inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
    imasac = 0;

    eqp->theta = thetex;

    if (stride == 3)
      cs_balance_vector(idtvar,
                        f_id,
                        imasac,
                        inc,
                        ivisep,
                        eqp,
                        nullptr, /* pvar == pvara */
                        (const cs_real_3_t *)pvara,
                        bc_coeffs,
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
                        nullptr, /* pvar == pvara */
                        (const cs_real_6_t *)pvara,
                        bc_coeffs,
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
    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      smbini[c_id][isou] = smbrp[c_id][isou];
      smbrp[c_id][isou] = 0.;
    }
  });

  /* pvar is initialized on n_cells_ext to avoid a synchronization */
  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t isou = 0; isou < stride; isou++)
      pvar[c_id][isou] = pvark[c_id][isou];
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

  if (CS_F_(vel) != NULL && CS_F_(vel)->id == f_id) {
    f = cs_field_by_name_try("velocity_explicit_balance");

    if (f != nullptr) {
      cs_real_3_t *cpro_cv_df_v = (cs_real_3_t *)f->val;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          cpro_cv_df_v[c_id][isou] = smbrp[c_id][isou];
      });
    }
  }

  /* Dynamic relaxation*/
  if (iswdyp >= 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        rhs0[c_id][isou] = smbrp[c_id][isou];

        cs_real_t diff = 0.;
        for (cs_lnum_t j = 0; j < stride; j++)
          diff += fimp[c_id][isou][j]*(pvar[c_id][j] - pvara[c_id][j]);

        smbini[c_id][isou] -= diff;
        smbrp[c_id][isou] += smbini[c_id][isou];

        adxkm1[c_id][isou] = 0.;
        adxk[c_id][isou] = 0.;
        dpvar[c_id][isou] = 0.;
      }
    });

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t isou = 0; isou < stride; isou++) {

        cs_real_t diff = 0.;
        for (cs_lnum_t j = 0; j < stride; j++)
          diff += fimp[c_id][isou][j]*(pvar[c_id][j] - pvara[c_id][j]);

        smbini[c_id][isou] -= diff;

        smbrp[c_id][isou] += smbini[c_id][isou];
      }
    });
  }

  ctx.wait();

  /* --- Right hand side residual */
  residu = sqrt(cs_gdot(stride*n_cells, (cs_real_t *)smbrp, (cs_real_t *)smbrp));

  /* Normalization residual
     (L2-norm of B.C. + source terms + non-orthogonality terms)

     Caution: when calling a matrix-vector product, here for a variable
     which is not "by increments" and is assumed initialized, including
     for ghost values. */

  /* Allocate a temporary array */
  var_t *w1, *w2;
  CS_MALLOC_HD(w1, n_cells_ext, var_t, amode);
  CS_MALLOC_HD(w2, n_cells_ext, var_t, amode);

  cs_real_t *pvar_i;
  CS_MALLOC_HD(pvar_i, n_cells_ext, cs_real_t, amode);

  /* Compute the L2 norm of the variable */
  for (cs_lnum_t i = 0; i < stride; i++) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      pvar_i[c_id] = pvar[c_id][i];

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar_i);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X_%d^n = %f\n", i, p_mean);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      w2[c_id][i] = (pvar[c_id][i] - p_mean);

  }
  CS_FREE_HD(pvar_i);

  cs_matrix_vector_native_multiply(symmetric,
                                   db_size,
                                   eb_size,
                                   f_id,
                                   (cs_real_t *)dam,
                                   xam,
                                   (cs_real_t *)w2,
                                   (cs_real_t *)w1);

  ctx.wait(); // matrix vector multiply uses the same stream as the ctx

  if (iwarnp >= 2) {
    const cs_real_t *_w1 = (cs_real_t *)w1, *_smbrp = (cs_real_t *)smbrp;
    bft_printf("L2 norm ||AX^n|| = %f\n",
               sqrt(cs_gdot(stride*n_cells, _w1, _w1)));
    bft_printf("L2 norm ||B^n|| = %f\n",
               sqrt(cs_gdot(stride*n_cells, _smbrp, _smbrp)));
  }

  if (has_dc == 1) {
    int *c_disable_flag = mq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
                                /* Remove contributions from penalized cells */
                                for (cs_lnum_t i = 0; i < stride; i++)
                                  w1[c_id][i] += smbrp[c_id][i];

                                if (c_disable_flag[c_id] != 0)
                                  for (cs_lnum_t i = 0; i < stride; i++)
                                    w1[c_id][i] = 0.;
                              });
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
                                for (cs_lnum_t i = 0; i < stride; i++)
                                  w1[c_id][i] += smbrp[c_id][i];
                              });
  }

  ctx.wait();

  cs_real_t rnorm2 = cs_gdot(stride*n_cells, (cs_real_t *)w1, (cs_real_t *)w1);
  cs_real_t rnorm = sqrt(rnorm2);

  CS_FREE_HD(w1);
  CS_FREE_HD(w2);

  sinfo.rhs_norm = rnorm;

  /* Warning: for Weight Matrix, one and only one sweep is done. */
  int nswmod = CS_MAX(eqp->nswrsm, 1);

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
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          dpvarm1[c_id][isou] = dpvar[c_id][isou];
          dpvar[c_id][isou] = 0.;
        }
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          dpvar[c_id][isou] = 0.;
      });
    }

    ctx.wait();

    /*  Solver residual */
    ressol = residu;

    if (conv_diff_mg)
      cs_sles_setup_native_conv_diff(f_id,
                                     var_name,
                                     db_size,
                                     eb_size,
                                     (cs_real_t *)dam,
                                     xam,
                                     true);

    cs_sles_solve_native(f_id,
                         var_name,
                         symmetric,
                         db_size,
                         eb_size,
                         (cs_real_t *)dam,
                         xam,
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
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          adxkm1[c_id][isou] = adxk[c_id][isou];
          adxk[c_id][isou] = - rhs0[c_id][isou];
        }
      });

      ctx.wait();

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

      /* ||E.dx^k-E.0||^2 */
      nadxk = cs_gdot(stride*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxk);

      /* < E.dx^k-E.0; r^k > */
      cs_real_t paxkrk = cs_gdot(stride*n_cells,
                                 (cs_real_t *)smbrp,
                                 (cs_real_t *)adxk);

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {
        /* < E.dx^(k-1)-E.0; r^k > */
        paxm1rk = cs_gdot(stride*n_cells, (cs_real_t *)smbrp, (cs_real_t *)adxkm1);

        /* < E.dx^(k-1)-E.0; E.dx^k-E.0 > */
        paxm1ax = cs_gdot(stride*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxkm1);

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
        alph = -paxkrk/cs_math_fmax(nadxk, 1.e-30*rnorm2);
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/cs_math_fmax(nadxk, 1.e-30*rnorm2);
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
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          pvar[c_id][isou] += dpvar[c_id][isou];

          /* smbini already contains unsteady terms and mass source terms
           * of the RHS updated at each sweep */

        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; j++)
            diff += fimp[c_id][isou][j]*dpvar[c_id][j];

          smbini[c_id][isou] -= diff;
          smbrp[c_id][isou] = smbini[c_id][isou];
        }

      });
    }
    else if (iswdyp == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t isou = 0; isou < stride; isou++)
          pvar[c_id][isou] += alph*dpvar[c_id][isou];

          /* smbini already contains unsteady terms and mass source terms
           * of the RHS updated at each sweep */

        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; isou++)
            diff += fimp[c_id][isou][j]*alph*dpvar[c_id][j];

          smbini[c_id][isou] -= diff;
          smbrp[c_id][isou] = smbini[c_id][isou];
        }
      });
    }
    else if (iswdyp >= 2) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          pvar[c_id][isou] +=   alph*dpvar[c_id][isou]
                              + beta*dpvarm1[c_id][isou];
        }

        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          cs_real_t diff = 0.;
          for (cs_lnum_t j = 0; j < stride; j++)
            diff += fimp[c_id][isou][j]*(  alph*dpvar[c_id][j]
                                         + beta*dpvarm1[c_id][j]);

          smbini[c_id][isou] -= diff;
          smbrp[c_id][isou] = smbini[c_id][isou];
        }
      });
    }

    /* --> Handle parallelism and periodicity */

    if (cs_glob_rank_id >= 0 || cs_glob_mesh->n_init_perio > 0) {
      if (stride == 3)
        cs_mesh_sync_var_vect((cs_real_t *)pvar);
      else if (stride == 6)
        cs_mesh_sync_var_sym_tens((cs_real_6_t *)pvar);
    }

    ctx.wait();

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
    residu = sqrt(cs_gdot(stride*n_cells,
                          (cs_real_t *)smbrp,
                          (cs_real_t *)smbrp));

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
  if (cs_math_fabs(rnorm)/sqrt(cs_real_t(stride)) > cs_math_epzero)
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

    /* smbini already contains unsteady terms and mass source terms
       of the RHS updated at each sweep */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        cs_real_t diff = 0.;
        for (cs_lnum_t jsou = 0; jsou < stride; jsou++)
          diff += fimp[c_id][isou][jsou]*dpvar[c_id][jsou];

        smbrp[c_id][isou] = smbini[c_id][isou] - diff;
      }
    });

    ctx.wait();

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
      for (cs_lnum_t isou = 0; isou < stride; isou++)
        eswork[c_id][isou] = cs_math_pow2(smbrp[c_id][isou] / cell_vol[c_id]);
    });
    ctx.wait();

  }

  /*==========================================================================
   * Free solver setup
   *==========================================================================*/

  cs_sles_free_native(f_id, var_name);

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
  CS_FREE_HD(xam);
  CS_FREE_HD(smbini);
  CS_FREE_HD(dpvar);
  if (iswdyp >= 1) {
    CS_FREE_HD(adxk);
    CS_FREE_HD(adxkm1);
    CS_FREE_HD(dpvarm1);
    CS_FREE_HD(rhs0);
  }

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
  cs_real_t relaxp = eqp->relaxv;
  cs_real_t thetap = eqp->theta;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
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
  cs_real_t *dam = nullptr, *xam = nullptr, *smbini = nullptr;
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

  CS_MALLOC_HD(dam, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xam, isym*n_i_faces, cs_real_t, cs_alloc_mode);

  cs_matrix_wrapper_scalar(iconvp,
                           idiffp,
                           ndircp,
                           isym,
                           thetap,
                           imucpp,
                           bc_coeffs,
                           rovsdt,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm,
                           xcpp,
                           dam,
                           xam);

  /* Precaution if diagonal is 0, which may happen is all surrounding cells
   * are disabled
   * If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      if (CS_ABS(dam[cell_id]) < DBL_MIN)
        dam[cell_id] += 1.;
    });
  }

  /* For steady computations, the diagonal is relaxed */
  if (idtvar < 0) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      dam[cell_id] /= relaxp;
    });
  }

  /*==========================================================================
   * 2. Iterative process to handle non orthogonalities (starting from the
   *    second iteration).
   *==========================================================================*/

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
                      smbrp);

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
                    smbrp);

  if (iswdyp >= 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      rhs0[cell_id] = smbrp[cell_id];
      smbini[cell_id] -= rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      smbrp[cell_id]  += smbini[cell_id];

      adxkm1[cell_id] = 0.;
      adxk[cell_id] = 0.;
      dpvar[cell_id] = 0.;
    });

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      smbini[cell_id] -= rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      smbrp[cell_id]  += smbini[cell_id];
    });
  }

  ctx.wait();

  /* --- Right hand side residual */
  residu = sqrt(cs_gdot(n_cells, smbrp, smbrp));

  if (normp > 0.)
    rnorm = normp;

  else {
    /* Normalization residual
       (L2-norm of B.C. + source terms + non-orthogonality terms)

       Caution: when calling a matrix-vector product, here for a variable
       which is not "by increments" and is assumed initialized, including
       for ghost values. */

    /* Allocate a temporary array */
    CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_real_t *w2;
    CS_MALLOC_HD(w2, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X^n = %f\n", p_mean);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      w2[cell_id] = (pvar[cell_id]-p_mean);
    });

    ctx.wait();

    cs_matrix_vector_native_multiply(symmetric,
                                     1,  /* db_size */
                                     1,  /* eb_size */
                                     f_id,
                                     dam,
                                     xam,
                                     w2,
                                     w1);

    CS_FREE_HD(w2);

    if (iwarnp >= 2) {
      bft_printf("L2 norm ||AX^n|| = %f\n", sqrt(cs_gdot(n_cells, w1, w1)));
      bft_printf("L2 norm ||B^n|| = %f\n", sqrt(cs_gdot(n_cells, smbrp, smbrp)));
    }
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      w1[cell_id] += smbrp[cell_id];
      /* Remove contributions from penalized cells */
      if (has_dc * c_disable_flag[has_dc * cell_id] != 0)
        w1[cell_id] = 0.;
    });

    ctx.wait();

    rnorm2 = cs_gdot(n_cells, w1, w1);
    rnorm = sqrt(rnorm2);
  }

  sinfo.rhs_norm = rnorm;

  /* Free memory */
  CS_FREE_HD(w1);

  /* Warning: for weight matrix, one and only one sweep is done. */
  int nswmod = CS_MAX(eqp->nswrsm, 1);

  /* Reconstruction loop (beginning) */
  if (iterns <= 1)
    sinfo.n_it = 0;

  int isweep = 1;

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {

    /* Solving on the increment: dpvar */

    if (iswdyp >= 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        dpvarm1[cell_id] = dpvar[cell_id];
        dpvar[cell_id] = 0.;
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        dpvar[cell_id] = 0.;
      });
    }

    /* Solver residual */
    ressol = residu;

    if (conv_diff_mg)
      cs_sles_setup_native_conv_diff(f_id,
                                     var_name,
                                     1,  /* db_size */
                                     1,  /* eb_size */
                                     dam,
                                     xam,
                                     true);

    ctx.wait(); /* We now need dpvar, computed by ctx */

    cs_sles_solve_native(f_id,
                         var_name,
                         symmetric,
                         1,  /* db_size */
                         1,  /* eb_size */
                         dam,
                         xam,
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
                        adxk);

      /* ||E.dx^(k-1)-E.0||^2 */
      nadxkm1 = nadxk;

      /* ||E.dx^k-E.0||^2 */
      nadxk = cs_gdot(n_cells, adxk, adxk);

      /* < E.dx^k-E.0; r^k > */
      paxkrk = cs_gdot(n_cells, smbrp, adxk);

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {

        /* < E.dx^(k-1)-E.0; r^k > */
        paxm1rk = cs_gdot(n_cells, smbrp, adxkm1);

        /* < E.dx^(k-1)-E.0; E.dx^k-E.0 > */
        paxm1ax = cs_gdot(n_cells, adxk, adxkm1);

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
        alph = -paxkrk/CS_MAX(nadxk, 1e-30*rnorm2);
        if (isnan(alph))
          alph = 0;
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/CS_MAX(nadxk, 1e-30*rnorm2);
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
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[cell_id] -= rovsdt[cell_id]*dpvar[cell_id];
        smbrp[cell_id]   = smbini[cell_id];
      });
    }
    else if (iswdyp == 1) {
      if (alph < 0.) break;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        pvar[cell_id] += alph*dpvar[cell_id];
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[cell_id] -= rovsdt[cell_id]*alph*dpvar[cell_id];
        smbrp[cell_id]   = smbini[cell_id];
      });
    }
    else if (iswdyp >= 2) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        pvar[cell_id] += alph*dpvar[cell_id] + beta*dpvarm1[cell_id];
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[cell_id]
          -= rovsdt[cell_id]*(alph*dpvar[cell_id]+beta*dpvarm1[cell_id]);
        smbrp[cell_id] = smbini[cell_id];
      });
    }

    ctx.wait();

    /*  ---> Handle parallelism and periodicity */
    if (cs_glob_rank_id >= 0 || m->n_init_perio > 0) {
      cs_mesh_sync_var_scal(pvar);
    }

    /* Compute the beta (min/max) limiter */
    if (f_id > -1)
      cs_beta_limiter_building(f_id, inc, rovsdt);

    /* The added convective scalar mass flux is:
       (thetex*Y_\face-imasac*Y_\celli)*mf.
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
                      smbrp);

    /* --- Convergence test */
    residu = sqrt(cs_gdot(n_cells, smbrp, smbrp));

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

  /* Rebuild final flux at faces */

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

      cs_real_t *i_flux = cs_field_by_id(i_flux_id)->val;
      cs_real_t *b_flux = cs_field_by_id(b_flux_id)->val;

      /* rebuild before-last value of variable */
      cs_real_t *prev_s_pvar;
      CS_MALLOC_HD(prev_s_pvar, n_cells_ext, cs_real_t, cs_alloc_mode);
      ctx_c.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        prev_s_pvar[cell_id] = pvar[cell_id]-dpvar[cell_id];
      });

      cs_real_2_t *i_flux2;
      CS_MALLOC_HD(i_flux2, n_i_faces, cs_real_2_t, cs_alloc_mode);
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
        i_flux2[face_id][0] = 0.;
        i_flux2[face_id][1] = 0.;
      });

      inc  = 1;
      imasac = 0; /* mass accumulation not taken into account */

      ctx.wait();

      /* If thetex = 0, no need to do more */
      if (fabs(thetex) > cs_math_epzero) {
        thetap = eqp->theta;
        eqp->theta = thetex;

        cs_face_convection_scalar(idtvar,
                                  f_id,
                                  *eqp,
                                  icvflb,
                                  inc,
                                  imasac,
                                  nullptr, /* pvar == pvara */
                                  pvara,
                                  icvfli,
                                  bc_coeffs,
                                  i_massflux,
                                  b_massflux,
                                  i_flux2,
                                  b_flux);

        eqp->theta = thetap;
      }

      ctx_c.wait();

      cs_face_convection_scalar(idtvar,
                                f_id,
                                *eqp,
                                icvflb,
                                inc,
                                imasac,
                                prev_s_pvar,
                                pvara,
                                icvfli,
                                bc_coeffs,
                                i_massflux,
                                b_massflux,
                                i_flux2,
                                b_flux);
      CS_FREE_HD(prev_s_pvar);

      /* last increment in upwind to fulfill exactly the considered
         balance equation */

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
                                i_flux2,
                                b_flux);

      /* FIXME diffusion part */

      /* Store the convectif flux,
       * Note that the two sides are equal if imasac=0 */
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        i_flux[face_id] += i_flux2[face_id][0];
      });
      ctx.wait();

      CS_FREE_HD(i_flux2);
    }
  }

  /*==========================================================================
   * 3. After having computed the new value, an estimator is computed for the
   * prediction step of the velocity.
   *==========================================================================*/

  if (iescap > 0) {
    /*  Computation of the estimator of the current component */

    /* smbini already contains unsteady terms and mass source terms
       of the RHS updated at each sweep */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      smbrp[cell_id] = smbini[cell_id] - rovsdt[cell_id]*dpvar[cell_id];
    });

    inc = 1;

    ctx.wait();

    /* Without relaxation even for a stationnary computation */

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
                      smbrp);

    /* Contribution of the current component to the L2 norm stored in eswork */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      eswork[cell_id] = cs_math_pow2(smbrp[cell_id] / cell_vol[cell_id]);
    });
    ctx.wait();

  }

  /*==========================================================================
   * 4. Free solver setup
   *==========================================================================*/

  cs_sles_free_native(f_id, var_name);

  ctx.wait();

  /*  Free memory */
  CS_FREE_HD(dam);
  CS_FREE_HD(xam);
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
                                   const cs_field_bc_coeffs_t *bc_coeffs_v,
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
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
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
                                   const cs_field_bc_coeffs_t *bc_coeffs_ts,
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
