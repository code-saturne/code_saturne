/*============================================================================
 * Convection diffusion equation solver.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

BEGIN_C_DECLS

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_equation_iterative_solve.c
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
 * \param[in]     name          associated name if f_id < 0, or NULL
 * \param[in]     iescap        compute the predictor indicator if 1
 * \param[in]     imucpp        indicator
 *                               - 0 do not multiply the convectiv term by Cp
 *                               - 1 do multiply the convectiv term by Cp
 * \param[in]     normp         Reference norm to solve the system (optional)
 *                              if negative: recomputed here
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ a^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ a^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in]     coefap        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]     coefbp        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]     cofafp        boundary condition array for the diffusion
 *                               of the variable (explicit part)
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
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
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_t       pvara[],
                                   const cs_real_t       pvark[],
                                   const cs_real_t       coefap[],
                                   const cs_real_t       coefbp[],
                                   const cs_real_t       cofafp[],
                                   const cs_real_t       cofbfp[],
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

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int iconvp = var_cal_opt->iconv;
  int idiffp = var_cal_opt->idiff;
  int iwarnp = var_cal_opt->verbosity;
  int iswdyp = var_cal_opt->iswdyn;
  int ndircp = var_cal_opt->ndircl;
  double epsrsp = var_cal_opt->epsrsm;
  double epsilp = var_cal_opt->epsilo;
  double relaxp = var_cal_opt->relaxv;
  double thetap = var_cal_opt->thetav;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  int isym, inc, isweep, niterf, nswmod;
  int lvar, imasac, key_sinfo_id;
  double residu, rnorm, ressol;
  double thetex, nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk;

  double rnorm2 = 0, alph = 0, beta = 0;

  cs_solving_info_t sinfo;

  int coupling_id = -1;
  cs_field_t *f = NULL;
  cs_real_t *dam = NULL, *xam = NULL, *smbini = NULL;
  cs_real_t *w1 = NULL;

  bool conv_diff_mg = false;

  /*============================================================================
   * 0.  Initialization
   *==========================================================================*/

  /* Name */
  const char *var_name = cs_sles_name(f_id, name);

  if (iwarnp >= 1)
    bft_printf("Equation iterative solve of: %s\n", var_name);

  /* solving info */
  key_sinfo_id = cs_field_key_id("solving_info");
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

  cs_real_t *adxk = NULL, *adxkm1 = NULL, *dpvarm1 = NULL, *rhs0 = NULL;

  if (iswdyp >= 1) {
    BFT_MALLOC(adxk, n_cells_ext, cs_real_t);
    BFT_MALLOC(adxkm1, n_cells_ext, cs_real_t);
    BFT_MALLOC(dpvarm1, n_cells_ext, cs_real_t);
    BFT_MALLOC(rhs0, n_cells_ext, cs_real_t);
  }

  /* Symmetric matrix, except if advection */
  isym = 1;
  if (iconvp > 0) isym = 2;

  bool symmetric = (isym == 1) ? true : false;

  /* Periodicity has to be taken into account */

  /* Initialization for test before matrix vector product
     for computation of initial residual */

  /*============================================================================
   * 1.  Building of the "simplified" matrix
   *==========================================================================*/

  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(xam, isym*n_i_faces, cs_real_t);

  cs_matrix_wrapper_scalar(iconvp,
                           idiffp,
                           ndircp,
                           isym,
                           thetap,
                           imucpp,
                           coefbp,
                           cofbfp,
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
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      if (CS_ABS(dam[cell_id]) < DBL_MIN) {
        dam[cell_id] += 1.;
      }
    }
  }

  /* For steady computations, the diagonal is relaxed */
  if (idtvar < 0) {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      dam[cell_id] /= relaxp;
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
    inc    = 1;

    /* The added convective scalar mass flux is:
       (thetex*Y_\face-imasac*Y_\celli)*mf.
       When building the explicit part of the rhs, one
       has to impose 0 on mass accumulation. */
    imasac = 0;

    var_cal_opt->thetav = thetex;

    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      var_cal_opt,
                      NULL, /* pvar == pvara */
                      pvara,
                      coefap,
                      coefbp,
                      cofafp,
                      cofbfp,
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

    var_cal_opt->thetav = thetap;

  }

  /* Before looping, the RHS without reconstruction is stored in smbini */

  BFT_MALLOC(smbini, n_cells_ext, cs_real_t);

  cs_lnum_t has_dc = mq->has_disable_flag;
# pragma omp parallel if(n_cells > CS_THR_MIN)
  {
#   pragma omp for nowait
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      smbini[cell_id] = smbrp[cell_id];
      smbrp[cell_id] = 0.;
    }

#   pragma omp for nowait
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
      pvar[cell_id] = pvark[cell_id];
  }

  /* In the following, cs_balance_scalar is called with inc=1,
     except for Weight Matrix (nswrsp=-1) */
  inc = 1;

  if (var_cal_opt->nswrsm == -1) {
    var_cal_opt->nswrsm = 1;
    inc = 0;
  }

  /* Incrementation and rebuild of right hand side */

  /* The added convective scalar mass flux is:
     (thetap*Y_\face-imasac*Y_\celli)*mf.
     When building the implicit part of the rhs, one
     has to impose 1 on mass accumulation. */
  imasac = 1;

  cs_balance_scalar(idtvar,
                    f_id,
                    imucpp,
                    imasac,
                    inc,
                    var_cal_opt,
                    pvar,
                    pvara,
                    coefap,
                    coefbp,
                    cofafp,
                    cofbfp,
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
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      rhs0[cell_id] = smbrp[cell_id];
      smbini[cell_id] -= rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      smbrp[cell_id]  += smbini[cell_id];

      adxkm1[cell_id] = 0.;
      adxk[cell_id] = 0.;
      dpvar[cell_id] = 0.;
    }

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      smbini[cell_id] -= rovsdt[cell_id]*(pvar[cell_id] - pvara[cell_id]);
      smbrp[cell_id]  += smbini[cell_id];
    }
  }

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
    BFT_MALLOC(w1, n_cells_ext, cs_real_t);

    cs_real_t *w2;
    BFT_MALLOC(w2, n_cells_ext, cs_real_t);

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X^n = %f\n", p_mean);
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      w2[cell_id] = (pvar[cell_id]-p_mean);

    cs_matrix_vector_native_multiply(symmetric,
                                     1,  /* db_size */
                                     1,  /* eb_size */
                                     f_id,
                                     dam,
                                     xam,
                                     w2,
                                     w1);

    BFT_FREE(w2);

    if (iwarnp >= 2) {
      bft_printf("L2 norm ||AX^n|| = %f\n", sqrt(cs_gdot(n_cells, w1, w1)));
      bft_printf("L2 norm ||B^n|| = %f\n", sqrt(cs_gdot(n_cells, smbrp, smbrp)));
    }
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      w1[cell_id] += smbrp[cell_id];
      /* Remove contributions from penalized cells */
      if (has_dc * mq->c_disable_flag[has_dc * cell_id] != 0)
        w1[cell_id] = 0.;
    }

    rnorm2 = cs_gdot(n_cells, w1, w1);
    rnorm = sqrt(rnorm2);
  }

  sinfo.rhs_norm = rnorm;

  /* Free memory */
  BFT_FREE(w1);

  /* Warning: for weight matrix, one and only one sweep is done. */
  nswmod = CS_MAX(var_cal_opt->nswrsm, 1);

  /* Reconstruction loop (beginning) */
  if (iterns <= 1)
    sinfo.n_it = 0;
  isweep = 1;

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {

    /* Solving on the increment: dpvar */

    if (iswdyp >= 1) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        dpvarm1[cell_id] = dpvar[cell_id];
        dpvar[cell_id] = 0.;
      }
    }
    else {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        dpvar[cell_id] = 0.;
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

#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        adxkm1[cell_id] = adxk[cell_id];
        adxk[cell_id] = - rhs0[cell_id];
      }

      cs_balance_scalar(idtvar,
                        lvar,
                        imucpp,
                        imasac,
                        inc,
                        var_cal_opt,
                        dpvar,
                        NULL, /* dpvara == dpvar */
                        coefap,
                        coefbp,
                        cofafp,
                        cofbfp,
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

    /* --- Update the solution with the increment */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        pvar[cell_id] += dpvar[cell_id];
    }
    else if (iswdyp == 1) {
      if (alph < 0.) break;
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        pvar[cell_id] += alph*dpvar[cell_id];
    }
    else if (iswdyp >= 2) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        pvar[cell_id] += alph*dpvar[cell_id] + beta*dpvarm1[cell_id];
    }

    /*  ---> Handle parallelism and periodicity */
    if (cs_glob_rank_id >= 0 || cs_glob_mesh->n_init_perio > 0) {
      cs_mesh_sync_var_scal(pvar);
    }

    /* --- Update the right hand side And compute the new residual */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[iel] -= rovsdt[iel]*dpvar[iel];
        smbrp[iel]   = smbini[iel];
      }
    }
    else if (iswdyp == 1) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[iel] -= rovsdt[iel]*alph*dpvar[iel];
        smbrp[iel]   = smbini[iel];
      }
    }
    else if (iswdyp >= 2) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
           of the RHS updated at each sweep */
        smbini[iel] -= rovsdt[iel]*(alph*dpvar[iel]+beta*dpvarm1[iel]);
        smbrp[iel]   = smbini[iel];
      }
    }

    /* Compute the beta (min/max) limiter */
    if (f_id > -1)
      cs_beta_limiter_building(f_id, inc, rovsdt);

    /* The added convective scalar mass flux is:
       (thetex*Y_\face-imasac*Y_\celli)*mf.
       When building the implicit part of the rhs, one
       has to impose 1 on mass accumulation. */
    imasac = 1;

    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      var_cal_opt,
                      pvar,
                      pvara,
                      coefap,
                      coefbp,
                      cofafp,
                      cofbfp,
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

      cs_field_t *i_flux = cs_field_by_id(i_flux_id);
      cs_field_t *b_flux = cs_field_by_id(b_flux_id);

      /* rebuild before-last value of variable */
      cs_real_t *prev_s_pvar;
      BFT_MALLOC(prev_s_pvar, n_cells_ext, cs_real_t);
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        prev_s_pvar[iel] = pvar[iel]-dpvar[iel];
      }

      cs_real_2_t *i_flux2;
      BFT_MALLOC(i_flux2, n_i_faces, cs_real_2_t);
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        i_flux2[face_id][0] = 0.;
        i_flux2[face_id][1] = 0.;
      }

      inc  = 1;
      imasac = 0; /* mass accumluation not taken into account */

      /* If thetex = 0, no need to do more */
      if (fabs(thetex) > cs_math_epzero) {
        thetap = var_cal_opt->thetav;
        var_cal_opt->thetav = thetex;

        cs_face_convection_scalar(idtvar,
                                  f_id,
                                  *var_cal_opt,
                                  icvflb,
                                  inc,
                                  imasac,
                                  NULL, /* pvar == pvara */
                                  pvara,
                                  icvfli,
                                  coefap,
                                  coefbp,
                                  i_massflux,
                                  b_massflux,
                                  i_flux2,
                                  b_flux->val);

        var_cal_opt->thetav = thetap;
      }

      cs_face_convection_scalar(idtvar,
                                f_id,
                                *var_cal_opt,
                                icvflb,
                                inc,
                                imasac,
                                prev_s_pvar,
                                pvara,
                                icvfli,
                                coefap,
                                coefbp,
                                i_massflux,
                                b_massflux,
                                i_flux2,
                                b_flux->val);
      BFT_FREE(prev_s_pvar);

      /* last increment in upwind to fulfill exactly the considered
         balance equation */

      inc  = 0;

      cs_var_cal_opt_t var_cal_opt_loc;
      int k_id = cs_field_key_id("var_cal_opt");
      cs_field_get_key_struct(f, k_id, &var_cal_opt_loc);

      var_cal_opt_loc.nswrgr = 0;
      var_cal_opt_loc.blencv = 0.;

      cs_face_convection_scalar(idtvar,
                                f_id,
                                var_cal_opt_loc,
                                icvflb,
                                inc,
                                imasac,
                                dpvar,
                                pvara,
                                icvfli,
                                coefap,
                                coefbp,
                                i_massflux,
                                b_massflux,
                                i_flux2,
                                b_flux->val);

      /* FIXME diffusion part */

      /* Store the convectif flux,
       * Note that the two sides are equal if imasac=0 */
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
        i_flux->val[face_id] += i_flux2[face_id][0];
      BFT_FREE(i_flux2);
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

#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      smbrp[iel] = smbini[iel] - rovsdt[iel]*dpvar[iel];

    inc    = 1;

    /* Without relaxation even for a stationnary computation */

    cs_balance_scalar(idtvar,
                      f_id,
                      imucpp,
                      imasac,
                      inc,
                      var_cal_opt,
                      pvar,
                      pvara,
                      coefap,
                      coefbp,
                      cofafp,
                      cofbfp,
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

#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      eswork[iel] = pow(smbrp[iel] / cell_vol[iel],2);
  }

  /*==========================================================================
   * 4. Free solver setup
   *==========================================================================*/

  cs_sles_free_native(f_id, var_name);

  /*  Free memory */
  BFT_FREE(dam);
  BFT_FREE(xam);

  BFT_FREE(smbini);
  if (iswdyp >= 1) {
    BFT_FREE(adxk);
    BFT_FREE(adxkm1);
    BFT_FREE(dpvarm1);
    BFT_FREE(rhs0);
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
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]      iescap        compute the predictor indicator if >= 1
 * \param[in]      var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]      pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]      pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               \c pvara (usually \c pvar= \c pvara)
 * \param[in]      coefav        boundary condition array for the variable
 *                               (explicit part)
 * \param[in]      coefbv        boundary condition array for the variable
 *                               (implicit part)
 * \param[in]      cofafv        boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]      cofbfv        boundary condition array for the diffusion
 *                               of the variable (Implicit part)
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
 * \param[in]      smbrp         Right hand side \f$ \vect{Rhs}^k \f$
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
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_3_t     pvara[],
                                   const cs_real_3_t     pvark[],
                                   const cs_real_3_t     coefav[],
                                   const cs_real_33_t    coefbv[],
                                   const cs_real_3_t     cofafv[],
                                   const cs_real_33_t    cofbfv[],
                                   const cs_real_t       i_massflux[],
                                   const cs_real_t       b_massflux[],
                                   cs_real_t             i_viscm[],
                                   const cs_real_t       b_viscm[],
                                   const cs_real_t       i_visc[],
                                   const cs_real_t       b_visc[],
                                   const cs_real_t       i_secvis[],
                                   const cs_real_t       b_secvis[],
                                   cs_real_6_t           viscel[],
                                   const cs_real_2_t     weighf[],
                                   const cs_real_t       weighb[],
                                   int                   icvflb,
                                   const int             icvfli[],
                                   cs_real_33_t          fimp[],
                                   cs_real_3_t           smbrp[],
                                   cs_real_3_t           pvar[],
                                   cs_real_3_t           eswork[])
{
  /* Local variables */

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int iconvp = var_cal_opt->iconv;
  int idiffp = var_cal_opt->idiff;
  int iwarnp = var_cal_opt->verbosity;
  int iswdyp = var_cal_opt->iswdyn;
  int idftnp = var_cal_opt->idften;
  int ndircp = var_cal_opt->ndircl;
  double epsrsp = var_cal_opt->epsrsm;
  double epsilp = var_cal_opt->epsilo;
  double relaxp = var_cal_opt->relaxv;
  double thetap = var_cal_opt->thetav;

  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  int isym, inc, isweep, niterf, nswmod;
  int key_sinfo_id;
  int lvar, imasac;
  double residu, rnorm, ressol, thetex;
  double paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax;

  double alph = 0., beta = 0.;

  cs_solving_info_t sinfo;

  cs_field_t *f = NULL;;

  cs_real_t    *xam = NULL;
  cs_real_33_t *dam = NULL;;
  cs_real_3_t  *dpvar = NULL, *smbini = NULL, *w1 = NULL, *w2 = NULL;
  cs_real_3_t  *adxk = NULL, *adxkm1 = NULL, *dpvarm1 = NULL, *rhs0 = NULL;

  /*============================================================================
   * 0.  Initialization
   *==========================================================================*/

  /* Name */
  const char *var_name = cs_sles_name(f_id, name);

  if (iwarnp >= 1)
    bft_printf("Equation iterative solve of: %s\n", var_name);

  /* Matrix block size */
  cs_lnum_t db_size = 3;
  cs_lnum_t eb_size = 1; /* CS_ISOTROPIC_DIFFUSION
                            or CS_ANISOTROPIC_RIGHT_DIFFUSION */
  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION)
    eb_size = 3;

  if (cs_glob_porous_model == 3) { //FIXME iphydr + other option?
    if (iconvp > 0)
      eb_size = 3;
  }

  /* Allocate temporary arrays */
  BFT_MALLOC(dam, n_cells_ext, cs_real_33_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(smbini, n_cells_ext, cs_real_3_t);

  if (iswdyp >= 1) {
    BFT_MALLOC(adxk, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(adxkm1, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(dpvarm1, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(rhs0, n_cells_ext, cs_real_3_t);
  }

  cs_real_3_t *i_pvar = NULL;
  cs_real_3_t *b_pvar = NULL;
  cs_field_t *i_vf = NULL;
  cs_field_t *b_vf = NULL;

  /* Storing face values for kinetic energy balance and initialize them */
  if (CS_F_(vel)->id == f_id) {

    i_vf = cs_field_by_name_try("inner_face_velocity");
    if (i_vf != NULL) {
      cs_array_real_fill_zero(3*n_i_faces, i_vf->val);
      i_pvar = (cs_real_3_t *)i_vf->val;
    }

    b_vf = cs_field_by_name_try("boundary_face_velocity");
    if (b_vf != NULL) {
      cs_array_real_fill_zero(3*n_b_faces, b_vf->val);
      b_pvar = (cs_real_3_t *)b_vf->val;
    }

  }

  /* solving info */
  key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);
  }

  /* Symmetric matrix, except if advection */
  isym = 1;
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
  BFT_MALLOC(xam, eb_stride*isym*n_faces, cs_real_t);

  /*==========================================================================
   * 1.  Building of the "simplified" matrix
   *==========================================================================*/

  int tensorial_diffusion = 1;

  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION)
    tensorial_diffusion = 2;

  cs_matrix_wrapper_vector(iconvp,
                           idiffp,
                           tensorial_diffusion,
                           ndircp,
                           isym,
                           eb_size,
                           thetap,
                           coefbv,
                           cofbfv,
                           (const cs_real_33_t *)fimp,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm,
                           dam,
                           xam);

  /* Precaution if diagonal is 0, which may happen is all surrounding cells
   * are disabled
   * If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        if (CS_ABS(dam[cell_id][i][i]) < DBL_MIN) {
          dam[cell_id][i][i] += 1.;
        }
      }
    }
  }

  /*  For steady computations, the diagonal is relaxed */
  if (idtvar < 0) {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          dam[iel][isou][jsou] /= relaxp;
      }
    }
  }

  /*===========================================================================
   * 2. Iterative process to handle non orthogonlaities (starting from the
   * second iteration).
   *===========================================================================*/

  /* Application du theta schema */

  /* On calcule le bilan explicite total */
  thetex = 1. - thetap;

  /* Si THETEX=0, ce n'est pas la peine d'en rajouter */
  if (fabs(thetex) > cs_math_epzero) {
    inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
    imasac = 0;

    var_cal_opt->thetav = thetex;

    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,
                      ivisep,
                      var_cal_opt,
                      NULL, /* pvar == pvara */
                      pvara,
                      coefav,
                      coefbv,
                      cofafv,
                      cofbfv,
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
                      i_pvar,
                      b_pvar,
                      smbrp);

    /* Save (1-theta)* face_value at previous time step if needed */
    if (i_vf != NULL && b_vf != NULL) {
      cs_field_current_to_previous(i_vf);
      cs_field_current_to_previous(b_vf);
    }

    var_cal_opt->thetav = thetap;
  }

  /* Before looping, the RHS without reconstruction is stored in smbini */

  cs_lnum_t has_dc = mq->has_disable_flag;
# pragma omp parallel if(n_cells > CS_THR_MIN)
  {
#   pragma omp for nowait
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        smbini[cell_id][isou] = smbrp[cell_id][isou];
        smbrp[cell_id][isou] = 0.;
      }
    }

    /* pvar is initialized on n_cells_ext to avoid a synchronization */
#   pragma omp for nowait
    for (cs_lnum_t iel = 0; iel < n_cells_ext; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        pvar[iel][isou] = pvark[iel][isou];
    }
  }

  /* In the following, cs_balance_vector is called with inc=1,
   * except for Weight Matrix (nswrsp=-1) */
  inc = 1;

  if (var_cal_opt->nswrsm == -1) {
    var_cal_opt->nswrsm = 1;
    inc = 0;
  }

  /*  ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE */

  /*  On est entre avec un smb explicite base sur PVARA.
   *  si on initialise avec PVAR avec autre chose que PVARA
   *  on doit donc corriger SMBR (c'est le cas lorsqu'on itere sur navsto) */

  /* The added convective scalar mass flux is:
   *      (thetap*Y_\face-imasac*Y_\celli)*mf.
   * When building the implicit part of the rhs, one
   * has to impose 1 on mass accumulation. */
  imasac = 1;

  cs_balance_vector(idtvar,
                    f_id,
                    imasac,
                    inc,
                    ivisep,
                    var_cal_opt,
                    pvar,
                    pvara,
                    coefav,
                    coefbv,
                    cofafv,
                    cofbfv,
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
                    NULL,
                    NULL,
                    smbrp);

  if (CS_F_(vel)->id == f_id) {
    f = cs_field_by_name_try("velocity_explicit_balance");

    if (f != NULL) {
      cs_real_3_t *cpro_cv_df_v = (cs_real_3_t *)f->val;
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          cpro_cv_df_v[iel][isou] = smbrp[iel][isou];
      }
    }
  }

  /* Dynamic relaxation*/
  if (iswdyp >= 1) {
#   pragma omp parallel for  if(n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        rhs0[iel][isou] = smbrp[iel][isou];
        smbini[iel][isou] = smbini[iel][isou]
                          -fimp[iel][isou][0]*(pvar[iel][0] - pvara[iel][0])
                          -fimp[iel][isou][1]*(pvar[iel][1] - pvara[iel][1])
                          -fimp[iel][isou][2]*(pvar[iel][2] - pvara[iel][2]);
        smbrp[iel][isou] += smbini[iel][isou];

        adxkm1[iel][isou] = 0.;
        adxk[iel][isou] = 0.;
        dpvar[iel][isou] = 0.;
      }
    }

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
#   pragma omp parallel for  if(n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        smbini[iel][isou] = smbini[iel][isou]
                          -fimp[iel][isou][0]*(pvar[iel][0] - pvara[iel][0])
                          -fimp[iel][isou][1]*(pvar[iel][1] - pvara[iel][1])
                          -fimp[iel][isou][2]*(pvar[iel][2] - pvara[iel][2]);
        smbrp[iel][isou] += smbini[iel][isou];
      }
    }
  }

  /* --- Convergence test */
  residu = sqrt(cs_gdot(3*n_cells, (cs_real_t *)smbrp, (cs_real_t *)smbrp));

  /* ---> RESIDU DE NORMALISATION
   *    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE) */

  /* Allocate a temporary array */
  BFT_MALLOC(w1, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_3_t);

  cs_real_t *pvar_i;
  BFT_MALLOC(pvar_i, n_cells_ext, cs_real_t);

  /* Compute the L2 norm of the variable */
  for (cs_lnum_t i = 0; i < 3; i++) {

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      pvar_i[cell_id] = pvar[cell_id][i];

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar_i);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X_%d^n = %f\n", i, p_mean);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      w2[cell_id][i] = (pvar[cell_id][i] - p_mean);

  }
  BFT_FREE(pvar_i);

  cs_matrix_vector_native_multiply(symmetric,
                                   db_size,
                                   eb_size,
                                   f_id,
                                   (cs_real_t *)dam,
                                   xam,
                                   (cs_real_t *)w2,
                                   (cs_real_t *)w1);

  BFT_FREE(w2);

  if (iwarnp >= 2) {
    const cs_real_t *_w1 = (cs_real_t *)w1, *_smbrp = (cs_real_t *)smbrp;
    bft_printf("L2 norm ||AX^n|| = %f\n", sqrt(cs_gdot(3*n_cells, _w1, _w1)));
    bft_printf("L2 norm ||B^n|| = %f\n",
               sqrt(cs_gdot(3*n_cells, _smbrp, _smbrp)));
  }

# pragma omp parallel for
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (int i = 0; i < 3; i++)
      w1[cell_id][i] += smbrp[cell_id][i];
    /* Remove contributions from penalized cells */
    if (has_dc * mq->c_disable_flag[has_dc * cell_id] != 0)
      for (int i = 0; i < 3; i++)
        w1[cell_id][i] = 0.;
  }

  double  rnorm2 = cs_gdot(3*n_cells, (cs_real_t *)w1, (cs_real_t *)w1);
  rnorm = sqrt(rnorm2);
  sinfo.rhs_norm = rnorm;

  /* Free memory */
  BFT_FREE(w1);

  /* Warning: for Weight Matrix, one and only one sweep is done. */
  nswmod = CS_MAX(var_cal_opt->nswrsm, 1);

  isweep = 1;

  /* Reconstruction loop (beginning)
   *-------------------------------- */
  if (iterns <= 1)
    sinfo.n_it = 0;

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {
    /* --- Solving on the increment dpvar */

    /*  Dynamic relaxation of the system */
    if (iswdyp >= 1) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          dpvarm1[iel][isou] = dpvar[iel][isou];
          dpvar[iel][isou] = 0.;
        }
      }
    }
    else {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          dpvar[iel][isou] = 0.;
      }
    }

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

#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          adxkm1[iel][isou] = adxk[iel][isou];
          adxk[iel][isou] = - rhs0[iel][isou];
        }
      }

      cs_balance_vector(idtvar,
                        lvar,
                        imasac,
                        inc,
                        ivisep,
                        var_cal_opt,
                        dpvar,
                        NULL, /* dpvar */
                        coefav,
                        coefbv,
                        cofafv,
                        cofbfv,
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
                        NULL,
                        NULL,
                        adxk);

      /* ||E.dx^(k-1)-E.0||^2 */
      nadxkm1 = nadxk;

      /* ||E.dx^k-E.0||^2 */
      nadxk = cs_gdot(3*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxk);

      /* < E.dx^k-E.0; r^k > */
      paxkrk = cs_gdot(3*n_cells, (cs_real_t *)smbrp, (cs_real_t *)adxk);

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {
        /* < E.dx^(k-1)-E.0; r^k > */
        paxm1rk = cs_gdot(3*n_cells, (cs_real_t *)smbrp, (cs_real_t *)adxkm1);

        /* < E.dx^(k-1)-E.0; E.dx^k-E.0 > */
        paxm1ax = cs_gdot(3*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxkm1);

        if (   (nadxkm1 > 1.e-30*rnorm2)
            && (nadxk*nadxkm1-pow(paxm1ax,2)) > 1.e-30*rnorm2)
          beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-pow(paxm1ax,2));
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
        alph = -paxkrk/CS_MAX(nadxk, 1.e-30*rnorm2);
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/CS_MAX(nadxk, 1.e-30*rnorm2);
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

    /* --- Update the solution with the increment */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          pvar[iel][isou] += dpvar[iel][isou];
      }
    }
    else if (iswdyp == 1) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          pvar[iel][isou] += alph*dpvar[iel][isou];
      }
    }
    else if (iswdyp >= 2) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          pvar[iel][isou] += alph*dpvar[iel][isou] + beta*dpvarm1[iel][isou];
      }
    }

    /* --> Handle parallelism and periodicity */

    if (cs_glob_rank_id  >=0 || cs_glob_mesh->n_init_perio > 0)
      cs_mesh_sync_var_vect((cs_real_t *)pvar);

    /* --- Update the right hand and compute the new residual */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*dpvar[iel][0]
                            - fimp[iel][isou][1]*dpvar[iel][1]
                            - fimp[iel][isou][2]*dpvar[iel][2];
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }
    else if (iswdyp == 1) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*alph*dpvar[iel][0]
                            - fimp[iel][isou][1]*alph*dpvar[iel][1]
                            - fimp[iel][isou][2]*alph*dpvar[iel][2];
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }
    else if (iswdyp == 2) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*(alph*dpvar[iel][0]
                                                + beta*dpvarm1[iel][0])
                            - fimp[iel][isou][1]*(alph*dpvar[iel][1]
                                                + beta*dpvarm1[iel][1])
                            - fimp[iel][isou][2]*(alph*dpvar[iel][2]
                                                + beta*dpvarm1[iel][2]);
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }

    /* Increment face value with theta * face_value at current time step
     * if needed
     * Reinit the previous value before */
    if (i_vf != NULL && b_vf != NULL) {
      cs_array_real_copy(3 * n_i_faces, i_vf->val_pre, i_vf->val);
      cs_array_real_copy(3 * n_b_faces, b_vf->val_pre, b_vf->val);
    }

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the implicit part of the rhs, one
     * has to impose 1 on mass accumulation. */
    imasac = 1;

    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,
                      ivisep,
                      var_cal_opt,
                      pvar,
                      pvara,
                      coefav,
                      coefbv,
                      cofafv,
                      cofbfv,
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
                      i_pvar,
                      b_pvar,
                      smbrp);

    /* --- Convergence test */
    residu = sqrt(cs_gdot(3*n_cells, (cs_real_t *)smbrp, (cs_real_t *)smbrp));

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
  if (fabs(rnorm)/sqrt(3.) > cs_math_epzero)
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
   * 3. After having computed the new value, an estimator is computed for the
   * prediction step of the velocity.
   *==========================================================================*/

  if (iescap > 0) {
    /* ---> Computation of the estimator of the current component */

    /* smbini already contains unsteady terms and mass source terms
       of the RHS updated at each sweep */

#   pragma omp parallel for  if(n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        smbrp[iel][isou] = smbini[iel][isou] - fimp[iel][isou][0]*dpvar[iel][0]
                                             - fimp[iel][isou][1]*dpvar[iel][1]
                                             - fimp[iel][isou][2]*dpvar[iel][2];
    }

    inc = 1;

    /* Without relaxation even for a stationnary computation */

    cs_balance_vector(idtvar,
                      f_id,
                      imasac,
                      inc,
                      ivisep,
                      var_cal_opt,
                      pvar,
                      pvara,
                      coefav,
                      coefbv,
                      cofafv,
                      cofbfv,
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
                      NULL,
                      NULL,
                      smbrp);

    /* Contribution of the current component to the L2 norm stored in eswork */

#   pragma omp parallel for  if(n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        eswork[iel][isou] = pow(smbrp[iel][isou] / cell_vol[iel],2);
    }
  }

  /*==========================================================================
   * 4. Free solver setup
   *==========================================================================*/

  cs_sles_free_native(f_id, var_name);

  /* Save diagonal in case we want to use it */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
        fimp[cell_id][i][j] = dam[cell_id][i][j];
  }

  /* Free memory */
  BFT_FREE(dam);
  BFT_FREE(xam);
  BFT_FREE(smbini);
  BFT_FREE(dpvar);
  if (iswdyp >= 1) {
    BFT_FREE(adxk);
    BFT_FREE(adxkm1);
    BFT_FREE(dpvarm1);
    BFT_FREE(rhs0);
  }
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
 * \param[in]     name          associated name if f_id < 0, or NULL
 * \param[in]     var_cal_opt   pointer to a cs_var_cal_opt_t structure which
 *                              contains variable calculation options
 * \param[in]     pvara         variable at the previous time step
 *                               \f$ \vect{a}^n \f$
 * \param[in]     pvark         variable at the previous sub-iteration
 *                               \f$ \vect{a}^k \f$.
 *                               If you sub-iter on Navier-Stokes, then
 *                               it allows to initialize by something else than
 *                               pvara (usually pvar=pvara)
 * \param[in]     coefats       boundary condition array for the variable
 *                               (Explicit part)
 * \param[in]     coefbts       boundary condition array for the variable
 *                               (Impplicit part)
 * \param[in]     cofafts       boundary condition array for the diffusion
 *                               of the variable (Explicit part)
 * \param[in]     cofbfts       boundary condition array for the diffusion
 *                               of the variable (Implicit part)
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
 * \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
 * \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
 * \param[in,out] pvar          current variable
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_iterative_solve_tensor(int                   idtvar,
                                   int                   f_id,
                                   const char           *name,
                                   cs_var_cal_opt_t     *var_cal_opt,
                                   const cs_real_6_t     pvara[],
                                   const cs_real_6_t     pvark[],
                                   const cs_real_6_t     coefats[],
                                   const cs_real_66_t    coefbts[],
                                   const cs_real_6_t     cofafts[],
                                   const cs_real_66_t    cofbfts[],
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
                                   const cs_real_66_t    fimp[],
                                   cs_real_6_t           smbrp[],
                                   cs_real_6_t           pvar[])
{
  /* Local variables */

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int iconvp = var_cal_opt->iconv;
  int idiffp = var_cal_opt->idiff;
  int iwarnp = var_cal_opt->verbosity;
  int iswdyp = var_cal_opt->iswdyn;
  int idftnp = var_cal_opt->idften;
  int ndircp = var_cal_opt->ndircl;
  double epsrsp = var_cal_opt->epsrsm;
  double epsilp = var_cal_opt->epsilo;
  double relaxp = var_cal_opt->relaxv;
  double thetap = var_cal_opt->thetav;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  int isym, inc, isweep, niterf, nswmod;
  int key_sinfo_id;
  int lvar, imasac;
  double residu, rnorm, ressol, thetex;
  double paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax;

  double rnorm2 = 0, alph = 0, beta = 0;

  cs_solving_info_t sinfo;

  cs_field_t *f = NULL;

  cs_real_t    *xam = NULL;
  cs_real_66_t *dam = NULL;
  cs_real_6_t  *dpvar = NULL, *smbini = NULL, *w1 = NULL, *w2 = NULL;
  cs_real_6_t  *adxk = NULL, *adxkm1 = NULL, *dpvarm1 = NULL, *rhs0 = NULL;

  /*==========================================================================
   * 0.  Initialization
   *==========================================================================*/

  /* Name */
  const char *var_name = cs_sles_name(f_id, name);

  if (iwarnp >= 1)
    bft_printf("Equation iterative solve of: %s\n", var_name);

  /* Matrix block size */
  cs_lnum_t db_size = 6;
  cs_lnum_t eb_size = 1; /* CS_ISOTROPIC_DIFFUSION
                            or CS_ANISOTROPIC_RIGHT_DIFFUSION */
  if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) eb_size = 6;

  /* Allocate temporary arrays */
  BFT_MALLOC(dam, n_cells_ext, cs_real_66_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(smbini, n_cells_ext, cs_real_6_t);

  if (iswdyp >= 1) {
    BFT_MALLOC(adxk, n_cells_ext, cs_real_6_t);
    BFT_MALLOC(adxkm1, n_cells_ext, cs_real_6_t);
    BFT_MALLOC(dpvarm1, n_cells_ext, cs_real_6_t);
    BFT_MALLOC(rhs0, n_cells_ext, cs_real_6_t);
  }

  /* solving info */
  key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);
  }

  /* Determine if we are in a case with special requirements */

  bool conv_diff_mg = false;
  if (iconvp > 0) {
    cs_sles_t *sc = cs_sles_find_or_add(f_id, name);
    const char *sles_type = cs_sles_get_type(sc);
    if (strcmp(sles_type, "cs_multigrid_t") == 0)
      conv_diff_mg = true;
  }

  /* Symmetric matrix, except if advection */
  isym = 1;
  if (iconvp > 0) isym = 2;

  bool symmetric = (isym == 1) ? true : false;

  /*  be carefull here, xam is interleaved*/
  cs_lnum_t eb_stride = eb_size*eb_size;
  BFT_MALLOC(xam, eb_stride*isym*n_faces, cs_real_t);

  /*==========================================================================
   * 1.  Building of the "simplified" matrix
   *==========================================================================*/

  int tensorial_diffusion = 1;
  if (eb_size == 6)
    tensorial_diffusion = 2;

  cs_matrix_wrapper_tensor(iconvp,
                           idiffp,
                           tensorial_diffusion,
                           ndircp,
                           isym,
                           thetap,
                           coefbts,
                           cofbfts,
                           fimp,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm,
                           dam,
                           xam);

  /* Precaution if diagonal is 0, which may happen is all surrounding cells
   * are disabled
   * If a whole line of the matrix is 0, the diagonal is set to 1 */
  if (mq->has_disable_flag == 1) {
# pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t i = 0; i < 6; i++) {
        if (CS_ABS(dam[cell_id][i][i]) < DBL_MIN) {
          dam[cell_id][i][i] += 1.;
        }
      }
    }
  }

  /*  For steady computations, the diagonal is relaxed */
  if (idtvar < 0) {
#   pragma omp parallel for  if(n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 6; jsou++)
          dam[iel][isou][jsou] /= relaxp;
      }
    }
  }

  /*===========================================================================
   * 2. Iterative process to handle non orthogonlaities (starting from the
   * second iteration).
   *===========================================================================*/

  /* Application du theta schema */

  /* On calcule le bilan explicite total */
  thetex = 1. - thetap;

  /* Si THETEX=0, ce n'est pas la peine d'en rajouter */
  if (fabs(thetex) > cs_math_epzero) {
    inc = 1;

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the explicit part of the rhs, one
     * has to impose 0 on mass accumulation. */
    imasac = 0;

    var_cal_opt->thetav = thetex;

    cs_balance_tensor(idtvar,
                      f_id,
                      imasac,
                      inc,
                      var_cal_opt,
                      NULL, /* pvar == pvara */
                      pvara,
                      coefats,
                      coefbts,
                      cofafts,
                      cofbfts,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      smbrp);

    var_cal_opt->thetav = thetap;
  }

  /* Before looping, the RHS without reconstruction is stored in smbini */

  cs_lnum_t has_dc = mq->has_disable_flag;
# pragma omp parallel  if(n_cells > CS_THR_MIN)
  {
#   pragma omp for nowait
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        smbini[cell_id][isou] = smbrp[cell_id][isou];
        smbrp[cell_id][isou] = 0.;
      }
    }

  /* pvar is initialized on n_cells_ext to avoid a synchronization */
#   pragma omp for nowait
    for (cs_lnum_t iel = 0; iel < n_cells_ext; iel++) {
      for (cs_lnum_t isou = 0; isou < 6; isou++)
        pvar[iel][isou] = pvark[iel][isou];
    }
  }

  /* In the following, cs_balance_vector is called with inc=1,
   * except for Weight Matrix (nswrsp=-1) */
  inc = 1;

  if (var_cal_opt->nswrsm == -1) {
    var_cal_opt->nswrsm = 1;
    inc = 0;
  }

  /*  ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE */

  /*  On est entre avec un smb explicite base sur PVARA.
   *  si on initialise avec PVAR avec autre chose que PVARA
   *  on doit donc corriger SMBR (c'est le cas lorsqu'on itere sur navsto) */

  /* The added convective scalar mass flux is:
   *      (thetap*Y_\face-imasac*Y_\celli)*mf.
   * When building the implicit part of the rhs, one
   * has to impose 1 on mass accumulation. */
  imasac = 1;

  cs_balance_tensor(idtvar,
                    f_id,
                    imasac,
                    inc,
                    var_cal_opt,
                    pvar,
                    pvara,
                    coefats,
                    coefbts,
                    cofafts,
                    cofbfts,
                    i_massflux,
                    b_massflux,
                    i_visc,
                    b_visc,
                    viscel,
                    weighf,
                    weighb,
                    icvflb,
                    icvfli,
                    smbrp);

  /* Dynamic relaxation*/
  if (iswdyp >= 1) {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        rhs0[iel][isou] = smbrp[iel][isou];
        smbini[iel][isou] = smbini[iel][isou]
                          -fimp[iel][isou][0]*(pvar[iel][0] - pvara[iel][0])
                          -fimp[iel][isou][1]*(pvar[iel][1] - pvara[iel][1])
                          -fimp[iel][isou][2]*(pvar[iel][2] - pvara[iel][2])
                          -fimp[iel][isou][3]*(pvar[iel][3] - pvara[iel][3])
                          -fimp[iel][isou][4]*(pvar[iel][4] - pvara[iel][4])
                          -fimp[iel][isou][5]*(pvar[iel][5] - pvara[iel][5]);
        smbrp[iel][isou] += smbini[iel][isou];

        adxkm1[iel][isou] = 0.;
        adxk[iel][isou] = 0.;
        dpvar[iel][isou] = 0.;
    }

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
  }
  else {
#   pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        smbini[iel][isou] =   smbini[iel][isou]
                            - fimp[iel][isou][0]*(pvar[iel][0] - pvara[iel][0])
                            - fimp[iel][isou][1]*(pvar[iel][1] - pvara[iel][1])
                            - fimp[iel][isou][2]*(pvar[iel][2] - pvara[iel][2])
                            - fimp[iel][isou][3]*(pvar[iel][3] - pvara[iel][3])
                            - fimp[iel][isou][4]*(pvar[iel][4] - pvara[iel][4])
                            - fimp[iel][isou][5]*(pvar[iel][5] - pvara[iel][5]);
        smbrp[iel][isou] += smbini[iel][isou];
      }
    }
  }

  /* --- Convergence test */
  residu = sqrt(cs_gdot(6*n_cells, (cs_real_t *)smbrp, (cs_real_t *)smbrp));

  /* ---> RESIDU DE NORMALISATION
   *    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE) */

  /* Allocate a temporary array */
  BFT_MALLOC(w1, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);

  cs_real_t *pvar_i;
  BFT_MALLOC(pvar_i, n_cells_ext, cs_real_t);

  /* Compute the L2 norm of the variable */
  for (cs_lnum_t i = 0; i < 6; i++) {

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      pvar_i[cell_id] = pvar[cell_id][i];

    cs_real_t p_mean = cs_gmean(n_cells, mq->cell_vol, pvar_i);

    if (iwarnp >= 2)
      bft_printf("Spatial average of X_%d^n = %e\n", i, p_mean);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      w2[cell_id][i] = (pvar[cell_id][i] - p_mean);

  }
  BFT_FREE(pvar_i);

  cs_matrix_vector_native_multiply(symmetric,
                                   db_size,
                                   eb_size,
                                   f_id,
                                   (cs_real_t *)dam,
                                   xam,
                                   (cs_real_t *)w2,
                                   (cs_real_t *)w1);

  BFT_FREE(w2);

  if (iwarnp >= 2) {
    const cs_real_t *_w1 = (cs_real_t *)w1, *_smbrp = (cs_real_t *)smbrp;
    bft_printf("L2 norm ||AX^n|| = %e\n", sqrt(cs_gdot(6*n_cells, _w1, _w1)));
    bft_printf("L2 norm ||B^n|| = %e\n",
               sqrt(cs_gdot(6*n_cells, _smbrp, _smbrp)));
  }


# pragma omp parallel for  if (n_cells > CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    for (cs_lnum_t i = 0; i < 6; i++)
      w1[cell_id][i] += smbrp[cell_id][i];
    /* Remove contributions from penalized cells */
    if (has_dc * mq->c_disable_flag[has_dc * cell_id] != 0) {
      for (cs_lnum_t i = 0; i < 6; i++)
        w1[cell_id][i] = 0.;
    }
  }

  rnorm2 = cs_gdot(6*n_cells, (cs_real_t *)w1, (cs_real_t *)w1);
  rnorm = sqrt(rnorm2);
  sinfo.rhs_norm = rnorm;

  /* Free memory */
  BFT_FREE(w1);

  /* Warning: for Weight Matrix, one and only one sweep is done. */
  nswmod = CS_MAX(var_cal_opt->nswrsm, 1);

  isweep = 1;

  /* Reconstruction loop (beginning)
   *-------------------------------- */
  sinfo.n_it = 0;

  while ((isweep <= nswmod && residu > epsrsp*rnorm) || isweep == 1) {
    /* Solve on the increment dpvar */

    /*  Dynamic relaxation of the system */
    if (iswdyp >= 1) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          dpvarm1[iel][isou] = dpvar[iel][isou];
          dpvar[iel][isou] = 0.;
        }
      }
    }
    else {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 6; isou++)
          dpvar[iel][isou] = 0.;
      }
    }

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

#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          adxkm1[iel][isou] = adxk[iel][isou];
          adxk[iel][isou] = - rhs0[iel][isou];
        }
      }

      cs_balance_tensor(idtvar,
                        lvar,
                        imasac,
                        inc,
                        var_cal_opt,
                        dpvar,
                        NULL, /* dpvar */
                        coefats,
                        coefbts,
                        cofafts,
                        cofbfts,
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        viscel,
                        weighf,
                        weighb,
                        icvflb,
                        icvfli,
                        adxk);

      /* ||E.dx^(k-1)-E.0||^2 */
      nadxkm1 = nadxk;

      /* ||E.dx^k-E.0||^2 */
      nadxk = cs_gdot(6*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxk);

      /* < E.dx^k-E.0; r^k > */
      paxkrk = cs_gdot(6*n_cells, (cs_real_t *)smbrp, (cs_real_t *)adxk);

      /* Relaxation with respect to dx^k and dx^(k-1) */
      if (iswdyp >= 2) {
        /* < E.dx^(k-1)-E.0; r^k > */
        paxm1rk = cs_gdot(6*n_cells, (cs_real_t *)smbrp, (cs_real_t *)adxkm1);

        /* < E.dx^(k-1)-E.0; E.dx^k-E.0 > */
        paxm1ax = cs_gdot(6*n_cells, (cs_real_t *)adxk, (cs_real_t *)adxkm1);

        if ((nadxkm1 > 1.e-30*rnorm2)
         && (nadxk*nadxkm1-pow(paxm1ax,2)) > 1.e-30*rnorm2)
          beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-pow(paxm1ax,2));
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
        alph = -paxkrk/CS_MAX(nadxk, 1.e-30*rnorm2);
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/CS_MAX(nadxk, 1.e-30*rnorm2);
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

    /* --- Update the solution with the increment */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        for (cs_lnum_t isou = 0; isou < 6; isou++)
          pvar[iel][isou] += dpvar[iel][isou];
    }
    else if (iswdyp == 1) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        for (cs_lnum_t isou = 0; isou < 6; isou++)
           pvar[iel][isou] += alph*dpvar[iel][isou];
    }
    else if (iswdyp >= 2) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        for (cs_lnum_t isou = 0; isou < 6; isou++)
          pvar[iel][isou] += alph*dpvar[iel][isou] + beta*dpvarm1[iel][isou];
    }

    /* --> Handle parallelism and periodicity */

    if (cs_glob_rank_id  >=0 || cs_glob_mesh->n_init_perio > 0)
      cs_mesh_sync_var_sym_tens(pvar);

    /* --- Update the right hand and compute the new residual */

    if (iswdyp <= 0) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*dpvar[iel][0]
                            - fimp[iel][isou][1]*dpvar[iel][1]
                            - fimp[iel][isou][2]*dpvar[iel][2]
                            - fimp[iel][isou][3]*dpvar[iel][3]
                            - fimp[iel][isou][4]*dpvar[iel][4]
                            - fimp[iel][isou][5]*dpvar[iel][5];
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }
    else if (iswdyp == 1) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*alph*dpvar[iel][0]
                            - fimp[iel][isou][1]*alph*dpvar[iel][1]
                            - fimp[iel][isou][2]*alph*dpvar[iel][2]
                            - fimp[iel][isou][3]*alph*dpvar[iel][3]
                            - fimp[iel][isou][4]*alph*dpvar[iel][4]
                            - fimp[iel][isou][5]*alph*dpvar[iel][5];
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }
    else if (iswdyp == 2) {
#     pragma omp parallel for  if (n_cells > CS_THR_MIN)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        /* smbini already contains unsteady terms and mass source terms
         * of the RHS updated at each sweep */
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          smbini[iel][isou] = smbini[iel][isou]
                            - fimp[iel][isou][0]*(  alph*dpvar[iel][0]
                                                  + beta*dpvarm1[iel][0])
                            - fimp[iel][isou][1]*(  alph*dpvar[iel][1]
                                                  + beta*dpvarm1[iel][1])
                            - fimp[iel][isou][2]*(  alph*dpvar[iel][2]
                                                  + beta*dpvarm1[iel][2])
                            - fimp[iel][isou][3]*(  alph*dpvar[iel][3]
                                                  + beta*dpvarm1[iel][3])
                            - fimp[iel][isou][4]*(  alph*dpvar[iel][4]
                                                  + beta*dpvarm1[iel][4])
                            - fimp[iel][isou][5]*(  alph*dpvar[iel][5]
                                                  + beta*dpvarm1[iel][5]);
          smbrp[iel][isou] = smbini[iel][isou];
        }
      }
    }

    /* The added convective scalar mass flux is:
     *      (thetex*Y_\face-imasac*Y_\celli)*mf.
     * When building the implicit part of the rhs, one
     * has to impose 1 on mass accumulation. */
    imasac = 1;

    cs_balance_tensor(idtvar,
                      f_id,
                      imasac,
                      inc,
                      var_cal_opt,
                      pvar,
                      pvara,
                      coefats,
                      coefbts,
                      cofafts,
                      cofbfts,
                      i_massflux,
                      b_massflux,
                      i_visc,
                      b_visc,
                      viscel,
                      weighf,
                      weighb,
                      icvflb,
                      icvfli,
                      smbrp);

    /* --- Convergence test */
    residu = sqrt(cs_gdot(6*n_cells, (cs_real_t *)smbrp, (cs_real_t *)smbrp));

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
  if (fabs(rnorm)/sqrt(6.) > cs_math_epzero)
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
   * 3. Free solver setup
   *==========================================================================*/

  cs_sles_free_native(f_id, var_name);

  /* Free memory */
  BFT_FREE(dam);
  BFT_FREE(xam);
  BFT_FREE(smbini);
  BFT_FREE(dpvar);
  if (iswdyp >= 1) {
    BFT_FREE(adxk);
    BFT_FREE(adxkm1);
    BFT_FREE(dpvarm1);
    BFT_FREE(rhs0);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
