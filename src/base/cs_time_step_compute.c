/*============================================================================
 * Local time step and CFL computation.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"

#include "cs_array.h"
#include "cs_cf_compute.h"
#include "cs_coupling.h"
#include "cs_face_viscosity.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_log_iteration.h"
#include "cs_matrix_building.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_vof.h"

#include "cs_cf_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_step_compute.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_step_compute.c
        Local time step and CFL computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the local time step.
 *
 * \param[in]     itrale        ALE iteration number
 * \param[in]     iwarnp        verbosity
 * \param[in]     dt            time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_local_time_step_compute(int  itrale)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_real_3_t *cell_f_cen = (const cs_real_3_t *)fvq->cell_f_cen;

  cs_real_t *dt = CS_F_(dt)->val;
  const int nt_cur = cs_glob_time_step->nt_cur;
  int nt_max = cs_glob_time_step->nt_max;

  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;

  const int idtvar = cs_glob_time_step_options->idtvar;

  /* Initialization
     -------------- */

  bool log_is_active = cs_log_default_is_active();

  cs_field_t *vel = CS_F_(vel);
  cs_field_t *p = CS_F_(p);

  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);
  cs_equation_param_t *eqp_p = cs_field_get_equation_param(p);

  cs_field_t *f_courant_number = cs_field_by_name_try("courant_number");
  cs_field_t *f_fourier_number = cs_field_by_name_try("fourier_number");

  if (   !(eqp_vel->iconv >= 1 && f_courant_number != NULL)
      && !(eqp_vel->idiff >= 1  && f_fourier_number != NULL)
      && !(   (eqp_vel->iconv >= 1 || eqp_vel->idiff >= 1)
           && (eqp_vel->verbosity >= 2 || log_is_active))
      && !(   cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
           && (eqp_vel->verbosity >= 2 || log_is_active))
      && !(   idtvar == -1 || idtvar == 1 || idtvar == 2
           || (   (eqp_vel->verbosity >= 2 || log_is_active)
               && (   eqp_vel->idiff >= 1
                   || eqp_vel->iconv >= 1
                   || cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0))))
    return;

  /* Pointers to the mass fluxes */

  int iflmas_v
    = cs_field_get_key_int(vel, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux_vel = cs_field_by_id(iflmas_v)->val;

  int iflmab_v
    = cs_field_get_key_int(vel, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux_vel = cs_field_by_id(iflmab_v)->val;

  const cs_real_t *i_mass_flux_volf = NULL, *b_mass_flux_volf = NULL;
  if (cs_glob_vof_parameters->vof_model > 0) {
    cs_field_t *volf2 = CS_F_(void_f);

    int iflmas
      = cs_field_get_key_int(volf2, cs_field_key_id("inner_mass_flux_id"));
    i_mass_flux_volf = cs_field_by_id(iflmas)->val;

    int iflmab
      = cs_field_get_key_int(volf2, cs_field_key_id("boundary_mass_flux_id"));
    b_mass_flux_volf = cs_field_by_id(iflmab)->val;
  }

  /* Allocate temporary arrays for the time-step resolution */

  cs_real_t *i_visc, *b_visc, *dam, *cofbft, *coefbt;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  BFT_MALLOC(cofbft, n_b_faces, cs_real_t);
  BFT_MALLOC(coefbt, n_b_faces, cs_real_t);

  /* Allocate other arrays, depending on user options */
  cs_real_t *wcf = NULL;
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0)
    BFT_MALLOC(wcf, n_cells_ext, cs_real_t);

  /* Allocate work arrays */
  cs_real_t *w1, *w2, *w3;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(w3, n_cells_ext, cs_real_t);

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  /* Compute CFL like condition on the time step for
     positive density for the compressible module
     ------------------------------------------------ */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0)
    cs_cf_cfl_compute(wcf);

  /* Compute the diffusivity at the faces
     ------------------------------------ */

  if (eqp_vel->idiff >= 1) {
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w1[c_id] = viscl[c_id] + eqp_vel->idifft*visct[c_id];
    }

    cs_face_viscosity(mesh,
                      fvq,
                      eqp_vel->imvisf,
                      w1,
                      i_visc,
                      b_visc);
  }
  else {
    cs_array_real_set_scalar(n_i_faces, 0.0, i_visc);
    cs_array_real_set_scalar(n_b_faces, 0.0, b_visc);
  }

  /* Boundary condition for matrdt
     ----------------------------- */

  if (idtvar >= 0) {

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      if (b_mass_flux_vel[f_id] < 0.0) {
        const cs_lnum_t c_id = b_face_cells[f_id];
        const cs_real_t hint
          = eqp_vel->idiff*(  viscl[c_id]
                            + eqp_vel->idifft*visct[c_id])/b_dist[f_id];
        coefbt[f_id] = 0.0;
        cofbft[f_id] = hint;
      }
      else {
        coefbt[f_id] = 1.0;
        cofbft[f_id] = 0.0;
      }
    }

  }
  else {

    /* TODO for steady algorithm, check if using the third of the trace
       is appropriate, or if a better solution is available
       (algorithm was probably broken since velocity components are
       coupled) */

    const cs_real_t mult = 1.0 / 3.0;

    cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
    cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefbt[f_id] = cs_math_33_trace(coefb_vel[f_id]) * mult;
      cofbft[f_id] = cs_math_33_trace(cofbf_vel[f_id]) * mult;
    }

  }

  /* Steady algorithm
     ---------------- */

  if (idtvar >= 0) {

    /* Variable time step from imposed Courant and Fourie */

    /* We compute the max thermal time step
       (also when IDTVAR=0, for display)
       dttmax = 1/sqrt(max(0+,gradRO.g/RO) -> w3 */

    if (cs_glob_time_step_options->iptlro == 1) {

      /* Allocate a temporary array for the gradient calculation */
      cs_real_3_t *grad;
      BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

      cs_real_t *coefbr;
      BFT_MALLOC(coefbr, n_b_faces, cs_real_t);
      cs_array_real_set_scalar(n_b_faces, 0.0, coefbr);

      int imrgrp = eqp_p->imrgra;
      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
      cs_gradient_type_by_imrgra(imrgrp,
                                 &gradient_type,
                                 &halo_type);

      /* Compute gradient */

      cs_gradient_scalar("Work array",
                         gradient_type,
                         halo_type,
                         1,             /* inc */
                         eqp_p->nswrgr,
                         0,             /* iphydp */
                         1,             /* w_stride */
                         eqp_p->verbosity,
                         eqp_p->imligr,
                         eqp_p->epsrgr,
                         eqp_p->climgr,
                         NULL,          /* f_ext */
                         brom,          /* coefa */
                         coefbr,        /* coefb */
                         crom,          /* pvar */
                         NULL,          /* c_weight */
                         NULL,          /* cpl */
                         grad);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w3[c_id] = cs_math_3_dot_product(grad[c_id], gxyz) / crom[c_id];
        w3[c_id] = 1.0 / sqrt(cs_math_fmax(cs_math_epzero, w3[c_id]));
      }

      /* Free memory */
      BFT_FREE(grad);
      BFT_FREE(coefbr);

    }

    if (idtvar == 1 || idtvar == 2) {

      int icou = 0;
      int ifou = 0;

      /* Courant limitation
         ------------------ */

      if (cs_glob_time_step_options->coumax > 0.0 && eqp_vel->iconv >= 1) {

        /* icou = 1 indicates the existence of a limitation by the courant */
        icou = 1;

        cs_matrix_time_step(mesh,
                            eqp_vel->iconv,
                            0, /* Construction of U/DX (COURANT) = W1 */
                            2, /* Non symmetric matrix */
                            coefbt,
                            cofbft,
                            i_mass_flux_vel,
                            b_mass_flux_vel,
                            i_visc,
                            b_visc,
                            dam);

        /* Compute w1 = time step verifying CFL constraint given by the user */

        const cs_real_t coumax = cs_glob_time_step_options->coumax;

        /* When using VoF we also compute the volume Courant
           number (without rho). It replaces the mass Courant
           number constraint */

        if (cs_glob_vof_parameters->vof_model > 0) {

          cs_matrix_time_step(mesh,
                              eqp_vel->iconv,
                              0,
                              2,
                              coefbt,
                              cofbft,
                              i_mass_flux_volf,
                              b_mass_flux_volf,
                              i_visc,
                              b_visc,
                              dam);

          /* Compute w1 = time step verifying CFL constraint
             given by the user */

#         pragma omp parallel for if (n_cells > CS_THR_MIN)
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            cs_real_t d_vol = (cs_mesh_quantities_cell_is_active(fvq, c_id)) ?
              1.0 / cell_f_vol[c_id] : 0;
            w1[c_id] =   coumax
                       / cs_math_fmax(dam[c_id] * d_vol, cs_math_epzero);
          }

        }
        else {

#         pragma omp parallel for if (n_cells > CS_THR_MIN)
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            cs_real_t d_vol = (cs_mesh_quantities_cell_is_active(fvq, c_id)) ?
              1.0 / cell_f_vol[c_id] : 0;
            w1[c_id] =   coumax
                       / cs_math_fmax(dam[c_id] * d_vol / crom[c_id],
                                      cs_math_epzero);
          }

        }

        /* Uniform time step: we take the minimum of the constraint */
        if (idtvar == 1)  {
          cs_real_t w1min = cs_math_big_r;
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            w1min = cs_math_fmin(w1min, w1[c_id]);

          cs_parall_min(1, CS_REAL_TYPE, &w1min);

          cs_array_real_set_scalar(n_cells, w1min, w1);
        }

      }

      /* Fourier limitation
         ------------------ */

      if (cs_glob_time_step_options->foumax > 0.0 && eqp_vel->idiff >= 1) {

        /* ifou = 1 indicates the existence of a limitation by the fourier */
        ifou = 1;

        /* Construction of +2.NU/DX (FOURIER) = W2 */

        cs_matrix_time_step(mesh,
                            0, /*iconv0 */
                            eqp_vel->idiff,
                            1, /* Symmetric matrix */
                            coefbt,
                            cofbft,
                            i_mass_flux_vel,
                            b_mass_flux_vel,
                            i_visc,
                            b_visc,
                            dam);

        const cs_real_t foumax = cs_glob_time_step_options->foumax;

        /* Compute W2 = variable time step verifying the maximum
           Fourier number specified by the user */

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t d_vol = (cs_mesh_quantities_cell_is_active(fvq, c_id)) ?
            1.0 / cell_f_vol[c_id] : 0;
          cs_real_t w2_l = dam[c_id] * d_vol / crom[c_id];
          w2[c_id] = foumax / cs_math_fmax(w2_l, cs_math_epzero);
        }

        /* Uniform time step: we take the minimum of the constraint */

        if (idtvar == 1) {
          cs_real_t w2min = cs_math_big_r;
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            w2min = cs_math_fmin(w2min, w2[c_id]);

          cs_parall_min(1, CS_REAL_TYPE, &w2min);

          cs_array_real_set_scalar(n_cells, w2min, w2);
        }

      }

      /* Limitation for the compressible algorithm
         ----------------------------------------- */

      /* It's important to keep WCF intact: we reuse it
         below for display */

      int icoucf = 0;
      if (   cs_glob_time_step_options->coumax > 0.0
          && cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {

        icoucf = 1;

        /* Compute DAM = variable time step verifying the maximum
           CFL constraint specified by the user */

        const cs_real_t cflmmx = cs_glob_time_step_options->cflmmx;
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          dam[c_id] =   cflmmx
                      / cs_math_fmax(wcf[c_id], cs_math_epzero);
        }

        /* Uniform time step: we take the minimum of ther constraint */

        if (idtvar == 1) {
          cs_real_t w3min = cs_math_big_r;
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            w3min = cs_math_fmin(w3min, dam[c_id]);

          cs_parall_min(1, CS_REAL_TYPE, &w3min);

          cs_array_real_set_scalar(n_cells, w3min, dam);
        }

      }

      /* We take the most restrictive limitation
         --------------------------------------- */

      /* The minimum of the two if they exist, and
         the existing one if there is only one */

#     pragma omp parallel if (n_cells > CS_THR_MIN)
      {
        cs_lnum_t s_id, e_id;
        cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

        if (icou == 1 && ifou == 1) {
          for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
            w1[c_id] = cs_math_fmin(w1[c_id], w2[c_id]);
        }
        else if (icou == 0 && ifou == 1) {
          for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
            w1[c_id] = w2[c_id];
        }

        /* For compressible flows, the limitation associated with
           density must be taken into account. */

        if (icoucf == 1) {
          for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
            w1[c_id] = cs_math_fmin(w1[c_id], dam[c_id]);
        }

        /* Time step computation
           ----------------------- */

        /*  Progressive increase of the time step.
            Immediate time step down */

        for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {

          if (w1[c_id] >= dt[c_id]) {
            const cs_real_t unpvdt = 1.0 + cs_glob_time_step_options->varrdt;
            dt[c_id] = cs_math_fmin(unpvdt*dt[c_id], w1[c_id]);
          }
          else
            dt[c_id] = w1[c_id];

        }

      }  /* End of OpenMP section */

      /* We limit by the max "thermal" time step
         --------------------------------------- */

      /* DTTMAX = W3 = 1/SQRT(MAX(0+,gradRO.g/RO)
         we limit the time step to DTTMAX */

      if (cs_glob_time_step_options->iptlro == 1)  {

        /* Clip the time step to DTTMAX (display in ecrlis) */

        cs_lnum_t nclptr = 0;

        cs_real_t vmin = dt[0];
        cs_real_t vmax = dt[0];

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          vmin = cs_math_fmin(vmin, dt[c_id]);
          vmax = cs_math_fmax(vmax, dt[c_id]);
          if (dt[c_id] > w3[c_id])  {
            nclptr = nclptr + 1;
            dt[c_id] = w3[c_id];
          }

        }

        cs_log_iteration_clipping("dt (clip/dtrho)",
                                  1,
                                  0,
                                  nclptr,
                                  &vmin,
                                  &vmax);

        /* Uniform time step: we reuniform the time step */

        if (idtvar == 1) {
          cs_real_t w3min = cs_math_big_r;
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            w3min = cs_math_fmin(w3min, dt[c_id]);

          cs_parall_min(1, CS_REAL_TYPE, &w3min);

          cs_array_real_set_scalar(n_cells, w3min, dt);
        }

      }

      /* Clip the time step with respect to DTMIN and DTMAX
         --------------------------------------------------- */

      cs_lnum_t icfmin = 0;
      cs_lnum_t icfmax = 0;

      const cs_real_t dtmax = cs_glob_time_step_options->dtmax;
      const cs_real_t dtmin = cs_glob_time_step_options->dtmin;

      if (idtvar == 1) {

        cs_real_t dtloc = dt[0];
        if (dtloc > dtmax) {
          dtloc = dtmax;
          icfmax = n_cells;
        }
        if (dtloc < cs_glob_time_step_options->dtmin) {
          dtloc = cs_glob_time_step_options->dtmin;
          icfmin = n_cells;
        }

        int ntcam1 = nt_cur - 1;

        cs_coupling_sync_apps(0,      /* flags */
                              ntcam1,
                              &nt_max,
                              &dtloc);

        cs_log_iteration_clipping_field(CS_F_(dt)->id,
                                        icfmin,
                                        icfmax,
                                        dt,
                                        dt,
                                        &icfmin,
                                        &icfmax);

        int nt_prev = cs_glob_time_step->nt_prev;
        if (itrale > 0 && nt_max > nt_prev) {
          cs_time_step_update_dt(dtloc);

          if (log_is_active) {

            cs_log_printf
              (CS_LOG_DEFAULT,
               _("\n"
                 " INSTANT %18.9e   TIME STEP NUMBER %15d\n"
                 " ================================="
                 "============================\n\n\n"),
               cs_glob_time_step->t_cur, nt_cur);

          }
        }

        cs_array_real_set_scalar(n_cells, dtloc, dt);

      }
      else if (log_is_active || eqp_p->verbosity >= 2) {

        cs_real_t vmin = dt[0];
        cs_real_t vmax = dt[0];

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

          vmin = cs_math_fmin(vmin, dt[c_id]);
          vmax = cs_math_fmax(vmax, dt[c_id]);

          if (dt[c_id] > dtmax) {
            icfmax = icfmax + 1;
            dt[c_id] = dtmax;
          }
          if (dt[c_id] < dtmin) {
            icfmin = icfmin + 1;
            dt[c_id] = dtmin;
          }

        }

        cs_log_iteration_clipping_field(CS_F_(dt)->id,
                                        icfmin,
                                        icfmax,
                                        &vmin,
                                        &vmax,
                                        &icfmin,
                                        &icfmax);

      }

      if (eqp_p->verbosity >= 2) {
        cs_gnum_t cpt[2] = {icfmin, icfmax};
        cs_parall_counter(cpt, 2);

        cs_log_printf
          (CS_LOG_DEFAULT,
           _("\nDT CLIPPING: %llu at %11.4e, %llu at %11.4e\n"),
           (unsigned long long)cpt[0], dtmin,
           (unsigned long long)cpt[1], dtmax);
      }

  }

    /* Ratio DT/DTmax related to density effect (display in ecrlis) */

    if (   cs_glob_time_step_options->iptlro == 1
        && log_is_active) {

      cs_real_t *dtsdt0;
      BFT_MALLOC(dtsdt0, n_cells, cs_real_t);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        dtsdt0[c_id] = dt[c_id] / w3[c_id];

      cs_log_iteration_add_array("Dt/Dtrho max",
                                 "criterion",
                                 CS_MESH_LOCATION_CELLS,
                                 true,
                                 1,
                                 dtsdt0);

      BFT_FREE(dtsdt0);

    }

    /* Compute CFL constraint for compressible flow for log
       ---------------------------------------------------- */

    if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
        && (log_is_active || eqp_p->verbosity >= 2)) {

      cs_real_t cfmax = -HUGE_VAL;
      cs_real_t cfmin =  HUGE_VAL;
      cs_lnum_t icfmax = -1;
      cs_lnum_t icfmin = -1;

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w2[c_id] = wcf[c_id] * dt[c_id];

      cs_log_iteration_add_array("CFL / Mass",
                                 "criterion",
                                 CS_MESH_LOCATION_CELLS,
                                 true,
                                 1,
                                 w2);

      if (eqp_p->verbosity >= 2) {

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          if (w2[c_id] <= cfmin) {
            cfmin = w2[c_id];
            icfmin = c_id;
          }
          if (w2[c_id] >= cfmax) {
            cfmax = w2[c_id];
            icfmax = c_id;
          }
        }

        const cs_lnum_t min_c = cs_math_fmax(icfmin, 0);
        const cs_lnum_t max_c = cs_math_fmax(icfmax, 0);

        cs_real_t xyzmin[3];
        xyzmin[0] = cell_f_cen[min_c][0];
        xyzmin[1] = cell_f_cen[min_c][1];
        xyzmin[2] = cell_f_cen[min_c][2];

        cs_real_t xyzmax[3];
        xyzmax[0] = cell_f_cen[max_c][0];
        xyzmax[1] = cell_f_cen[max_c][1];
        xyzmax[2] = cell_f_cen[max_c][2];

        if (cs_glob_rank_id > -1) {
          cs_parall_min_loc_vals(3, &cfmin, xyzmin);
          cs_parall_max_loc_vals(3, &cfmax, xyzmax);
        }

        int icflag [2] = {icfmin, icfmax};
        cs_parall_max(2, CS_INT_TYPE, icflag);

        if (icflag[1] > -1)
          cs_log_printf
            (CS_LOG_DEFAULT,
             _("\n CFL/MAS MAX=%11.4e at (%11.4e %11.4e %11.4e)\n"),
             cfmax, xyzmax[0], xyzmax[1], xyzmax[2]);

        if (icflag[0] > -1)
          cs_log_printf
            (CS_LOG_DEFAULT,
             _("\n CFL/MAS MIN=%11.4e at (%11.4e %11.4e %11.4e)\n"),
             cfmin, xyzmin[0], xyzmin[1], xyzmin[2]);

        if (icflag[0] == -1 || icflag[1] == -1)
          cs_log_printf
            (CS_LOG_DEFAULT,
             _("\n CFL/MAS Too big to be displayed\n"));
      }

    }

    /* Steady algorithm
       ---------------- */
}
  else { /* if idtvar < 0 */

    int isym = 1;
    if (eqp_vel->iconv > 0)
      isym = 2;

    cs_matrix_time_step(mesh,
                        eqp_vel->iconv,
                        eqp_vel->idiff,
                        isym,
                        coefbt,
                        cofbft,
                        i_mass_flux_vel,
                        b_mass_flux_vel,
                        i_visc,
                        b_visc,
                        dt);

    const cs_real_t relaxv = eqp_vel->relaxv;
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      dt[c_id] =   relaxv * crom[c_id] * cell_f_vol[c_id]
                 / cs_math_fmax(dt[c_id], cs_math_epzero);
    }

  }

  /* Free memory */
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(dam);
  BFT_FREE(coefbt);
  BFT_FREE(cofbft);
  BFT_FREE(wcf);
  BFT_FREE(w1);
  BFT_FREE(w2);
  BFT_FREE(w3);
}

/*----------------------------------------------------------------------------*/
/*!
 *\brief Compute the local Courant and Fourier number to the log.
 *
 * This function has access to the boundary face type, except for the
 * first time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_courant_fourier_compute(void)
{

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_real_3_t *cell_f_cen = (const cs_real_3_t *)fvq->cell_f_cen;

  const int idtvar = cs_glob_time_step_options->idtvar;

  cs_real_t *dt = CS_F_(dt)->val;

  /* Initialization
     -------------- */

  bool log_is_active = cs_log_default_is_active();

  cs_field_t *vel = CS_F_(vel);
  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);
  cs_field_t *f_courant_number = cs_field_by_name_try("courant_number");
  cs_field_t *f_fourier_number = cs_field_by_name_try("fourier_number");

  if (   !(eqp_vel->iconv >= 1 && f_courant_number != NULL)
      && !(eqp_vel->idiff >= 1  && f_fourier_number != NULL)
      && !(   (eqp_vel->iconv >= 1 || eqp_vel->idiff >= 1)
           && (eqp_vel->verbosity >= 2 || log_is_active))
      && !(   cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
           && (eqp_vel->verbosity >= 2 || log_is_active))
      && !(   idtvar == -1 || idtvar == 1 || idtvar == 2
           || (   (eqp_vel->verbosity >= 2 || log_is_active)
               && (   eqp_vel->idiff >= 1
                   || eqp_vel->iconv >= 1
                   || cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0))))
    return;

  /* Pointers to the mass fluxes */

  int iflmas_v
    = cs_field_get_key_int(vel, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux_vel = cs_field_by_id(iflmas_v)->val;

  int iflmab_v
    = cs_field_get_key_int(vel, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux_vel = cs_field_by_id(iflmab_v)->val;

  const cs_real_t *i_mass_flux_volf = NULL, *b_mass_flux_volf = NULL;
  if (cs_glob_vof_parameters->vof_model > 0) {
    cs_field_t *volf2 = CS_F_(void_f);

    int iflmas
      = cs_field_get_key_int(volf2, cs_field_key_id("inner_mass_flux_id"));
    i_mass_flux_volf = cs_field_by_id(iflmas)->val;

    int iflmab
      = cs_field_get_key_int(volf2, cs_field_key_id("boundary_mass_flux_id"));
    b_mass_flux_volf = cs_field_by_id(iflmab)->val;
  }

  /* Allocate temporary arrays for the time-step resolution */
  cs_real_t *i_visc, *b_visc, *dam, *cofbft, *coefbt;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(cofbft, n_b_faces, cs_real_t);
  BFT_MALLOC(coefbt, n_b_faces, cs_real_t);

  /* Allocate work arrays */
  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  const cs_real_t *crom = CS_F_(rho)->val;

  cs_real_t *courant_number = f_courant_number->val;
  cs_real_t *fourier_number = f_fourier_number->val;

  cs_real_t *vol_courant_number = NULL;
  if (cs_glob_vof_parameters->vof_model > 0)
    vol_courant_number = cs_field_by_name("volume_courant_number")->val;

  /* Compute the diffusivity at the faces
     ------------------------------------ */

  if (eqp_vel->idiff >= 1) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      w1[c_id] = viscl[c_id] + eqp_vel->idifft*visct[c_id];

    cs_face_viscosity(mesh,
                      fvq,
                      eqp_vel->imvisf,
                      w1,
                      i_visc,
                      b_visc);
  }
  else {
    cs_array_real_set_scalar(n_i_faces, 0.0, i_visc);
    cs_array_real_set_scalar(n_b_faces, 0.0, b_visc);
  }

  /* Boundary condition for matrdt
     ----------------------------- */

  if (idtvar >= 0) {

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      if (b_mass_flux_vel[f_id] < 0.0) {
        const cs_lnum_t c_id = b_face_cells[f_id];

        const cs_real_t hint
          =   eqp_vel->idiff*(viscl[c_id]
            + eqp_vel->idifft*visct[c_id])/b_dist[f_id];

        coefbt[f_id] = 0.0;
        cofbft[f_id] = hint;
      }
      else {
        coefbt[f_id] = 1.0;
        cofbft[f_id] = 0.0;
      }
    }

  }
  else {

    /* TODO for steady algorithm, check if using the third of the trace
       is appropriate, or if a better solution is available
       (algorithm was probably broken since velocity components are
       coupled) */

    const cs_real_t mult = 1.0 / 3.0;

    cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
    cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefbt[f_id] = cs_math_33_trace(coefb_vel[f_id]) * mult;
      cofbft[f_id] = cs_math_33_trace(cofbf_vel[f_id]) * mult;
    }

  }

  /* 1: Compute Courant number for log
     2: Compute VoF Courant number for log
     3: Compute Fourier number for log
     4: Compute Courant/Fourier ratio for log */

  bool is_courant[4] = {(eqp_vel->iconv >= 1 && f_courant_number != NULL),
                        (eqp_vel->iconv >= 1 && vol_courant_number != NULL),
                        (eqp_vel->idiff >= 1 && f_fourier_number != NULL),
                        ((   eqp_vel->idiff >= 1 || eqp_vel->iconv >= 1)
                         && cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)};

  int idiff[4] = {0, 0, eqp_vel->idiff, eqp_vel->idiff};
  int iconv[4] = {eqp_vel->iconv, eqp_vel->iconv, 0, eqp_vel->iconv};
  int isym[4] = {2, 2, 1, 1};
  if (eqp_vel->iconv > 0)
    isym[3] = 2;

  const cs_real_t *i_mass_flux_list[4] = {i_mass_flux_vel, i_mass_flux_volf,
                                          i_mass_flux_vel, i_mass_flux_vel};

  const cs_real_t *b_mass_flux_list[4] = {b_mass_flux_vel, b_mass_flux_volf,
                                          b_mass_flux_vel, b_mass_flux_vel};

  const char *info[] = {"COURANT", "VOLUME COURANT",
                        "FOURIER", "COURANT/FOURIER"};

  cs_real_t *cpro_tab_list[4] = {courant_number, vol_courant_number,
                                 fourier_number, w1};

  for (int i = 0; i < 4; i++) {

    if (!(is_courant[i]))
      continue;

    /* 1: Build matrix of  U/DX(Courant) = w1
       2: Build matrix of  U/DX(Courant) = w1
       3: Build matrix +2.NU/DX (Fourier)= w1
       4: Build matrix    U/DX +2.NU/DX  (COURANT +FOURIER) = w1 */

    const cs_real_t *i_mass_flux = i_mass_flux_list[i];
    const cs_real_t *b_mass_flux = b_mass_flux_list[i];

    cs_matrix_time_step(mesh,
                        iconv[i],
                        idiff[i],
                        isym[i], /* 1 -> symmetric matrix, 2 not */
                        coefbt,
                        cofbft,
                        i_mass_flux,
                        b_mass_flux,
                        i_visc,
                        b_visc,
                        dam);

    /* Compute min and max Courant numbers */
    cs_real_t cfmax = -HUGE_VAL;
    cs_real_t cfmin =  HUGE_VAL;
    cs_lnum_t icfmax = -1;
    cs_lnum_t icfmin = -1;

    cs_real_t *cpro_tab = cpro_tab_list[i];

    if (i != 1) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t d_vol = (cs_mesh_quantities_cell_is_active(fvq, c_id)) ?
          1.0 / cell_f_vol[c_id] : 0;
        cpro_tab[c_id]
          = dam[c_id] * d_vol * dt[c_id] / crom[c_id];
      }
    }
    else if (i == 1) { /* Volume courant */
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t d_vol = (cs_mesh_quantities_cell_is_active(fvq, c_id)) ?
          1.0 / cell_f_vol[c_id] : 0;
        cpro_tab[c_id] = dam[c_id] * d_vol * dt[c_id];
      }
    }

    if (i == 3)
      cs_log_iteration_add_array("Courant/Fourier",
                                 "criterion",
                                 CS_MESH_LOCATION_CELLS,
                                 true,
                                 1,
                                 cpro_tab);

    if (eqp_vel->verbosity >= 2) {

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (cpro_tab[c_id] <= cfmin) {
          cfmin = cpro_tab[c_id];
          icfmin = c_id;
        }
        if (cpro_tab[c_id] >= cfmax) {
          cfmax = cpro_tab[c_id];
          icfmax = c_id;
        }
      }

      const cs_lnum_t min_c = cs_math_fmax(icfmin, 0);
      const cs_lnum_t max_c = cs_math_fmax(icfmax, 0);

      cs_real_t xyzmin[3];
      xyzmin[0] = cell_f_cen[min_c][0];
      xyzmin[1] = cell_f_cen[min_c][1];
      xyzmin[2] = cell_f_cen[min_c][2];

      cs_real_t xyzmax[3];
      xyzmax[0] = cell_f_cen[max_c][0];
      xyzmax[1] = cell_f_cen[max_c][1];
      xyzmax[2] = cell_f_cen[max_c][2];

      if (cs_glob_rank_id > -1) {
        cs_parall_min_loc_vals(3, &cfmin, xyzmin);
        cs_parall_max_loc_vals(3, &cfmax, xyzmax);
      }

      int icflag [2] = {icfmin, icfmax};
      cs_parall_max(2, CS_INT_TYPE, icflag);

      if (icflag[1] > -1)
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("\n %-s MAX=%11.4e at (%11.4e %11.4e %11.4e)\n"),
           info[i], cfmax, xyzmax[0],xyzmax[1],xyzmax[2]);

      if (icflag[0] > -1)
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("\n %-s MIN=%11.4e at (%11.4e %11.4e %11.4e)\n"),
           info[i], cfmin, xyzmin[0],xyzmin[1],xyzmin[2]);

      if (icflag[0] == -1 || icflag[1] == -1)
        cs_log_printf(CS_LOG_DEFAULT,
                      _("\n %-s Too big to be displayed\n"),
                      info[i]);

    }
  }

  /* Free memory */
  BFT_FREE(w1);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(dam);
  BFT_FREE(cofbft);
  BFT_FREE(coefbt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
