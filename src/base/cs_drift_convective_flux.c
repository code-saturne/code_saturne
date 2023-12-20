/*============================================================================
 * Compute the modified convective flux for scalars with a drift.
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
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"

#include "cs_array.h"
#include "cs_assert.h"
#include "cs_balance.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_convection_diffusion.h"
#include "cs_divergence.h"
#include "cs_face_viscosity.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_drift_convective_flux.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  ! \file cs_drift_scalar_convective_flux_compute.c
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
 * \brief Compute the modified convective flux for scalars with a drift.
 *
 * \param[in]     f_sc          drift scalar field
 * \param[in,out] i_mass_flux   scalar mass flux at interior face centers
 * \param[in,out] b_mass_flux   scalar mass flux at boundary face centers
 * \param[in,out] divflu        divergence of drift flux
 */
/*----------------------------------------------------------------------------*/

void
cs_drift_convective_flux(cs_field_t  *f_sc,
                         cs_real_t    i_mass_flux[],
                         cs_real_t    b_mass_flux[],
                         cs_real_t    divflu[])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)mesh->b_face_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;
  const cs_real_t *cell_vol = fvq->cell_f_vol;

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int keyccl = cs_field_key_id("scalar_class");
  const int keydri = cs_field_key_id("drift_scalar_model");
  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

  const int iscdri = cs_field_get_key_int(f_sc, keydri);
  const int icla = cs_field_get_key_int(f_sc, keyccl);

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_turb_model_type_t iturb  = cs_glob_turb_model->iturb;
  const int itytur = cs_glob_turb_model->itytur;
  const int n_fields = cs_field_n_fields();
  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  const int idtvar = cs_glob_time_step_options->idtvar;
  const int *bc_type = cs_glob_bc_type;

  /* Mass fraction of gas */

  cs_field_t *f_xc = cs_field_by_name_try("x_c");
  cs_real_t *x1 = NULL, *b_x1 = NULL;
  cs_real_t *i_mass_flux_gas = NULL;
  cs_real_t *b_mass_flux_gas = NULL;

  if (f_xc != NULL) {
    x1 = f_xc->val;

    /* Mass fraction of the gas at the boundary */
    cs_field_t *f_b_xc = cs_field_by_name("b_x_c");
    b_x1 = f_b_xc->val;

    /* Get the mass flux of the continuous phase (gas)
       that is the negative scalar class */

    int iflmas = -1;
    int iflmab = -1;

    for (int i = 0; i < n_fields; i++) {
      cs_field_t *f = cs_field_by_id(i);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;

      const int jcla = cs_field_get_key_int(f, keyccl);
      if (jcla == -1) {
        iflmas = cs_field_get_key_int(f, kimasf);
        iflmab = cs_field_get_key_int(f, kbmasf);
        break;
      }
    }

    assert(iflmab > -1);
    i_mass_flux_gas = cs_field_by_id(iflmas)->val;
    /* Pointer to the Boundary mass flux */
    b_mass_flux_gas = cs_field_by_id(iflmab)->val;
  }

  cs_field_t *f_vel = CS_F_(vel);
  cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_sc);
  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(f_vel);

  /* Pointers to the mass fluxes of the mix (based on mix velocity) */

  const int iflmas_v = cs_field_get_key_int(f_vel, kimasf);
  const int iflmab_v = cs_field_get_key_int(f_vel, kbmasf);
  cs_real_t *i_mass_flux_mix = cs_field_by_id(iflmas_v)->val;
  cs_real_t *b_mass_flux_mix = cs_field_by_id(iflmab_v)->val;

  /* Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;
  cs_real_3_t *vel_pre = (cs_real_3_t *)f_vel->val_pre;

  /* Initialization
     -------------- */

  /* Physical properties */

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_k = CS_F_(k);

  cs_real_6_t *rij = NULL;
  cs_real_t *k = NULL;

  if (f_rij != NULL)
    rij = (cs_real_6_t *)f_rij->val;
  if (f_k != NULL)
    k = f_k->val;

  /* Brownian diffusivity */
  cs_real_t *cpro_viscls = NULL;
  int ifcvsl = cs_field_get_key_int(f_sc, kivisl);
  if (ifcvsl >= 0)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  /* Vector containing all the additional convective terms */

  cs_real_t *w1, *viscce;
  cs_real_3_t *dudt;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscce, n_cells_ext, cs_real_t);
  BFT_MALLOC(dudt, n_cells_ext, cs_real_3_t);

  cs_real_t *coefap, *coefbp, *cofafp, *cofbfp;
  cs_real_3_t *coefa1;
  cs_real_33_t *coefb1;

  BFT_MALLOC(coefap, n_b_faces, cs_real_t);
  BFT_MALLOC(coefbp, n_b_faces, cs_real_t);
  BFT_MALLOC(cofafp, n_b_faces, cs_real_t);
  BFT_MALLOC(cofbfp, n_b_faces, cs_real_t);
  BFT_MALLOC(coefa1, n_b_faces, cs_real_3_t);
  BFT_MALLOC(coefb1, n_b_faces, cs_real_33_t);

  cs_real_t *i_visc, *flumas;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(flumas, n_i_faces, cs_real_t);

  cs_real_t *b_visc, *flumab;
  BFT_MALLOC(flumab, n_b_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  if (iscdri & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX) {

    /* Index of the corresponding relaxation time */
    cs_real_t *cpro_taup = NULL;
    {
      cs_field_t *f_tau
        = cs_field_by_composite_name_try("drift_tau", f_sc->name);

      if (f_tau != NULL)
        cpro_taup = f_tau->val;
    }

    /* Index of the corresponding relaxation time (cpro_taup) */
    cs_real_3_t *cpro_drift = NULL;
    {
      cs_field_t *f_drift_vel
        = cs_field_by_composite_name_try("drift_vel", f_sc->name);

      if (f_drift_vel != NULL)
        cpro_drift = (cs_real_3_t *)f_drift_vel->val;
    }

    /* Index of the corresponding interaction time
       particle--eddies (drift_turb_tau) */

    cs_real_t *cpro_taufpt = NULL;
    if (iscdri & CS_DRIFT_SCALAR_TURBOPHORESIS) {
      cs_field_t *f_drift_turb_tau
        = cs_field_by_composite_name("drift_turb_tau", f_sc->name);

      cpro_taufpt = f_drift_turb_tau->val;
    }

    /* Initialization of the convection flux for the current particle class */

    cs_array_real_fill_zero(n_i_faces, i_visc);
    cs_array_real_fill_zero(n_i_faces, flumas);

    cs_array_real_fill_zero(n_b_faces, b_visc);
    cs_array_real_fill_zero(n_b_faces, flumab);

    /* Initialization of the gas "class" convective flux by the
       first particle "class":
       it is initialized by the mass flux of the bulk */

    if (icla == 1 && f_xc != NULL) {
      cs_array_real_copy(n_i_faces, i_mass_flux_mix, i_mass_flux_gas);
      cs_array_real_copy(n_b_faces, b_mass_flux_mix, b_mass_flux_gas);
    }

    /* Initialize the additional convective flux with the gravity term
       --------------------------------------------------------------- */

    /* Test if a deviation velocity of particles class exists */

    if (icla >= 1) {

      char var_name[15];
      snprintf(var_name, 14, "vd_p_%02d", icla);
      var_name[14] = '\0';

      cs_field_t *f_vdp_i = cs_field_by_name_try(var_name);
      cs_real_3_t *vdp_i = NULL;

      if (f_vdp_i != NULL) {
        vdp_i = (cs_real_3_t *)f_vdp_i->val;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

          const cs_real_t rho = crom[c_id];
          cpro_drift[c_id][0] = rho * vdp_i[c_id][0];
          cpro_drift[c_id][1] = rho * vdp_i[c_id][1];
          cpro_drift[c_id][2] = rho * vdp_i[c_id][2];

        }
      }
    }

    else if (icla >= 0 && cpro_taup != NULL && cpro_drift != NULL) {

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t rho = crom[c_id];
        cpro_drift[c_id][0] = rho * cpro_taup[c_id] * gxyz[0];
        cpro_drift[c_id][1] = rho * cpro_taup[c_id] * gxyz[1];
        cpro_drift[c_id][2] = rho * cpro_taup[c_id] * gxyz[2];
      }

    }

    /* Computation of the turbophoresis and the thermophoresis terms
       ------------------------------------------------------------- */

    /* Initialized to 0 */
    cs_array_real_fill_zero(n_cells, viscce);

    if ((iscdri & CS_DRIFT_SCALAR_TURBOPHORESIS) && iturb != CS_TURB_NONE) {

      /* The diagonal part is easy to implicit (Grad (K) . n = (K_j - K_i)/IJ)
         Compute the K=1/3*trace(R) coefficient (diffusion of Zaichik) */

      if (itytur == 3) {

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t rtrace = cs_math_6_trace(rij[c_id]);

          /* Correction by Omega */
          const cs_real_t omega = cpro_taup[c_id] / cpro_taufpt[c_id];
          /* FIXME: use idifft or not? */
          viscce[c_id] = 1.0/3.0 * cpro_taup[c_id] / (1.0 + omega) * rtrace;
        }

      }
      else if (itytur == 2 || itytur == 5 || iturb == CS_TURB_K_OMEGA) {

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          /* Correction by Omega */
          const cs_real_t omega = cpro_taup[c_id] / cpro_taufpt[c_id];
          viscce[c_id] = 2.0/3.0 * cpro_taup[c_id] / (1.0 + omega) * k[c_id];
        }

      }

    } /* End turbophoresis */

    if (iscdri & CS_DRIFT_SCALAR_THERMOPHORESIS) {

      /* cpro_viscls[c_id]: contains the Brownian motion
         ------------------------------------------------ */

      if (ifcvsl >= 0) {

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          viscce[c_id] += cpro_viscls[c_id] / crom[c_id];

      }
      else {

        const int kvisl0 = cs_field_key_id("diffusivity_ref");
        const cs_real_t visls_0 = cs_field_get_key_double(f_sc, kvisl0);

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          viscce[c_id] += viscce[c_id] + visls_0 / crom[c_id];

      }

    } /* End thermophoresis */

    if (   (iscdri & CS_DRIFT_SCALAR_TURBOPHORESIS)
        || (iscdri & CS_DRIFT_SCALAR_THERMOPHORESIS)) {

      /* Face diffusivity of rho to compute rho*(Grad K . n)_face */
      cs_array_real_copy(n_cells, crom, w1);

      if (mesh->halo != NULL)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, w1);

      cs_face_viscosity(mesh,
                        fvq,
                        eqp_sc->imvisf,
                        w1,
                        i_visc,
                        b_visc);

      /* Homogeneous Neumann BC */
#     pragma omp parallel for if (n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        cs_boundary_conditions_set_neumann_scalar_hmg(&coefap[face_id],
                                                      &cofafp[face_id],
                                                      &coefbp[face_id],
                                                      &cofbfp[face_id]);

      /* The computed convective flux has the dimension of rho*velocity */

      cs_face_diffusion_potential(-1,
                                  mesh,
                                  fvq,
                                  0, /* init */
                                  1, /* inc */
                                  eqp_sc->imrgra,
                                  eqp_sc->nswrgr,
                                  eqp_sc->imligr,
                                  0, /* iphydr */
                                  0, /* iwgrp */
                                  eqp_sc->verbosity,
                                  eqp_sc->epsrgr,
                                  eqp_sc->climgr,
                                  NULL, /* frcxt */
                                  viscce,
                                  coefap, coefbp,
                                  cofafp, cofbfp,
                                  i_visc, b_visc,
                                  w1,
                                  flumas, flumab);

      /* TODO add extradiagonal part */

    } /* End turbophoresis or thermophoresis */

    /* Centrifugal force (particular derivative Du/Dt)
       ----------------------------------------------- */

    if (iscdri & CS_DRIFT_SCALAR_CENTRIFUGALFORCE) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t rhovdt = crom[c_id] * cell_vol[c_id] / dt[c_id];

        dudt[c_id][0] = - rhovdt * (vel[c_id][0]-vel_pre[c_id][0]);
        dudt[c_id][1] = - rhovdt * (vel[c_id][1]-vel_pre[c_id][1]);
        dudt[c_id][2] = - rhovdt * (vel[c_id][2]-vel_pre[c_id][2]);
      }

      /* Reset i_visc and b_visc */
      cs_array_real_fill_zero(n_i_faces, i_visc);
      cs_array_real_fill_zero(n_b_faces, b_visc);

      /* Get Boundary conditions of the velocity */
      cs_real_3_t  *coefa_vel = (cs_real_3_t  *)f_vel->bc_coeffs->a;
      cs_real_33_t *coefb_vel = (cs_real_33_t *)f_vel->bc_coeffs->b;
      cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)f_vel->bc_coeffs->af;
      cs_real_33_t *cofbf_vel = (cs_real_33_t *)f_vel->bc_coeffs->bf;

      /* The added convective scalar mass flux is:
         (thetap*Y_\face-imasac*Y_\celli)*mf.
         When building the implicit part of the rhs, one
         has to impose 1 on mass accumulation. */

      cs_equation_param_t eqp_loc = *eqp_vel;

      eqp_loc.iconv  = 1;
      eqp_loc.istat  = -1;
      eqp_loc.idiff  = 0;
      eqp_loc.idifft = -1;
      eqp_loc.iswdyn = -1;
      eqp_loc.nswrsm = -1;
      eqp_loc.iwgrec = 0;
      eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */
      eqp_loc.epsilo = -1;
      eqp_loc.epsrsm = -1;

      cs_balance_vector(idtvar,
                        CS_F_(vel)->id,
                        1, /* imasac */
                        1, /* inc */
                        0, /* ivisep */
                        &eqp_loc,
                        vel, vel,
                        coefa_vel, coefb_vel,
                        cofaf_vel, cofbf_vel,
                        i_mass_flux_mix, b_mass_flux_mix,
                        i_visc, b_visc,
                        NULL, NULL, /* secvif, secvib */
                        NULL, NULL, NULL,
                        0, NULL, /* icvflb, icvfli */
                        NULL, NULL,
                        dudt);

      /* Warning: cs_balance_vector adds "-( grad(u) . rho u)" */

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cpro_drift[c_id][0] =   cpro_drift[c_id][0]
                              + cpro_taup[c_id]*dudt[c_id][0]/cell_vol[c_id];

        cpro_drift[c_id][1] =   cpro_drift[c_id][1]
                              + cpro_taup[c_id]*dudt[c_id][1]/cell_vol[c_id];

        cpro_drift[c_id][2] =   cpro_drift[c_id][2]
                              + cpro_taup[c_id]*dudt[c_id][2]/cell_vol[c_id];
      }

    } /* End centrifugalforce */

    /* Electrophoresis term
       -------------------- */

    if (iscdri & CS_DRIFT_SCALAR_ELECTROPHORESIS) {

      /* TODO */
      bft_error(__FILE__, __LINE__, 0,
                _("The drift scalar electrophoresis "
                  "functionality is not yet available"));
    }

    /* Finalization of the mass flux of the current class
       -------------------------------------------------- */

    /* For all scalar with a drift excepted the gas phase which is deduced
       And for those whom mass flux is imposed elsewhere */

    if (icla >= 0 && !(iscdri & CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX)) {

      /* Homogeneous Neumann at the boundary */
      if (iscdri & CS_DRIFT_SCALAR_ZERO_BNDY_FLUX) {

#       pragma omp parallel for if (n_b_faces > CS_THR_MIN)
        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

          for (cs_lnum_t i = 0; i < 3; i++) {
            coefa1[face_id][i] = 0.0;

            for (cs_lnum_t j = 0; j < 3; j++)
              coefb1[face_id][i][j] = 0.0;
          }
        }

      }
      else if (iscdri & CS_DRIFT_SCALAR_ZERO_BNDY_FLUX_AT_WALLS) {

#       pragma omp parallel for if (n_b_faces > CS_THR_MIN)
        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

          for (cs_lnum_t i = 0; i < 3; i++) {
            coefa1[face_id][i] = 0.0;

            for (cs_lnum_t j = 0; j < 3; j++)
              coefb1[face_id][i][j] = 0.0;

            if (   bc_type[face_id] != CS_SMOOTHWALL
                && bc_type[face_id] != CS_ROUGHWALL)
              coefb1[face_id][i][i] = 1.0;

          }
        }
      }
      else {

#       pragma omp parallel for if (n_b_faces > CS_THR_MIN)
        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
          for (cs_lnum_t i = 0; i < 3; i++) {
            coefa1[face_id][i] = 0.0;

            for (cs_lnum_t j = 0; j < 3; j++)
              coefb1[face_id][i][j] = 0.0;

            coefb1[face_id][i][i] = 1.0;
          }
        }

      }

      cs_mass_flux(mesh,
                   fvq,
                   -1,
                   0, /* itypfl: drift has already been multiplied by rho */
                   0, /* iflmb0 */
                   0, /* init */
                   1, /* inc */
                   eqp_sc->imrgra,
                   eqp_sc->nswrgr,
                   eqp_sc->imligr,
                   eqp_sc->verbosity,
                   eqp_sc->epsrgr,
                   eqp_sc->climgr,
                   crom, brom,
                   cpro_drift,
                   coefa1, coefb1,
                   flumas, flumab);

      /* Update the convective flux, exception for the Gas "class" */
#     pragma omp parallel for if (n_i_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
        i_mass_flux[face_id] = i_mass_flux_mix[face_id] + flumas[face_id];

#     pragma omp parallel for if (n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        b_mass_flux[face_id] = b_mass_flux_mix[face_id] + flumab[face_id];

      /* Deduce the convective flux of the continuous "class" by removing
         the flux of the current particle "class":
         (rho x1 V1)_ij = (rho Vs)_ij - sum_classes (rho x2 V2)_ij
         ---------------------------------------------------------------- */

      if (icla >= 1) {

        char var_name[15];
        snprintf(var_name, 14, "x_p_%02d", icla);
        var_name[14] = '\0';

        cs_field_t *f_x_p_i = cs_field_by_name_try(var_name);
        cs_real_t *x2 = NULL;

        if (f_x_p_i != NULL) {
          x2 = f_x_p_i->val;

#         pragma omp parallel for if (n_i_faces > CS_THR_MIN)
          for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

            /* Upwind value of x2 at the face, consistent with the
               other transport equations */

            cs_lnum_t c_id_up = i_face_cells[face_id][1];

            if (i_mass_flux[face_id] >= 0.0)
              c_id_up = i_face_cells[face_id][0];

            i_mass_flux_gas[face_id] += -x2[c_id_up] * i_mass_flux[face_id];

          }

#         pragma omp parallel for if (n_b_faces > CS_THR_MIN)
          for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

            /* TODO Upwind value of x2 at the face, consistent with the
               other transport equations
               !if (bmasfl[face_id]>=0.d0) */
            cs_lnum_t c_id_up = b_face_cells[face_id];
            b_mass_flux_gas[face_id] += -x2[c_id_up] * b_mass_flux[face_id];

          }
        }
      }
    } /* End drift scalar imposed mass flux */

    /* Finalize the convective flux of the gas "class" by scaling by x1
       (rho x1 V1)_ij = (rho Vs)_ij - sum_classes (rho x2 V2)_ij
       Warning, x1 at the face must be computed so that it is consistent
       with an upwind scheme on (rho V1) */

    else if (icla == -1 && f_xc != NULL) {

#     pragma omp parallel for if (n_i_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

        /* Upwind value of x2 at the face, consistent with the
           other transport equations */
        cs_lnum_t c_id_up = i_face_cells[face_id][1];

        if (i_mass_flux_gas[face_id] >= 0.0)
          c_id_up = i_face_cells[face_id][0];

        i_mass_flux_gas[face_id] /= x1[c_id_up];

      }

#     pragma omp parallel for if (n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

        /* Upwind value of x1 at the face, consistent with the
           other transport equations */
        const cs_lnum_t c_id_up = b_face_cells[face_id];

        if (b_mass_flux_gas[face_id] < 0.0)
          b_mass_flux_gas[face_id] /= b_x1[face_id];
        else
          b_mass_flux_gas[face_id] /= x1[c_id_up];

      }

    }

  } /* End drift scalar add drift flux */

  /* Mass aggregation term of the additional part "div(rho(u_p-u_f))"
     ---------------------------------------------------------------- */

  /* Recompute the difference between mixture and the class */
  if (iscdri & CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX) {
#   pragma omp parallel for if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      flumas[face_id] = - i_mass_flux_mix[face_id];

#   pragma omp parallel for if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      flumab[face_id] = - b_mass_flux_mix[face_id];
  }
  else {
#   pragma omp parallel for if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      flumas[face_id] = i_mass_flux[face_id] - i_mass_flux_mix[face_id];

#   pragma omp parallel for if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      flumab[face_id] = b_mass_flux[face_id] - b_mass_flux_mix[face_id];
  }

  cs_divergence(mesh,
                1, /* init */
                flumas,
                flumab,
                divflu);

  /* Free memory */
  BFT_FREE(viscce);
  BFT_FREE(dudt);
  BFT_FREE(w1);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(flumas);
  BFT_FREE(flumab);
  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(cofafp);
  BFT_FREE(cofbfp);
  BFT_FREE(coefa1);
  BFT_FREE(coefb1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
