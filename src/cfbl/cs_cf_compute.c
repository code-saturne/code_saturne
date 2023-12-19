/*============================================================================
 * Operators for compressible flows
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
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_cf_thermo.h"
#include "cs_convection_diffusion.h"
#include "cs_divergence.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_volume_mass_injection.h"
#include "cs_mass_source_terms.h"
#include "cs_matrix_building.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"
//#include "cs_vof.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_compute.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 *! \file cs_cf_compute.c
 *
 * \brief Computations for compressible flows.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the "mass flux" at the faces for the CFL restriction
 *        computation and the solving of the pressure
 *
 * \param[in]   iterns            Navier-Stokes iteration number
 * \param[out]  i_mass_flux       Internal faces mass flux
 * \param[out]  b_mass_flux       Boundary faces mass flux
 */
/*----------------------------------------------------------------------------*/

static void
_compressible_pressure_mass_flux(int iterns, // cfmsfp en fortran
                                 cs_real_t i_mass_flux[],
                                 cs_real_t b_mass_flux[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_real_t *dt = CS_F_(dt)->val;

  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  const int itytur = cs_glob_turb_model->itytur;
  const cs_velocity_pressure_model_t *vp_model = cs_glob_velocity_pressure_model;
  const int idtvar = cs_glob_time_step_options->idtvar;

  /* Initialization
     -------------- */

  cs_field_t *vel = CS_F_(vel);
  cs_real_3_t *vela = (cs_real_3_t *)vel->val_pre;

  cs_real_3_t  *coefa_vel = (cs_real_3_t  *)vel->bc_coeffs->a;
  cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
  cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)vel->bc_coeffs->af;
  cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

  /* Allocate work arrays */
  cs_real_t *w1;
  cs_real_3_t *tsexp, *gavinj, *vel0;
  cs_real_33_t *coefbv, *tsimp;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(tsexp, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(gavinj, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(vel0, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(coefbv, n_b_faces, cs_real_33_t);
  BFT_MALLOC(tsimp, n_cells_ext, cs_real_33_t);

  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);
  cs_equation_param_t *eqp_p = cs_field_get_equation_param(CS_F_(p));

  cs_real_t *i_visc = NULL, *b_visc = NULL;
  cs_real_6_t *viscce = NULL;
  if (eqp_vel->idften & CS_ISOTROPIC_DIFFUSION) {
    BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
    BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  }
  else if (eqp_vel->idften & CS_ANISOTROPIC_LEFT_DIFFUSION) {
    BFT_MALLOC(i_visc, 9*n_i_faces, cs_real_t);
    BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
    BFT_MALLOC(viscce, n_cells_ext, cs_real_6_t);
  }

  cs_real_t *secvib = NULL, *secvif = NULL;
  if (vp_model->ivisse == 1) {
    BFT_MALLOC(secvif, n_i_faces, cs_real_t);
    BFT_MALLOC(secvib, n_b_faces, cs_real_t);
  }

  /* Density */

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  /* Mass flux at the faces
     ---------------------- */

  /* Source terms of the momentum equations
     -------------------------------------- */

  /* Some first tests (double expansion waves in a shock tube)
     has shown that taking into account all the
     momentum equation terms in the mass equation seems to be a
     bad idea (in particular the convective term, but the diffusive
     term, the transposed gradient, the mass and user source terms
     were all null in the considered tests).
     However, it may be due to a bug at that early stage of implementation
     of the algorithm (but we didn't find it).
     We thus recommand not to take into account the momentum source terms,
     except the gravity term (because it is in balance with the pressure
     gradient and because its effect is visible at equilibrium).
     However, we keep here the implementation of the preliminary tests
     (1.1.0.h version) with an overall test so that the correction is not
     active (thus, there is no user question and there is always the
     possibility to perform other tests in the future).
     Note that, with these terms, the thoeretical analysis is harder
     (Without these terms we are in the configuration Euler + gravity) */

  /* Initialization */
  cs_array_real_set_scalar(3*n_cells, 0.0, (cs_real_t *)tsexp);
  cs_array_real_set_scalar(9*n_cells, 0.0, (cs_real_t *)tsimp);

  /* Test on momentum source terms */
  int itsqdm = 0;

  if (itsqdm != 0) { /* we never enter here since itsqdm is fixed above */

    if (vp_model->ivisse == 1)
      cs_face_viscosity_secondary(secvif, secvib);

    cs_user_source_terms(cs_glob_domain,
                         vel->id,
                         (cs_real_t *)tsexp,
                         (cs_real_t *)tsimp);

    /* Mass flux computation */

    cs_mass_flux(mesh,
                 fvq,
                 vel->id,
                 1,       /* itypfl */
                 1,       /* iflum0 */
                 1,       /* init */
                 1,       /* inc */
                 eqp_vel->imrgra,
                 eqp_vel->nswrgr,
                 eqp_vel->imligr,
                 eqp_vel->verbosity,
                 eqp_vel->epsrgr,
                 eqp_vel->climgr,
                 crom,
                 brom,
                 vela,
                 coefa_vel,
                 coefb_vel,
                 i_mass_flux,
                 b_mass_flux);

    /* Face diffusivity for the velocity */
    int imvisp = eqp_vel->imvisf;

    if (eqp_vel->idiff >= 1) {

      const cs_real_t *viscl = CS_F_(mu)->val;
      const cs_real_t *visct = CS_F_(mu_t)->val;

      if (itytur == 3) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = viscl[c_id];
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          w1[c_id] = viscl[c_id] + eqp_vel->idifft * visct[c_id];
      }

      /* Scalar diffusivity (Default) */
      if (eqp_vel->idften & CS_ISOTROPIC_DIFFUSION) {

        cs_face_viscosity(mesh,
                          fvq,
                          imvisp, // a voir avec Thomas
                          w1,
                          i_visc,
                          b_visc);
      }
      /* Tensorial diffusion of the velocity (in case of tensorial porosity) */
      else if (eqp_vel->idften & CS_ANISOTROPIC_LEFT_DIFFUSION) {

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

          for (cs_lnum_t i = 0; i < 3; i++)
            viscce[c_id][i] = w1[c_id];

          for (cs_lnum_t i = 3; i < 6; i++)
            viscce[c_id][i] = 0.0;

        }

        cs_face_anisotropic_viscosity_vector(mesh,
                                             fvq,
                                             imvisp,
                                             viscce,
                                             (cs_real_33_t *)i_visc,
                                             b_visc);
      }

    }

    /* If no diffusion, viscosity is set to 0. */
    else {
      cs_array_real_set_scalar(n_i_faces, 0.0, i_visc);
      cs_array_real_set_scalar(n_b_faces, 0.0, b_visc);
    }

    /* The added convective scalar mass flux is:
       (thetap*Y_\face-imasac*Y_\celli)*mf.
       When building the implicit part of the rhs, one
       has to impose 1 on mass accumulation. */

    cs_equation_param_t eqp_vel_loc = *eqp_vel;

    eqp_vel_loc.istat  = -1;
    eqp_vel_loc.idifft = -1;
    eqp_vel_loc.iswdyn = -1;
    eqp_vel_loc.nswrsm = -1;
    eqp_vel_loc.iwgrec = 0;
    eqp_vel_loc.blend_st = 0; // Warning, may be overwritten if a field
    eqp_vel_loc.epsilo = -1;
    eqp_vel_loc.epsrsm = -1;

    int *icvfli = cs_cf_get_icvfli();

    cs_balance_vector(idtvar,
                      vel->id,
                      1, /* imasac */
                      1, /* inc */
                      vp_model->ivisse,
                      &eqp_vel_loc,
                      vela,
                      (const cs_real_3_t *)vela,
                      coefa_vel,
                      coefb_vel,
                      cofaf_vel,
                      cofbf_vel,
                      i_mass_flux,
                      b_mass_flux,
                      i_visc,
                      b_visc,
                      secvif,
                      secvib,
                      NULL,   /* viscel */
                      NULL,   /* weighf */ // rvoid a voir
                      NULL,   /* weighb */
                      0,      /* icvflb */
                      icvfli,
                      NULL, /* i_pvar */
                      NULL, /* b_pvar */
                      tsexp);

  }

  /* End of the test on momentum source terms */

  cs_lnum_t ncetsm
    = cs_volume_zone_n_type_cells(CS_VOLUME_ZONE_MASS_SOURCE_TERM);

  /* Mass source term */
  if (ncetsm > 0) {

    /* The momentum balance is used in its conservative form here
       so the mass source term is only composed of gamma*uinj
       => array of previous velocity has to be set to zero */

    cs_array_real_set_scalar(3*n_cells, 0.0, (cs_real_t *)vel0);

    int *itypsm = NULL;
    cs_lnum_t ncesmp = 0;
    cs_lnum_t *icetsm = NULL;
    cs_real_t *smacel_vel, *smacel_p = NULL;

    cs_volume_mass_injection_get_arrays(vel,
                                        &ncesmp,
                                        &icetsm,
                                        &itypsm,
                                        &smacel_vel,
                                        &smacel_p);

    cs_mass_source_terms(iterns,
                         3,
                         ncetsm,
                         icetsm,
                         itypsm,
                         cell_f_vol,
                         (cs_real_t *)vel0,
                         smacel_vel,
                         smacel_p,
                         (cs_real_t *)tsexp,
                         (cs_real_t *)tsimp,
                         (cs_real_t *)gavinj);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 3; i++)
        tsexp[c_id][i] += gavinj[c_id][i];
    }

  }


# pragma omp parallel if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Volumic forces term (gravity) */
    const cs_real_t rom = crom[c_id];
    for (cs_lnum_t i = 0; i < 3; i++)
      tsexp[c_id][i] = gxyz[i] + tsexp[c_id][i]/rom;

    /* Calculation of the convective "velocities at the cell centers
       (Calculation of u^n+dt*f^n) */
    for (cs_lnum_t i = 0; i < 3; i++)
      tsexp[c_id][i] *= dt[c_id];

  }

  /* Computation of the flux
     Volumic flux part based on dt*f^n */

  /* No contribution of f to the boundary mass flux */
  cs_array_real_set_scalar(9*n_b_faces, 0.0, (cs_real_t *)coefbv);

  cs_mass_flux(mesh,
               fvq,
               -1,      /* field_id */
               0,       /* itypfl, Velocity flux (crom, brom not used) */
               1,       /* iflum0 */
               1,       /* init */
               0,       /* inc */
               eqp_p->imrgra,
               0,       /* nswrgp */
               eqp_p->imligr,
               eqp_p->verbosity,
               eqp_p->epsrgr,
               eqp_p->climgr,
               crom,
               brom,
               tsexp,
               coefa_vel,
               coefbv,
               i_mass_flux,
               b_mass_flux);

  /* Volumic flux part based on velocity u^n
     take into account Dirichlet velocity boundary conditions */

  cs_mass_flux(mesh,
               fvq,
               vel->id,
               0,       /* itypfl, Velocity flux (crom, brom not used) */
               1,       /* iflum0 */
               0,       /* init */
               1,       /* inc */
               eqp_p->imrgra,
               0,       /* nswrgp */
               eqp_p->imligr,
               eqp_p->verbosity,
               eqp_p->epsrgr,
               eqp_p->climgr,
               crom,
               brom,
               vela,
               coefa_vel,
               coefb_vel,
               i_mass_flux,
               b_mass_flux);

  /* Free memory */
  BFT_FREE(w1);
  BFT_FREE(tsexp);
  BFT_FREE(gavinj);
  BFT_FREE(tsimp);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(secvif);
  BFT_FREE(secvib);
  BFT_FREE(viscce);
  BFT_FREE(coefbv);
  BFT_FREE(vel0);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the convective mass flux before the velocity prediction step.
 *        It is the first step of the compressible algorithm at each time
          iteration.
 *
 * This function solves the continuity equation in pressure formulation and then
 * updates the density and the mass flux.
 *
 * \param[in]     iterns        Navier-Stokes iteration number
*/
/*----------------------------------------------------------------------------*/

void
cs_cf_convective_mass_flux(int  iterns)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_real_t *dt = CS_F_(dt)->val;
  int idtvar = cs_glob_time_step_options->idtvar;

  /* Initialization
     -------------- */

  /* Allocate temporary arrays for the mass resolution */

  cs_real_t *i_visc, *wflmas, *ivolfl;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(wflmas, n_i_faces, cs_real_t);
  BFT_MALLOC(ivolfl, n_i_faces, cs_real_t);

  cs_real_t *b_visc, *wflmab, *bvolfl, *wbfa, *wbfb;
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  BFT_MALLOC(wflmab, n_b_faces, cs_real_t);
  BFT_MALLOC(bvolfl, n_b_faces, cs_real_t);
  BFT_MALLOC(wbfa, n_b_faces, cs_real_t);
  BFT_MALLOC(wbfb, n_b_faces, cs_real_t);

  cs_real_t *smbrs, *rovsdt;
  BFT_MALLOC(smbrs, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  /* Allocate work arrays */

  cs_real_t *w1, *w7, *w8, *w9, *w10, *dpvar;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(w7, n_cells_ext, cs_real_t);
  BFT_MALLOC(w8, n_cells_ext, cs_real_t);
  BFT_MALLOC(w9, n_cells_ext, cs_real_t);
  BFT_MALLOC(w10, n_cells_ext, cs_real_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_field_t *f_p = CS_F_(p);
  cs_field_t *e_tot = CS_F_(e_tot);

  /* Mass flux associated to energy */

  int iflmas_e
    = cs_field_get_key_int(e_tot, cs_field_key_id("inner_mass_flux_id"));
  cs_real_t *i_mass_flux_e = cs_field_by_id(iflmas_e)->val;

  int iflmab_e
    = cs_field_get_key_int(e_tot, cs_field_key_id("boundary_mass_flux_id"));
  cs_real_t *b_mass_flux_e = cs_field_by_id(iflmab_e)->val;

  cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  cs_real_t *cvar_pr = f_p->val;
  cs_real_t *cvar_pr_pre = f_p->val_pre;

  cs_equation_param_t *eqp_p = cs_field_get_equation_param(f_p);

  cs_real_t *cvar_fracv, *cvar_fracm, *cvar_frace;
  cvar_fracv = NULL;
  cvar_fracm = NULL;
  cvar_frace = NULL;

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > 1) {
    cvar_fracv = CS_F_(volume_f)->val;
    cvar_fracm = CS_F_(mass_f)->val;
    cvar_frace = CS_F_(energy_f)->val;
  }

  if (eqp_p->verbosity >= 1)
    cs_log_printf
      (CS_LOG_DEFAULT,
       _(" ** RESOLUTION FOR THE PRESSURE VARIABLE\n"
         "    ------------------------------------\n"));

  cs_real_t *cofaf_p = f_p->bc_coeffs->af;
  cs_real_t *cofbf_p = f_p->bc_coeffs->bf;

  const int icp = fluid_props->icp;
  const int icv = fluid_props->icv;
  cs_real_t *cpro_cp = NULL, *cpro_cv = NULL;

  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  if (icv >= 0)
    cpro_cv = cs_field_by_id(icv)->val;

  /* Computation of the boundary coefficients for the pressure
     gradient recontruction in accordance with the diffusion
     boundary coefficients (coefaf_p, coefbf_p) */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    const cs_lnum_t c_id = b_face_cells[f_id];
    const cs_real_t hint = dt[c_id] / b_dist[f_id];

    /* TODO check this: cofaf_p both input and output here...
       alternative:

       const cs_real_t qimp = coefaf_p[f_id];
       wbfa[f_id] = -qimp / cs_math_fmax(hint, 1e-300);
       wbfb[f_id] = 1.0;
    */

    cs_real_t bf_unused[1];

    cs_boundary_conditions_set_neumann_scalar(&wbfa[f_id],
                                              &cofaf_p[f_id],
                                              &wbfb[f_id],
                                              bf_unused,
                                              cofaf_p[f_id],
                                              hint);
  }

  /* Source terms
     ------------ */

  /* Initialization */

  cs_array_real_set_scalar(n_cells, 0.0, smbrs);
  cs_array_real_set_scalar(n_cells, 0.0, rovsdt);

  /* Mass source term
     ---------------- */

  cs_lnum_t ncetsm
    = cs_volume_zone_n_type_cells(CS_VOLUME_ZONE_MASS_SOURCE_TERM);

  if (ncetsm > 0) {

    cs_lnum_t ncesmp = 0;
    cs_lnum_t *icetsm = NULL;
    int *itpsm_p = NULL;
    cs_real_t *smcel_p = NULL; //, *gamma = NULL;

    cs_volume_mass_injection_get_arrays(f_p,
                                        &ncesmp,
                                        &icetsm,
                                        &itpsm_p,
                                        &smcel_p,
                                        NULL); // sinon mettre gamma

    for (cs_lnum_t ii = 0; ii < ncesmp; ii++) {
      const cs_lnum_t c_id = icetsm[ii];
      smbrs[c_id] = smbrs[c_id] + smcel_p[ii]*cell_f_vol[c_id];
    }
  }

  /* Unsteady term
     ------------- */

  /* Computation of the square of sound velocity c2.
     Pressure is an unsteady variable in this algorithm
     Varpos has been modified for that. */

  cs_real_t *c2;
  BFT_MALLOC(c2, n_cells_ext, cs_real_t);

  cs_cf_thermo_c_square(cpro_cp,
                        cpro_cv,
                        cvar_pr,
                        crom,
                        cvar_fracv,
                        cvar_fracm,
                        cvar_frace,
                        c2,
                        n_cells);

# pragma omp parallel if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rovsdt[c_id] =   rovsdt[c_id]
                   + eqp_p->istat*(cell_f_vol[c_id]/(dt[c_id]*c2[c_id]));
  }

  /* "Mass flux" and face "viscosity" computation
     -------------------------------------------- */

  /* Computation of the "convective flux" for the density */

  /* Volumic flux (u + dt f) */
  _compressible_pressure_mass_flux(iterns,
                                   ivolfl,
                                   bvolfl);

  /* Mass flux at internal faces (upwind scheme for the density)
     (negative because added to RHS) */

# pragma omp parallel if(n_i_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    const cs_lnum_t c_id1 = i_face_cells[f_id][0];
    const cs_lnum_t c_id2 = i_face_cells[f_id][1];

    wflmas[f_id] = -0.5 *
      (  crom[c_id1]*(ivolfl[f_id]+fabs(ivolfl[f_id]))
       + crom[c_id2]*(ivolfl[f_id]-fabs(ivolfl[f_id])));
  }

  /* Mass flux at boundary faces
     (negative because added to RHS) */

# pragma omp parallel if(n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    wflmab[f_id] = -brom[f_id] * bvolfl[f_id];

  cs_divergence(mesh,
                0, /* init */
                wflmas,
                wflmab,
                smbrs);

  cs_field_t *f_divu = cs_field_by_name_try("predicted_vel_divergence");

  if (f_divu != NULL) {
    cs_real_t * cpro_divu = f_divu->val;
    cs_array_real_copy(n_cells, smbrs, cpro_divu);
  }

  /* (Delta t)_ij is calculated as the "viscosity" associated to the pressure */

  cs_face_viscosity(mesh,
                    fvq,
                    1, /* harmonic visc_mean_type */
                    dt,
                    i_visc,
                    b_visc);

  cs_equation_param_t eqp_p_loc = *eqp_p;

  eqp_p_loc.istat  = -1;
  eqp_p_loc.icoupl = -1;
  eqp_p_loc.idifft = -1;
  eqp_p_loc.iwgrec = 0;   /* Warning, may be overwritten if a field */
  eqp_p_loc.blend_st = 0; /* Warning, may be overwritten if a field */

  cs_equation_iterative_solve_scalar(idtvar,
                                     0, /* init */
                                     f_p->id,
                                     NULL,
                                     0,      /* iescap */
                                     0,      /* imucpp */
                                     -1.0,   /* normp */
                                     &eqp_p_loc,
                                     cvar_pr_pre, cvar_pr_pre,
                                     wbfa, wbfb,
                                     cofaf_p, cofbf_p,
                                     wflmas, wflmab,
                                     i_visc, b_visc,
                                     i_visc, b_visc,
                                     NULL,   /* viscel */
                                     NULL, NULL, /* weighf, weighb */
                                     0,      /* icvflb (upwind conv. flux) */
                                     NULL,   /* icvfli */
                                     rovsdt,
                                     smbrs,
                                     cvar_pr, dpvar,
                                     NULL,   /* xcpp */
                                     NULL);  /* eswork */

  /* Printings and clippings
     ----------------------- */

  /* User intervention for a finer management of the bounds and possible
     corrective processing. */

  cs_cf_check_pressure(cvar_pr, n_cells);

  /* Explicit balance (see cs_equation_iterative_solve_scalar:
     the increment is removed) */

  if (eqp_p->verbosity >= 2) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrs[c_id] =    smbrs[c_id]
                    -  eqp_p->istat * (cell_f_vol[c_id] / dt[c_id])
                    * (cvar_pr[c_id] - cvar_pr_pre[c_id])
                    * cs_math_fmax(0, cs_math_fmin(eqp_p->nswrsm-2, 1));
    }
    const cs_real_t sclnor = sqrt(cs_gdot(n_cells, smbrs, smbrs));
    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n PRESSURE: EXPLICIT BALANCE = %14.5e\n"),
                  sclnor);
  }

  /* Communication of P
     ------------------ */

  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, cvar_pr);

  /* Acoustic mass flux computation at the faces
     ------------------------------------------- */

  /* Mass flux = [dt (grad P).n] + [rho (u + dt f)] */

  /* Computation of [dt (grad P).n] by cs_face_diffusion_potential */

  /* festx,y,z    = rvoid
     i_visc, b_visc = arithmetic mean at faces
     viscelx,y,z  = dt
     This flux is stored as the mass flux of the energy */

  cs_face_diffusion_potential(-1,
                              mesh,
                              fvq,
                              1,  /* init */
                              1,  /* inc */
                              eqp_p->imrgra,
                              eqp_p->nswrgr,
                              eqp_p->imligr,
                              0, /* iphydp */
                              0, /* iwgrec */
                              eqp_p->verbosity,
                              eqp_p->epsrgr,
                              eqp_p->climgr,
                              NULL, /* frcxt */
                              cvar_pr,
                              wbfa, wbfb,
                              cofaf_p, cofbf_p,
                              i_visc, b_visc,
                              dt,
                              i_mass_flux_e, b_mass_flux_e);

  /* Incrementation of the flux with [rho (u + dt f)].n = wflmas
     (added with a negative sign since wflmas,wflmab was used above
     in the right hand side). */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    i_mass_flux_e[f_id] = i_mass_flux_e[f_id] - wflmas[f_id];

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    b_mass_flux_e[f_id] = b_mass_flux_e[f_id] - wflmab[f_id];

  /* Updating of the density
     ----------------------- */

  if (cs_glob_velocity_pressure_param->igrdpp > 0) {

    cs_real_t *crom_pre = CS_F_(rho)->val_pre;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Backup of the current density values
         FIXMA: should be done only if iterns = 1 */
      crom_pre[c_id] = crom[c_id];

      /* Update of density values */
      crom[c_id] += (cvar_pr[c_id] - cvar_pr_pre[c_id]) / c2[c_id];
    }

    /* Density communication */

    if (mesh->halo != NULL) {
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, crom);
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, crom_pre);
    }

  }

  BFT_FREE(c2);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(smbrs);
  BFT_FREE(rovsdt);
  BFT_FREE(wflmas);
  BFT_FREE(wflmab);
  BFT_FREE(ivolfl);
  BFT_FREE(bvolfl);
  BFT_FREE(w1);
  BFT_FREE(w7);
  BFT_FREE(w8);
  BFT_FREE(w9);
  BFT_FREE(w10);
  BFT_FREE(dpvar);
  BFT_FREE(wbfa);
  BFT_FREE(wbfb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the constraint for the CFL (compressible algorithm)
 *
 * \param[in]  wcf  compressible constraint
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_cfl_compute(cs_real_t wcf[]) // before : cfdttv
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const int *c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag;

  const int icp = fluid_props->icp;
  const int icv = fluid_props->icv;

  /* Map field arrays */

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;

  cs_real_t *cvar_fracv = NULL;
  cs_real_t *cvar_fracm = NULL;
  cs_real_t *cvar_frace = NULL;

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2){
    cvar_fracv = (cs_real_t *)CS_F_(volume_f)->val;
    cvar_fracm = (cs_real_t *)CS_F_(mass_f)->val;
    cvar_frace = (cs_real_t *)CS_F_(energy_f)->val;
  }

  /* Initialization
     -------------- */

  /* Allocate temporary arrays */
  cs_real_t *i_visc, *i_mass_flux;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(i_mass_flux, n_i_faces, cs_real_t);

  cs_real_t *coefbt, *cofbft, *b_mass_flux, *b_visc;
  BFT_MALLOC(coefbt, n_b_faces, cs_real_t);
  BFT_MALLOC(cofbft, n_b_faces, cs_real_t);
  BFT_MALLOC(b_mass_flux, n_b_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  /* Allocate work arrays */
  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  /* Compute CFL condition associated to the pressure equation
     --------------------------------------------------------- */

  /* Map specific heats fields for sound celerity computation */

  cs_real_t *cpro_cp = NULL, *cpro_cv = NULL;
  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  if (icv >= 0)
    cpro_cv = cs_field_by_id(icv)->val;

  /* Computation of the convective flux associated to the density */
  cs_array_real_set_scalar(n_i_faces, 0.0, i_mass_flux);
  cs_array_real_set_scalar(n_b_faces, 0.0, b_mass_flux);

  int iterns = 1;
  _compressible_pressure_mass_flux(iterns, i_mass_flux, b_mass_flux);

  /* Summation at each cell taking only outward flux */

  int iconvp = 1;
  int idiffp = 0;
  int isym   = 2;

  cs_array_real_set_scalar(n_i_faces, 0.0, i_visc);
  cs_array_real_set_scalar(n_b_faces, 0.0, coefbt);
  cs_array_real_set_scalar(n_b_faces, 0.0, cofbft);
  cs_array_real_set_scalar(n_b_faces, 0.0, b_visc);

  cs_matrix_time_step(mesh,
                      iconvp,
                      idiffp,
                      isym,
                      coefbt,
                      cofbft,
                      i_mass_flux,
                      b_mass_flux,
                      i_visc,
                      b_visc,
                      w1);

  /* Compute the square of the sound celerity */
  cs_real_t *c2;
  BFT_MALLOC(c2, n_cells_ext, cs_real_t);

  cs_cf_thermo_c_square(cpro_cp,
                        cpro_cv,
                        cvar_pr,
                        crom,
                        cvar_fracv,
                        cvar_fracm,
                        cvar_frace,
                        c2,
                        n_cells);

  /* Compute the coefficient CFL/dt */
  cs_real_t psginf = cs_glob_cf_model->psginf;
  if (cs_glob_porous_model >= 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (has_dc * c_disable_flag[has_dc * c_id] != 0)
        wcf[c_id] = cs_math_epzero;
      else
        wcf[c_id] =   w1[c_id] * c2[c_id] * crom[c_id]
                    / ((cvar_pr[c_id] + psginf)*cell_f_vol[c_id]);
    }
  }
  else {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      wcf[c_id] =   w1[c_id] * c2[c_id] * crom[c_id]
                  / ((cvar_pr[c_id] + psginf)*cell_f_vol[c_id]);
  }

  /* Free memory */
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(w1);
  BFT_FREE(c2);
  BFT_FREE(coefbt);
  BFT_FREE(cofbft);
  BFT_FREE(i_mass_flux);
  BFT_FREE(b_mass_flux);
}

/*---------------------------------------------------------------------------- */

END_C_DECLS
