/*============================================================================
 * Weakly compressible algorithm (semi-analytic):
   Computation of scalar diffusion terms.
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
#include "cs_face_viscosity.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_physical_constants.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_dilatable_scalar_diff_st.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_dilatable_scalar_diff_st.c
        Weakly compressible algorithm (semi-analytic):
        Computation of scalar diffusion terms.
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
 * \brief Weakly compressible algorithm (semi-analytic):
          Computation of scalar diffusion terms.
 *
 * \param[in]  iterns      Navier-Stokes iteration number
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_dilatable_scalar_diff_st(int iterns)

{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const int icp = fluid_props->icp;
  const cs_real_t rair = fluid_props->r_pg_cnst;
  const cs_real_t cp0 = fluid_props->cp0;

  const int n_fields = cs_field_n_fields();
  const int keysca  = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int kscacp  = cs_field_key_id("is_temperature");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");

  /* Memory allocation */
  cs_real_t *vistot, *i_visc, *b_visc, *xcpp;
  BFT_MALLOC(vistot, n_cells_ext, cs_real_t);

  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xcpp, n_cells_ext, cs_real_t, cs_alloc_mode);

  const cs_real_t *cpro_cp = nullptr;
  if (icp >= 0)
    cpro_cp = CS_F_(cp)->val;

  /* Index for turbulent diffusivity */
  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Loop on scalars */
  for (int ii = 0; ii < n_fields; ii++) {

    cs_field_t *f_scal = cs_field_by_id(ii);

    if (!(f_scal->type & CS_FIELD_VARIABLE))
      continue;
    if (cs_field_get_key_int(f_scal, keysca) <= 0)
      continue;

    cs_real_t *cvar_scal = f_scal->val;

    /* Key id for buoyant field (inside the Navier Stokes loop) */
    const int key_coupled_with_vel_p = cs_field_key_id_try("coupled_with_vel_p");
    const int coupled_with_vel_p_fld = cs_field_get_key_int(f_scal, key_coupled_with_vel_p);

    if (   (coupled_with_vel_p_fld == 1 && iterns == -1)
        || (coupled_with_vel_p_fld == 0 && iterns != -1))
      continue;

    const int iscavr = cs_field_get_key_int(f_scal, kscavr);

    const int iscacp = (iscavr > 0) ?
      cs_field_get_key_int(cs_field_by_id(iscavr), kscacp):
      cs_field_get_key_int(f_scal, kscacp);

    const int imucpp = (iscacp == 1) ? 1 : 0;

    if (imucpp == 0)
      cs_array_real_set_scalar(n_cells, 1.0, xcpp);
    else if (imucpp == 1) { // TODO: humid air

      if (icp >= 0) {
        if (iscacp == 1)
          cs_array_real_copy(n_cells, cpro_cp, xcpp);
        else if (iscacp == 2) {
#         pragma omp parallel for if (n_cells > CS_THR_MIN)
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            xcpp[c_id] = cpro_cp[c_id] - rair;
        }
      }
      else {
        if (iscacp == 1)
          cs_array_real_set_scalar(n_cells, cp0, xcpp);
        else if (iscacp == 2)
          cs_array_real_set_scalar(n_cells, cp0-rair, xcpp);
      }
    }

    /* Handle parallelism and periodicity */
    if (mesh->halo != nullptr) {
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, xcpp);
    }

    cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_scal);

    /* Pointers to the mass fluxes */
    const int iflmas
      = cs_field_get_key_int(CS_F_(vel), cs_field_key_id("inner_mass_flux_id"));
    const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

    const int iflmab
      = cs_field_get_key_int(CS_F_(vel),
                             cs_field_key_id("boundary_mass_flux_id"));
    const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

    /* Diffusion velocity */

    /* Index for molecular diffusivity */
    const int ifcvsl = cs_field_get_key_int(f_scal, kivisl);
    cs_real_t *viscls = nullptr;
    if (ifcvsl > -1)
      viscls = cs_field_by_id(ifcvsl)->val;

    if (eqp_sc->idiff >= 1) {

      /* Only the positive part of mu_t is considered (MAX(mu_t,0)),
         Dynamic LES can cause negative mu_t (clipping on (mu+mu_t))
         The positive part of (K+K_t) would have been considered but
         should allow negative K_t that is considered non physical here */

      const cs_real_t turb_schmidt = cs_field_get_key_double(f_scal, ksigmas);

      if (ifcvsl < 0) {
        const cs_real_t visls_0 = cs_field_get_key_double(f_scal, kvisl0);
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          vistot[c_id] = visls_0
            + eqp_sc->idifft*xcpp[c_id]*cs_math_fmax(visct[c_id], 0.)
            / turb_schmidt;
      }
      else {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          vistot[c_id] = viscls[c_id]
            + eqp_sc->idifft*xcpp[c_id]*cs_math_fmax(visct[c_id], 0.)
            / turb_schmidt;
      }

      cs_face_viscosity(mesh,
                        fvq,
                        eqp_sc->imvisf,
                        vistot,
                        i_visc,
                        b_visc);

    }
    else {
      cs_array_real_fill_zero(n_i_faces, i_visc);
      cs_array_real_fill_zero(n_b_faces, b_visc);
      cs_array_real_fill_zero(n_cells, vistot);
    }

    /* Source term */
    char scalar_st_name[50];
    sprintf(scalar_st_name, "%s_%s", f_scal->name, "dila_st");
    cs_real_t *cpro_tsscal = cs_field_by_name(scalar_st_name)->val;

    /* Diffusion term calculation */
    const cs_field_bc_coeffs_t *bc_coeffs_sc = f_scal->bc_coeffs;

    cs_equation_param_t eqp_loc = *eqp_sc;

    /*  All boundary convective flux with upwind */

    eqp_loc.iconv  = 0; /* diffusion term only */
    eqp_loc.istat  = -1;
    eqp_loc.idiff  = 1;
    eqp_loc.idifft = -1;
    eqp_loc.idften = CS_ISOTROPIC_DIFFUSION; /* FIXME when activating GGDH */
    eqp_loc.iswdyn = -1;
    eqp_loc.ischcv = 1;
    eqp_loc.isstpc = 1;
    eqp_loc.nswrsm = -1;
    eqp_loc.iwgrec = 0;
    eqp_loc.theta = 1.0;
    eqp_loc.blencv = 0.0;
    eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */
    eqp_loc.epsilo = -1;
    eqp_loc.epsrsm = -1;

    cs_balance_scalar(cs_glob_time_step_options->idtvar,
                      -1,    /* f_id */
                      imucpp,
                      0,     /* imasac */
                      1,     /* inc */
                      &eqp_loc,
                      cvar_scal, cvar_scal,
                      bc_coeffs_sc,
                      i_mass_flux, b_mass_flux,
                      i_visc, b_visc,
                      nullptr,  /* viscel */
                      xcpp,
                      nullptr,  /* weighf */
                      nullptr,  /* weighb */
                      0,     /* icvflb; upwind scheme */
                      nullptr,
                      cpro_tsscal);
  } /* end loop on fields */

  /* Free memory */
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);
  CS_FREE_HD(xcpp);

  BFT_FREE(vistot);

}

/*---------------------------------------------------------------------------- */

END_C_DECLS
