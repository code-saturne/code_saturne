/*============================================================================
 * Compute distance to wall.
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
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_convection_diffusion.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_wall_distance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_distance.c
        Compute distance to wall.
*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool _initialized = false;

static cs_lnum_t n_wall = 0;

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
 * \brief Compute distance to wall by solving a 3d diffusion equation.
 * Solve
 *   \f[ -\divs ( \grad \varia ) = 1 \f]
 * with:
 *  - \f$ \varia_|b = 0 \f$  at the wall
 *  - \f$ \grad \varia \cdot \vect{n} = 0 \f$ elsewhere
 * The wall distance is then equal to:
 *  \f[
 *  d \simeq -|\grad \varia |
 *  + \sqrt{ \grad \varia \cdot \grad \varia +2 \varia }
 *  \f]
 *
 * \param[in]     iterns        iteration number on Navier-Stokes equations
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance(int iterns)
{

  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const int *bc_type = cs_glob_bc_type;
  const cs_halo_t *halo = mesh->halo;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  /* Initialization
     -------------- */

  cs_field_t *f_w_dist = cs_field_by_name("wall_distance");

  cs_equation_param_t *eqp_wd = cs_field_get_equation_param(f_w_dist);

  cs_real_t *wall_dist = f_w_dist->val;

  /* Previous value is stored in a specific field beacause
     the solved field is not directly the wall distance */

  cs_field_t *f_w_dist_aux_pre = cs_field_by_name_try("wall_distance_aux_pre");

  cs_real_t *wall_dist_pre = (f_w_dist_aux_pre != NULL) ?
                              f_w_dist_aux_pre->val:
                              f_w_dist->val_pre;

  /* Boundary conditions
     ------------------- */

  /* Boundary conditions for the resolved scalar T
     Dirichlet to 0 at paroi
     Neumann hmg elsewhere
     We also test for the presence of a Dirichlet */

  int ndircp = 0;

  cs_real_t *coefa_wd = f_w_dist->bc_coeffs->a;
  cs_real_t *coefb_wd = f_w_dist->bc_coeffs->b;
  cs_real_t *cofaf_wd = f_w_dist->bc_coeffs->af;
  cs_real_t *cofbf_wd = f_w_dist->bc_coeffs->bf;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {

      /* Dirichlet boundary condition
         ---------------------------- */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t pimp = 0.0;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_wd[f_id],
                                                  &cofaf_wd[f_id],
                                                  &coefb_wd[f_id],
                                                  &cofbf_wd[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);


      ndircp = ndircp + 1;
    }
    else {

      /* Neumann Boundary Conditions
         --------------------------- */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t qimp = 0.0;

      cs_boundary_conditions_set_neumann_scalar(&coefa_wd[f_id],
                                                &cofaf_wd[f_id],
                                                &coefb_wd[f_id],
                                                &cofbf_wd[f_id],
                                                qimp,
                                                hint);
    }
  }

  if (cs_glob_rank_id > -1)
    cs_parall_sum(1, CS_LNUM_TYPE, &ndircp);

  /* If no wall initialization to a big value */
  if (ndircp == 0) {
    cs_array_real_set_scalar(n_cells, cs_math_big_r, wall_dist);
    return;
  }

  /* Prepare system to solve
     ----------------------- */

  /* Allocate temporary arrays for the species resolution */
  cs_real_t *i_visc, *b_visc, *i_mass_flux, *b_mass_flux;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  BFT_MALLOC(i_mass_flux, n_i_faces, cs_real_t);
  BFT_MALLOC(b_mass_flux, n_b_faces, cs_real_t);

  cs_real_t *dpvar, *smbrp, *rovsdt;
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);
  BFT_MALLOC(smbrp, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  /* Allocate work arrays */
  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  /* Initialize variables to avoid compiler warnings */

  cs_array_real_fill_zero(n_i_faces, i_mass_flux);
  cs_array_real_fill_zero(n_b_faces, b_mass_flux);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rovsdt[c_id] = 0.0; /* Diagonal */
    w1[c_id] = 1.0; /* Diffusion at faces */
  }

  cs_face_viscosity(mesh,
                    fvq,
                    eqp_wd->imvisf,
                    w1,
                    i_visc,
                    b_visc);

  /* Solve system
     ------------ */

  /* Distance to wall is initialized to 0 for reconstruction */

  int nswrsp = eqp_wd->nswrsm;
  int ircflp = eqp_wd->ircflu;

  /* All boundary convective flux with upwind */

  int icvflb = 0;
  cs_real_t normp = -1.0;

  cs_equation_param_t eqp_loc = *eqp_wd;
  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.ndircl = ndircp;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;   /* Warning, may be overwritten if a field */
  eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */

  cs_array_real_fill_zero(n_cells_ext, dpvar);

  /* RHS */
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rovsdt[c_id] = 0.0;
    smbrp[c_id] = cell_f_vol[c_id];
  }

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     f_w_dist->id,
                                     NULL, /* name */
                                     0, /* iescap */
                                     0, /* imucpp */
                                     normp,
                                     &eqp_loc,
                                     wall_dist_pre, wall_dist_pre,
                                     coefa_wd, coefb_wd,
                                     cofaf_wd, cofbf_wd,
                                     i_mass_flux, b_mass_flux,
                                     i_visc, b_visc,
                                     i_visc, b_visc,
                                     NULL, /* viscel */
                                     NULL, /* weightf */
                                     NULL, /* weighb */
                                     icvflb,
                                     NULL, /* icvfli */
                                     rovsdt,
                                     smbrp,
                                     wall_dist, dpvar,
                                     NULL, /* xcpp */
                                     NULL); /* eswork */

  /* Count clippings */
  cs_lnum_t mmprpl = 0;
  cs_real_t dismin = cs_math_big_r;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (wall_dist[c_id] < 0.0) {
      mmprpl = mmprpl + 1;
      dismin = cs_math_fmin(wall_dist[c_id], dismin);
      wall_dist[c_id] = cs_math_epzero * pow(cell_f_vol[c_id], 1.0/3.0);
    }
  }

  if (cs_glob_rank_id > -1) {
    cs_parall_sum(1, CS_LNUM_TYPE, &mmprpl);
    cs_parall_min(1, CS_REAL_TYPE, &dismin);
  }

  /* Recompute wall distance without reconstruction
     (that ensure that it is positive) */

  if (mmprpl >= 1) {
    if (nswrsp > 0) {
      nswrsp = 0;
      ircflp = 0;
      /* Reset also in var_cal_opt structure because some basic routines
         directly use the field options in vcopt... */

      eqp_wd->nswrsm = nswrsp;
      eqp_wd->ircflu = ircflp;

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ @@ WARNING: Wall distance calculation\n"
           "@    =========\n"
           "@  The laplacian solution does not respect the maximum\n"
           "@  principle in %10d cells. We recompute the laplacien\n"
           "@  without reconstructions.\n"), mmprpl);

      /* Reset wall distance */
      cs_array_real_fill_zero(n_cells_ext, wall_dist);

      int n_iter = 0;
      do {

        mmprpl = 0;
        cs_array_real_fill_zero(n_cells_ext, dpvar);

        /* RHS */
# pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          rovsdt[c_id] = 0.0;
          smbrp[c_id] = cell_f_vol[c_id];
        }

        cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                           iterns,
                                           f_w_dist->id,
                                           NULL, /* name */
                                           0, /* iescap */
                                           0, /* imucpp */
                                           normp,
                                           &eqp_loc,
                                           wall_dist_pre, wall_dist_pre,
                                           coefa_wd, coefb_wd,
                                           cofaf_wd, cofbf_wd,
                                           i_mass_flux, b_mass_flux,
                                           i_visc, b_visc,
                                           i_visc, b_visc,
                                           NULL, /* viscel */
                                           NULL, /* weightf */
                                           NULL, /* weighb */
                                           icvflb,
                                           NULL, /* icvfli */
                                           rovsdt,
                                           smbrp,
                                           wall_dist, dpvar,
                                           NULL, /* xcpp */
                                           NULL); /* eswork */

        /* Count clippings */

        cs_real_t _dismin = cs_math_big_r;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          if (wall_dist[c_id] < 0.0) {
            mmprpl += 1;
            _dismin = cs_math_fmin(wall_dist[c_id], _dismin);
            wall_dist[c_id] = cs_math_epzero * pow(cell_f_vol[c_id], 1.0/3.0);
          }
        }

        if (cs_glob_rank_id > -1) {
          cs_parall_sum(1, CS_LNUM_TYPE, &mmprpl);
          cs_parall_min(1, CS_REAL_TYPE, &_dismin);
        }

        n_iter++;
        if (n_iter > 10)
          bft_error(__FILE__, __LINE__, 0,
                    _("Problem for the positivity of wall distance"));

      } while (mmprpl != 0);

    }
    else {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ @@ WARNING: Wall distance calculation\n"
           "@    =========\n"
           "@  The laplacian solution does not respect the maximum\n"
           "@  principle. (laplacian solution is negative : %14.6e)\n"),
         dismin);
    }
  }

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    dpvar[c_id] = cs_math_fmax(wall_dist[c_id], 0.0);

    /* Save working field for the next time step */
    if (f_w_dist_aux_pre != NULL)
      wall_dist_pre[c_id] = wall_dist[c_id];

  }

  if (f_w_dist_aux_pre != NULL)
    cs_halo_sync_var(halo, CS_HALO_STANDARD, wall_dist_pre);

  /* Compute distance to wall
     ------------------------ */

  /* Allocate a temporary array for the gradient calculation */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Compute current gradient */

  cs_field_gradient_scalar(f_w_dist,
                           false,
                           1, /* inc */
                           grad);

  cs_lnum_t counter = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t norm_grad = cs_math_3_dot_product(grad[c_id], grad[c_id]);

    if (norm_grad + 2.0 * dpvar[c_id] >= 0.0)
      wall_dist[c_id] = sqrt(norm_grad + 2.0*dpvar[c_id]) - sqrt(norm_grad);
    else
      counter += 1;

  }

  if (cs_glob_rank_id > -1)
    cs_parall_sum(1, CS_LNUM_TYPE, &counter);

  if (counter > 0)
    cs_log_printf
    (CS_LOG_DEFAULT,
     _("@\n"
       "@ @@ WARNING: Wall distance computation\n"
       "@    =========\n"
       "@  The associated variable does not converge in %10d cells.\n"),
     counter);

  /* Free memory */
  BFT_FREE(grad);

  /* Compute bounds and print info
     ----------------------------- */

  if (cs_glob_rank_id > -1 || mesh->periodicity != NULL)
    cs_halo_sync_var(halo, CS_HALO_EXTENDED, wall_dist);

  cs_real_t _dismax = -cs_math_big_r;
  cs_real_t _dismin =  cs_math_big_r;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    _dismin = cs_math_fmin(wall_dist[c_id], _dismin);
    _dismax = cs_math_fmax(wall_dist[c_id], _dismax);
  }

  if (cs_glob_rank_id > -1)  {
    cs_parall_min(1, CS_REAL_TYPE, &_dismin);
    cs_parall_max(1, CS_REAL_TYPE, &_dismax);
  }

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       " ** WALL DISTANCE\n"
       "    -------------\n\n"
       "  Min distance = %14.5e, Max distance = %14.5e.\n"),
     _dismin, _dismax);

  /* Free memory */
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(dpvar);
  BFT_FREE(smbrp);
  BFT_FREE(rovsdt);
  BFT_FREE(i_mass_flux);
  BFT_FREE(b_mass_flux);
  BFT_FREE(w1);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This subroutine computes the dimensionless distance to the wall
 *        solving a steady transport equation.
 *
 * This function solves the following steady pure convection equation on
 * \f$ \varia \f$:
 * \f[
 * \divs \left( \varia \vect{V} \right)
 *     - \divs \left( \vect{V} \right) \varia = 0
 * \f]
 * where the vector field \f$ \vect{V} \f$ is defined by:
 * \f[
 *  \vect{V} = \dfrac{ \grad y }{\norm{\grad y} }
 * \f]
 * The boundary conditions on \f$ \varia \f$ read:
 * \f[
 *  \varia = \dfrac{u_\star}{\nu} \textrm{ on walls}
 * \f]
 * \f[
 *  \dfrac{\partial \varia}{\partial n} = 0 \textrm{ elsewhere}
 * \f]
 *
 * Then the dimensionless distance is deduced by:
 * \f[
 *  y^+ = y \varia
 * \f]
 *
 * Then, Imposition of an amortization of Van Driest type for the LES.
 *        \f$ \nu_T \f$ is absorbed by \f$ (1-\exp(\dfrac{-y^+}{d^+}))^2 \f$
 *        where \f$ d^+ \f$ is set at 26.
 *
 * \param[in]     visvdr        dynamic viscosity in edge cells after
 *                              driest velocity amortization
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance_yplus(cs_real_t visvdr[])
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_real_t *b_dist = fvq->b_dist;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_real_t *b_face_surf = fvq->b_face_surf;

  const cs_halo_t *halo = mesh->halo;
  const int *bc_type = cs_glob_bc_type;
  cs_lnum_t nt_cur = cs_get_glob_time_step()->nt_cur;

  /* Initialization
     -------------- */

  cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;

  cs_real_t *b_uet = NULL;
  cs_field_t *boundary_ustar = cs_field_by_name_try("boundary_ustar");
  if (boundary_ustar != NULL)
    b_uet = boundary_ustar->val;

  cs_field_t *f_wall_dist = cs_field_by_name("wall_distance");
  cs_real_t *w_dist = f_wall_dist->val;

  cs_real_t *coefa_wd = f_wall_dist->bc_coeffs->a;
  cs_real_t *coefb_wd = f_wall_dist->bc_coeffs->b;
  cs_real_t *cofaf_wd = f_wall_dist->bc_coeffs->af;
  cs_real_t *cofbf_wd = f_wall_dist->bc_coeffs->bf;

  cs_field_t *f_yplus = cs_field_by_name("wall_yplus");
  cs_equation_param_t *eqp_yp = cs_field_get_equation_param(f_yplus);
  cs_real_t *yplus = f_yplus->val;

  cs_real_t *coefa_yp = f_yplus->bc_coeffs->a;
  cs_real_t *coefb_yp = f_yplus->bc_coeffs->b;
  cs_real_t *cofaf_yp = f_yplus->bc_coeffs->af;
  cs_real_t *cofbf_yp = f_yplus->bc_coeffs->bf;

  int iflmas = cs_field_get_key_int(f_yplus,
                                    cs_field_key_id("inner_mass_flux_id"));

  int iflmab = cs_field_get_key_int(f_yplus,
                                    cs_field_key_id("boundary_mass_flux_id"));

  /* Get pointer to the convective mass flux */
  cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;
  cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Number of wall faces */
  if (_initialized == false) {

    _initialized = true;

    n_wall = 0;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL)
        n_wall += 1;
    }

    if (cs_glob_rank_id > -1)
      cs_parall_sum(1, CS_LNUM_TYPE, &n_wall);

  }

  /* If no wall, no wall distance */
  if (n_wall == 0) {
    cs_array_real_set_scalar(n_cells_ext, cs_math_big_r, yplus);
    return;
  }

  /* At the first time step
     ---------------------- */

  /* At the fist time step, in general we have u* = 0 (or false)
     y+ is not computed (this is time-consuming, especially when u* is small,
     because we have to compute y+ up to a large distance from the walls. */

  if (nt_cur == 1) {
    cs_array_real_set_scalar(n_cells, cs_math_big_r, yplus);

    if (eqp_yp->verbosity >= 1) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           " ** DIMENSIONLESS WALL DISTANCE\n"
           "    ---------------------------\n\n"
           "  It is not computed at the first time step\n"));
    }
    return;
  }

  /* Allocate temporary arrays for the distance resolution */
  cs_real_t *dvarp, *smbdp, *rovsdp, *dpvar, *viscap;
  BFT_MALLOC(dvarp, n_cells_ext, cs_real_t);
  BFT_MALLOC(smbdp, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdp, n_cells_ext, cs_real_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscap, n_cells_ext, cs_real_t);

  cs_real_t *i_visc, *b_visc;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  /* Boundary conditions
     ------------------- */

  /* Dirichlet u*./nu at walls, homogeneous Neumann elsewhere */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {
      const cs_lnum_t c_id = b_face_cells[f_id];

      /* Dirichlet Boundary Condition */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t pimp = b_uet[f_id] * crom[c_id] / viscl[c_id];

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_yp[f_id],
                                                  &cofaf_yp[f_id],
                                                  &coefb_yp[f_id],
                                                  &cofbf_yp[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet Boundary Condition */

      const cs_real_t pimp_wd = 0.0;

      cs_boundary_conditions_set_dirichlet_scalar(&coefa_wd[f_id],
                                                  &cofaf_wd[f_id],
                                                  &coefb_wd[f_id],
                                                  &cofbf_wd[f_id],
                                                  pimp_wd,
                                                  hint,
                                                  cs_math_infinite_r);

    }
    else {

      /* Neumann boundary conditions */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t qimp = 0.0;

      cs_boundary_conditions_set_neumann_scalar(&coefa_yp[f_id],
                                                &cofaf_yp[f_id],
                                                &coefb_yp[f_id],
                                                &cofbf_yp[f_id],
                                                qimp,
                                                hint);

      /* Neumann Boundary Conditions */

      cs_boundary_conditions_set_neumann_scalar(&coefa_wd[f_id],
                                                &cofaf_wd[f_id],
                                                &coefb_wd[f_id],
                                                &cofbf_wd[f_id],
                                                qimp,
                                                hint);

    }

  } /* End loop on boundary faces */

  /* Compute the mass flux due to V = Grad(y)
     ---------------------------------------- */

  /* Take Dirichlet into account */
  const int inc    = 1;
  const int iphydp = 0;

  /* Pseudo viscosity, to compute the convective flux "1 grad(y). Sij" */
  cs_array_real_set_scalar(n_cells_ext, 1.0, viscap);

  cs_face_viscosity(mesh,
                    fvq,
                    eqp_yp->imvisf,
                    viscap,
                    i_visc,
                    b_visc);

  /* If the equation on the wall distance has no
     flux-reconstruction (ircflu=0) then no reconstruction
     on the mass-flux (nswrgr) */

  const int nswrgp = (eqp_yp->ircflu == 0) ? 0 : eqp_yp->nswrgr;

  /* Compute convective mass flux
     here -div(1 grad(y)) */

  cs_face_diffusion_potential(f_wall_dist->id,
                              mesh,
                              fvq,
                              1, /* Default initilization at 0 */
                              inc,
                              eqp_yp->imrgra,
                              nswrgp,
                              eqp_yp->imligr,
                              iphydp,
                              0, /* iwgrp */
                              eqp_yp->verbosity,
                              eqp_yp->epsrgr,
                              eqp_yp->climgr,
                              NULL,
                              w_dist,
                              coefa_wd, coefb_wd,
                              cofaf_wd, cofbf_wd,
                              i_visc, b_visc,
                              viscap,
                              i_mass_flux, b_mass_flux);

  /* Now take the opposite */
# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    i_mass_flux[f_id] = - i_mass_flux[f_id];

# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    b_mass_flux[f_id] = - b_mass_flux[f_id];

  /* Diagonal part of the matrix
     --------------------------- */

  cs_array_real_fill_zero(n_cells_ext, rovsdp);

  /* Reinforce diagonal */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_lnum_t c_id0 = i_face_cells[f_id][0];
    const cs_lnum_t c_id1 = i_face_cells[f_id][1];

    rovsdp[c_id0] +=   i_mass_flux[f_id];
    rovsdp[c_id1] += - i_mass_flux[f_id];
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    const cs_lnum_t c_id = b_face_cells[f_id];

    rovsdp[c_id] += b_mass_flux[f_id];
  }

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    rovsdp[c_id] = 1.e-6 * fabs(rovsdp[c_id]);

  if (halo != NULL)
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rovsdp);

  /* Time loop
     --------- */

  /* Initialization */

  /* Unknown
     In the case where the stationary state is not completely reached,
     yplus must not be zero in the zone where the boundary conditions
     have not been convected. Instead, we want Y to be maximum.
     If you use zero or a negative value as initialization,
     you risk ending up with values close to zero resulting from
     diffusion due to the upwind scheme in the vicinity of the convected
     front, and therefore with yplus values close to zero anywhere.
     We will therefore use the maximum value of u*./nu. */

     /* From the second time step, we also have the yplus of
        the previous time step */

  /* Compute the min and max */

  cs_real_t xusnmx = -cs_math_big_r;
  cs_real_t xusnmn =  cs_math_big_r;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    if (   bc_type[f_id] == CS_SMOOTHWALL
        || bc_type[f_id] == CS_ROUGHWALL) {

      xusnmx = cs_math_fmax(xusnmx, coefa_yp[f_id]);
      xusnmn = cs_math_fmin(xusnmn, coefa_yp[f_id]);
    }
  }

  if (cs_glob_rank_id > -1) {
    cs_parall_max(1, CS_REAL_TYPE, &xusnmx);
    cs_parall_min(1, CS_REAL_TYPE, &xusnmn);
  }

  if (nt_cur == 1)
    cs_array_real_set_scalar(n_cells_ext, xusnmx, dvarp);
  else {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t usna = yplus[c_id] / cs_math_fmax(w_dist[c_id],
                                                  cs_math_epzero);
      usna = cs_math_fmax(usna, xusnmn);
      usna = cs_math_fmin(usna, xusnmx);

      dvarp[c_id] = usna;
    }
  }

  /* L2 norm of (u* / nu) over wall boundary faces */

  cs_real_t xnorm0 = 0.0;
  cs_real_t wall_surf = 0.0;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    if (   bc_type[f_id] == CS_SMOOTHWALL
        || bc_type[f_id] == CS_ROUGHWALL) {

      wall_surf += b_face_surf[f_id];
      xnorm0 += cs_math_pow2(coefa_yp[f_id]) * b_face_surf[f_id];
    }
  }

  if (cs_glob_rank_id > -1) {
    cs_real_t sum[2] = {xnorm0, wall_surf};
    cs_parall_sum(2, CS_REAL_TYPE, sum);
    xnorm0 = sum[0];
    wall_surf = sum[1];
  }

  xnorm0 = sqrt(xnorm0 / wall_surf) * fvq->tot_vol;

  if (cs_glob_rank_id > -1 || mesh->periodicity != NULL)
    cs_halo_sync_var(halo, CS_HALO_STANDARD, dvarp);

  /* Right hand side
     --------------- */

  cs_array_real_fill_zero(n_cells, smbdp);

  /* Solving
     ------- */

  /* All boundary convective flux with upwind */
  int icvflb = 0;

  /* There are some Dirichlet BCs */
  int ndircp = 1;

  /* Warning: no diffusion so no need of other diffusive
     boundary coefficient */

  cs_equation_param_t eqp_loc = *eqp_yp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.ndircl = ndircp;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0; /* Warning, may be overwritten if a field */
  eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */

  cs_equation_iterative_solve_scalar(0, /* No steady state algo */
                                     -1, /* No over loops */
                                     f_yplus->id,
                                     NULL, /* name */
                                     0, /* No error estimate */
                                     0, /* No error estimate */
                                     xnorm0,
                                     &eqp_loc,
                                     dvarp,
                                     dvarp,
                                     coefa_yp, coefb_yp,
                                     cofaf_yp, cofbf_yp,
                                     i_mass_flux, b_mass_flux,
                                     i_mass_flux, b_mass_flux,
                                     i_mass_flux, b_mass_flux,
                                     NULL, NULL,
                                     NULL,
                                     icvflb,
                                     NULL,
                                     rovsdp,
                                     smbdp,
                                     dvarp, dpvar,
                                     NULL, NULL);

  /* Warning: no diffusion so no need of other diffusive
     boundary coefficient */

  /* Finalization and printing
     ------------------------- */

  cs_real_t dismax = - cs_math_big_r;
  cs_real_t dismin =   cs_math_big_r;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Clipping: essential if you initialize with (u * /nu) of
       the previous time step */

    dvarp[c_id] = cs_math_fmax(dvarp[c_id], xusnmn);
    dvarp[c_id] = cs_math_fmin(dvarp[c_id], xusnmx);

    yplus[c_id] = dvarp[c_id] * w_dist[c_id];

    dismin = cs_math_fmin(yplus[c_id], dismin);
    dismax = cs_math_fmax(yplus[c_id], dismax);

  }

  if (cs_glob_rank_id > -1) {
    cs_parall_min(1, CS_REAL_TYPE, &dismin);
    cs_parall_max(1, CS_REAL_TYPE, &dismax);
  }

  if (eqp_yp->verbosity >= 1) {

    cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       " ** DIMENSIONLESS WALL DISTANCE\n"
       "    ---------------------------\n\n"
       "  Min distance+ = %14.5e, Max distance+ = %14.5e.\n"),
     dismin, dismax);
  }

  /* Van Driest amortization
     ----------------------- */

  cs_real_t *visct = CS_F_(mu_t)->val;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    visct[c_id] *= cs_math_pow2(1.0 - exp(-yplus[c_id] / cs_turb_cdries));

    /* For the wall cells we add the turbulent viscosity which was absorbed
       in clptur and which has served to calculate the boundary conditions */

    if (visvdr[c_id] > -900.0)
      visct[c_id] = visvdr[c_id];
  }

  /* Free memory */
  BFT_FREE(dvarp);
  BFT_FREE(smbdp);
  BFT_FREE(rovsdp);
  BFT_FREE(dpvar);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
  BFT_FREE(viscap);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
