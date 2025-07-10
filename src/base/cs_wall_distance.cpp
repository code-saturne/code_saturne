/*============================================================================
 * Compute distance to wall.
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

#include <float.h>
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

#include "bft/bft_printf.h"
#include "base/cs_array.h"
#include "base/cs_dispatch.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "alge/cs_convection_diffusion.h"
#include "base/cs_dispatch.h"
#include "base/cs_equation_iterative_solve.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_reducers.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_wall_distance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_distance.cpp
        Compute distance to wall.
*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool _initialized = false;

static cs_lnum_t n_wall = 0;

/* Fortran mapping int values for wall distance. These should eventually
 * be removed in the future.
 */

static cs_wall_distance_options_t _wall_distance_options = {
  .need_compute = 0,
  .is_up_to_date = 0,
  .method = 1
};

const cs_wall_distance_options_t *cs_glob_wall_distance_options
= &_wall_distance_options;

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
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const int *bc_type = cs_glob_bc_type;
  const cs_halo_t *halo = mesh->halo;
  const cs_real_t *b_dist = mq->b_dist;
  const cs_real_t *cell_f_vol = mq->cell_vol;

  cs_dispatch_context ctx;

  /* Initialization
     -------------- */

  cs_field_t *f_w_dist = cs_field_by_name("wall_distance");

  cs_equation_param_t *eqp_wd = cs_field_get_equation_param(f_w_dist);

  cs_real_t *wall_dist = f_w_dist->val;

  /* Previous value is stored in a specific field because
     the solved field is not directly the wall distance */

  cs_field_t *f_w_dist_aux_pre = cs_field_by_name_try("wall_distance_aux_pre");

  cs_real_t *wall_dist_pre = (f_w_dist_aux_pre != nullptr) ?
                              f_w_dist_aux_pre->val:
                              f_w_dist->val_pre;

  /* Boundary conditions
     ------------------- */

  /* Boundary conditions for the resolved scalar T
     Dirichlet to 0 at wall
     Neumann hmg elsewhere
     We also test for the presence of a Dirichlet */

  int *info;
  CS_MALLOC_HD(info, 2, int, cs_alloc_mode);
  info[0] = 0;
  info[1] = 1;

  cs_field_bc_coeffs_t *bc_coeffs_wd = f_w_dist->bc_coeffs;
  cs_real_t *coefa_wd = bc_coeffs_wd->a;
  cs_real_t *coefb_wd = bc_coeffs_wd->b;
  cs_real_t *cofaf_wd = bc_coeffs_wd->af;
  cs_real_t *cofbf_wd = bc_coeffs_wd->bf;

  /* Fixed mesh: update only if BC's have changed (i.e. in restart) */

  if (mesh->time_dep == CS_MESH_FIXED) {

    /* Check if they are difference between coefa_pre and coefa for wall_dist */
    info[1] = 0;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_real_t a_prev = coefa_wd[f_id];
      cs_real_t b_prev = coefb_wd[f_id];

      if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {

        const cs_real_t hint = 1.0 / b_dist[f_id];
        const cs_real_t pimp = 0.0;

        cs_boundary_conditions_set_dirichlet_scalar(coefa_wd[f_id],
                                                    cofaf_wd[f_id],
                                                    coefb_wd[f_id],
                                                    cofbf_wd[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);

        info[0] = 1;
      }
      else
        cs_boundary_conditions_set_neumann_scalar_hmg(coefa_wd[f_id],
                                                      cofaf_wd[f_id],
                                                      coefb_wd[f_id],
                                                      cofbf_wd[f_id]);

      cs_real_t d =   cs::abs(a_prev - coefa_wd[f_id])
                    + cs::abs(b_prev - coefb_wd[f_id]);

      if (d > cs_math_epzero)
        info[1] = 1;
    });

    ctx.wait();

    cs_parall_max(1, CS_INT_TYPE, &info[1]);
  }

  /* Time evolving mesh: always update */

  else {
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {
        const cs_real_t hint = 1.0 / b_dist[f_id];
        const cs_real_t pimp = 0.0;
        cs_boundary_conditions_set_dirichlet_scalar(coefa_wd[f_id],
                                                    cofaf_wd[f_id],
                                                    coefb_wd[f_id],
                                                    cofbf_wd[f_id],
                                                    pimp,
                                                    hint,
                                                    cs_math_infinite_r);
        info[0] = 1;
      }
      else
        cs_boundary_conditions_set_neumann_scalar_hmg(coefa_wd[f_id],
                                                      cofaf_wd[f_id],
                                                      coefb_wd[f_id],
                                                      cofbf_wd[f_id]);
    });
  }

  ctx.wait();

  /* Check that wall distance if initialized/read if no diff in BC's (for
     the strange case where BC's would be read but nut the wall distance) */
  if (info[1] == 0) {
    double d = cs_dot_xx(n_cells, wall_dist);
    if (d <= 0)
      info[1] = 1;
  }

  cs_lnum_t c[2] = {info[0], info[1]};
  cs_parall_sum(2, CS_LNUM_TYPE, c);
  int ndircp = c[0];
  int have_diff = c[1];

  /* Without wall or if value already computed (i.e. same BC's), return */

  if (ndircp == 0 || have_diff == 0) {

    /* If no wall, initialization to a big value */
    if (ndircp == 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        wall_dist[c_id] = cs_math_big_r;
      });
    }

    ctx.wait();
    CS_FREE_HD(info);
    return;
  }

  /* Prepare system to solve
     ----------------------- */

  cs_real_t *smbrp, *rovsdt;
  CS_MALLOC_HD(rovsdt, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbrp, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Allocate temporary arrays for the species resolution */
  cs_real_t *dpvar, *i_visc, *b_visc, *i_mass_flux, *b_mass_flux;
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(i_mass_flux, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_mass_flux, n_b_faces, cs_real_t, cs_alloc_mode);

  /* Allocate work arrays */
  cs_real_t *w1;
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Initialize variables to avoid compiler warnings */

  ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    i_mass_flux[f_id] = 0.;
  });

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    b_mass_flux[f_id] = 0.;
  });

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    w1[c_id] = 1.0; /* Diffusion at faces */
    dpvar[c_id] = 0.;

     /* RHS */
    rovsdt[c_id] = 0.0;
    smbrp[c_id] = cell_f_vol[c_id];
  });

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    dpvar[c_id] = 0.;
  });

  ctx.wait();

  cs_face_viscosity(mesh,
                    mq,
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

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     f_w_dist->id,
                                     nullptr, /* name */
                                     0, /* iescap */
                                     0, /* imucpp */
                                     normp,
                                     &eqp_loc,
                                     wall_dist_pre, wall_dist_pre,
                                     bc_coeffs_wd,
                                     i_mass_flux, b_mass_flux,
                                     i_visc, b_visc,
                                     i_visc, b_visc,
                                     nullptr, /* viscel */
                                     nullptr, /* weightf */
                                     nullptr, /* weighb */
                                     icvflb,
                                     nullptr, /* icvfli */
                                     rovsdt,
                                     smbrp,
                                     wall_dist, dpvar,
                                     nullptr, /* xcpp */
                                     nullptr); /* eswork */

  /* Count clippings */

  struct cs_data_1int_1float rd;
  struct cs_reduce_sum1i_min1float reducer;

  ctx.parallel_for_reduce
    (n_cells, rd, reducer,
     [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_1int_1float &res) {

    res.i[0] = 0;
    res.r[0] = cs_math_infinite_r;

    if (wall_dist[c_id] < 0.0) {
      res.i[0] = 1;
      res.r[0] = wall_dist[c_id];
      wall_dist[c_id] = cs_math_epzero * pow(cell_f_vol[c_id], 1.0/3.0);
    }
  });

  ctx.wait();

  cs_lnum_t mmprpl = rd.i[0];
  cs_real_t dismin = rd.r[0];

  cs_parall_sum(1, CS_LNUM_TYPE, &mmprpl);
  cs_parall_min(1, CS_REAL_TYPE, &dismin);

  /* Recompute wall distance without reconstruction
     (that ensure that it is positive) */

  if (mmprpl >= 1) {
    if (nswrsp > 0) {
      nswrsp = 0;
      ircflp = 0;
      /* Reset also in equation param structure because some basic routines
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
      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        wall_dist[c_id] = 0.0;
      });

      int n_iter = 0;
      do {

        mmprpl = 0;
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          dpvar[c_id] = 0.0;

          /* RHS */
          rovsdt[c_id] = 0.0;
          smbrp[c_id] = cell_f_vol[c_id];
        });

        ctx.wait();

        cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                           iterns,
                                           f_w_dist->id,
                                           nullptr, /* name */
                                           0, /* iescap */
                                           0, /* imucpp */
                                           normp,
                                           &eqp_loc,
                                           wall_dist_pre, wall_dist_pre,
                                           bc_coeffs_wd,
                                           i_mass_flux, b_mass_flux,
                                           i_visc, b_visc,
                                           i_visc, b_visc,
                                           nullptr, /* viscel */
                                           nullptr, /* weightf */
                                           nullptr, /* weighb */
                                           icvflb,
                                           nullptr, /* icvfli */
                                           rovsdt,
                                           smbrp,
                                           wall_dist, dpvar,
                                           nullptr, /* xcpp */
                                           nullptr); /* eswork */

        /* Count clippings */

        struct cs_data_1int_1float rd_1i_1r;
        struct cs_reduce_sum1i_min1float reducer_1i_1r;

        ctx.parallel_for_reduce
          (n_cells, rd_1i_1r, reducer_1i_1r,
           [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_1int_1float &res) {

          res.i[0] = 0;
          res.r[0] = cs_math_infinite_r;

          if (wall_dist[c_id] < 0.0) {
            res.i[0] = 1;
            res.r[0] = wall_dist[c_id];
            wall_dist[c_id] = cs_math_epzero * pow(cell_f_vol[c_id], 1.0/3.0);
          }
        });

        ctx.wait();

        mmprpl = rd_1i_1r.i[0];
        cs_parall_sum(1, CS_LNUM_TYPE, &mmprpl);
        cs_parall_min(1, CS_REAL_TYPE, &rd_1i_1r.r[0]);

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
           "@    =======\n"
           "@  The laplacian solution does not respect the maximum\n"
           "@  principle. (laplacian solution is negative: %14.6e)\n"),
         dismin);
    }
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    dpvar[c_id] = cs::max(wall_dist[c_id], 0.0);

    /* Save working field for the next time step */
    if (f_w_dist_aux_pre != nullptr)
      wall_dist_pre[c_id] = wall_dist[c_id];
  });

  ctx.wait();

  if (f_w_dist_aux_pre != nullptr)
    cs_halo_sync_var(halo, CS_HALO_STANDARD, wall_dist_pre);

  /* Compute distance to wall
     ------------------------ */

  /* Allocate a temporary array for the gradient calculation */
  cs_real_3_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Compute current gradient */

  cs_field_gradient_scalar(f_w_dist,
                           false,
                           1, /* inc */
                           grad);

  struct cs_data_1int_2float rd_1i_2r;
  struct cs_reduce_min1float_max1float_sum1int reducer_1i_2r;

  ctx.parallel_for_reduce
    (n_cells, rd_1i_2r, reducer_1i_2r,
     [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_1int_2float &res) {

    res.i[0] = 0;

    const cs_real_t norm_grad = cs_math_3_dot_product(grad[c_id], grad[c_id]);

    if (norm_grad + 2.0 * dpvar[c_id] >= 0.0)
      wall_dist[c_id] = sqrt(norm_grad + 2.0*dpvar[c_id]) - sqrt(norm_grad);
    else
      res.i[0] = 1;

    res.r[0] = wall_dist[c_id];
    res.r[1] = wall_dist[c_id];
  });

  ctx.wait();

  /* Print info
     ---------- */

  cs_halo_sync_var(halo, CS_HALO_EXTENDED, wall_dist);

  cs_parall_sum(1, CS_LNUM_TYPE, &rd_1i_2r.i[0]);

  if (rd_1i_2r.i[0] > 0)
    cs_log_printf
    (CS_LOG_DEFAULT,
     _("@\n"
       "@ @@ WARNING: Wall distance computation\n"
       "@    =========\n"
       "@  The associated variable does not converge in %10d cells.\n"),
     rd_1i_2r.i[0]);

  cs_real_t _dismin = rd_1i_2r.r[0];
  cs_real_t _dismax = rd_1i_2r.r[1];
  cs_parall_min(1, CS_REAL_TYPE, &_dismin);
  cs_parall_max(1, CS_REAL_TYPE, &_dismax);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       " ** WALL DISTANCE\n"
       "    -------------\n\n"
       "  Min distance = %14.5e, Max distance = %14.5e.\n"),
     _dismin, _dismax);

  /* Free memory */
  CS_FREE_HD(grad);
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);
  CS_FREE_HD(dpvar);
  CS_FREE_HD(smbrp);
  CS_FREE_HD(i_mass_flux);
  CS_FREE_HD(b_mass_flux);
  CS_FREE_HD(rovsdt);
  CS_FREE_HD(w1);
  CS_FREE_HD(info);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dimensionless distance to the wall
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
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_real_t *b_dist = mq->b_dist;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_real_t *b_face_surf = mq->b_face_surf;

  const cs_halo_t *halo = mesh->halo;
  const int *bc_type = cs_glob_bc_type;
  cs_lnum_t nt_cur = cs_get_glob_time_step()->nt_cur;

  /* Initialization
     -------------- */

  cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;

  cs_real_t *b_uet = nullptr;
  cs_field_t *boundary_ustar = cs_field_by_name_try("boundary_ustar");
  if (boundary_ustar != nullptr)
    b_uet = boundary_ustar->val;

  cs_field_t *f_wall_dist = cs_field_by_name("wall_distance");
  cs_real_t *w_dist = f_wall_dist->val;
  cs_field_bc_coeffs_t *bc_coeffs_wd = f_wall_dist->bc_coeffs;
  cs_real_t *coefa_wd = bc_coeffs_wd->a;
  cs_real_t *coefb_wd = bc_coeffs_wd->b;
  cs_real_t *cofaf_wd = bc_coeffs_wd->af;
  cs_real_t *cofbf_wd = bc_coeffs_wd->bf;

  cs_field_t *f_yplus = cs_field_by_name("wall_yplus");
  cs_equation_param_t *eqp_yp = cs_field_get_equation_param(f_yplus);
  cs_real_t *yplus = f_yplus->val;

  cs_field_bc_coeffs_t *bc_coeffs_yp = f_yplus->bc_coeffs;
  cs_real_t *coefa_yp = bc_coeffs_yp->a;
  cs_real_t *coefb_yp = bc_coeffs_yp->b;
  cs_real_t *cofaf_yp = bc_coeffs_yp->af;
  cs_real_t *cofbf_yp = bc_coeffs_yp->bf;

  int iflmas = cs_field_get_key_int(f_yplus,
                                    cs_field_key_id("inner_mass_flux_id"));

  int iflmab = cs_field_get_key_int(f_yplus,
                                    cs_field_key_id("boundary_mass_flux_id"));

  /* Get pointer to the convective mass flux */
  cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;
  cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  cs_dispatch_context ctx;

  /* Number of wall faces */
  if (_initialized == false) {
    _initialized = true;

    n_wall = 0;
    ctx.parallel_for_reduce_sum
      (n_b_faces, n_wall, [=] CS_F_HOST_DEVICE
       (cs_lnum_t f_id, CS_DISPATCH_REDUCER_TYPE(cs_lnum_t) &sum) {
      if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL)
        sum += 1;
    });
    ctx.wait();
    cs_parall_sum(1, CS_LNUM_TYPE, &n_wall);
  }

  /* If no wall, no wall distance */
  if (n_wall == 0) {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      yplus[c_id] = cs_math_big_r;
    });
    ctx.wait();
    return;
  }

  /* At the first time step
     ---------------------- */

  /* At the fist time step, in general we have u* = 0 (or false)
     y+ is not computed (this is time-consuming, especially when u* is small,
     because we have to compute y+ up to a large distance from the walls. */

  if (nt_cur == 1) {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      yplus[c_id] = cs_math_big_r;
    });

    ctx.wait();

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
  CS_MALLOC_HD(dvarp, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbdp, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(rovsdp, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscap, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *i_visc, *b_visc;
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  /* Boundary conditions
     ------------------- */

  /* Dirichlet u*./nu at walls, homogeneous Neumann elsewhere */

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {
      const cs_lnum_t c_id = b_face_cells[f_id];

      /* Dirichlet boundary condition */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t pimp = b_uet[f_id] * crom[c_id] / viscl[c_id];

      cs_boundary_conditions_set_dirichlet_scalar(coefa_yp[f_id],
                                                  cofaf_yp[f_id],
                                                  coefb_yp[f_id],
                                                  cofbf_yp[f_id],
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);

      /* Dirichlet boundary condition */

      const cs_real_t pimp_wd = 0.0;

      cs_boundary_conditions_set_dirichlet_scalar(coefa_wd[f_id],
                                                  cofaf_wd[f_id],
                                                  coefb_wd[f_id],
                                                  cofbf_wd[f_id],
                                                  pimp_wd,
                                                  hint,
                                                  cs_math_infinite_r);

    }
    else {

      /* Neumann boundary conditions */

      const cs_real_t hint = 1.0 / b_dist[f_id];
      const cs_real_t qimp = 0.0;

      cs_boundary_conditions_set_neumann_scalar(coefa_yp[f_id],
                                                cofaf_yp[f_id],
                                                coefb_yp[f_id],
                                                cofbf_yp[f_id],
                                                qimp,
                                                hint);

      /* Neumann Boundary Conditions */

      cs_boundary_conditions_set_neumann_scalar(coefa_wd[f_id],
                                                cofaf_wd[f_id],
                                                coefb_wd[f_id],
                                                cofbf_wd[f_id],
                                                qimp,
                                                hint);

    }

  }); /* End loop on boundary faces */

  /* Compute the mass flux due to V = Grad(y)
     ---------------------------------------- */

  /* Take Dirichlet into account */
  const int inc = 1;
  const int iphydp = 0;

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    /* Pseudo viscosity, to compute the convective flux "1 grad(y). Sij" */
    viscap[c_id] = 1.0;

    rovsdp[c_id] = 0.;
  });

  ctx.wait();

  cs_face_viscosity(mesh,
                    mq,
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
                              mq,
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
                              nullptr,
                              w_dist,
                              f_wall_dist->bc_coeffs,
                              i_visc, b_visc,
                              viscap,
                              i_mass_flux, b_mass_flux);

  /* Diagonal part of the matrix
     --------------------------- */

  /* Reinforce diagonal */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    /* Now take the opposite */
    i_mass_flux[f_id] = - i_mass_flux[f_id];

    const cs_lnum_t c_id0 = i_face_cells[f_id][0];
    const cs_lnum_t c_id1 = i_face_cells[f_id][1];

    rovsdp[c_id0] +=   i_mass_flux[f_id];
    rovsdp[c_id1] += - i_mass_flux[f_id];
  }

  /* Unknown
     In the case where the steady state is not completely reached,
     yplus must not be zero in the zone where the boundary conditions
     have not been convected. Instead, we want Y to be maximum.
     If you use zero or a negative value as initialization,
     you risk ending up with values close to zero resulting from
     diffusion due to the upwind scheme in the vicinity of the convected
     front, and therefore with yplus values close to zero anywhere.
     We will therefore use the maximum value of u*./nu. */

  /* From the second time step, we also have the yplus of
     the previous time step */

  struct cs_data_double_n<4> rd;
  struct cs_reduce_min1r_max1r_sum2r reducer;

  ctx.parallel_for_reduce
    (n_b_faces, rd, reducer,
     [=] CS_F_HOST_DEVICE (cs_lnum_t f_id, cs_data_double_n<4> &res) {
    /* Now take the opposite */
    b_mass_flux[f_id] = - b_mass_flux[f_id];

    const cs_lnum_t c_id = b_face_cells[f_id];

    rovsdp[c_id] += b_mass_flux[f_id];

    res.r[0] =  cs_math_infinite_r;
    res.r[1] = -cs_math_infinite_r;
    res.r[2] = 0.;
    res.r[3] = 0.;

    if (   bc_type[f_id] == CS_SMOOTHWALL
        || bc_type[f_id] == CS_ROUGHWALL) {

      /* Compute the min and max */

      res.r[0] = coefa_yp[f_id];
      res.r[1] = coefa_yp[f_id];

      /* L2 norm of (u* / nu) over wall boundary faces */

      res.r[2] = cs_math_pow2(coefa_yp[f_id]) * b_face_surf[f_id];
      res.r[3] = b_face_surf[f_id];
    }
  });

  ctx.wait();

  cs_parall_min(1, CS_DOUBLE, &rd.r[0]);
  cs_parall_max(1, CS_DOUBLE, &rd.r[1]);

  const cs_real_t xusnmn = rd.r[0];
  const cs_real_t xusnmx = rd.r[1];

  cs_real_t sum[2] = {rd.r[2], rd.r[3]};
  cs_parall_sum(2, CS_DOUBLE, sum);
  cs_real_t xnorm0 = sum[0];
  cs_real_t wall_surf = sum[1];
  xnorm0 = sqrt(xnorm0 / wall_surf) * mq->tot_vol;

  /* Time loop
     --------- */

  /* Initialization */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    rovsdp[c_id] = 1.e-6 * cs::abs(rovsdp[c_id]);

    /* Right hand side */
    smbdp[c_id] = 0.;

    if (nt_cur == 1) {
      dvarp[c_id] = xusnmx;
    }
    else {
      cs_real_t usna = yplus[c_id] / cs::max(w_dist[c_id],
                                             cs_math_epzero);
      usna = cs::max(usna, xusnmn);
      usna = cs::min(usna, xusnmx);

      dvarp[c_id] = usna;
    }
  });

  ctx.wait();

  cs_halo_sync_var(halo, CS_HALO_STANDARD, rovsdp);
  cs_halo_sync_var(halo, CS_HALO_STANDARD, dvarp);

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
                                     nullptr, /* name */
                                     0, /* No error estimate */
                                     0, /* No error estimate */
                                     xnorm0,
                                     &eqp_loc,
                                     dvarp,
                                     dvarp,
                                     bc_coeffs_yp,
                                     i_mass_flux, b_mass_flux,
                                     i_mass_flux, b_mass_flux,
                                     i_mass_flux, b_mass_flux,
                                     nullptr, nullptr,
                                     nullptr,
                                     icvflb,
                                     nullptr,
                                     rovsdp,
                                     smbdp,
                                     dvarp, dpvar,
                                     nullptr, nullptr);

  /* Warning: no diffusion so no need of other diffusive
     boundary coefficient */

  /* Finalization and printing
     ------------------------- */

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t cdries = cs_turb_cdries;

  struct cs_double_n<2> rd_2r;
  struct cs_reduce_min_max_nr<1> reducer_2r;

  ctx.parallel_for_reduce
    (n_cells, rd_2r, reducer_2r,
     [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_double_n<2> &res) {

    /* Clipping: essential if you initialize with (u * /nu) of
       the previous time step */

    dvarp[c_id] = cs::max(dvarp[c_id], xusnmn);
    dvarp[c_id] = cs::min(dvarp[c_id], xusnmx);

    yplus[c_id] = dvarp[c_id] * w_dist[c_id];
    res.r[0] = yplus[c_id];
    res.r[1] = yplus[c_id];

    /* Van Driest amortization */
    visct[c_id] *= cs_math_pow2(1.0 - exp(-yplus[c_id] / cdries));
    if (visvdr[c_id] > -900.0)
      visct[c_id] = visvdr[c_id];
  });

  ctx.wait();
  cs_real_t dismin = rd_2r.r[0];
  cs_real_t dismax = rd_2r.r[1];
  cs_parall_min(1, CS_REAL_TYPE, &dismin);
  cs_parall_max(1, CS_REAL_TYPE, &dismax);

  if (eqp_yp->verbosity >= 1) {

    cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       " ** DIMENSIONLESS WALL DISTANCE\n"
       "    ---------------------------\n\n"
       "  Min distance+ = %14.5e, Max distance+ = %14.5e.\n"),
     dismin, dismax);
  }

  /* Free memory */
  CS_FREE_HD(dvarp);
  CS_FREE_HD(smbdp);
  CS_FREE_HD(dpvar);
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);
  CS_FREE_HD(rovsdp);
  CS_FREE_HD(viscap);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes distance to wall by a brute force geometric approach
 *        (serial only)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_distance_geometric(void)
{
  const cs_mesh_t  *mesh          = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const int *bc_type = cs_glob_bc_type;

  const cs_lnum_t n_cells       = mesh->n_cells;
  const cs_lnum_t n_b_faces     = mesh->n_b_faces;
  const cs_real_3_t *b_face_cog = fvq->b_face_cog;
  const cs_real_3_t *cell_cen   = fvq->cell_cen;

  // Usually one would not use MPI here but just in case...
  if (mesh->halo != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error function cannot be used in parallel"
                " or with periodic mesh"));

  cs_field_t *f_w_dist = cs_field_by_name("wall_distance");
  cs_real_t *wall_dist = f_w_dist->val;

  /* Deprecated model to compute wall distance */

  /* One must be careful in parallel or on periodic domains
     (one wall can be closer when crossing one edge...)*/

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    wall_dist[c_id] = cs_math_big_r * cs_math_big_r;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    if (bc_type[f_id] == CS_SMOOTHWALL || bc_type[f_id] == CS_ROUGHWALL) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t xdis = cs_math_3_square_distance(b_face_cog[f_id],
                                                   cell_cen[c_id]);
        if (wall_dist[c_id] > xdis)
          wall_dist[c_id] = xdis;

      }
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    wall_dist[c_id] = sqrt(wall_dist[c_id]);
  }

  /* Compute bounds and print infoc--------------------------------------------*/

  cs_real_t dismax = -cs_math_big_r;
  cs_real_t dismin =  cs_math_big_r;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    dismin = cs::min(wall_dist[c_id], dismin);
    dismax = cs::max(wall_dist[c_id], dismax);
  }

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       " ** WALL DISTANCE (brute force algorithm) \n"
       "    ------------- \n"
       " \n"
       " Min distance = %14.5f, Max distance = %14.5f \n"),
     dismin, dismax);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide read/write access to cs_glob_wall_distance
 *
 * \return pointer to global wall distance structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_distance_options_t *
cs_get_glob_wall_distance_options(void)
{
  return &_wall_distance_options;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
