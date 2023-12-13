/*============================================================================
 * Hybrid Temporal LES turbulence model.
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

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_array.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_equation.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gradient.h"
#include "cs_lagr.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_mass_source_terms.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rotation.h"
#include "cs_velocity_pressure.h"

#include "bft_printf.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_htles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_htles.c

  Solve the HTLES method for incompressible flows
  for one time step.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the HTLES method.
 *
 * Solve the HTLES for incompressible flows
 * for one time step.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_htles(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_real_t *dt = CS_F_(dt)->val;

  const cs_turb_model_t *turb_model = cs_get_glob_turb_model();

  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_real_t *crom       = (cs_real_t *)CS_F_(rho)->val;
  cs_real_t *cpro_pcvlo = (cs_real_t *)CS_F_(mu)->val;

  cs_real_3_t *vel    = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *cvar_k   = (cs_real_t *)CS_F_(k)->val;
  cs_real_t *cvar_omg = NULL;
  cs_real_t *cvar_eps = NULL;
  if (turb_model->iturb == CS_TURB_K_OMEGA) {
    cvar_omg = (cs_real_t *)CS_F_(omg)->val;
  }
  else {
    cvar_eps = (cs_real_t *)CS_F_(eps)->val;
  }

  const cs_real_t *w_dist =  cs_field_by_name("wall_distance")->val;

  const double d2s3 = 2./3.;

  cs_real_t *mean_omg = NULL;
  cs_real_t *kwsst_f1 = NULL;

  // TODO use standard time moments...
  cs_real_t *mean_u   = cs_field_by_name("vel_mag_mean")->val;
  cs_real_t *mean_k   = cs_field_by_name("k_tot")->val;
  cs_real_t *mean_km  = cs_field_by_name("k_mod")->val;
  cs_real_t *mean_kr  = cs_field_by_name("k_res")->val;
  cs_real_t *mean_eps = cs_field_by_name("eps_mod")->val;
  if (turb_model->iturb == CS_TURB_K_OMEGA) {
    mean_omg = cs_field_by_name("omg_mod")->val;
    kwsst_f1 = cs_field_by_name("f1_kwsst")->val;
  }

  cs_real_t *hyb_psi = cs_field_by_name("htles_psi")->val;
  cs_real_t *hyb_r   = cs_field_by_name("htles_r")->val;
  cs_real_t *hyb_t   = cs_field_by_name("htles_t")->val;
  cs_real_t *hyb_icc = cs_field_by_name("htles_icc")->val;
  cs_real_t *hyb_fs  = cs_field_by_name("htles_fs")->val;
  cs_real_t *dlt_max = cs_field_by_name("Delta_max")->val;
  cs_real_t *hybrid_fd_coeff = cs_field_by_name("hybrid_blend")->val;

  /* TEMP - TIME AVERAGED */
  cs_real_t time_mean = CS_MAX(cs_glob_turb_hybrid_model->n_iter_mean
                               * cs_glob_time_step->dt_ref,
                               cs_glob_turb_hybrid_model->time_mean);

  cs_real_3_t *mean_vel = (cs_real_3_t *)cs_field_by_name("velocity_mean")->val;
  cs_real_3_t *mean_ui2 = (cs_real_3_t *)cs_field_by_name("ui2_mean")->val;

# pragma omp for nowait
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Variables */
    cs_real_t xmu  = cpro_pcvlo[c_id];
    cs_real_t xro  = crom[c_id];
    cs_real_t xnu  = xmu/xro;
    cs_real_t xcmu  = cs_turb_cmu;
    cs_real_t xdmax = dlt_max[c_id];

    /* Dt divided by the exponential time filter width */
    cs_real_t factor = dt[c_id] / CS_MIN(time_mean, cs_glob_time_step->t_cur);

    /* Time averaged velocity magnitude */
    mean_u[c_id] += factor * (cs_math_3_norm(vel[c_id]) - mean_u[c_id]);

    mean_vel[c_id][0] += factor * (vel[c_id][0] - mean_vel[c_id][0]);
    mean_vel[c_id][1] += factor * (vel[c_id][1] - mean_vel[c_id][1]);
    mean_vel[c_id][2] += factor * (vel[c_id][2] - mean_vel[c_id][2]);
    /* Time averaged second order moments */
    mean_ui2[c_id][0] += factor * (vel[c_id][0]*vel[c_id][0] - mean_ui2[c_id][0]);
    mean_ui2[c_id][1] += factor * (vel[c_id][1]*vel[c_id][1] - mean_ui2[c_id][1]);
    mean_ui2[c_id][2] += factor * (vel[c_id][2]*vel[c_id][2] - mean_ui2[c_id][2]);

    mean_km[c_id] += factor * (cvar_k[c_id] - mean_km[c_id]);
    mean_kr[c_id] = 0.5*( mean_ui2[c_id][0] - cs_math_pow2(mean_vel[c_id][0])
                        + mean_ui2[c_id][1] - cs_math_pow2(mean_vel[c_id][1])
                        + mean_ui2[c_id][2] - cs_math_pow2(mean_vel[c_id][2])
                        );
    mean_k[c_id] = mean_km[c_id] + mean_kr[c_id];

    /* Time averaged turbulent dissipation */
    if (turb_model->iturb == CS_TURB_K_OMEGA) {
      mean_omg[c_id] += factor * (cvar_omg[c_id] - mean_omg[c_id]);
      mean_eps[c_id] = xcmu*mean_km[c_id]*mean_omg[c_id];
    }
    else {
      mean_eps[c_id] += factor * (cvar_eps[c_id] - mean_eps[c_id]);
    }

    /* Temporally averaged variables */
    cs_real_t xpsi0 = hyb_psi[c_id];
    cs_real_t xum   = mean_u[c_id];
    cs_real_t xkt   = mean_k[c_id];
    cs_real_t xkm   = mean_km[c_id];
    cs_real_t xkr   = mean_kr[c_id];
    cs_real_t xepsm = mean_eps[c_id];

    /* Shielding function */
    cs_real_t xfs = 1.0;
    if (cs_glob_turb_hybrid_model->ishield == 1) {
      cs_real_t xdist  = CS_MAX(w_dist[c_id], cs_math_epzero);
      cs_real_t xsik   = 45.0 * pow(xnu, 0.75)/(pow(xpsi0*xepsm, 0.25)*xdist);
      cs_real_t xsid   = pow(3.0, 1.0/6.0)*xdmax/xdist;
      xfs = 1.0 - tanh(CS_MAX(pow(xsik, 8.0),pow(xsid, 6.0)));
    }

    /* Analytic energy ratio r */
    cs_real_t xbt0   = cs_turb_chtles_bt0;
    cs_real_t xdelta = pow(cell_f_vol[c_id], 1./3.);
    cs_real_t xus = xum + sqrt(d2s3)*sqrt(xkt);
    cs_real_t xwc = CS_MIN(cs_math_pi/dt[c_id], xus*cs_math_pi/xdelta);
    cs_real_t xrk = 1.0/xbt0 * pow(xus/sqrt(xkt), d2s3)
      *pow(xwc*xkt/(xpsi0*xepsm), -d2s3);

    /* Energy ratio */
    cs_real_t xr = 1.0;
    if (cs_glob_turb_hybrid_model->ishield == 1) {
      xr  = (1.0 - xfs) + xfs*CS_MIN(1.0, xrk);
    }
    else if (cs_glob_turb_hybrid_model->ishield == 0) {
      xr  = CS_MIN(1.0, xrk);
    }

    /* Hybrid scheme (same for ICC) */
    cs_real_t xrc = 1.0;
    if (xr > 0.9999) {
      xrc = 0.0;
    }
    else {
      if (cs_glob_turb_hybrid_model->ishield == 1) {
        xrc = xfs;
      }
      else {
        xrc = 1.0;
      }
    }

    /* Internal Consistency Constraint (ICC) */
    cs_real_t xicc = 1.0;
    if (cs_glob_turb_hybrid_model->iicc == 1)
      xicc = xrc;

    /* Hybridation function psi */
    cs_real_t xpsi = 0;
    if (turb_model->iturb == CS_TURB_K_OMEGA) {
      cs_real_t xxf1   = kwsst_f1[c_id];
      cs_real_t xgamma = xxf1 * cs_turb_ckwgm1 + (1. - xxf1)*cs_turb_ckwgm2;
      cs_real_t xbeta  = xxf1 * cs_turb_ckwbt1 + (1. - xxf1)*cs_turb_ckwbt2;
      xpsi   = xbeta/(xcmu*xgamma + xr*(xbeta - xcmu*xgamma));
    }
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      cs_real_t xce1 = cs_turb_cpale1;
      cs_real_t xce2 = cs_turb_cpale2;
      xpsi   = xce2/(xce1 + xr*(xce2 - xce1));
    }

    /* Modeled time scale T */
    cs_real_t xt = 0;
    if (turb_model->iturb == CS_TURB_K_OMEGA)
      xt = xr/xpsi*(xkm+xicc*xkr)/(xcmu * mean_omg[c_id] * xkm);
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K)
      xt = xr/xpsi*(xkm+xicc*xkr)/(xepsm);

    /* Save fields */
    hybrid_fd_coeff[c_id] = xrc;
    hyb_psi[c_id] = xpsi;
    hyb_r[c_id]   = xr;
    hyb_t[c_id]   = xt;
    hyb_icc[c_id] = xicc;
    hyb_fs[c_id]  = xfs;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization of the fields for HTLES and Delta_max calculation
 */
/*----------------------------------------------------------------------------*/

void
cs_htles_initialization(void) {

  const cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_turb_model_t *turb_model = cs_get_glob_turb_model();

  cs_real_3_t *vel    = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *cvar_k   = (cs_real_t *)CS_F_(k)->val;
  cs_real_t *cvar_omg = NULL;
  cs_real_t *cvar_eps = NULL;
  if (turb_model->iturb == CS_TURB_K_OMEGA) {
    cvar_omg = (cs_real_t *)CS_F_(omg)->val;
  }
  else {
    cvar_eps = (cs_real_t *)CS_F_(eps)->val;
  }

  cs_real_t xcmu  = cs_turb_cmu;

  cs_real_t *mean_u   = cs_field_by_name("vel_mag_mean")->val;
  cs_real_t *mean_k   = cs_field_by_name("k_tot")->val;
  cs_real_t *mean_km  = cs_field_by_name("k_mod")->val;
  cs_real_t *mean_kr  = cs_field_by_name("k_res")->val;
  cs_real_t *mean_eps = cs_field_by_name("eps_mod")->val;

  cs_real_t *mean_omg = NULL;
  cs_real_t *kwsst_f1 = NULL;

  if (turb_model->iturb == CS_TURB_K_OMEGA) {
    mean_omg = cs_field_by_name("omg_mod")->val;
    kwsst_f1 = cs_field_by_name("f1_kwsst")->val;
  }

  cs_real_t *hyb_psi = cs_field_by_name("htles_psi")->val;
  cs_real_t *hyb_r   = cs_field_by_name("htles_r")->val;
  cs_real_t *hyb_t   = cs_field_by_name("htles_t")->val;
  cs_real_t *hyb_icc = cs_field_by_name("htles_icc")->val;
  cs_real_t *hyb_fs  = cs_field_by_name("htles_fs")->val;
  cs_real_t *dlt_max = cs_field_by_name("Delta_max")->val;

  /* Time averaged */
  cs_real_3_t *mean_vel = (cs_real_3_t *)cs_field_by_name("velocity_mean")->val;
  cs_real_3_t *mean_ui2 = (cs_real_3_t *)cs_field_by_name("ui2_mean")->val;

  /**************************************/
  /* Initialization of the HTLES fields */
  /**************************************/

# pragma omp for nowait
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Initialization of HTLES fields */
    hyb_psi[c_id] = 1.0;
    hyb_r[c_id]   = 1.0;
    hyb_icc[c_id] = 1.0;
    hyb_fs[c_id]  = 1.0;
    if (turb_model->iturb == CS_TURB_K_OMEGA) {
      hyb_t[c_id] = 1.0/(xcmu*cvar_omg[c_id]);
    }
    else {
      hyb_t[c_id] = cvar_k[c_id]/cvar_eps[c_id];
    }

    /* Time averaged velocity */
    mean_u[c_id] = cs_math_3_norm(vel[c_id]);
    mean_vel[c_id][0] = vel[c_id][0];
    mean_vel[c_id][1] = vel[c_id][1];
    mean_vel[c_id][2] = vel[c_id][2];
    mean_ui2[c_id][0] = vel[c_id][0]*vel[c_id][0];
    mean_ui2[c_id][1] = vel[c_id][1]*vel[c_id][1];
    mean_ui2[c_id][2] = vel[c_id][2]*vel[c_id][2];

    /* Time averaged turbulent energy */
    mean_km[c_id] = cvar_k[c_id];
    mean_kr[c_id] = 0.0;
    mean_k[c_id]  = mean_km[c_id] + mean_kr[c_id];

    /* Time averaged turbulent dissipation */
    if (turb_model->iturb == CS_TURB_K_OMEGA) {
      mean_omg[c_id] = cvar_omg[c_id];
      mean_eps[c_id] = xcmu*mean_km[c_id]*mean_omg[c_id];
      kwsst_f1[c_id] = 1.0;
    }
    else {
      mean_eps[c_id] = cvar_eps[c_id];
    }

  }

  /****************************/
  /* Calculation of Delta max */
  /****************************/

  /* Calculate the maximum number of edges per cell (with the redundancy)
     as a function of the maximum number of faces per cell
     and of the maximum number of edges per face (with diagonals):
     For hexaedral cells, it is 6 and 6, respectively,
     For tetrahedral cells, it is 4 and 3, respectively.  */
  cs_lnum_t max_faces_per_cells = 6;
  cs_lnum_t max_edg_per_faces = 6;
  cs_lnum_t max_edg_per_cells = max_faces_per_cells*max_edg_per_faces;

  cs_lnum_t c_id_1, c_id_2;
  cs_lnum_t vtx1_id, vtx2_id;
  cs_lnum_t vtx1_edg1_id, vtx2_edg1_id;
  cs_lnum_t vtx1_edg2_id, vtx2_edg2_id;
  cs_lnum_t vtx_start, vtx_end;
  cs_lnum_t cpt_edg_cell;

  cs_real_t cnx, cny, cnz;

  cs_lnum_t *vtx1_edg_per_cells = NULL;
  cs_lnum_t *vtx2_edg_per_cells = NULL;
  cs_lnum_t *cpt_edg_per_cells = NULL;

  BFT_MALLOC(vtx1_edg_per_cells, n_cells*max_edg_per_cells, cs_lnum_t);
  BFT_MALLOC(vtx2_edg_per_cells, n_cells*max_edg_per_cells, cs_lnum_t);
  BFT_MALLOC(cpt_edg_per_cells, n_cells, cs_lnum_t);

  /* The counter is set to zero */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cpt_edg_per_cells[c_id] = 0;
  }

  /* For all interior faces
     save the vertices of each edge */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    c_id_1 = m->i_face_cells[f_id][0];
    c_id_2 = m->i_face_cells[f_id][1];

    vtx_start = m->i_face_vtx_idx[f_id];
    vtx_end   = m->i_face_vtx_idx[f_id + 1];
    for (cs_lnum_t vtx1 = vtx_start; vtx1 < vtx_end-1; vtx1++) {
      for (cs_lnum_t vtx2 = vtx1+1; vtx2 < vtx_end; vtx2++) {

        vtx1_id = m->i_face_vtx_lst[vtx1];
        vtx2_id = m->i_face_vtx_lst[vtx2];

        /* Not take into account the ghost cells in case of periodic boundaries */
        if (c_id_1 < n_cells) {
          vtx1_edg_per_cells[c_id_1*max_edg_per_cells + cpt_edg_per_cells[c_id_1]] = vtx1_id;
          vtx2_edg_per_cells[c_id_1*max_edg_per_cells + cpt_edg_per_cells[c_id_1]] = vtx2_id;
          cpt_edg_per_cells[c_id_1] += 1;
        }
        if (c_id_2 < n_cells) {
          vtx1_edg_per_cells[c_id_2*max_edg_per_cells + cpt_edg_per_cells[c_id_2]] = vtx1_id;
          vtx2_edg_per_cells[c_id_2*max_edg_per_cells + cpt_edg_per_cells[c_id_2]] = vtx2_id;
          cpt_edg_per_cells[c_id_2] += 1;
        }
      }
    }
  }

  /* For all boundary faces
     save the vertices of each edge */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    c_id_1 = m->b_face_cells[f_id];

    vtx_start = m->b_face_vtx_idx[f_id];
    vtx_end   = m->b_face_vtx_idx[f_id + 1];
    for (cs_lnum_t vtx1 = vtx_start; vtx1 < vtx_end-1; vtx1++) {
      for (cs_lnum_t vtx2 = vtx1+1; vtx2 < vtx_end; vtx2++) {

        vtx1_id = m->b_face_vtx_lst[vtx1];
        vtx2_id = m->b_face_vtx_lst[vtx2];

        vtx1_edg_per_cells[c_id_1*max_edg_per_cells + cpt_edg_per_cells[c_id_1]] = vtx1_id;
        vtx2_edg_per_cells[c_id_1*max_edg_per_cells + cpt_edg_per_cells[c_id_1]] = vtx2_id;
        cpt_edg_per_cells[c_id_1] += 1;
      }
    }
  }

  /* Verification
   * TODO correct it in parallel
   * */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cpt_edg_per_cells[c_id] != 36) {
      bft_printf("HTLES: Error counter: %d \n", cpt_edg_per_cells[c_id]);
    }
  }

  /* Calculation of Delta max :
     Loop over all the cells, and all the edges, cutting the diagonal */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t dmax = 0.0;

    cpt_edg_cell = cpt_edg_per_cells[c_id];
    for (cs_lnum_t edg1_id = 0; edg1_id < cpt_edg_cell-1; edg1_id++) {
      for (cs_lnum_t edg2_id = edg1_id+1; edg2_id < cpt_edg_cell; edg2_id++) {

        /* if the edge appears twice, it is not a diagonal
           (i.e. the pair of vertices (vtx1, vtx2) is the same for edg1 and edg2)
           so the edge is taken into account to determine Delta_max */
        vtx1_edg1_id = vtx1_edg_per_cells[c_id*max_edg_per_cells + edg1_id];
        vtx2_edg1_id = vtx2_edg_per_cells[c_id*max_edg_per_cells + edg1_id];
        vtx1_edg2_id = vtx1_edg_per_cells[c_id*max_edg_per_cells + edg2_id];
        vtx2_edg2_id = vtx2_edg_per_cells[c_id*max_edg_per_cells + edg2_id];
        if (((vtx1_edg1_id == vtx1_edg2_id) && (vtx2_edg1_id == vtx2_edg2_id))
            || ((vtx1_edg1_id == vtx2_edg2_id) && (vtx2_edg1_id == vtx1_edg2_id))) {
          vtx1_id = vtx1_edg1_id;
          vtx2_id = vtx2_edg1_id;

          /* Coordinates and length of the edge between vtx1 and vtx2 */
          cnx = m->vtx_coord[vtx1_id*3 + 0] - m->vtx_coord[vtx2_id*3 + 0];
          cny = m->vtx_coord[vtx1_id*3 + 1] - m->vtx_coord[vtx2_id*3 + 1];
          cnz = m->vtx_coord[vtx1_id*3 + 2] - m->vtx_coord[vtx2_id*3 + 2];
          cs_real_t dxyz = sqrt(pow(cnx,2) + pow(cny, 2) + pow(cnz, 2));

          /* Delta_max of the cell */
          dmax = CS_MAX(dmax, dxyz);
        }
      }
    }
    dlt_max[c_id] = dmax;

  }

  BFT_FREE(vtx1_edg_per_cells);
  BFT_FREE(vtx2_edg_per_cells);
  BFT_FREE(cpt_edg_per_cells);

}

END_C_DECLS
