/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <math.h>
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-turbomachinery.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function).
 *
 * This example is a part of the \subpage turbomachinery example.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * User function to find the closest cell and the associated rank of a point
 * in a given rotation zone.
 *
 * parameters:
 *   domain  <-- pointer to a cs_domain_t structure
 *   r       <-- id of the rotation zone
 *   coords  <-- point coordinates
 *   node    --> cell id
 *   rank    --> rank id
 *----------------------------------------------------------------------------*/

static void
_findpt_r(cs_domain_t          *domain,
          const cs_rotation_t  *r,
          const cs_real_3_t     coords,
          cs_lnum_t            *node,
          cs_lnum_t            *rank)
{
  cs_real_t d[3];

  cs_lnum_t n_cells = domain->mesh-> n_cells;
  cs_real_3_t *cell_cen = (cs_real_3_t *)domain->mesh_quantities->cell_cen;

  *node = (int)(n_cells + 1)/2 - 1;

  for (int i = 0; i < 3; i++)
    d[i] = coords[i] - cell_cen[*node][i];
  cs_real_t dis2mn = cs_math_3_square_norm(d);

  /*! [extra_tbm_get_rotor] */

  const cs_lnum_t *rotor_num = NULL;
  rotor_num = cs_turbomachinery_get_cell_rotor_num();

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_rotation_t *rot = cs_glob_rotation + rotor_num[cell_id];
    /* ... */

  /*! [extra_tbm_get_rotor] */

    if (r == rot) {
      for (int i = 0; i < 3; i++)
        d[i] = coords[i] - cell_cen[*node][i];
      cs_real_t dis2 = cs_math_3_square_norm(d);
      if (dis2 < dis2mn) {
        *node = cell_id;
        dis2mn = dis2;
      }
    }
  }

  if (cs_glob_rank_id >= 0)
    cs_parall_min_id_rank_r(node, rank, dis2mn);
  else
    *rank = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Example of extra operations for turbomachinery studies.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{

  /* Mesh-related variables */

  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)domain->mesh_quantities->cell_cen;

  /* 0. Initialization
     ================= */

  /* Local constants defintion */

  const cs_real_t _PI = acos(-1.);
  const cs_real_t _G = 9.81;

  /* Flow properties */

  double ro0 = cs_glob_fluid_properties->ro0;
  cs_real_t flowrate = 0.238;

  /* Assert there is a rotation zone */

  if (cs_glob_rotation == NULL)  return;

  /*! [extra_tbm_get_rotor_info] */

  /* Access to a specific rotation zone (the first one) */

  cs_lnum_t rotor_id = 1;
  const cs_rotation_t *ref_rot = cs_glob_rotation + rotor_id;

  /* Rotation zone parameters */

  double omega = ref_rot->omega;         /* angular velocity [rad/s] */
  cs_real_3_t axis = {ref_rot->axis[0],  /* rotation axis (normalized vector) */
                      ref_rot->axis[1], ref_rot->axis[2]};

  /*! [extra_tbm_get_rotor_info] */

  /* Example 1: calculation of the machinery characteristics
     ======================================================= */

  /*! [extra_tbm_post_util] */

  /* Postprocessing of the couple: axial moment resulting from the stresses */
  /* on the rotor blades */

  cs_lnum_t n_elts;
  cs_lnum_t *elt_list;

  BFT_MALLOC(elt_list, n_b_faces, cs_lnum_t);
  cs_selector_get_b_face_list("rotor_blades", &n_elts, elt_list);

  cs_real_t c = cs_post_moment_of_force(n_elts, elt_list, axis);

  BFT_FREE(elt_list);

  /* Postprocessing of the head: total pressure increase through the machinery */

  cs_real_t turbomachinery_head
    = cs_post_turbomachinery_head
        ("inlet",                          /* selection criteria at suction */
         CS_MESH_LOCATION_BOUNDARY_FACES,  /* associated mesh location */
         "z > 2.725 and z < 2.775",        /* selection criteria at discharge */
         CS_MESH_LOCATION_CELLS);          /* associated mesh location */

  /*! [extra_tbm_post_util] */

  cs_real_t manometric_head = turbomachinery_head/(ro0*_G);

  cs_real_t power = c*omega;

  cs_real_t efficiency = turbomachinery_head*flowrate/power;

  /* Print in the log */

  bft_printf("Turbomachinery characteristics:\n\n"
             "  %17s%17s%17s%17s\n"
             "  %17.9e%17.9e%17.9e%17.9e\n",
             "Flowrate [m3/s]","Head [m]", "Power [W]", "Efficiency [1]",
             flowrate, manometric_head, power, efficiency);

  /* Example 2: extraction of a velocity profile in cylindrical corordinates
     ======================================================================= */

  if (domain->time_step->nt_cur == domain->time_step->nt_max){

    cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

    cs_field_t *f_mom_wr = NULL, *f_mom_wt = NULL;

    /* In case of transient rotor/stator computation, the moments for
       the velocity in cylindrical coordinates are assumed to be defined
       in the dedicated user file (cs_user_parameters.c, see the user example
       for time moments definition) */

    cs_turbomachinery_model_t tbm_model = cs_turbomachinery_get_model();

    if (tbm_model == CS_TURBOMACHINERY_TRANSIENT) {
      f_mom_wr = cs_time_moment_get_field(1);
      f_mom_wt = cs_time_moment_get_field(2);
    }

    FILE *f1 = NULL, *f2 = NULL;

    /* Only process of rank 0 (parallel) or -1 (scalar) writes to this file. */

    if (cs_glob_rank_id <= 0) {

      f1 = fopen("Wr_mean.dat","w");
      fprintf(f1, "# %17s%17s\n", "angle", "Wr");

      f2 = fopen("Wt_mean.dat","w");
      fprintf(f2, "# %17s%17s\n", "angle", "Wt");
    }

    cs_real_t radius = 0.214;
    cs_real_t axicoo = 2.e-2;

    cs_lnum_t npoint = 360;
    cs_lnum_t cell_id1 = -999;
    cs_lnum_t rank_id1 = -999;

    for (cs_lnum_t point_id = 0; point_id < npoint; point_id++) {

      cs_real_t xth0 = -_PI/npoint;
      cs_real_3_t xyz = {radius*cos(xth0 - (double)point_id/npoint*2.*_PI),
                         radius*sin(xth0 - (double)point_id/npoint*2.*_PI),
                         axicoo};

      cs_lnum_t cell_id, rank_id;

      /* Find the closest cell in this rotor */
      _findpt_r(domain, ref_rot, xyz, &cell_id, &rank_id);

      if ((cell_id != cell_id1) || (rank_id != rank_id1)) {
        cell_id1 = cell_id;
        rank_id1 = rank_id;

        cs_real_t xr, xtheta, xvr, xvt;

        /* Set temporary variables for the process containing
         * the point and then send it to other processes. */
        if (cs_glob_rank_id == rank_id) {

          /* Radius (first component of the x vector in cylindrical coords) */
          cs_real_t xc[3];
          cs_rotation_cyl_v(ref_rot, cell_cen[cell_id], cell_cen[cell_id], xc);
          xr = xc[0];

          /* Angle in [0, 2pi] */
          double xthet1, xthet2;
          xthet1 = acos((cell_cen[cell_id][0] - ref_rot->invariant[0])/xr);
          xthet2 = asin((cell_cen[cell_id][1] - ref_rot->invariant[1])/xr);
          if (xthet2 > 0)
            xtheta = xthet1;
          else if (xthet1 < _PI/2.)
            xtheta = 2.*_PI + xthet2;
          else
            xtheta = _PI - xthet2;

          if (tbm_model == CS_TURBOMACHINERY_FROZEN) {
            cs_real_3_t xvc;
            cs_rotation_cyl_v(ref_rot, cell_cen[cell_id], vel[cell_id], xvc);
            xvr = xvc[0];
            /* Relative tangential velocity */
            xvt = xvc[1] - omega*xr;
          }
          else {
            xvr = f_mom_wr->val[cell_id];
            xvt = f_mom_wt->val[cell_id];
          }

        } else {
          xtheta = 0.;
          xvr = 0.;
          xvt = 0.;
        }

        /* Broadcast to other ranks in parallel */
        cs_parall_bcast(rank_id, 1, CS_DOUBLE, &xtheta);
        cs_parall_bcast(rank_id, 1, CS_DOUBLE, &xvr);
        cs_parall_bcast(rank_id, 1, CS_DOUBLE, &xvt);

        if (cs_glob_rank_id <= 0) {
          fprintf(f1,"  %17.9e%17.9e\n", xtheta, xvr);
          fprintf(f2,"  %17.9e%17.9e\n", xtheta, xvt);
        }
      }

    } /* end of loop on point_id */
    if (cs_glob_rank_id <= 0) {
      fclose(f1);
      fclose(f2);
    }

  } /* end of if statement on current time step */

}

END_C_DECLS
