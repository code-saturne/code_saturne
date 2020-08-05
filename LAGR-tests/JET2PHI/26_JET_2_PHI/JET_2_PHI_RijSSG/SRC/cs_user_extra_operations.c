/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user subroutine)
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

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

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_parall.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"
#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user subroutine)
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Local user function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Sort an array based on coordinates
 *----------------------------------------------------------------------------*/

static void
_sortc2(cs_real_t  coo[],
        cs_real_t  val[],
        int        n)
{
  int ns = 0;
  for (int ii = 0; ii < n - 1; ii++) {
    cs_real_t c_min = 1e+20;
    int kk = 0;
    for (int jj = ns; jj < n; jj++) {
      if (c_min > coo[jj]) {
        c_min = coo[jj];
        kk     = jj;
      }
    }
    cs_real_t xx = val[kk];
    coo[kk] = coo[ns];
    val[kk] = val[ns];
    coo[ns] = c_min;
    val[ns] = xx;
    ns++;
  }
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  cs_lnum_t  n_cells     = cs_glob_mesh->n_cells;
  cs_lnum_t  n_b_faces   = cs_glob_mesh->n_b_faces;
  cs_lnum_t  n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)cs_glob_mesh->i_face_cells;
  const cs_real_t *restrict f_weight = cs_glob_mesh_quantities->weight;

  /* Initialization */
  /* -------------- */

  if (cs_glob_time_step->nt_cur < cs_glob_time_step->nt_max)
    return;

  cs_real_t xlist[100][8];

  const cs_real_3_t *i_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->i_face_cog;

  /* Allocate a temporary array for cells or interior/boundary faces selection */

  cs_lnum_t  *elts = NULL;
  BFT_MALLOC(elts, CS_MAX(CS_MAX(n_cells, n_i_faces), n_b_faces), cs_lnum_t);

  const cs_real_t zz[] = {0.005, 0.025, 0.050, 0.075,
                          0.100, 0.150, 0.200, 0.250};

  const char *name[] = {"xl01", "xl05", "xl10", "xl15",
                        "xl20", "xl30", "xl40", "xl50"};

  const  cs_lagr_stat_type_t s_types[]
    = {cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY),
       cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY),
       CS_LAGR_STAT_VOLUME_FRACTION,
       CS_LAGR_STAT_CUMULATIVE_WEIGHT};
  const  cs_lagr_stat_moment_t m_types[]
    = {CS_LAGR_MOMENT_MEAN,
       CS_LAGR_MOMENT_VARIANCE,
       CS_LAGR_MOMENT_MEAN,
       CS_LAGR_MOMENT_MEAN};

  for (int plane_id = 0; plane_id < 8; plane_id++) {

    cs_field_t *f_p[4];
    for (int s_id = 0; s_id < 3; s_id++)
      f_p[s_id] = cs_lagr_stat_get_moment(s_types[s_id],
                                          CS_LAGR_STAT_GROUP_PARTICLE,
                                          m_types[s_id],
                                          0, -1);
    f_p[3] = cs_lagr_stat_get_stat_weight(0);

    int id = 0;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      if (   i_face_cog[face_id][0] <  zz[plane_id] + 1.e-5
          && i_face_cog[face_id][0] >= zz[plane_id] - 1.e-5) {

        cs_lnum_t ic0 = i_face_cells[face_id][0];
        cs_lnum_t ic1 = i_face_cells[face_id][1];

        cs_real_t fw = f_weight[face_id];

        xlist[id][0] = i_face_cog[face_id][0];
        xlist[id][1] = i_face_cog[face_id][2] * (1.e3 / 5.);

        int col_id = 2;

        for (int s_id = 0; s_id < 4; s_id++) {

          const cs_field_t *f = f_p[s_id];

          if (f->dim == 3) {
            xlist[id][col_id] =   f->val[ic0*3 + 0] * fw
                                + f->val[ic1*3 + 0] * (1 -fw);
            col_id++;
            xlist[id][col_id] =   f->val[ic0*3 + 2] * fw
                                + f->val[ic1*3 + 2] * (1 -fw);
            col_id++;
          }
          else if (f->dim == 6) {
            xlist[id][col_id] =   f->val[ic0*6 + 0] * fw
                                + f->val[ic1*6 + 0] * (1 -fw);
            col_id++;
            xlist[id][col_id] =   f->val[ic0*6 + 2] * fw
                                + f->val[ic1*6 + 2] * (1 -fw);
            col_id++;
          }
          else if (f->dim == 1) {
            xlist[id][col_id] =   f->val[ic0] * fw
                                + f->val[ic1] * (1 -fw);
            col_id++;
          }

        }

        id++;

      }

    }

    cs_gnum_t n_g = id;
    cs_parall_counter(&n_g, 1);

    cs_real_t *gxlist;
    BFT_MALLOC(gxlist, n_g*8, cs_real_t);

    if (cs_glob_n_ranks > 1)
      cs_parall_allgather_r(id*8, n_g*8, (cs_real_t *)xlist, gxlist);
    else {
      for (int i = 0; i < id; i++) {
        for (int j = 0; j < 8; j++)
          gxlist[i*8 + j] = xlist[i][j];
      }
    }

    if (cs_glob_rank_id < 1) {

      cs_real_t *coo_aux, *val_aux;
      BFT_MALLOC(coo_aux, n_g, cs_real_t);
      BFT_MALLOC(val_aux, n_g*8, cs_real_t);

      for (int col_id = 0; col_id < 8; col_id ++) {

        for (int i = 0; i < (int)n_g; i++) {
          coo_aux[i] = gxlist[i*8 + 1];
          val_aux[n_g*col_id + i] = gxlist[i*8 + col_id];
        }
        _sortc2(coo_aux, val_aux+n_g*col_id, n_g);

      }

      FILE *f = fopen(name[plane_id], "w");

      for (int i = 0; i < (int)n_g; i++) {

        for (int col_id = 4; col_id < 6; col_id ++)
          val_aux[n_g*col_id + i] = sqrt(fmax(0., val_aux[n_g*col_id + i]));

        for (int col_id = 0; col_id < 8; col_id ++)
          fprintf(f, " %13.5e", val_aux[n_g*col_id + i]);
          fprintf(f, "\n");
      }

      fclose(f);

      BFT_FREE(val_aux);
      BFT_FREE(coo_aux);

    }

    BFT_FREE(gxlist);
  }

  BFT_FREE(elts);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
