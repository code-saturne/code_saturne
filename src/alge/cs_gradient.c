/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/* Basic per gradient compuation options and logging */
/*---------------------------------------------------*/

typedef struct _cs_gradient_info_t {

  char                *name;               /* System name */
  cs_gradient_type_t   type;               /* Gradient type */

  unsigned             n_calls;            /* Number of times system solved */

  cs_timer_counter_t   t_tot;              /* Total time used */

} cs_gradient_info_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int cs_glob_gradient_n_systems = 0;      /* Current number of systems */
static int cs_glob_gradient_n_max_systems = 0;  /* Max. number of sytems for
                                                   cs_glob_gradient_systems. */

/* System info array */
static cs_gradient_info_t **cs_glob_gradient_systems = NULL;

/* Short names for gradient computation types */

const char *cs_gradient_type_name[] = {N_("Iterative reconstruction"),
                                       N_("Least-squares (standard)"),
                                       N_("Least-squares (extended)"),
                                       N_("Least-squares (partially extended)"),
                                       N_("Least-squares then iterative")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to new gradient computation info structure.
 *
 * parameters:
 *   name --> system name
 *   type --> resolution method
 *
 * returns:
 *   pointer to newly created linear system info structure
 *----------------------------------------------------------------------------*/

static cs_gradient_info_t *
_gradient_info_create(const char          *name,
                      cs_gradient_type_t   type)
{
  cs_gradient_info_t *new_info = NULL;

  BFT_MALLOC(new_info, 1, cs_gradient_info_t);
  BFT_MALLOC(new_info->name, strlen(name) + 1, char);

  strcpy(new_info->name, name);
  new_info->type = type;

  new_info->n_calls = 0;

  CS_TIMER_COUNTER_INIT(new_info->t_tot);

  return new_info;
}

/*----------------------------------------------------------------------------
 * Destroy gradient computationw info structure.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure pointer
 *----------------------------------------------------------------------------*/

static void
_gradient_info_destroy(cs_gradient_info_t  **this_info)
{
  if (*this_info != NULL) {
    BFT_FREE((*this_info)->name);
    BFT_FREE(*this_info);
  }
}

/*----------------------------------------------------------------------------
 * Output information regarding gradient computation.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure
 *----------------------------------------------------------------------------*/

static void
_gradient_info_dump(cs_gradient_info_t *this_info)
{
  int n_calls = this_info->n_calls;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Summary of gradient computations pour \"%s\" (%s):\n\n"
                  "  Number of calls:     %12d\n"
                  "  Total elapsed time:  %12.3f\n"),
                this_info->name, cs_gradient_type_name[this_info->type],
                n_calls, this_info->t_tot.wall_nsec*1e-9);
}

/*----------------------------------------------------------------------------
 * Return pointer to gradient computation info.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name --> system name
 *   type --> resolution method
 *----------------------------------------------------------------------------*/

static cs_gradient_info_t *
_find_or_add_system(const char          *name,
                    cs_gradient_type_t   type)
{
  int ii, start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find system */

  start_id = 0;
  end_id = cs_glob_gradient_n_systems - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((cs_glob_gradient_systems[mid_id])->name, name);
    if (cmp_ret == 0)
      cmp_ret = (cs_glob_gradient_systems[mid_id])->type - type;
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If found, return */

  if (cmp_ret == 0)
    return cs_glob_gradient_systems[mid_id];

  /* Reallocate global array if necessary */

  if (cs_glob_gradient_n_systems >= cs_glob_gradient_n_max_systems) {

    if (cs_glob_gradient_n_max_systems == 0)
      cs_glob_gradient_n_max_systems = 10;
    else
      cs_glob_gradient_n_max_systems *= 2;
    BFT_REALLOC(cs_glob_gradient_systems,
                cs_glob_gradient_n_max_systems,
                cs_gradient_info_t*);

  }

  /* Insert in sorted list */

  for (ii = cs_glob_gradient_n_systems; ii > mid_id; ii--)
    cs_glob_gradient_systems[ii] = cs_glob_gradient_systems[ii - 1];

  cs_glob_gradient_systems[mid_id] = _gradient_info_create(name,
                                                           type);
  cs_glob_gradient_n_systems += 1;

  return cs_glob_gradient_systems[mid_id];
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a scalar if necessary. This function deals with
 * the standard or extended neighborhood.
 *
 * parameters:
 *   imrgra         <-- type of computation for the gradient
 *   imligp         <-- type of clipping for the computation of the gradient
 *   iwarnp         <-- output level
 *   idimtr         <-- 0 for scalars or without rotational periodicity,
 *                      1 or 2 for vectors or tensors in case of rotational
 *                      periodicity
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   var            <-- variable
 *   dpdx           <-> X component of the pressure gradient
 *   dpdy           <-> Y component of the pressure gradient
 *   dpdz           <-> Z component of the pressure gradient
 *----------------------------------------------------------------------------*/

static void
_scalar_gradient_clipping(const cs_int_t   *imrgra,
                          const cs_int_t   *imligp,
                          const cs_int_t   *iwarnp,
                          const cs_int_t   *idimtr,
                          const cs_real_t  *climgp,
                          cs_real_t         var[],
                          cs_real_t         dpdx[],
                          cs_real_t         dpdy[],
                          cs_real_t         dpdz[])
{
  cs_lnum_t  i, i1, i2, j, k;
  cs_real_t  dist[3];
  cs_real_t  dvar, dist1, dist2, dpdxf, dpdyf, dpdzf;
  cs_real_t  global_min_factor, global_max_factor, factor1, factor2;

  cs_gnum_t  n_clip = 0, n_g_clip =0;
  cs_real_t  min_factor = 1;
  cs_real_t  max_factor = 0;
  cs_real_t  *restrict buf = NULL, *restrict clip_factor = NULL;
  cs_real_t  *restrict denom = NULL, *restrict denum = NULL;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t  *cell_cen = cs_glob_mesh_quantities->cell_cen;

  const cs_halo_t *halo = mesh->halo;

  if (*imligp < 0)
    return;

  if (*imrgra == 2 || *imrgra ==  3)
    halo_type = CS_HALO_EXTENDED;

  /* Synchronize variable */

  if (halo != NULL) {

    cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, var);

    /* Exchange for the gradients. Not useful for working array */

    if (*imligp == 1) {

      if (*idimtr > 0){
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdx);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdy);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdz);
      }
      else {
        cs_halo_sync_var(halo, halo_type, dpdx);
        cs_halo_sync_var(halo, halo_type, dpdy);
        cs_halo_sync_var(halo, halo_type, dpdz);
        cs_halo_perio_sync_var_vect_ni(halo, halo_type,
                                       dpdx, dpdy, dpdz);
      }

    } /* End if imligp == 1 */

  } /* End if halo */

  /* Allocate and initialize working buffers */

  if (*imligp == 1)
    BFT_MALLOC(buf, 3*n_cells_wghosts, cs_real_t);
  else
    BFT_MALLOC(buf, 2*n_cells_wghosts, cs_real_t);

  denum = buf;
  denom = buf + n_cells_wghosts;

  if (*imligp == 1)
    clip_factor = buf + 2*n_cells_wghosts;

  for (i = 0; i < n_cells_wghosts; i++) {
    denum[i] = 0;
    denom[i] = 0;
  }

  /* First computations:
      denum holds the maximum variation of the gradient
      denom holds the maximum variation of the variable */

  if (*imligp == 0) {

    for (i = 0; i < n_i_faces; i++) {

      i1 = mesh->i_face_cells[2*i] - 1;
      i2 = mesh->i_face_cells[2*i + 1] - 1;

      for (j = 0; j < 3; j++)
        dist[j] = cell_cen[3*i1 + j] - cell_cen[3*i2 + j];

      dist1 = CS_ABS(dist[0]*dpdx[i1] + dist[1]*dpdy[i1] + dist[2]*dpdz[i1]);
      dist2 = CS_ABS(dist[0]*dpdx[i2] + dist[1]*dpdy[i2] + dist[2]*dpdz[i2]);

      dvar = CS_ABS(var[i1] - var[i2]);

      denum[i1] = CS_MAX(denum[i1], dist1);
      denum[i2] = CS_MAX(denum[i2], dist2);
      denom[i1] = CS_MAX(denom[i1], dvar);
      denom[i2] = CS_MAX(denom[i2], dvar);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (i1 = 0; i1 < n_cells; i1++) {
        for (j = cell_cells_idx[i1] - 1; j < cell_cells_idx[i1+1] - 1; j++) {

          i2 = cell_cells_lst[j] - 1;

          for (k = 0; k < 3; k++)
            dist[k] = cell_cen[3*i1 + k] - cell_cen[3*i2 + k];

          dist1 = CS_ABS(  dist[0]*dpdx[i1]
                         + dist[1]*dpdy[i1]
                         + dist[2]*dpdz[i1]);
          dvar = CS_ABS(var[i1] - var[i2]);

          denum[i1] = CS_MAX(denum[i1], dist1);
          denom[i1] = CS_MAX(denom[i1], dvar);

        }
      }

    } /* End for extended halo */

  }
  else if (*imligp == 1) {

    for (i = 0; i < n_i_faces; i++) {

      i1 = mesh->i_face_cells[2*i] - 1;
      i2 = mesh->i_face_cells[2*i + 1] - 1;

      for (j = 0; j < 3; j++)
        dist[j] = cell_cen[3*i1 + j] - cell_cen[3*i2 + j];

      dpdxf = 0.5 * (dpdx[i1] + dpdx[i2]);
      dpdyf = 0.5 * (dpdy[i1] + dpdy[i2]);
      dpdzf = 0.5 * (dpdz[i1] + dpdz[i2]);

      dist1 = CS_ABS(dist[0]*dpdxf + dist[1]*dpdyf + dist[2]*dpdzf);
      dvar = CS_ABS(var[i1] - var[i2]);

      denum[i1] = CS_MAX(denum[i1], dist1);
      denum[i2] = CS_MAX(denum[i2], dist1);
      denom[i1] = CS_MAX(denom[i1], dvar);
      denom[i2] = CS_MAX(denom[i2], dvar);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (i1 = 0; i1 < n_cells; i1++) {
        for (j = cell_cells_idx[i1] - 1; j < cell_cells_idx[i1+1] - 1; j++) {

          i2 = cell_cells_lst[j] - 1;

          for (k = 0; k < 3; k++)
            dist[k] = cell_cen[3*i1 + k] - cell_cen[3*i2 + k];

          dpdxf = 0.5 * (dpdx[i1] + dpdx[i2]);
          dpdyf = 0.5 * (dpdy[i1] + dpdy[i2]);
          dpdzf = 0.5 * (dpdz[i1] + dpdz[i2]);

          dist1 = CS_ABS(dist[0]*dpdxf + dist[1]*dpdyf + dist[2]*dpdzf);
          dvar = CS_ABS(var[i1] - var[i2]);

          denum[i1] = CS_MAX(denum[i1], dist1);
          denom[i1] = CS_MAX(denom[i1], dvar);

        }
      }

    } /* End for extended neighborhood */

  } /* End if *imligp == 1 */

  /* Clipping of the gradient if denum/denom > climgp */

  if (*imligp == 0) {

    for (i = 0; i < n_cells; i++) {

      if (denum[i] > *climgp * denom[i]) {

        factor1 = *climgp * denom[i]/denum[i];
        dpdx[i] *= factor1;
        dpdy[i] *= factor1;
        dpdz[i] *= factor1;

        min_factor = CS_MIN( factor1, min_factor);
        max_factor = CS_MAX( factor1, max_factor);
        n_clip++;

      } /* If clipping */

    } /* End of loop on cells */

  }
  else if (*imligp == 1) {

    for (i = 0; i < n_cells_wghosts; i++)
      clip_factor[i] = (cs_real_t)DBL_MAX;


    /* Synchronize variable */

    if (halo != NULL) {
      if (*idimtr > 0) {
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, denom);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, denum);
      }
      else {
        cs_halo_sync_var(halo, halo_type, denom);
        cs_halo_sync_var(halo, halo_type, denum);
      }
    }

    for (i = 0; i < n_i_faces; i++) {

      i1 = mesh->i_face_cells[2*i] - 1;
      i2 = mesh->i_face_cells[2*i + 1] - 1;

      factor1 = 1.0;
      if (denum[i1] > *climgp * denom[i1])
        factor1 = *climgp * denom[i1]/denum[i1];

      factor2 = 1.0;
      if (denum[i2] > *climgp * denom[i2])
        factor2 = *climgp * denom[i2]/denum[i2];

      min_factor = CS_MIN(factor1, factor2);

      clip_factor[i1] = CS_MIN( clip_factor[i1], min_factor);
      clip_factor[i2] = CS_MIN( clip_factor[i2], min_factor);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (i1 = 0; i1 < n_cells; i1++) {

        factor1 = 1.0;

        for (j = cell_cells_idx[i1] - 1; j < cell_cells_idx[i1+1] - 1; j++) {

          i2 = cell_cells_lst[j] - 1;
          factor2 = 1.0;

          if (denum[i2] > *climgp * denom[i2])
            factor2 = *climgp * denom[i2]/denum[i2];

          factor1 = CS_MIN(factor1, factor2);

        }

        clip_factor[i1] = CS_MIN(clip_factor[i1], factor1);

      } /* End of loop on cells */

    } /* End for extended neighborhood */

    for (i = 0; i < n_cells; i++) {

      dpdx[i] *= clip_factor[i];
      dpdy[i] *= clip_factor[i];
      dpdz[i] *= clip_factor[i];

      if (clip_factor[i] < 0.99) {

        max_factor = CS_MAX(max_factor, clip_factor[i]);
        min_factor = CS_MIN(min_factor, clip_factor[i]);
        n_clip++;

      }

    } /* End of loop on cells */

  } /* End if *imligp == 1 */

  /* Update min/max and n_clip in case of parallelism */

#if defined(HAVE_MPI)

  if (mesh->n_domains > 1) {

    assert(sizeof(cs_real_t) == sizeof(double));

    /* Global Max */

    MPI_Allreduce(&max_factor, &global_max_factor, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    max_factor = global_max_factor;

    /* Global min */

    MPI_Allreduce(&min_factor, &global_min_factor, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    min_factor = global_min_factor;

    /* Sum number of clippings */

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_GNUM,
                  MPI_SUM, cs_glob_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(HAVE_MPI) */

  /* Output warning if necessary */

  if (*iwarnp > 1)
    bft_printf(_(" Gradient limitation in %llu cells\n"
                 "   minimum factor = %14.5e; maximum factor = %14.5e\n"),
               (unsigned long long)n_clip, min_factor, max_factor);

  /* Synchronize dpdx, dpdy, dpdz */

  if (halo != NULL) {

    if (*idimtr > 0) {

      /* If the gradient is not treated as a "true" vector */

      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdx);
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdy);
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdz);

    }
    else {

      cs_halo_sync_var(halo, halo_type, dpdx);
      cs_halo_sync_var(halo, halo_type, dpdy);
      cs_halo_sync_var(halo, halo_type, dpdz);

      cs_halo_perio_sync_var_vect_ni(halo,
                                     halo_type,
                                     dpdx, dpdy, dpdz);

    }

  }

  BFT_FREE(buf);
}

/*----------------------------------------------------------------------------
 * Clip the gradient of a vector if necessary. This function deals with the
 * standard or extended neighborhood.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   clipping_type  <-- type of clipping for the computation of the gradient
 *   verbosity      <-- output level
 *   climgp         <-- clipping coefficient for the computation of the gradient
 *   pvar           <-- variable
 *   gradv          <-> gradient of pvar
 *   pvar           <-- variable
 *----------------------------------------------------------------------------*/

static void
_vector_gradient_clipping(const cs_mesh_t              *m,
                          const cs_mesh_quantities_t   *fvq,
                          const cs_halo_type_t          halo_type,
                          const cs_int_t                clipping_type,
                          const cs_int_t                verbosity,
                          const cs_real_t               climgp,
                          const cs_real_3_t   *restrict pvar,
                          cs_real_33_t        *restrict gradv)
{
  /* Local variables */

  cs_lnum_t  cell_id, cell_id1, cell_id2, cidx, face_id, i, j;
  cs_real_3_t dist, grad_dist1, grad_dist2;
  cs_real_t  dvar_sq, dist_sq1, dist_sq2;
  cs_real_t  global_min_factor, global_max_factor, factor1, factor2;

  cs_gnum_t  n_clip = 0, n_g_clip =0;
  cs_real_t  min_factor = 1;
  cs_real_t  max_factor = 0;
  cs_real_t  clipp_coef_sq = climgp*climgp;
  cs_real_t  *restrict buf = NULL, *restrict clip_factor = NULL;
  cs_real_t  *restrict denom = NULL, *restrict denum = NULL;

  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  const cs_halo_t *halo = m->halo;

  if (clipping_type < 0)
    return;

  /* The gradient and the variable must be already synchronized */

  /* Allocate and initialize working buffers */

  if (clipping_type == 1)
    BFT_MALLOC(buf, 3*n_cells_ext, cs_real_t);
  else
    BFT_MALLOC(buf, 2*n_cells_ext, cs_real_t);

  denum = buf;
  denom = buf + n_cells_ext;

  if (clipping_type == 1)
    clip_factor = buf + 2*n_cells_ext;

  /* Initialization */

  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    denum[cell_id] = 0;
    denom[cell_id] = 0;
    if (clipping_type == 1)
      clip_factor[cell_id] = (cs_real_t)DBL_MAX;
  }

  /* Remark:
     denum: holds the maximum l2 norm of the variation of the gradient squared
     denom: holds the maximum l2 norm of the variation of the variable squared */

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clipping_type == 0) {

    for (face_id = 0; face_id < n_i_faces; face_id++) {

      cell_id1 = i_face_cells[face_id][0] - 1;
      cell_id2 = i_face_cells[face_id][1] - 1;

      for (i = 0; i < 3; i++)
        dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

      for (i = 0; i < 3; i++) {

        grad_dist1[i] = gradv[cell_id1][0][i] * dist[0]
                      + gradv[cell_id1][1][i] * dist[1]
                      + gradv[cell_id1][2][i] * dist[2];

        grad_dist2[i] = gradv[cell_id2][0][i] * dist[0]
                      + gradv[cell_id2][1][i] * dist[1]
                      + gradv[cell_id2][2][i] * dist[2];

      }

      dist_sq1 = grad_dist1[0]*grad_dist1[0]
               + grad_dist1[1]*grad_dist1[1]
               + grad_dist1[2]*grad_dist1[2];

      dist_sq2 = grad_dist2[0]*grad_dist2[0]
               + grad_dist2[1]*grad_dist2[1]
               + grad_dist2[2]*grad_dist2[2];

      dvar_sq = (pvar[cell_id1][0]-pvar[cell_id2][0])
               *(pvar[cell_id1][0]-pvar[cell_id2][0])
              + (pvar[cell_id1][1]-pvar[cell_id2][1])
               *(pvar[cell_id1][1]-pvar[cell_id2][1])
              + (pvar[cell_id1][2]-pvar[cell_id2][2])
               *(pvar[cell_id1][2]-pvar[cell_id2][2]);

      denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
      denum[cell_id2] = CS_MAX(denum[cell_id2], dist_sq2);
      denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);
      denom[cell_id2] = CS_MAX(denom[cell_id2], dvar_sq);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
        for (cidx = cell_cells_idx[cell_id1];
             cidx < cell_cells_idx[cell_id1+1];
             cidx++) {

          cell_id2 = cell_cells_lst[cidx-1] - 1;

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++)
            grad_dist1[i] = gradv[cell_id1][0][i] * dist[0]
                          + gradv[cell_id1][1][i] * dist[1]
                          + gradv[cell_id1][2][i] * dist[2];


          dist_sq1 = grad_dist1[0]*grad_dist1[0]
                   + grad_dist1[1]*grad_dist1[1]
                   + grad_dist1[2]*grad_dist1[2];

          dvar_sq = (pvar[cell_id1][0]-pvar[cell_id2][0])
                   *(pvar[cell_id1][0]-pvar[cell_id2][0])
                  + (pvar[cell_id1][1]-pvar[cell_id2][1])
                   *(pvar[cell_id1][1]-pvar[cell_id2][1])
                  + (pvar[cell_id1][2]-pvar[cell_id2][2])
                   *(pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);

        }
      }

    } /* End for extended halo */

  }

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clipping_type == 1) {

    for (face_id = 0; face_id < n_i_faces; face_id++) {

      cell_id1 = i_face_cells[face_id][0] - 1;
      cell_id2 = i_face_cells[face_id][1] - 1;

      for (i = 0; i < 3; i++)
        dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

      for (i = 0; i < 3; i++)
        grad_dist1[i] = 0.5*(
                       (gradv[cell_id1][0][i]+gradv[cell_id2][0][i])*dist[0]
                      +(gradv[cell_id1][1][i]+gradv[cell_id2][1][i])*dist[1]
                      +(gradv[cell_id1][2][i]+gradv[cell_id2][2][i])*dist[2]);

      dist_sq1 = grad_dist1[0]*grad_dist1[0]
               + grad_dist1[1]*grad_dist1[1]
               + grad_dist1[2]*grad_dist1[2];

      dvar_sq = (pvar[cell_id1][0]-pvar[cell_id2][0])
               *(pvar[cell_id1][0]-pvar[cell_id2][0])
              + (pvar[cell_id1][1]-pvar[cell_id2][1])
               *(pvar[cell_id1][1]-pvar[cell_id2][1])
              + (pvar[cell_id1][2]-pvar[cell_id2][2])
               *(pvar[cell_id1][2]-pvar[cell_id2][2]);

      denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
      denum[cell_id2] = CS_MAX(denum[cell_id2], dist_sq1);
      denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);
      denom[cell_id2] = CS_MAX(denom[cell_id2], dvar_sq);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
        for (cidx = cell_cells_idx[cell_id1];
             cidx < cell_cells_idx[cell_id1+1];
             cidx++) {

          cell_id2 = cell_cells_lst[cidx-1] - 1;

          for (i = 0; i < 3; i++)
            dist[i] = cell_cen[cell_id1][i] - cell_cen[cell_id2][i];

          for (i = 0; i < 3; i++)
            grad_dist1[i] = 0.5*(
                           (gradv[cell_id1][0][i]+gradv[cell_id2][0][i])*dist[0]
                          +(gradv[cell_id1][1][i]+gradv[cell_id2][1][i])*dist[1]
                          +(gradv[cell_id1][2][i]+gradv[cell_id2][2][i])*dist[2]);

          dist_sq1 = grad_dist1[0]*grad_dist1[0]
                   + grad_dist1[1]*grad_dist1[1]
                   + grad_dist1[2]*grad_dist1[2];

          dvar_sq = (pvar[cell_id1][0]-pvar[cell_id2][0])
                   *(pvar[cell_id1][0]-pvar[cell_id2][0])
                  + (pvar[cell_id1][1]-pvar[cell_id2][1])
                   *(pvar[cell_id1][1]-pvar[cell_id2][1])
                  + (pvar[cell_id1][2]-pvar[cell_id2][2])
                   *(pvar[cell_id1][2]-pvar[cell_id2][2]);

          denum[cell_id1] = CS_MAX(denum[cell_id1], dist_sq1);
          denom[cell_id1] = CS_MAX(denom[cell_id1], dvar_sq);

        }
      }

    } /* End for extended neighborhood */

    /* Synchronize variable */

    if (halo != NULL) {
      cs_mesh_sync_var_scal(denom);
      cs_mesh_sync_var_scal(denum);
    }

  } /* End if clipping_type == 1 */

  /* Clipping of the gradient if denum/denom > climgp**2 */

  /* First clipping Algorithm: based on the cell gradient */
  /*------------------------------------------------------*/

  if (clipping_type == 0) {

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      if (denum[cell_id] > clipp_coef_sq * denom[cell_id]) {

        factor1 = sqrt(clipp_coef_sq * denom[cell_id]/denum[cell_id]);


        for (i = 0; i < 3; i++)
          for (j = 0; j < 3; j++)
            gradv[cell_id][i][j] *= factor1;

        min_factor = CS_MIN( factor1, min_factor);
        max_factor = CS_MAX( factor1, max_factor);
        n_clip++;

      } /* If clipping */

    } /* End of loop on cells */

  }

  /* Second clipping Algorithm: based on the face gradient */
  /*-------------------------------------------------------*/

  else if (clipping_type == 1) {

    for (face_id = 0; face_id < n_i_faces; face_id++) {

      cell_id1 = i_face_cells[face_id][0] - 1;
      cell_id2 = i_face_cells[face_id][1] - 1;

      factor1 = 1.0;
      if (denum[cell_id1] > clipp_coef_sq * denom[cell_id1])
        factor1 = sqrt(clipp_coef_sq * denom[cell_id1]/denum[cell_id1]);

      factor2 = 1.0;
      if (denum[cell_id2] > clipp_coef_sq * denom[cell_id2])
        factor2 = sqrt(clipp_coef_sq * denom[cell_id2]/denum[cell_id2]);

      min_factor = CS_MIN(factor1, factor2);

      clip_factor[cell_id1] = CS_MIN( clip_factor[cell_id1], min_factor);
      clip_factor[cell_id2] = CS_MIN( clip_factor[cell_id2], min_factor);

    } /* End of loop on faces */

    /* Complement for extended neighborhood */

    if (cell_cells_idx != NULL && halo_type == CS_HALO_EXTENDED) {

      for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {

        min_factor = 1.0;

        for (cidx = cell_cells_idx[cell_id1]; cidx < cell_cells_idx[cell_id1+1]; cidx++) {

          cell_id2 = cell_cells_lst[cidx-1] - 1;
          factor2 = 1.0;

          if (denum[cell_id2] > clipp_coef_sq * denom[cell_id2])
            factor2 = sqrt(clipp_coef_sq * denom[cell_id2]/denum[cell_id2]);

          min_factor = CS_MIN(min_factor, factor2);

        }

        clip_factor[cell_id1] = CS_MIN(clip_factor[cell_id1], min_factor);

      } /* End of loop on cells */

    } /* End for extended neighborhood */

    for (cell_id = 0; cell_id < n_cells; cell_id++) {


      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          gradv[cell_id][i][j] *= clip_factor[cell_id];

      if (clip_factor[cell_id] < 0.99) {

        max_factor = CS_MAX(max_factor, clip_factor[cell_id]);
        min_factor = CS_MIN(min_factor, clip_factor[cell_id]);
        n_clip++;

      }

    } /* End of loop on cells */

  } /* End if clipping_type == 1 */

  /* Update min/max and n_clip in case of parallelism */
  /*--------------------------------------------------*/

#if defined(HAVE_MPI)

  if (m->n_domains > 1) {

    assert(sizeof(cs_real_t) == sizeof(double));

    /* Global Max */

    MPI_Allreduce(&max_factor, &global_max_factor, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    max_factor = global_max_factor;

    /* Global min */

    MPI_Allreduce(&min_factor, &global_min_factor, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    min_factor = global_min_factor;

    /* Sum number of clippings */

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_GNUM,
                  MPI_SUM, cs_glob_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(HAVE_MPI) */

  /* Output warning if necessary */

  if (verbosity > 1)
    bft_printf(_(" Gradient of a vector limitation in %llu cells\n"
                 "   minimum factor = %14.5e; maximum factor = %14.5e\n"),
               (unsigned long long)n_clip, min_factor, max_factor);

  /* Synchronize the updated Gradient */

  if (halo != NULL)
    cs_mesh_sync_var_tens((cs_real_t *)gradv);

  BFT_FREE(buf);
}

/*----------------------------------------------------------------------------
 * Initialize the gradient of a vector for gradient reconstruction.
 *
 * A non-reconstructed gradient is computed at this stage.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *----------------------------------------------------------------------------*/

static void
_initialize_vector_gradient(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            const cs_int_t                inc,
                            const cs_real_3_t   *restrict coefav,
                            const cs_real_33_t  *restrict coefbv,
                            const cs_real_3_t   *restrict pvar,
                            cs_real_33_t        *restrict gradv)
{
  /* Local variables */

  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;
  const int n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;

  cs_lnum_t  cell_id, face_id, i, j, cell_id1, cell_id2;
  cs_real_t  pfac, pond, epzero;
  cs_real_t  dvol, dvol1, dvol2;

  epzero = 1.e-12;

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (m->halo != NULL)
    cs_mesh_sync_var_vect((cs_real_t *)pvar);

  /* Computation without reconstruction */
  /*------------------------------------*/

  /* Initialization */

  for (cell_id = 0; cell_id < n_cells_ext; cell_id++)
    for (j = 0; j < 3; j++)
      for (i = 0; i < 3; i++)
        gradv[cell_id][j][i] = 0.0;

  /* Interior face treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {
    cell_id1 = i_face_cells[face_id][0] - 1;
    cell_id2 = i_face_cells[face_id][1] - 1;

    pond = weight[face_id];
    dvol1 = 1./cell_vol[cell_id1];
    dvol2 = 1./cell_vol[cell_id2];

    for (i = 0; i < 3; i++) {
      pfac   = pond * pvar[cell_id1][i] + (1.0-pond) * pvar[cell_id2][i];
      for (j = 0; j < 3; j++) {
        gradv[cell_id1][j][i] += pfac * i_face_normal[face_id][j] * dvol1;
        gradv[cell_id2][j][i] -= pfac * i_face_normal[face_id][j] * dvol2;
      }
    }
  }

  /* Boundar face treatment */

  for (face_id = 0; face_id < n_b_faces; face_id++) {
    cell_id = b_face_cells[face_id] - 1;

    dvol = 1./cell_vol[cell_id];

    for (i = 0; i < 3; i++) {
      pfac = inc*coefav[face_id][i] + coefbv[face_id][0][i] * pvar[cell_id][0]
                                    + coefbv[face_id][1][i] * pvar[cell_id][1]
                                    + coefbv[face_id][2][i] * pvar[cell_id][2];
      for (j = 0; j < 3; j++)
        gradv[cell_id][j][i] += pfac * b_face_normal[face_id][j]*dvol;

    }
  }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL)
    cs_mesh_sync_var_tens((cs_real_t *)gradv);

}
/*----------------------------------------------------------------------------
 * Compute the gradient of a vector with an iterative technique in order to
 * handle non-orthoganalities (nswrgp > 1).
 *
 * We do not take into account any volumic force here.
 *
 * Compute cocg at the first call and if needed.
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   var_num        <-- variable's number (0 if non-solved variable)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   nswrgp         --> >1: with reconstruction
 *   verbosity      --> verbosity level
 *   epsrgp         --> precision for iterative gradient calculation
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   gradv          <-> gradient of pvar
 *----------------------------------------------------------------------------*/

static void
_iterative_vector_gradient(const cs_mesh_t              *m,
                           const cs_mesh_quantities_t   *fvq,
                           int                           var_num,
                           int                           inc,
                           int                           nswrgp,
                           int                           verbosity,
                           double                        epsrgp,
                           const cs_real_3_t   *restrict coefav,
                           const cs_real_33_t  *restrict coefbv,
                           const cs_real_3_t   *restrict pvar,
                           cs_real_33_t        *restrict gradv)
{
  /* Local variables */

  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;
  const int n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;
  cs_real_33_t *restrict cocg = fvq->cocg_it;

  cs_real_33_t *rhs;

  cs_lnum_t  cell_id, face_id, i, j, k, cell_id1, cell_id2;
  cs_real_t  pfac, l2_norm, l2_residual, vecfac, pond, epzero;
  cs_real_t  dvol, dvol1, dvol2;

  int isweep;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_33_t);

  epzero = 1.e-12;

  /* Gradient reconstruction to handle non-orthogonal meshes */
  /*---------------------------------------------------------*/

  /* L2 norm */

  l2_norm = sqrt(cs_dot(9*n_cells, (cs_real_t *)gradv, (cs_real_t *)gradv));
  l2_residual = l2_norm;

  if (l2_norm > epzero) {

    /* Iterative process */
    /*-------------------*/

    for (isweep = 1; isweep < nswrgp && l2_residual > epsrgp*l2_norm; isweep++) {

      /* Computation of the Right Hand Side*/

      for (cell_id = 0; cell_id < n_cells_ext; cell_id++)
        for (j = 0; j < 3; j++)
          for (i = 0; i < 3; i++)
            rhs[cell_id][j][i] = -gradv[cell_id][j][i];


      /* Interior face treatment */

      for (face_id = 0; face_id < n_i_faces; face_id++) {

        cell_id1 = i_face_cells[face_id][0] - 1;
        cell_id2 = i_face_cells[face_id][1] - 1;
        pond = weight[face_id];

        dvol1 = 1./cell_vol[cell_id1];
        dvol2 = 1./cell_vol[cell_id2];

        for (i = 0; i < 3; i++) {
          pfac = pond*pvar[cell_id1][i] + (1.0-pond)*pvar[cell_id2][i];

          for (k = 0; k < 3; k++)
            pfac += 0.5*(gradv[cell_id1][k][i] + gradv[cell_id2][k][i])
                  * dofij[face_id][k];

          for (j = 0; j < 3; j++) {

            rhs[cell_id1][j][i] += pfac * i_face_normal[face_id][j] * dvol1;
            rhs[cell_id2][j][i] -= pfac * i_face_normal[face_id][j] * dvol2;

          }
        }
      }

      /* Boundar face treatment */

      for (face_id = 0; face_id < n_b_faces; face_id++) {

        cell_id = b_face_cells[face_id] - 1;
        dvol = 1./cell_vol[cell_id];

        for (i = 0; i < 3; i++) {

          pfac = inc*coefav[face_id][i];

          for (k = 0; k < 3; k++) {

            vecfac =  pvar[cell_id][k]
                   + gradv[cell_id][0][k] * diipb[face_id][0]
                   + gradv[cell_id][1][k] * diipb[face_id][1]
                   + gradv[cell_id][2][k] * diipb[face_id][2];
            pfac += coefbv[face_id][k][i] * vecfac;

          }

          for (j = 0; j < 3; j++)
            rhs[cell_id][j][i] += pfac * b_face_normal[face_id][j] * dvol;

        }
      }

      /* Increment of the gradient */

      for (cell_id = 0; cell_id < n_cells; cell_id++)
        for (j = 0; j < 3; j++)
          for (i = 0; i < 3; i++)
            for (k = 0; k < 3; k++)
              gradv[cell_id][j][i] += rhs[cell_id][k][i] * cocg[cell_id][k][j];

      /* Periodicity and parallelism treatment */

      if (m->halo != NULL)
        cs_mesh_sync_var_tens((cs_real_t *)gradv);

      /* Convergence test (L2 norm) */

      l2_residual = sqrt(cs_dot(9*n_cells, (cs_real_t *)rhs, (cs_real_t *)rhs));

    } /* End of the iterative process */

    /* Printing */

    if (l2_residual < epsrgp*l2_norm) {
      if (verbosity >= 2) {
        bft_printf
          (_(" %s: isweep = %d, residue norm: %e, norm: %e, var_num = %d\n"),
           __func__, isweep, l2_residual/l2_norm, l2_norm, var_num);
      }
    }
    else if (isweep >= nswrgp) {
      if (verbosity >= 0) {
        bft_printf(" %s: isweep = %d, residu norm: %e, norm: %e, var_num = %d\n",
                   __func__, isweep, l2_residual/l2_norm, l2_norm, var_num);
        bft_printf("@ @@ warning: non convergence of grdvec\n");
      }
    }
  }

  BFT_FREE(rhs);

}

/*----------------------------------------------------------------------------
 * Compute cell gradient of a vector using least-squares reconstruction for
 * non-orthogonal meshes (nswrgp > 1).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   gradv          <-> gradient of pvar
 *----------------------------------------------------------------------------*/

static void
_lsq_vector_gradient(const cs_mesh_t              *m,
                     const cs_mesh_quantities_t   *fvq,
                     const cs_halo_type_t          halo_type,
                     const cs_int_t                inc,
                     const cs_real_3_t   *restrict coefav,
                     const cs_real_33_t  *restrict coefbv,
                     const cs_real_3_t   *restrict pvar,
                     cs_real_33_t        *restrict gradv)
{
  /* Local variables */

  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)m->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)m->cell_cells_lst;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  cs_real_33_t *restrict cocg = fvq->cocg_lsq;

  cs_lnum_t  cell_id, cidx, face_id, cell_id1, cell_id2, i, j, k;
  int        g_id, t_id;
  cs_real_t  pfac, ddc;
  cs_real_3_t  dc;
  cs_real_4_t  fctb;

  cs_real_33_t *rhs;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_33_t);

  /* By default, handle the gradient as a tensor
     (i.e. we assume it is the gradient of a vector field) */

  if (m->halo != NULL)
    cs_mesh_sync_var_vect((cs_real_t *)pvar);

  /* Compute Right-Hand Side */
  /*-------------------------*/

    for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          rhs[cell_id][i][j] = 0.0;
    }

  /* Contribution from interior faces */

  for (g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for private(face_id, cell_id1, cell_id2,\
                                    i, j, pfac, dc, fctb, ddc)
    for (t_id = 0; t_id < n_i_threads; t_id++) {

      for (face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cell_id1 = i_face_cells[face_id][0] - 1;
        cell_id2 = i_face_cells[face_id][1] - 1;

        for (i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (i = 0; i < 3; i++) {
          pfac =  (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

          for (j = 0; j < 3; j++) {
            fctb[j] = dc[j] * pfac;
            rhs[cell_id1][j][i] += fctb[j];
            rhs[cell_id2][j][i] += fctb[j];
          }
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */

  /* Contribution from extended neighborhood */

  if (halo_type == CS_HALO_EXTENDED) {

#   pragma omp parallel for private(cidx, cell_id2, dc, pfac, ddc, i, j)
    for (cell_id1 = 0; cell_id1 < n_cells; cell_id1++) {
      for (cidx = cell_cells_idx[cell_id1];
           cidx < cell_cells_idx[cell_id1+1];
           cidx++) {

        cell_id2 = cell_cells_lst[cidx - 1] - 1;

        for (i = 0; i < 3; i++)
          dc[i] = cell_cen[cell_id2][i] - cell_cen[cell_id1][i];

        ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

        for (i = 0; i < 3; i++) {

          pfac = (pvar[cell_id2][i] - pvar[cell_id1][i]) * ddc;

          for (j = 0; j < 3; j++) {
            rhs[cell_id1][j][i] += dc[j] * pfac;
          }
        }
      }
    }

  } /* End for extended neighborhood */

  /* Contribution from boundary faces */

  for (g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for private(face_id, cell_id1, i, j, pfac, dc, ddc)
    for (t_id = 0; t_id < n_b_threads; t_id++) {

      for (face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cell_id1 = b_face_cells[face_id] - 1;

        for (i = 0; i < 3; i++)
          dc[i] = b_face_cog[face_id][i] - cell_cen[cell_id1][i];

        ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);


        for (i = 0; i < 3; i++) {
          pfac = (coefav[face_id][i]*inc
               + ( coefbv[face_id][0][i] * pvar[cell_id1][0]
                 + coefbv[face_id][1][i] * pvar[cell_id1][1]
                 + coefbv[face_id][2][i] * pvar[cell_id1][2]
                 -                         pvar[cell_id1][i])) * ddc;

          for (j = 0; j < 3; j++)
            rhs[cell_id1][j][i] += dc[j] * pfac;
        }

      } /* loop on faces */

    } /* loop on threads */

  } /* loop on thread groups */


  /* Compute gradient */
  /*------------------*/

  for (cell_id = 0; cell_id < n_cells; cell_id++)
    for (j = 0; j < 3; j++)
      for (i = 0; i < 3; i++) {

        gradv[cell_id][j][i] = 0.0;

        for (k = 0; k < 3; k++)
          gradv[cell_id][j][i] += rhs[cell_id][k][i] * cocg[cell_id][k][j];

      }

  /* Periodicity and parallelism treatment */

  if (m->halo != NULL)
    cs_mesh_sync_var_tens((cs_real_t *)gradv);

  BFT_FREE(rhs);
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Fortran Interface :
 *
 * SUBROUTINE CGDCEL
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdcel, CGDCEL)
(
 const cs_int_t   *const ncelet,      /* --> number of extended cells         */
 const cs_int_t   *const ncel,        /* --> number of cells                  */
 const cs_int_t   *const nfac,        /* --> number of internal faces         */
 const cs_int_t   *const nfabor,      /* --> number of boundary faces         */
 const cs_int_t   *const ncelbr,      /* --> number of cells on boundary      */
 const cs_int_t   *const ivar,
 const cs_int_t   *const imrgra,      /* --> gradient computation mode        */
 const cs_int_t   *const inc,         /* --> 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* --> 1 or 0: recompute COCG or not    */
 const cs_int_t   *const nswrgp,      /* --> >1: with reconstruction          */
 const cs_int_t   *const idimtr,      /* --> 0, 1, 2: scalar, vector, tensor
                                             in case of rotation              */
 const cs_int_t   *const iphydp,      /* --> use hydrosatatic pressure        */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
 const cs_int_t   *const nfecra,      /* --> standard output unit             */
 const cs_int_t   *const imligp,      /* --> type of clipping                 */
 const cs_real_t  *const epsrgp,      /* --> precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* --> extrapolate gradient at boundary */
 const cs_real_t  *const climgp,      /* --> clipping coefficient             */
 const cs_int_t          ifacel[],
 const cs_int_t          ifabor[],
 const cs_int_t          icelbr[],    /* --> list of cells on boundary        */
 const cs_int_t          isympa[],    /* --> indicator for symmetry faces     */
 const cs_real_t         volume[],
 const cs_real_t         surfac[],
 const cs_real_t         surfbo[],
 const cs_real_t         surfbn[],
 const cs_real_t         pond[],      /* --> interior faces geometric weight  */
 const cs_real_t         dist[],      /* --> interior faces I' to J' distance */
 const cs_real_t         distbr[],    /* --> boundary faces I' to J' distance */
 const cs_real_t         dijpf[],     /* --> interior faces I'J' vector       */
 const cs_real_t         diipb[],     /* --> boundary faces II' vector        */
 const cs_real_t         dofij[],
       cs_real_t         fextx[],     /* --> components of the exterior force */
       cs_real_t         fexty[],     /*     generating the hydrostatic       */
       cs_real_t         fextz[],     /*     pressure                         */
 const cs_real_t         xyzcen[],
 const cs_real_t         cdgfac[],
 const cs_real_t         cdgfbo[],
 const cs_real_t         coefap[],    /* --> boundary condition term          */
 const cs_real_t         coefbp[],    /* --> boundary condition term          */
       cs_real_t         pvar[],      /* --> gradient's base variable         */
       cs_real_t         cocgb[],     /* <-> contribution to COCG of cells
                                             on boundary's internal faces     */
       cs_real_t         cocg[],      /* <-> contribution to COCG of cells
                                             on boundary's boundary faces     */
       cs_real_t         cocib[],     /* <-> contribution to COCG of cells
                                             on boundary's internal faces     */
       cs_real_t         coci[],      /* <-> contribution to COCG of cells
                                             on boundary's boundary faces     */
       cs_real_t         grad[]       /* <-- gradient x component             */
)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_halo_t  *halo = mesh->halo;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  cs_int_t  *ipcvse = NULL;
  cs_int_t  *ielvse = NULL;

  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict bx, *restrict by, *restrict bz;
  cs_real_t  *restrict dpdx, *restrict dpdy, *restrict dpdz;

  bool update_stats = true;
  cs_gradient_type_t gradient_type = CS_GRADIENT_N_TYPES;

  /* Allocate work arrays */

  BFT_MALLOC(_aux_vectors, (*ncelet)*3, cs_real_t);

  bx = _aux_vectors;
  by = _aux_vectors +  *ncelet;
  bz = _aux_vectors + (*ncelet)*2;

  /* Associate gradient component arrays */

  dpdx = grad;
  dpdy = grad +  *ncelet;
  dpdz = grad + (*ncelet)*2;

  /* Choose gradient type */

  switch (*imrgra) {
  case 0: gradient_type = CS_GRADIENT_ITER; break;
  case 1: gradient_type = CS_GRADIENT_LSQ_STD; break;
  case 2: gradient_type = CS_GRADIENT_LSQ_EXT; break;
  case 3: gradient_type = CS_GRADIENT_LSQ_EXT_RED; break;
  case 4: gradient_type = CS_GRADIENT_LSQ_ITER; break;
  default: break;
  }

  if (update_stats == true) {
    char var_name[32];
    snprintf(var_name, 31, "Var. %2d", *ivar); var_name[31] = '\0';
    t0 = cs_timer_time();
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  if (*imrgra == 2 || *imrgra ==  3)
    halo_type = CS_HALO_EXTENDED;

  /* Synchronize variable */

  if (*imrgra != 0 && halo != NULL) {

    if (*idimtr > 0)
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, pvar);
    else
      cs_halo_sync_var(halo, halo_type, pvar);

    /* TODO: check if fext* components are all up to date, in which
     *       case we need no special treatment for *idimtr > 0 */

    if (*iphydp != 0) {

      if (*idimtr > 0){
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fextx);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fexty);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fextz);
      }
      else {
        cs_halo_sync_var(halo, halo_type, fextx);
        cs_halo_sync_var(halo, halo_type, fexty);
        cs_halo_sync_var(halo, halo_type, fextz);
        cs_halo_perio_sync_var_vect_ni(halo, halo_type, fextx, fexty, fextz);
      }
    }

  }

  /* "cell -> cells" connectivity for the extended neighborhood */

  ipcvse = mesh->cell_cells_idx;
  ielvse = mesh->cell_cells_lst;

  /* Compute gradient */

  if (*imrgra == 0) {

    CS_PROCF (gradrc, GRADRC)
      (ncelet, ncel  , nfac  , nfabor, ncelbr,
       imrgra, inc   , iccocg, nswrgp, idimtr, iphydp,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ivar  ,
       volume, surfac, surfbo, pond  , xyzcen, cdgfac, cdgfbo,
       dijpf , diipb , dofij , fextx , fexty , fextz ,
       coefap, coefbp, pvar  ,
       cocgb , cocg  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

  }
  else if (*imrgra == 1 || *imrgra == 2 || *imrgra == 3) {

    CS_PROCF(gradmc, GRADMC)
      (ncelet, ncel  , nfac  , nfabor, ncelbr,
       inc   , iccocg, nswrgp, idimtr, iphydp, imrgra,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ipcvse, ielvse, isympa,
       volume, surfac, surfbo, surfbn, pond  ,
       dist  , distbr, dijpf , diipb ,
       fextx , fexty , fextz , xyzcen,
       cdgfac, cdgfbo, coefap, coefbp, pvar  ,
       cocgb , cocg  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

  }
  else if (*imrgra == 4) {

    const cs_int_t  _imlini = 1;
    const cs_real_t _climin = 1.5;

    CS_PROCF(gradmc, GRADMC)
      (ncelet, ncel  , nfac  , nfabor, ncelbr,
       inc   , iccocg, nswrgp, idimtr, iphydp, imrgra,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ipcvse, ielvse, isympa,
       volume, surfac, surfbo, surfbn, pond  ,
       dist  , distbr, dijpf , diipb ,
       fextx , fexty , fextz , xyzcen,
       cdgfac, cdgfbo, coefap, coefbp, pvar  ,
       cocgb , cocg  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

    _scalar_gradient_clipping(imrgra,
                             &_imlini,
                              iwarnp,
                              idimtr,
                             &_climin,
                              pvar,
                              dpdx,
                              dpdy,
                              dpdz);

    CS_PROCF (gradrc, GRADRC)
      (ncelet, ncel  , nfac  , nfabor, ncelbr,
       imrgra, inc   , iccocg, nswrgp, idimtr, iphydp,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ivar  ,
       volume, surfac, surfbo, pond  , xyzcen, cdgfac, cdgfbo,
       dijpf , diipb , dofij , fextx , fexty , fextz ,
       coefap, coefbp, pvar  ,
       cocib , coci  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

  }

  _scalar_gradient_clipping(imrgra,
                            imligp,
                            iwarnp,
                            idimtr,
                            climgp,
                            pvar,
                            dpdx,
                            dpdy,
                            dpdz);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }

  BFT_FREE(_aux_vectors);
}

/*----------------------------------------------------------------------------
 *
 * Fortran Interface :
 *
 * SUBROUTINE CGDVEC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgdvec, CGDVEC)
(
 const cs_int_t         *const ivar,
 const cs_int_t         *const imrgra,  /* --> gradient computation mode      */
 const cs_int_t         *const inc,     /* --> 0 or 1: increment or not       */
 const cs_int_t         *const nswrgp,  /* --> >1: with reconstruction        */
 const cs_int_t         *const iwarnp,  /* --> verbosity level                */
 const cs_int_t         *const imligp,  /* --> type of clipping               */
 const cs_real_t        *const epsrgp,  /* --> precision for iterative gradient
                                               calculation                    */
 const cs_real_t        *const climgp,  /* --> clipping coefficient           */
 const cs_real_3_t   *restrict coefav,  /* --> boundary condition term        */
 const cs_real_33_t  *restrict coefbv,  /* --> boundary condition term        */
 const cs_real_3_t   *restrict pvar,    /* --> gradient's base variable       */
       cs_real_33_t  *restrict gradv    /* <-- gradient of the variable       */
)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_gradient_info_t *gradient_info = NULL;
  cs_timer_t t0, t1;

  cs_int_t  *ipcvse = NULL;
  cs_int_t  *ielvse = NULL;

  bool update_stats = true;
  cs_gradient_type_t gradient_type = CS_GRADIENT_N_TYPES;


  switch (*imrgra) {
  case 0: gradient_type = CS_GRADIENT_ITER; break;
  case 1: gradient_type = CS_GRADIENT_LSQ_STD; break;
  case 2: gradient_type = CS_GRADIENT_LSQ_EXT; break;
  case 3: gradient_type = CS_GRADIENT_LSQ_EXT_RED; break;
  case 4: gradient_type = CS_GRADIENT_LSQ_ITER; break;
  default: break;
  }

  if (update_stats == true) {
    char var_name[32];
    snprintf(var_name, 31, "Var. %2d", *ivar); var_name[31] = '\0';
    t0 = cs_timer_time();
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  if (*imrgra == 2 || *imrgra ==  3)
    halo_type = CS_HALO_EXTENDED;

  /* "cell -> cells" connectivity for the extended neighborhood */

  ipcvse = mesh->cell_cells_idx;
  ielvse = mesh->cell_cells_lst;

  /* Compute gradient */

  if (*imrgra == 0) {

    _initialize_vector_gradient(mesh,
                                fvq,
                               *inc,
                                coefav,
                                coefbv,
                                pvar,
                                gradv);

    /* If reconstructions are required */

    if (*nswrgp > 1)
      _iterative_vector_gradient(mesh,
                                 fvq,
                                *ivar,
                                *inc,
                                *nswrgp,
                                *iwarnp,
                                *epsrgp,
                                 coefav,
                                 coefbv,
                                 pvar,
                                 gradv);

  }
  else if (*imrgra == 1 || *imrgra == 2 || *imrgra == 3) {

    /* If NO reconstruction are required */

    if (*nswrgp <= 1)
      _initialize_vector_gradient(mesh,
                                  fvq,
                                 *inc,
                                  coefav,
                                  coefbv,
                                  pvar,
                                  gradv);

    /* Reconstruction with Least square method */

    else
      _lsq_vector_gradient(mesh,
                           fvq,
                           halo_type,
                          *inc,
                           coefav,
                           coefbv,
                           pvar,
                           gradv);

  }
  else if (*imrgra == 4) {

    /* Clipping algorithm and clipping factor */

    const cs_int_t  _imlini = 1;
    const cs_real_t _climin = 1.5;

    /* Initialization by the Least square method */

    _lsq_vector_gradient(mesh,
                         fvq,
                         halo_type,
                        *inc,
                         coefav,
                         coefbv,
                         pvar,
                         gradv);

    _vector_gradient_clipping(mesh,
                              fvq,
                              halo_type,
                             _imlini,
                             *iwarnp,
                             _climin,
                              pvar,
                              gradv);

    _iterative_vector_gradient(mesh,
                               fvq,
                              *ivar,
                              *inc,
                              *nswrgp,
                              *iwarnp,
                              *epsrgp,
                               coefav,
                               coefbv,
                               pvar,
                               gradv);

  }

   _vector_gradient_clipping(mesh,
                             fvq,
                             halo_type,
                            *imligp,
                            *iwarnp,
                            *climgp,
                             pvar,
                             gradv);

  if (update_stats == true) {
    gradient_info->n_calls += 1;
    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(gradient_info->t_tot), &t0, &t1);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_initialize(void)
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  assert(mesh != NULL);
}

/*----------------------------------------------------------------------------
 * Finalize gradient computation API.
 *----------------------------------------------------------------------------*/

void
cs_gradient_finalize(void)
{
  int ii;

  /* Free system info */

  for (ii = 0; ii < cs_glob_gradient_n_systems; ii++) {
    _gradient_info_dump(cs_glob_gradient_systems[ii]);
    _gradient_info_destroy(&(cs_glob_gradient_systems[ii]));
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  BFT_FREE(cs_glob_gradient_systems);

  cs_glob_gradient_n_systems = 0;
  cs_glob_gradient_n_max_systems = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
