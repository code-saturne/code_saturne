/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_perio.h"
#include "cs_prototypes.h"

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

  double               wt_tot;             /* Total wall-clock time used */
  double               cpu_tot;            /* Total (local) CPU used */

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

  new_info->wt_tot = 0.0;
  new_info->cpu_tot = 0.0;

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

  bft_printf(_("\n"
               "Summary of gradient computations pour \"%s\" (%s):\n\n"
               "  Number of calls:                  %d\n"
               "  Total elapsed time:               %12.3f\n"),
             this_info->name, cs_gradient_type_name[this_info->type],
             n_calls, this_info->wt_tot);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double cpu_min, cpu_max;
    double cpu_loc = this_info->cpu_tot;

    MPI_Allreduce(&cpu_loc, &cpu_min, 1, MPI_DOUBLE, MPI_MIN, cs_glob_mpi_comm);
    MPI_Allreduce(&cpu_loc, &cpu_max, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    bft_printf(_("  Min local total CPU time:         %12.3f\n"
                 "  Max local total CPU time:         %12.3f\n"),
               cpu_min, cpu_max);

  }

#endif

  if (cs_glob_n_ranks == 1)
    bft_printf(_("  Total CPU time:                   %12.3f\n"),
               this_info->cpu_tot);
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
 * Clip the gradient if necessary. This function deals with the standard or
 * extended neighborhood.
 *
 * parameters:
 *   imrgra         --> type of computation for the gradient
 *   imligp         --> type of clipping for the computation of the gradient
 *   iwarnp         --> output level
 *   itenso         --> for rotational periodicity
 *   climgp         --> clipping coefficient for the computation of the gradient
 *   var            --> variable
 *   dpdx           --> X component of the pressure gradient
 *   dpdy           --> Y component of the pressure gradient
 *   dpdz           --> Z component of the pressure gradient
 *----------------------------------------------------------------------------*/

static void
_gradient_clipping(const cs_int_t   *imrgra,
                   const cs_int_t   *imligp,
                   const cs_int_t   *iwarnp,
                   const cs_int_t   *itenso,
                   const cs_real_t  *climgp,
                   cs_real_t         var[],
                   cs_real_t         dpdx[],
                   cs_real_t         dpdy[],
                   cs_real_t         dpdz[])
{
  cs_int_t  i, i1, i2, j, k;
  cs_real_t  dist[3];
  cs_real_t  dvar, dist1, dist2, dpdxf, dpdyf, dpdzf;
  cs_real_t  global_min_factor, global_max_factor, factor1, factor2;

  cs_int_t  n_clip = 0, n_g_clip =0;
  cs_real_t  min_factor = 1;
  cs_real_t  max_factor = 0;
  cs_real_t  *restrict buf = NULL, *restrict clip_factor = NULL;
  cs_real_t  *restrict denom = NULL, *restrict denum = NULL;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_cells_wghosts = mesh->n_cells_with_ghosts;
  const cs_int_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_int_t  *cell_cells_lst = mesh->cell_cells_lst;
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

      if (*itenso == 2){
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdx);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdy);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdz);
      }
      else {
        cs_halo_sync_var(halo, halo_type, dpdx);
        cs_halo_sync_var(halo, halo_type, dpdy);
        cs_halo_sync_var(halo, halo_type, dpdz);
        cs_perio_sync_var_vect(halo, halo_type, CS_PERIO_ROTA_COPY,
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
      if (*itenso == 2) {
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

    assert(sizeof(cs_int_t) == sizeof(int));

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_INT,
                  MPI_SUM, cs_glob_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(HAVE_MPI) */

  /* Output warning if necessary */

  if (*iwarnp > 1)
    bft_printf(_(" GRADIENT LIMITATION in %10d cells\n"
                 "    MINIMUM FACTOR = %14.5e; MAXIMUM FACTOR = %14.5e\n"),
               n_clip, min_factor, max_factor);

  /* Synchronize dpdx, dpdy, dpdz */

  if (halo != NULL) {

    if (*itenso == 2) {

      /* If the gradient is not treated as a "true" vector */

      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdx);
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdy);
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, dpdz);

    }
    else {

      cs_halo_sync_var(halo, halo_type, dpdx);
      cs_halo_sync_var(halo, halo_type, dpdy);
      cs_halo_sync_var(halo, halo_type, dpdz);

      cs_perio_sync_var_vect(halo,
                             halo_type,
                             CS_PERIO_ROTA_COPY,
                             dpdx, dpdy, dpdz);

    }

  }

  BFT_FREE(buf);
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
 const cs_int_t   *const idimte,      /* --> 0, 1, 2: scalar, vector, tensor  */
 const cs_int_t   *const itenso,      /* --> for rotational periodicity       */
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
  double  wt_start = 0.0, wt_stop = 0.0;
  double  cpu_start = 0.0, cpu_stop = 0.0;

  char *var_name = NULL;

  cs_int_t  *ipcvse = NULL;
  cs_int_t  *ielvse = NULL;

  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict bx, *restrict by, *restrict bz;
  cs_real_t  *restrict dpdx, *restrict dpdy, *restrict dpdz;

  cs_bool_t update_stats = true;
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

  BFT_MALLOC(var_name, 7, char);
  sprintf(var_name, "Var. %1d", *ivar);

  if (update_stats == true) {
    wt_start =bft_timer_wtime();
    cpu_start =bft_timer_cpu_time();
    gradient_info = _find_or_add_system(var_name, gradient_type);
  }

  if (*imrgra == 2 || *imrgra ==  3)
    halo_type = CS_HALO_EXTENDED;

  /* Synchronize variable */

  if (*imrgra != 0 && halo != NULL) {

    if (*itenso == 2)
      cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, pvar);
    else
      cs_halo_sync_var(halo, halo_type, pvar);

    /* TODO: check if fext* components are all up to date, in which
     *       case we need no special treatment for *itenso = 2 */

    if (*iphydp != 0) {

      if (*itenso == 2){
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fextx);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fexty);
        cs_halo_sync_component(halo, halo_type, CS_HALO_ROTATION_IGNORE, fextz);
      }
      else {
        cs_halo_sync_var(halo, halo_type, fextx);
        cs_halo_sync_var(halo, halo_type, fexty);
        cs_halo_sync_var(halo, halo_type, fextz);
        cs_perio_sync_var_vect(halo, halo_type, CS_PERIO_ROTA_COPY,
                               fextx, fexty, fextz);
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
       imrgra, inc   , iccocg, nswrgp, idimte, itenso, iphydp,
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
       inc   , iccocg, nswrgp, idimte, itenso, iphydp, imrgra,
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
       inc   , iccocg, nswrgp, idimte, itenso, iphydp, imrgra,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ipcvse, ielvse, isympa,
       volume, surfac, surfbo, surfbn, pond  ,
       dist  , distbr, dijpf , diipb ,
       fextx , fexty , fextz , xyzcen,
       cdgfac, cdgfbo, coefap, coefbp, pvar  ,
       cocgb , cocg  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

    _gradient_clipping(imrgra, &_imlini, iwarnp, itenso, &_climin,
                       pvar  , dpdx    , dpdy  , dpdz  );

    CS_PROCF (gradrc, GRADRC)
      (ncelet, ncel  , nfac  , nfabor, ncelbr,
       imrgra, inc   , iccocg, nswrgp, idimte, itenso, iphydp,
       iwarnp, nfecra, epsrgp, extrap,
       ifacel, ifabor, icelbr, ivar  ,
       volume, surfac, surfbo, pond  , xyzcen, cdgfac, cdgfbo,
       dijpf , diipb , dofij , fextx , fexty , fextz ,
       coefap, coefbp, pvar  ,
       cocib , coci  ,
       dpdx  , dpdy  , dpdz  ,
       bx    , by    , bz    );

  }

  _gradient_clipping(imrgra, imligp, iwarnp, itenso, climgp,
                     pvar  , dpdx  , dpdy  , dpdz  );

  if (update_stats == true) {

    wt_stop =bft_timer_wtime();
    cpu_stop =bft_timer_cpu_time();

    gradient_info->n_calls += 1;

    gradient_info->wt_tot += (wt_stop - wt_start);
    gradient_info->cpu_tot += (cpu_stop - cpu_start);

  }

  BFT_FREE(var_name);
  BFT_FREE(_aux_vectors);
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

  BFT_FREE(cs_glob_gradient_systems);

  cs_glob_gradient_n_systems = 0;
  cs_glob_gradient_n_max_systems = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
