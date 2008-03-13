/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#if defined(_CS_HAVE_MPI)
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

#include "cs_ext_neighborhood.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_parall.h"
#include "cs_perio.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Encapsulation of the call to GRADMC (Fortran routine for the computation of
 * gradients by the least squares method). Add data for taking into account
 * the extended neighborhood.
 *
 * Fortran Interface :
 *
 * SUBROUTINE CGRDMC
 * *****************
 *
 *   ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
 *     INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP , IMRGRA ,
 *     IWARNP , NFECRA , EPSRGP , EXTRAP ,
 *     IFACEL , IFABOR , IA(IICELB) , IA(IISYMP) ,
 *     VOLUME , SURFAC , SURFBO , RA(ISRFBN) , RA(IPOND)   ,
 *     RA(IDIST)   , RA(IDISTB) ,
 *                   RA(IDIJPF) , RA(IDIIPB)  ,
 *     FEXTX  , FEXTY  , FEXTZ  ,
 *     XYZCEN , CDGFAC , CDGFBO , COEFAP , COEFBP , PVAR   ,
 *     RA(ICOCGB)  , RA(ICOCG)  ,
 *     DPDX   , DPDY   , DPDZ   ,
 *     DPDXA  , DPDYA  , DPDZA  )
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (cgrdmc, CGRDMC)
(
 const cs_int_t   *const ncelet,      /* --> number of extended cells         */
 const cs_int_t   *const ncel,        /* --> number of cells                  */
 const cs_int_t   *const nfac,        /* --> number of internal faces         */
 const cs_int_t   *const nfabor,      /* --> number of boundary faces         */
 const cs_int_t   *const ncelbr,      /* --> number of cells on boundary      */
 const cs_int_t   *const inc,         /* --> 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* --> 1 or 0: recompute COCG or not    */
 const cs_int_t   *const nswrgp,      /* --> >1: with reconstruction          */
 const cs_int_t   *const idimte,      /* --> 0, 1, 2: scalar, vector, tensor  */
 const cs_int_t   *const itenso,      /* --> for rotational periodicity       */
 const cs_int_t   *const iphydp,      /* --> use hydrosatatic pressure        */
 const cs_int_t   *const imrgra,      /* --> gradient computation mode        */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
 const cs_int_t   *const nfecra,      /* --> standard output unit             */
 const cs_real_t  *const epsrgp,      /* --> precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* --> extrapolate gradient at boundary */
 const cs_int_t          ifacel[],    /* --> interior face->cell connectivity */
 const cs_int_t          ifabor[],    /* --> boundary face->cell connectivity */
 const cs_int_t          icelbr[],    /* --> list of cells on boundary        */
 const cs_int_t          isympa[],    /* --> indicator for symmetry faces     */
 const cs_real_t         volume[],    /* --> cell volumes                     */
 const cs_real_t         surfac[],    /* --> surfaces of internal faces       */
 const cs_real_t         surfbo[],    /* --> surfaces of boundary faces       */
 const cs_real_t         surfbn[],    /* --> norm of surfbo                   */
 const cs_real_t         pond[],      /* --> interior faces geometric weight  */
 const cs_real_t         dist[],      /* --> interior faces I' to J' distance */
 const cs_real_t         distbr[],    /* --> boundary faces I' to J' distance */
 const cs_real_t         dijpf[],     /* --> interior faces I'J' vector       */
 const cs_real_t         diipb[],     /* --> boundary faces II' vector        */
       cs_real_t         fextx[],     /* --> components of the exterior force */
       cs_real_t         fexty[],     /*     generating the hydrostatic       */
       cs_real_t         fextz[],     /*     pressure                         */
 const cs_real_t         xyzcen[],    /* --> cell centers                     */
 const cs_real_t         cdgfac[],    /* --> interior face centers of gravity */
 const cs_real_t         cdgfbo[],    /* --> boundary face centers of gravity */
 const cs_real_t         coefap[],    /* --> boundary condition term          */
 const cs_real_t         coefbp[],    /* --> boundary condition term          */
       cs_real_t         pvar[],      /* --> gradient's base variable         */
       cs_real_t         cocgb[],     /* <-> contribution to COCG of cells
                                             on boundary's internal faces     */
       cs_real_t         cocg[],      /* <-> contribution to COCG of cells
                                             on boundary's boundary faces     */
       cs_real_t         dpdx[],      /* <-- gradient x component             */
       cs_real_t         dpdy[],      /* <-- gradient y component             */
       cs_real_t         dpdz[],      /* <-- gradient z component             */
       cs_real_t         bx[],        /* --- local work array                 */
       cs_real_t         by[],        /* --- local work array                 */
       cs_real_t         bz[]         /* --- local work array                 */
)
{
  cs_int_t  *ipcvse = NULL;
  cs_int_t  *ielvse = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  /* Sync data in case of parallelism or periodicity */

  if (mesh->n_domains > 1)
    cs_parall_sync_cells(pvar, CS_MESH_HALO_EXTENDED, 1);

  if (mesh->n_init_perio > 0)
    cs_perio_sync_var_scal(pvar,
                           CS_PERIO_ROTA_IGNORE,
                           CS_MESH_HALO_EXTENDED);

  if (*iphydp != 0) {

    if (mesh->n_domains > 1) {
      cs_parall_sync_cells(fextx, CS_MESH_HALO_EXTENDED, 1);
      cs_parall_sync_cells(fexty, CS_MESH_HALO_EXTENDED, 1);
      cs_parall_sync_cells(fextz, CS_MESH_HALO_EXTENDED, 1);
    }

    if (mesh->n_init_perio > 0)
      cs_perio_sync_var_vect(fextx, fexty, fextz,
                             CS_PERIO_ROTA_IGNORE,
                             CS_MESH_HALO_EXTENDED);

  }

  /* "cell -> cells" connectivity for the extended neighborhood */

  ipcvse = mesh->cell_cells_idx;
  ielvse = mesh->cell_cells_lst;

  /* Compute gradient */

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

/*----------------------------------------------------------------------------
 * Clip the gradient if necessary. This function deals with the standard or
 * extended neighborhood.
 *
 * Fortran Interface :
 *
 * SUBROUTINE CLMGRD
 * *****************
 *
 *    & ( IMRGRA , IMLIGP , IWARNP , CLIMGP ,
 *    &   VAR    , DPDX   , DPDY   , DPDZ   )
 *
 * parameters:
 *   imrgra         --> type of computation for the gradient
 *   imligp         --> type of clipping for the computation of the gradient
 *   idimte         --> dimension of the variable
 *                      0: scalar, 1: vector, 2: tensor
 *   itenso         --> only for periodicity when there is a rotation
 *   iwarnp         --> output level
 *   climgp         --> clipping coefficient for the computation of the gradient
 *   var            --> variable
 *   dpdx           --> X component of the pressure gradient
 *   dpdy           --> Y component of the pressure gradient
 *   dpdz           --> Z component of the pressure gradient
 *----------------------------------------------------------------------------*/

void
CS_PROCF (clmgrd, CLMGRD)(const cs_int_t   *imrgra,
                          const cs_int_t   *imligp,
                          const cs_int_t   *iwarnp,
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

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_cells_wghosts = mesh->n_cells_with_ghosts;
  const cs_int_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_int_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t  *cell_cen = cs_glob_mesh_quantities->cell_cen;

  if (*imligp < 0)
    return;

  /* Exchange for parallelism and periodicity */

  /* Caution: a processor may have no values to receive,
     but have values to send, so all ranks must call
     cs_parall_sync_cells() in parallel */

  if (mesh->n_domains > 1)
    cs_parall_sync_cells(var, CS_MESH_HALO_EXTENDED, 1);

  if (mesh->n_init_perio > 0)
    cs_perio_sync_var_scal(var,
                           CS_PERIO_ROTA_IGNORE,
                           CS_MESH_HALO_EXTENDED);

  /* Exchange for the gradients. Not useful for working array */

  if (*imligp == 1) {

    if (mesh->n_domains > 1) {
      cs_parall_sync_cells(dpdx, CS_MESH_HALO_EXTENDED, 1);
      cs_parall_sync_cells(dpdy, CS_MESH_HALO_EXTENDED, 1);
      cs_parall_sync_cells(dpdz, CS_MESH_HALO_EXTENDED, 1);
    }

    if (mesh->n_init_perio > 0)
      cs_perio_sync_var_vect(dpdx, dpdy, dpdz,
                             CS_PERIO_ROTA_IGNORE,
                             CS_MESH_HALO_EXTENDED);

  } /* End if imligp == 1 */

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

    if (cell_cells_idx != NULL) {

      if ( (*imrgra == 2) || (*imrgra == 3) ) {

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

      } /* End if imrgra == 2 or 3 */

    } /* End if there is an extended neighborhood */

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

    if (cell_cells_idx != NULL) {

      if ( (*imrgra == 2) || (*imrgra == 3) ) {

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

      } /* End if imrgra == 2 or 3 */

    } /* End if there is an extended neighborhood */

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

    if (mesh->n_domains > 1) {

      if (*imrgra == 2 || *imrgra ==  3) {
        cs_parall_sync_cells(denom, CS_MESH_HALO_EXTENDED, 1);
        cs_parall_sync_cells(denum, CS_MESH_HALO_EXTENDED, 1);
      }
      else {
        cs_parall_sync_cells(denom, CS_MESH_HALO_STANDARD, 1);
        cs_parall_sync_cells(denum, CS_MESH_HALO_STANDARD, 1);
      }

    }

    if (mesh->n_init_perio > 0) {

      if (*imrgra == 2 || *imrgra ==  3) {

        cs_perio_sync_var_scal(denom,
                               CS_PERIO_ROTA_IGNORE,
                               CS_MESH_HALO_EXTENDED);
        cs_perio_sync_var_scal(denum,
                               CS_PERIO_ROTA_IGNORE,
                               CS_MESH_HALO_EXTENDED);

      }
      else {

        cs_perio_sync_var_scal(denom,
                               CS_PERIO_ROTA_IGNORE,
                               CS_MESH_HALO_STANDARD);
        cs_perio_sync_var_scal(denum,
                               CS_PERIO_ROTA_IGNORE,
                               CS_MESH_HALO_STANDARD);

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

    if (cell_cells_idx != NULL) {

      if ( (*imrgra == 2) || (*imrgra == 3) ) {

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

      } /* End if imrgra == 2 or 3 */

    } /* End if there is an extended neighborhood */

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

#if defined(_CS_HAVE_MPI)

  if (mesh->n_domains > 1) {

    assert(sizeof(cs_real_t) == sizeof(double));

    /* Global Max */

    MPI_Allreduce(&max_factor, &global_max_factor, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_base_mpi_comm);

    max_factor = global_max_factor;

    /* Global min */

    MPI_Allreduce(&min_factor, &global_min_factor, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_base_mpi_comm);

    min_factor = global_min_factor;

    /* Sum number of clippings */

    assert(sizeof(cs_int_t) == sizeof(int));

    MPI_Allreduce(&n_clip, &n_g_clip, 1, CS_MPI_INT,
                  MPI_SUM, cs_glob_base_mpi_comm);

    n_clip = n_g_clip;

  } /* If n_domains > 1 */

#endif /* defined(_CS_HAVE_MPI) */

  /* Output warning if necessary */

  if (*iwarnp > 1)
    bft_printf(_(" LIMITATION DE GRADIENT en %10d CELLULES\n"
                 "    FACTEUR MINIMUM = %14.5e ; FACTEUR MAXIMUM = %14.5e\n"),
               n_clip, min_factor, max_factor);

  /* Sync dpdx, dpdy, dpdz in case of parallelism or periodicity */

  if (mesh->n_domains > 1) {

    cs_parall_sync_cells(dpdx, CS_MESH_HALO_STANDARD, 1);
    cs_parall_sync_cells(dpdy, CS_MESH_HALO_STANDARD, 1);
    cs_parall_sync_cells(dpdz, CS_MESH_HALO_STANDARD, 1);

  }

  if (mesh->n_init_perio > 0)
    cs_perio_sync_var_vect(dpdx, dpdy, dpdz,
                           CS_PERIO_ROTA_IGNORE,
                           CS_MESH_HALO_STANDARD);

  BFT_FREE(buf);

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

