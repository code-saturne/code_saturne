/*============================================================================
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

#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Fortran function/subroutine prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute coarsening array for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (autmgr, AUTMGR)
(
 const cs_int_t   *const igr,         /* --> new grid level (0 = base)        */
 const cs_int_t   *const isym,        /* --> 1: symmetric; 2 nonsymmteric     */
 const cs_int_t   *const iagmax,      /* --> max fine cells per coarse cell   */
 const cs_int_t   *const nagmax,      /* --> fine cells per coarse cell limit */
 const cs_int_t   *const ncelf,       /* --> number of cells in fine grid     */
 const cs_int_t   *const ncelfe,      /* --> n. of cells w. halo in fine grid */
 const cs_int_t   *const nfacf,       /* --> number of faces in fine grid     */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
 const cs_int_t          ifacef[],    /* --> fine grid face->cell connect.    */
 const cs_real_t         daf[],       /* --> diagonal terms of fine grid      */
 const cs_real_t         xaf[],       /* --> extradiagonal terms of fine grid */
 const cs_real_t         surfaf[],    /* --> fine grid face surface vectors   */
 const cs_real_t         volumf[],    /* --> fine grid cell volumes           */
 const cs_real_t         xyzfin[],    /* --> fine grid cell centers           */
       cs_int_t          irscel[],    /* <-- Fine -> coarse cell connectivity */
       cs_int_t          indic[],     /* --> work array of size ncelfe        */
       cs_int_t          inombr[],    /* --> work array of size ncelfe        */
       cs_int_t          irsfac[],    /* --> work array of size nfacf         */
       cs_int_t          indicf[],    /* --> work array of size nfacf         */
       cs_real_t         w1[],        /* --> work array of size ncelfe        */
       cs_real_t         w2[]         /* --> work array of size ncelfe        */
);

/*----------------------------------------------------------------------------
 * Compute coarsening grid values for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (crstgr, CRSTGR)
(
 const cs_int_t   *const iappel,      /* --> call number (0 or 1)             */
 const cs_int_t   *const isym,        /* --> 1: symmetric; 2 nonsymmteric     */
 const cs_int_t   *const igr,         /* --> new grid level (0 = base)        */
 const cs_int_t   *const ncelf,       /* --> number of cells in fine grid     */
 const cs_int_t   *const ncelg,       /* --> number of cells in coarse grid   */
 const cs_int_t   *const ncelfe,      /* --> n. of cells w. halo in fine grid */
 const cs_int_t   *const ncelge,      /* --> n. of cells w. halo coarse grid  */
 const cs_int_t   *const nfacf,       /* --> number of faces in fine grid     */
 const cs_int_t   *const nfacg,       /* --> number of faces in coarse grid   */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
 const cs_int_t          ifacef[],    /* --> fine grid face->cell connect.    */
 const cs_int_t          ifaceg[],    /* --> coarse grid face->cell connect.  */
 const cs_int_t          irscel[],    /* <-- Fine -> coarse cell connectivity */
 const cs_int_t          irsfac[],    /* <-- Fine -> coarse face connectivity */
 const cs_real_t         volumf[],    /* --> fine grid cell volumes           */
 const cs_real_t         xyzfin[],    /* --> fine grid cell centers           */
 const cs_real_t         surfaf[],    /* --> fine grid face surface vectors   */
 const cs_real_t         xaf0[],      /* --> symmetrized extradiagonal, fine  */
 const cs_real_t         xaf0ij[],    /* --> matrix coarsening term, fine     */
 const cs_real_t         daf[],       /* --> diagonal terms of fine grid      */
 const cs_real_t         xaf[],       /* --> extradiagonal terms of fine grid */
 const cs_real_t         volumg[],    /* --> coarse grid cell volumes         */
 const cs_real_t         xyzgro[],    /* --> coarse grid cell centers         */
 const cs_real_t         surfag[],    /* --> coarse grid face surface vectors */
 const cs_real_t         xag0[],      /* --> symmetrized extradiag., coarse   */
 const cs_real_t         xag0ij[],    /* --> matrix coarsening term, coarse   */
 const cs_real_t         dag[],       /* --> diagonal terms of coarse grid    */
 const cs_real_t         xag[],       /* --> extradiagonal terms, coarse grid */
       cs_real_t         rwc1[],      /* --> work array of size ncelfe        */
       cs_real_t         rwc2[],      /* --> work array of size ncelfe        */
       cs_real_t         rwc3[],      /* --> work array of size ncelfe        */
       cs_real_t         rwc4[]       /* --> work array of size ncelfe        */
);

/*----------------------------------------------------------------------------
 * Compute gradients using least squares method (standard or extended
 * neighborhood)
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (gradmc, GRADMC)
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
 const cs_int_t          ipcvse[],    /* --> cells -> extended neighborhood
                                             cells index                      */
 const cs_int_t          ielvse[],    /* --> cells -> extende neighborhood
                                             cells list                       */
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
 const cs_real_t         fextx[],     /* --> components of the exterior force */
 const cs_real_t         fexty[],     /*     generating the hydrostatic       */
 const cs_real_t         fextz[],     /*     pressure                         */
 const cs_real_t         xyzcen[],    /* --> cell centers                     */
 const cs_real_t         cdgfac[],    /* --> interior face centers of gravity */
 const cs_real_t         cdgfbo[],    /* --> boundary face centers of gravity */
 const cs_real_t         coefap[],    /* --> boundary condition term          */
 const cs_real_t         coefbp[],    /* --> boundary condition term          */
 const cs_real_t         pvar[],      /* --> gradient's base variable         */
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
);

/*----------------------------------------------------------------------------
 * Mesh renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (numvec, NUMVEC)
(
 const cs_int_t   *ncelet,        /* --> Number of cells, halo included */
 const cs_int_t   *ncel,          /* --> Number of local cells */
 const cs_int_t   *nfac,          /* --> Number of interior faces */
 const cs_int_t   *nfabor,        /* --> Number of boundary faces */
 const cs_int_t   *nsom,          /* --> Number of vertices */
 const cs_int_t   *lndfac,        /* --> Size of interior face connectivity */
 const cs_int_t   *lndfbr,        /* --> Size of boundary face connectivity */
 cs_int_t         *irveci,        /* <-- Interior face vectorization indic. */
 cs_int_t         *irvecb,        /* <-- Boundary face vectorization indic. */
 cs_int_t         *ifacel,        /* <-> Interior face->cell connectivity */
 cs_int_t         *ifabor,        /* <-> Boundary face->cell connectivity */
 cs_int_t         *ifmfbr,        /* <-> Boundary face group class number */
 cs_int_t         *ipnfac,        /* <-> Interior face->vertex index */
 cs_int_t         *nodfac,        /* <-> Interior face->vertex connectivity */
 cs_int_t         *ipnfbr,        /* <-> Boundary face->vertex index */
 cs_int_t         *nodfbr,        /* <-> Boundary face->vertex connectivity */
 cs_int_t         *inumfi,        /* <-> Interior faces renumbering array
                                         (size: nfac) */
 cs_int_t         *inumfb,        /* <-> Boundary faces renumbering array
                                         (size: nfac) */
 cs_int_t         *iworkf,        /* <-> Work array, size: max(nfac, nfabor) */
 cs_int_t         *ismbs,         /* <-> Work array, size: ncelet */
 cs_int_t         *ismbv,         /* <-> Work array, size: ncelet */
 cs_int_t         *ipnfaw,        /* <-> Work array, size: nfac+1 */
 cs_int_t         *nodfaw,        /* <-> Work array, size: lndfac */
 cs_int_t         *ipnfbw,        /* <-> Work array, size: nfabor+1 */
 cs_int_t         *nodfbw,        /* <-> Work array, size: lndfbr */
 cs_real_t        *rworkf,        /* <-> Work array, size: max(nfac, nfabor) */
 cs_real_t        *rsmbs,         /* <-> Work array, size: ncelet */
 cs_real_t        *rsmbv          /* <-> Work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * User function to compute coarsening array for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (ustmgr, USTMGR)
(
 const cs_int_t   *const iappel,      /* --> 1: initialization call
                                             2: computional call              */
 const cs_int_t   *const igr,         /* --> new grid level (0 = base)        */
 const cs_int_t   *const isym,        /* --> 1: symmetric; 2 nonsymmteric     */
 const cs_int_t   *const ncelf,       /* --> number of cells in fine grid     */
 const cs_int_t   *const ncelfe,      /* --> n. of cells w. halo in fine grid */
 const cs_int_t   *const nfacf,       /* --> number of faces in fine grid     */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
       cs_int_t   *const iusmgr,      /* <-  0: automatic method
                                             1: use this sub-routine          */
       cs_int_t   *const niw,         /* <-  size of iw for call 2            */
       cs_int_t   *const nrw,         /* <-  size of rw for call 2            */
 const cs_int_t          ifacef[],    /* --> fine grid face->cell connect.    */
 const cs_real_t         daf[],       /* --> diagonal terms of fine grid      */
 const cs_real_t         xaf[],       /* --> extradiagonal terms of fine grid */
 const cs_real_t         surfaf[],    /* --> fine grid face surface vectors   */
 const cs_real_t         volumf[],    /* --> fine grid cell volumes           */
 const cs_real_t         xyzfin[],    /* --> fine grid cell centers           */
       cs_int_t          irscel[],    /* <-- Fine -> coarse cell connectivity */
       cs_int_t          iw[],        /* --> work array of size niw (call 2)  */
       cs_real_t         rw[]         /* --> work array of size nrw (call 2)  */
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_PROTOTYPES_H__ */
