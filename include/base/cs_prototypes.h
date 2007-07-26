/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
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

void CS_PROCF (numvec, NUMVEC)
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
 * Matrix-vector product
 *----------------------------------------------------------------------------*/

void CS_PROCF (promav, PROMAV)
(
 const cs_int_t   *ncelet,        /* --> Number of cells, halo included */
 const cs_int_t   *ncel,          /* --> Number of local cells */
 const cs_int_t   *nfac,          /* --> Number of faces */
 const cs_int_t   *isym,          /* --> Symmetry indicator:
                                         1: symmetric; 2: not symmetric */
 const cs_int_t   *iinvpe,        /* --> Indicator to cancel increments
                                         in rotational periodicty (2) or
                                         to exchange them as scalars (1) */
 const cs_int_t   *ifacel,        /* --> Face -> cell connectivity  */
 const cs_real_t  *da,            /* --> Matrix diagonal */
 const cs_real_t  *xa,            /* --> Matrix extra-diagonal terms */
 const cs_real_t  *vx,            /* --> Vector to be multiplied */
 const cs_real_t  *vy             /* <-- Resulting vector */
);

/*----------------------------------------------------------------------------
 * Matrix-vector product (non-symmetric part only)
 *----------------------------------------------------------------------------*/

void CS_PROCF (proxav, PROXAV)
(
 const cs_int_t   *ncelet,        /* --> Number of cells, halo included */
 const cs_int_t   *ncel,          /* --> Number of local cells */
 const cs_int_t   *nfac,          /* --> Number of faces */
 const cs_int_t   *isym,          /* --> Symmetry indicator:
                                         1: symmetric; 2: not symmetric */
 const cs_int_t   *iinvpe,        /* --> Indicator to cancel increments
                                         in rotational periodicty (2) or
                                         to exchange them as scalars (1) */
 const cs_int_t   *ifacel,        /* --> Face -> cell connectivity  */
 const cs_real_t  *xa,            /* --> Matrix extra-diagonal terms */
 const cs_real_t  *vx,            /* --> Vector to be multiplied */
 const cs_real_t  *vy             /* <-- Resulting vector */
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_PROTOTYPES_H__ */
