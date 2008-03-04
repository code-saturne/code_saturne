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

#ifndef __CS_GRADIENT_H__
#define __CS_GRADIENT_H__

/*============================================================================
 * Gradient reconstruction.
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

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
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
);

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
                          cs_real_t         dpdz[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_GRADIENT__ */
