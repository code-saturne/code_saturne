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

BEGIN_C_DECLS

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
 const cs_int_t   *igr,         /* <-- new grid level (0 = base)              */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric           */
 const cs_int_t   *iagmax,      /* <-- max fine cells per coarse cell         */
 const cs_int_t   *nagmax,      /* <-- fine cells per coarse cell limit       */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid           */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid           */
 const cs_int_t   *iwarnp,      /* <-- verbosity level                        */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity      */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid            */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid       */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors         */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes                 */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers                 */
       cs_int_t    irscel[],    /* --> Fine -> coarse cell connectivity       */
       cs_int_t    indic[],     /* --- work array of size ncelfe              */
       cs_int_t    inombr[],    /* --- work array of size ncelfe              */
       cs_int_t    irsfac[],    /* --- work array of size nfacf               */
       cs_int_t    indicf[],    /* --- work array of size nfacf               */
       cs_real_t   w1[],        /* --- work array of size ncelfe              */
       cs_real_t   w2[]         /* --- work array of size ncelfe              */
);

/*----------------------------------------------------------------------------
 * Main Fortran subroutine
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (caltri, CALTRI)
(
 const cs_int_t   *iverif,   /* <-- activate elementary tests                 */
 const cs_int_t   *nideve,   /* <-- size of IDEVEL integer array              */
 const cs_int_t   *nrdeve,   /* <-- size of rdevel floating-point array       */
 const cs_int_t   *nituse,   /* <-- size of ITUSER integer array              */
 const cs_int_t   *nrtuse,   /* <-- size of RTUSER floating-point array       */
 const cs_int_t   *ifacel,   /* <-- interior faces -> cells connectivity      */
 const cs_int_t   *ifabor,   /* <-- boundary faces -> cells connectivity      */
 const cs_int_t   *ifmfbr,   /* <-- boundary face family (group class) number */
 const cs_int_t   *ifmcel,   /* <-- cell family (group class) number          */
 const cs_int_t   *iprfml,   /* <-- family (group class) properties           */
 const cs_int_t   *ipnfac,   /* <-- interior faces -> vertices connect. index */
 const cs_int_t   *nodfac,   /* <-- interior faces -> vertices connectivity   */
 const cs_int_t   *ipnfbr,   /* <-- boundary faces -> vertices connect. index */
 const cs_int_t   *nodfbr,   /* <-- boundary faces -> vertices connectivity   */
 cs_int_t         *idevel,   /* --> IDEVEL integer array                      */
 cs_int_t         *ituser,   /* --> ITUSER integer array                      */
 cs_int_t         *ia,       /* --> IA integer array                          */
 cs_real_t        *xyzcen,   /* <-> points associated with cell centers       */
 cs_real_t        *surfac,   /* <-> interior face surface vectors             */
 cs_real_t        *surfbo,   /* <-> boundary face surface vectors             */
 cs_real_t        *cdgfac,   /* <-> interior face centers                     */
 cs_real_t        *cdgfbr,   /* <-> boundary face vectors                     */
 cs_real_t        *xyznod,   /* <-> vertex coordinates (optional)             */
 cs_real_t        *volume,   /* <-> cell volumes                              */
 cs_real_t        *rdevel,   /* --> RDEVEL floating-point array               */
 cs_real_t        *rtuser,   /* --> RTUSER floating-point array               */
 cs_real_t        *ra        /* --> RA floating-point array                   */
);

/*----------------------------------------------------------------------------
 * Compute coarsening grid values for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (crstgr, CRSTGR)
(
 const cs_int_t   *iappel,      /* <-- call number (0 or 1)                   */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric           */
 const cs_int_t   *igr,         /* <-- new grid level (0 = base)              */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid           */
 const cs_int_t   *ncelg,       /* <-- number of cells in coarse grid         */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *ncelge,      /* <-- number of cells with halo coarse grid  */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid           */
 const cs_int_t   *nfacg,       /* <-- number of faces in coarse grid         */
 const cs_int_t   *iwarnp,      /* <-- verbosity level                        */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity      */
 const cs_int_t    ifaceg[],    /* <-- coarse grid face->cell connectivity    */
 const cs_int_t    irscel[],    /* <-- Fine -> coarse cell connectivity       */
 const cs_int_t    irsfac[],    /* <-- Fine -> coarse face connectivity       */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes                 */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers                 */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors         */
 const cs_real_t   xaf0[],      /* <-- symmetrized extradiagonal, fine        */
 const cs_real_t   xaf0ij[],    /* <-- matrix coarsening term, fine           */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid            */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid       */
 const cs_real_t   volumg[],    /* <-- coarse grid cell volumes               */
 const cs_real_t   xyzgro[],    /* <-- coarse grid cell centers               */
 const cs_real_t   surfag[],    /* <-- coarse grid face surface vectors       */
 cs_real_t         xag0[],      /* --> symmetrized extradiagonal, coarse      */
 cs_real_t         xag0ij[],    /* --> matrix coarsening term, coarse         */
 cs_real_t         dag[],       /* --> diagonal terms of coarse grid          */
 cs_real_t         xag[],       /* --> extradiagonal terms, coarse grid       */
 cs_real_t         rwc1[],      /* --- work array of size ncelfe              */
 cs_real_t         rwc2[],      /* --- work array of size ncelfe              */
 cs_real_t         rwc3[],      /* --- work array of size ncelfe              */
 cs_real_t         rwc4[]       /* --- work array of size ncelfe              */
);

/*----------------------------------------------------------------------------
 * Close log (listing) handled by Fortran: (CLose LIsting)
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csclli, CSCLLI)
(
 void
);

/*----------------------------------------------------------------------------
 * Flush standard output.
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csflsh, CSFLSH)
(
 void
);

/*----------------------------------------------------------------------------
 * Initialize Fortran log (listing) files
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csinit, CSINIT)
(
 const cs_int_t  *ifoenv,  /* <-- 0: SolCom mesh; 1: Preprocesor mesh         */
 const cs_int_t  *iparal,  /* <-- MPI Rank in parallel, -1 otherwise          */
 const cs_int_t  *nparal,  /* <-- Number of MPI processes, or 1               */
 const cs_int_t  *ilisr0,  /* <-- Output of main log (listing (rank 0):       */
                           /*     0: non redirected; 1: to 'listing' file     */
 const cs_int_t  *ilisrp   /* <-- Output of logs for ranks > 0:               */
                           /*     0: non redirected; 1: to 'listing_n*' files */
                           /*     2: to '/dev/null' (suppressed)              */
);

/*----------------------------------------------------------------------------
 * Print a message to standard output.
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csprnt, CSPRNT)
(
  char       *cs_buf_print,
  cs_int_t   *msgsize
);

/*----------------------------------------------------------------------------
 * Find the nearest cell's center from a node
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (findpt, FINDPT)
(
 const cs_int_t   *ncelet,      /* <-- number of extended cells               */
 const cs_int_t   *ncel,        /* <-- number of cells                        */
 const cs_real_t  *xyzcen,      /* <-- cell centers                           */
 const cs_real_t  *xx,          /* <-- node coordinate X                      */
 const cs_real_t  *yy,          /* <-- node coordinate Y                      */
 const cs_real_t  *zz,          /* <-- node coordinate Z                      */
       cs_int_t   *node,        /* --> node we are looking for, zero if error */
       cs_int_t   *ndrang       /* --> rank of associated process             */
);

/*----------------------------------------------------------------------------
 * Compute gradients using least squares method (standard or extended
 * neighborhood)
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (gradmc, GRADMC)
(
 const cs_int_t   *ncelet,      /* <-- number of extended cells               */
 const cs_int_t   *ncel,        /* <-- number of cells                        */
 const cs_int_t   *nfac,        /* <-- number of internal faces               */
 const cs_int_t   *nfabor,      /* <-- number of boundary faces               */
 const cs_int_t   *ncelbr,      /* <-- number of cells on boundary            */
 const cs_int_t   *inc,         /* <-- 0 or 1: increment or not               */
 const cs_int_t   *iccocg,      /* <-- 1 or 0: recompute COCG or not          */
 const cs_int_t   *nswrgp,      /* <-- >1: with reconstruction                */
 const cs_int_t   *idimte,      /* <-- 0, 1, 2: scalar, vector, tensor        */
 const cs_int_t   *itenso,      /* <-- for rotational periodicity             */
 const cs_int_t   *iphydp,      /* <-- use hydrosatatic pressure              */
 const cs_int_t   *imrgra,      /* <-- gradient computation mode              */
 const cs_int_t   *iwarnp,      /* <-- verbosity level                        */
 const cs_int_t   *nfecra,      /* <-- standard output unit                   */
 const cs_real_t  *epsrgp,      /* <-- precision for iterative gradient calc. */
 const cs_real_t  *extrap,      /* <-- extrapolate gradient at boundary       */
 const cs_int_t    ifacel[],    /* <-- interior face->cell connectivity       */
 const cs_int_t    ifabor[],    /* <-- boundary face->cell connectivity       */
 const cs_int_t    icelbr[],    /* <-- list of cells on boundary              */
 const cs_int_t    ipcvse[],    /* <-- cells -> ext. neighborhood cells index */
 const cs_int_t    ielvse[],    /* <-- cells -> ext. neighborhood cells list  */
 const cs_int_t    isympa[],    /* <-- indicator for symmetry faces           */
 const cs_real_t   volume[],    /* <-- cell volumes                           */
 const cs_real_t   surfac[],    /* <-- surfaces of internal faces             */
 const cs_real_t   surfbo[],    /* <-- surfaces of boundary faces             */
 const cs_real_t   surfbn[],    /* <-- norm of surfbo                         */
 const cs_real_t   pond[],      /* <-- interior faces geometric weight        */
 const cs_real_t   dist[],      /* <-- interior faces I' to J' distance       */
 const cs_real_t   distbr[],    /* <-- boundary faces I' to J' distance       */
 const cs_real_t   dijpf[],     /* <-- interior faces I'J' vector             */
 const cs_real_t   diipb[],     /* <-- boundary faces II' vector              */
 const cs_real_t   fextx[],     /* <-- components of the exterior force       */
 const cs_real_t   fexty[],     /*     generating the hydrostatic pressure    */
 const cs_real_t   fextz[],     /*                                            */
 const cs_real_t   xyzcen[],    /* <-- cell centers                           */
 const cs_real_t   cdgfac[],    /* <-- interior face centers of gravity       */
 const cs_real_t   cdgfbo[],    /* <-- boundary face centers of gravity       */
 const cs_real_t   coefap[],    /* <-- boundary condition term                */
 const cs_real_t   coefbp[],    /* <-- boundary condition term                */
 const cs_real_t   pvar[],      /* <-- gradient's base variable               */
       cs_real_t   cocgb[],     /* <-> contribution to COCG of cells on
                                       on boundary's internal faces           */
       cs_real_t   cocg[],      /* <-> contribution to COCG of cells on
                                       on boundary's boundary faces           */
       cs_real_t   dpdx[],      /* --> gradient x component                   */
       cs_real_t   dpdy[],      /* --> gradient y component                   */
       cs_real_t   dpdz[],      /* --> gradient z component                   */
       cs_real_t   bx[],        /* --- local work array                       */
       cs_real_t   by[],        /* --- local work array                       */
       cs_real_t   bz[]         /* --- local work array                       */
);

/*----------------------------------------------------------------------------
 * Main Fortran options initialization
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (initi1, INITI1)
(
 const cs_int_t  *iverif          /* <-- Activate elementary tests            */
);

/*----------------------------------------------------------------------------
 * Read geometric entities in "SolCom" format
 *----------------------------------------------------------------------------*/

void CS_PROCF (letgeo, LETGEO)
(
 const cs_int_t  *ndim,      /* <-- space dimension                           */
 const cs_int_t  *ncelet,    /* <-- number of extended cells                  */
 const cs_int_t  *ncel,      /* <-- number of cells                           */
 const cs_int_t  *nfac,      /* <-- number of interior faces                  */
 const cs_int_t  *nfabor,    /* <-- number of boundary faces                  */
 const cs_int_t  *nfml,      /* <-- number of familes                         */
 const cs_int_t  *nprfml,    /* <-- number of family properties               */
 const cs_int_t  *nnod,      /* <-- number of vertices                        */
 const cs_int_t  *lndfac,    /* <-- size of interior faces -> vertices array  */
 const cs_int_t  *lndfbr,    /* <-- size of boundary faces -> vertices array  */
 const cs_int_t  *ntetra,    /* <-- number of tetrahedra                      */
 const cs_int_t  *npyram,    /* <-- number of pyramids                        */
 const cs_int_t  *nprism,    /* <-- number of prisms                          */
 const cs_int_t  *nhexae,    /* <-- number of hexahedra                       */
       cs_int_t  *inodal,    /* --> indicates if we should read the nodal     */
                             /*     connectivity for post-processing          */
       cs_int_t   ifacel[],  /* --> interior faces -> cells connectivity      */
       cs_int_t   ifabor[],  /* --> boundary faces -> cells connectivity      */
       cs_int_t   ifmfbr[],  /* --> boundary face families                    */
       cs_int_t   ifmcel[],  /* --> cell families                             */
       cs_int_t   iprfml[],  /* --> list of family properties                 */
       cs_int_t   icotet[],  /* --> nodal connectivity for tetrahedra         */
       cs_int_t   icopyr[],  /* --> nodal connectivity for pyramids           */
       cs_int_t   icopri[],  /* --> nodal conncetivity for prisms             */
       cs_int_t   icohex[],  /* --> nodal connectivity for hexahedra          */
       cs_int_t   ipnfac[],  /* --> interior faces -> vertices index          */
       cs_int_t   nodfac[],  /* --> interior faces -> vertices connectivity   */
       cs_int_t   ipnfbr[],  /* --> boundary faces -> vertices index          */
       cs_int_t   nodfbr[],  /* --> boundary faces -> vertices connectivity   */
       cs_real_t  xyzcen[],  /* --> cell centers                              */
       cs_real_t  surfac[],  /* --> interior face surface vectors             */
       cs_real_t  surfbo[],  /* --> boundary face surface vectors             */
       cs_real_t  cdgfac[],  /* --> interior face centers                     */
       cs_real_t  cdgfbo[],  /* --> boundary face centers                     */
       cs_real_t  xyznod[]   /* --> vertex coordinates                        */
);

/*----------------------------------------------------------------------------
 * Update mesh dimensions in Fortran commons
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (majgeo, MAJGEO)
(
 const cs_int_t   *ncel,    /* <-- number of cells                            */
 const cs_int_t   *ncelet,  /* <-- number of halo cells                       */
 const cs_int_t   *nfac,    /* <-- number of internal faces                   */
 const cs_int_t   *nfabor,  /* <-- number of border faces                     */
 const cs_int_t   *nsom,    /* <-- number of vertices                         */
 const cs_int_t   *lndfac,  /* <-- internal face -> vertices array size       */
 const cs_int_t   *lndfbr,  /* <-- boundary face -> vertices array size       */
 const cs_int_t   *ncelgb,  /* <-- global number of cells                     */
 const cs_int_t   *nfacgb,  /* <-- global number of internal faces            */
 const cs_int_t   *nfbrgb,  /* <-- global number of boundary faces            */
 const cs_int_t   *nsomgb   /* <-- global number of vertices                  */
);

/*----------------------------------------------------------------------------
 * Initialize Fortran working array sizes
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (memini, MEMINI)
(
 cs_int_t  *iasize,   /* --- size of IA integer array                         */
 cs_int_t  *rasize,   /* --- size of RA integer array                         */
 cs_int_t  *nideve,   /* --- size of IDEVEL integer array                     */
 cs_int_t  *nrdeve,   /* --- size of rdevel floating-point array              */
 cs_int_t  *nituse,   /* --- size of ITUSER integer array                     */
 cs_int_t  *nrtuse    /* --- size of RTUSER floating-point array              */
);

/*----------------------------------------------------------------------------
 * Mesh renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (numvec, NUMVEC)
(
 const cs_int_t   *ncelet,        /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,          /* <-- Number of local cells */
 const cs_int_t   *nfac,          /* <-- Number of interior faces */
 const cs_int_t   *nfabor,        /* <-- Number of boundary faces */
 cs_int_t         *irveci,        /* <-> Interior face vectorization indic. */
 cs_int_t         *irvecb,        /* <-> Boundary face vectorization indic. */
 cs_int_t         *ifacel,        /* <-> Interior face->cell connectivity */
 cs_int_t         *ifabor,        /* <-> Boundary face->cell connectivity */
 cs_int_t         *inumfi,        /* <-> Interior faces renumbering array
                                         (size: nfac) */
 cs_int_t         *inumfb,        /* <-> Boundary faces renumbering array
                                         (size: nfac) */
 cs_int_t         *iworkf,        /* --> Work array, size: max(nfac, nfabor) */
 cs_int_t         *ismbs          /* --> Work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * Test renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (tstvec, TSTVEC)
(
 const cs_int_t   *ncelet,        /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,          /* <-- Number of local cells */
 const cs_int_t   *nfac,          /* <-- Number of interior faces */
 const cs_int_t   *nfabor,        /* <-- Number of boundary faces */
 cs_int_t         *ifacel,        /* <-> Interior face->cell connectivity */
 cs_int_t         *ifabor,        /* <-> Boundary face->cell connectivity */
 cs_int_t         *iworkf,        /* --> Work array, size: max(nfac, nfabor) */
 cs_int_t         *ismbs,         /* --> Work array, size: ncelet */
 cs_int_t         *ismbv,         /* --> Work array, size: ncelet */
 cs_real_t        *rworkf,        /* --> Work array, size: max(nfac, nfabor) */
 cs_real_t        *rsmbs,         /* --> Work array, size: ncelet */
 cs_real_t        *rsmbv          /* --> Work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * User subroutine for geometry modification
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (usmodg, USMODG)
(
 const cs_int_t  *ndim,      /* <-- spatial dimension                         */
 const cs_int_t  *ncelet,    /* <-- number of extended (real + ghost) cells   */
 const cs_int_t  *ncel,      /* <-- number of cells                           */
 const cs_int_t  *nfac,      /* <-- number of internal faces                  */
 const cs_int_t  *nfabor,    /* <-- number of boundary faces                  */
 const cs_int_t  *nfml,      /* <-- number of families (group classes)        */
 const cs_int_t  *nprfml,    /* <-- number of family (group class) properties */
 const cs_int_t  *nnod,      /* <-- number of vertices                        */
 const cs_int_t  *lndfac,    /* <-- size of nodfac                            */
 const cs_int_t  *lndfbr,    /* <-- size of nodfbr                            */
 const cs_int_t   ifacel[],  /* <-- interior faces / cells connectivity       */
 const cs_int_t   ifabor[],  /* <-- boundary faces / cell connectivity        */
 const cs_int_t   ifmfbr[],  /* <-- boundary face families                    */
 const cs_int_t   ifmcel[],  /* <-- cell families                             */
 const cs_int_t   iprfml[],  /* <-- list of family properties                 */
 const cs_int_t   ipnfac[],  /* <-- interior faces -> vertices connect. index */
 const cs_int_t   nodfac[],  /* <-- interior faces -> vertices connectivity   */
 const cs_int_t   ipnfbr[],  /* <-- boundary faces -> vertices connect. index */
 const cs_int_t   nodfbr[],  /* <-- boundary faces -> vertices connectivity   */
       cs_real_t  xyznod[]   /* --> vertex coordinates                        */
);

/*----------------------------------------------------------------------------
 * User function to compute coarsening array for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (ustmgr, USTMGR)
(
 const cs_int_t   *iappel,      /* <-- 1: initialization call
                                       2: computional call                    */
 const cs_int_t   *igr,         /* <-- new grid level (0 = base)              */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric           */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid           */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid           */
 const cs_int_t   *iwarnp,      /* <-- verbosity level                        */
       cs_int_t   *iusmgr,      /* --> 0: automatic method
                                       1: use this sub-routine                */
       cs_int_t   *niw,         /* --> size of iw for call 2                  */
       cs_int_t   *nrw,         /* --> size of rw for call 2                  */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity      */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid            */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid       */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors         */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes                 */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers                 */
       cs_int_t    irscel[],    /* --> Fine -> coarse cell connectivity       */
       cs_int_t    iw[],        /* --> work array of size niw (call 2)        */
       cs_real_t   rw[]         /* --> work array of size nrw (call 2)        */
);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
