/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
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
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

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
 const cs_int_t   *igr,         /* <-- new grid level (0 = base) */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric */
 const cs_int_t   *iagmax,      /* <-- max fine cells per coarse cell */
 const cs_int_t   *nagmax,      /* <-- fine cells per coarse cell limit */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid */
 const cs_int_t   *iwarnp,      /* <-- verbosity level */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers */
       cs_int_t    irscel[],    /* --> Fine -> coarse cell connectivity */
       cs_int_t    indic[],     /* --- work array of size ncelfe */
       cs_int_t    inombr[],    /* --- work array of size ncelfe */
       cs_int_t    irsfac[],    /* --- work array of size nfacf */
       cs_int_t    indicf[],    /* --- work array of size nfacf */
       cs_real_t   w1[],        /* --- work array of size ncelfe */
       cs_real_t   w2[]         /* --- work array of size ncelfe */
);

/*----------------------------------------------------------------------------
 * Main Fortran subroutine
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (caltri, CALTRI)
(
 const cs_int_t   *iverif   /* <-- activate elementary tests */
);

/*----------------------------------------------------------------------------
 * Compute coarsening grid values for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (crstgr, CRSTGR)
(
 const cs_int_t   *iappel,      /* <-- call number (0 or 1) */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric */
 const cs_int_t   *igr,         /* <-- new grid level (0 = base) */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid */
 const cs_int_t   *ncelg,       /* <-- number of cells in coarse grid */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *ncelge,      /* <-- number of cells with halo coarse grid */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid */
 const cs_int_t   *nfacg,       /* <-- number of faces in coarse grid */
 const cs_int_t   *iwarnp,      /* <-- verbosity level */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity */
 const cs_int_t    ifaceg[],    /* <-- coarse grid face->cell connectivity */
 const cs_int_t    irscel[],    /* <-- Fine -> coarse cell connectivity */
 const cs_int_t    irsfac[],    /* <-- Fine -> coarse face connectivity */
 const cs_real_t  *rlxp1,       /* <-- P0/P1 relaxation parameter */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors */
 const cs_real_t   xaf0[],      /* <-- symmetrized extradiagonal, fine */
 const cs_real_t   xaf0ij[],    /* <-- matrix coarsening term, fine */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid */
 const cs_real_t   volumg[],    /* <-- coarse grid cell volumes */
 const cs_real_t   xyzgro[],    /* <-- coarse grid cell centers */
 const cs_real_t   surfag[],    /* <-- coarse grid face surface vectors */
 cs_real_t         xag0[],      /* --> symmetrized extradiagonal, coarse */
 cs_real_t         xag0ij[],    /* --> matrix coarsening term, coarse */
 cs_real_t         dag[],       /* --> diagonal terms of coarse grid */
 cs_real_t         xag[],       /* --> extradiagonal terms, coarse grid */
 cs_real_t         rwc1[],      /* --- work array of size ncelfe */
 cs_real_t         rwc2[],      /* --- work array of size ncelfe */
 cs_real_t         rwc3[],      /* --- work array of size ncelfe */
 cs_real_t         rwc4[]       /* --- work array of size ncelfe */
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
 * Initialize Fortran base common bloc values
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csinit, CSINIT)
(
 const cs_int_t  *irgpar,  /* <-- MPI Rank in parallel, -1 otherwise */
 const cs_int_t  *nrgpar,  /* <-- Number of MPI processes, or 1 */
 const cs_int_t  *nthpar   /* <-- Number of threads */
);

/*----------------------------------------------------------------------------
 * Initialize Fortran log (listing) files
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csopli, CSOPLI)
(
 const cs_int_t  *irkpar,  /* <-- MPI Rank in parallel, -1 otherwise */
 const cs_int_t  *nrkpar,  /* <-- Number of MPI processes, or 1 */
 const cs_int_t  *ilogr0,  /* <-- Output of main log (listing (rank 0): */
                           /*     0: non redirected; 1: to 'listing' file */
 const cs_int_t  *ilogrp   /* <-- Output of logs for ranks > 0: */
                           /*     0: non redirected; 1: to 'listing_n*' files */
                           /*     2: to '/dev/null' (suppressed) */
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
 * Developper function for output of variables on a post-processing mesh
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (dvvpst, DVVPST)
(
 const cs_int_t  *nummai,    /* <-- number or post-processing mesh */
 const cs_int_t  *numtyp,    /* <-- number or post-processing type
                              *     (-1 as volume, -2 as boundary, or nummai) */
 const cs_int_t  *nvar,      /* <-- number of variables */
 const cs_int_t  *nscal,     /* <-- number of scalars */
 const cs_int_t  *nvlsta,    /* <-- number of statistical variables (lagr) */
 const cs_int_t  *nvisbr,    /* <-- number of boundary stat. variables (lagr) */
 const cs_int_t  *ncelps,    /* <-- number of post-processed cells */
 const cs_int_t  *nfacps,    /* <-- number of post processed interior faces */
 const cs_int_t  *nfbrps,    /* <-- number of post processed boundary faces */
 const cs_int_t   itypps[3], /* <-- flag (0 or 1) for presence of cells, */
 const cs_int_t   lstcel[],  /* <-- list of post-processed cells */
 const cs_int_t   lstfac[],  /* <-- list of post-processed interior faces */
 const cs_int_t   lstfbr[],  /* <-- list of post-processed boundary faces */
 const cs_real_t  dt[],      /* <-- local time step */
 const cs_real_t  rtpa[],    /* <-- cell variables at previous time step */
 const cs_real_t  rtp[],     /* <-- cell variables */
 const cs_real_t  propce[],  /* <-- cell physical properties */
 const cs_real_t  propfa[],  /* <-- interior face physical properties */
 const cs_real_t  propfb[],  /* <-- boundary face physical properties */
 const cs_real_t  coefa[],   /* <-- boundary conditions array */
 const cs_real_t  coefb[],   /* <-- boundary conditions array */
 const cs_real_t  statce[],  /* <-- cell statistics (Lagrangian) */
 const cs_real_t  stativ[],  /* <-- cell variance statistics (Lagrangian) */
 const cs_real_t  statfb[],  /* <-- boundary face statistics (Lagrangian) */
 cs_real_t        tracel[],  /* --- work array for output cells */
 cs_real_t        trafac[],  /* --- work array for output interior faces */
 cs_real_t        trafbr[],  /* --- work array for output boundary faces */
 cs_real_t        ra[]       /* <-- RA floating-point array */
);

/*----------------------------------------------------------------------------
 * Find the nearest cell's center from a node
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (findpt, FINDPT)
(
 const cs_int_t   *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t   *ncel,     /* <-- number of cells */
 const cs_real_t  *xyzcen,   /* <-- cell centers */
 const cs_real_t  *xx,       /* <-- node coordinate X */
 const cs_real_t  *yy,       /* <-- node coordinate Y */
 const cs_real_t  *zz,       /* <-- node coordinate Z */
       cs_int_t   *node,     /* --> node we are looking for, zero if error */
       cs_int_t   *ndrang    /* --> rank of associated process */
);

/*----------------------------------------------------------------------------
 * Compute gradients using least squares method (standard or extended
 * neighborhood)
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (gradmc, GRADMC)
(
 const cs_int_t   *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t   *ncel,     /* <-- number of cells */
 const cs_int_t   *nfac,     /* <-- number of interior faces */
 const cs_int_t   *nfabor,   /* <-- number of boundary faces */
 const cs_int_t   *ncelbr,   /* <-- number of cells on boundary */
 const cs_int_t   *inc,      /* <-- 0 or 1: increment or not */
 const cs_int_t   *iccocg,   /* <-- 1 or 0: recompute COCG or not */
 const cs_int_t   *nswrgp,   /* <-- >1: with reconstruction */
 const cs_int_t   *idimte,   /* <-- 0, 1, 2: scalar, vector, tensor */
 const cs_int_t   *itenso,   /* <-- for rotational periodicity */
 const cs_int_t   *iphydp,   /* <-- use hydrosatatic pressure */
 const cs_int_t   *imrgra,   /* <-- gradient computation mode */
 const cs_int_t   *iwarnp,   /* <-- verbosity level */
 const cs_int_t   *nfecra,   /* <-- standard output unit */
 const cs_real_t  *epsrgp,   /* <-- precision for iterative gradient calc. */
 const cs_real_t  *extrap,   /* <-- extrapolate gradient at boundary */
 const cs_int_t    ifacel[], /* <-- interior face->cell connectivity */
 const cs_int_t    ifabor[], /* <-- boundary face->cell connectivity */
 const cs_int_t    icelbr[], /* <-- list of cells on boundary */
 const cs_int_t    ipcvse[], /* <-- cells -> ext. neighborhood cells index */
 const cs_int_t    ielvse[], /* <-- cells -> ext. neighborhood cells list */
 const cs_int_t    isympa[], /* <-- indicator for symmetry faces */
 const cs_real_t   volume[], /* <-- cell volumes */
 const cs_real_t   surfac[], /* <-- surfaces of interior faces */
 const cs_real_t   surfbo[], /* <-- surfaces of boundary faces */
 const cs_real_t   surfbn[], /* <-- norm of surfbo */
 const cs_real_t   pond[],   /* <-- interior faces geometric weight */
 const cs_real_t   dist[],   /* <-- interior faces I' to J' distance */
 const cs_real_t   distbr[], /* <-- boundary faces I' to J' distance */
 const cs_real_t   dijpf[],  /* <-- interior faces I'J' vector */
 const cs_real_t   diipb[],  /* <-- boundary faces II' vector */
 const cs_real_t   fextx[],  /* <-- components of the exterior force */
 const cs_real_t   fexty[],  /*     generating the hydrostatic pressure */
 const cs_real_t   fextz[],
 const cs_real_t   xyzcen[], /* <-- cell centers */
 const cs_real_t   cdgfac[], /* <-- interior face centers of gravity */
 const cs_real_t   cdgfbo[], /* <-- boundary face centers of gravity */
 const cs_real_t   coefap[], /* <-- boundary condition term */
 const cs_real_t   coefbp[], /* <-- boundary condition term */
 const cs_real_t   pvar[],   /* <-- gradient's base variable */
       cs_real_t   cocgb[],  /* <-> contribution to COCG of cells on
                                    on boundary's interior faces */
       cs_real_t   cocg[],   /* <-> contribution to COCG of cells on
                                    on boundary's boundary faces */
       cs_real_t   dpdx[],   /* --> gradient x component */
       cs_real_t   dpdy[],   /* --> gradient y component */
       cs_real_t   dpdz[],   /* --> gradient z component */
       cs_real_t   bx[],     /* --- local work array */
       cs_real_t   by[],     /* --- local work array */
       cs_real_t   bz[]      /* --- local work array */
);

/*----------------------------------------------------------------------------
 * Compute gradients using reconstruction method
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (gradrc, GRADRC)
(
 const cs_int_t   *const ncelet,      /* --> number of extended cells         */
 const cs_int_t   *const ncel,        /* --> number of cells                  */
 const cs_int_t   *const nfac,        /* --> number of internal faces         */
 const cs_int_t   *const nfabor,      /* --> number of boundary faces         */
 const cs_int_t   *const ncelbr,      /* --> number of cells on boundary      */
 const cs_int_t   *const imrgra,      /* --> gradient computation mode        */
 const cs_int_t   *const inc,         /* --> 0 or 1: increment or not         */
 const cs_int_t   *const iccocg,      /* --> 1 or 0: recompute COCG or not    */
 const cs_int_t   *const nswrgp,      /* --> >1: with reconstruction          */
 const cs_int_t   *const idimte,      /* --> 0, 1, 2: scalar, vector, tensor  */
 const cs_int_t   *const itenso,      /* --> for rotational periodicity       */
 const cs_int_t   *const iphydp,      /* --> use hydrosatatic pressure        */
 const cs_int_t   *const iwarnp,      /* --> verbosity level                  */
 const cs_int_t   *const nfecra,      /* --> standard output unit             */
 const cs_real_t  *const epsrgp,      /* --> precision for iterative gradient
                                             calculation                      */
 const cs_real_t  *const extrap,      /* --> extrapolate gradient at boundary */
 const cs_int_t          ifacel[],    /* --> interior face->cell connectivity */
 const cs_int_t          ifabor[],    /* --> boundary face->cell connectivity */
 const cs_int_t          icelbr[],    /* --> list of cells on boundary        */
 const cs_int_t   *const ivar,        /* --> variable number                  */
 const cs_real_t         volume[],    /* --> cell volumes                     */
 const cs_real_t         surfac[],    /* --> surfaces of internal faces       */
 const cs_real_t         surfbo[],    /* --> surfaces of boundary faces       */
 const cs_real_t         pond[],      /* --> interior faces geometric weight  */
 const cs_real_t         xyzcen[],    /* --> cell centers                     */
 const cs_real_t         cdgfac[],    /* --> interior face centers of gravity */
 const cs_real_t         cdgfbo[],    /* --> boundary face centers of gravity */
 const cs_real_t         dijpf[],     /* --> interior faces I'J' vector       */
 const cs_real_t         diipb[],     /* --> boundary faces II' vector        */
 const cs_real_t         dofij[],     /* --> interior faces OF vector         */
 const cs_real_t         fextx[],     /* --> components of the exterior force */
 const cs_real_t         fexty[],     /*     generating the hydrostatic       */
 const cs_real_t         fextz[],     /*     pressure                         */
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
 * Main Fortran options initialization
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (initi1, INITI1)
(
 const cs_int_t  *iverif          /* <-- Activate elementary tests */
);

/*----------------------------------------------------------------------------
 * Get parameters necessary for post-processing initialization from
 * the Fortran API.
 *----------------------------------------------------------------------------*/

void CS_PROCF (inipst, INIPST)
(
 const cs_int_t  *ichrvl,    /* <-- fluid volume post processing indicator */
 const cs_int_t  *ichrbo,    /* <-- boundary faces post processing indicator */
 const cs_int_t  *ipstmd,    /* <-- deformable mesh flag
                              *     0: no time dependency
                              *     1: deformable post processing meshes
                              *     2: output of a displacement field */
 const cs_int_t  *ntchr,     /* <-- frequency of post processing output */
 const cs_real_t *frchr,     /* <-- frequency of post processing output */
 const char      *fmtchr,    /* <-- name of main output format */
 const char      *optchr     /* <-- options for main output format */
);

/*----------------------------------------------------------------------------
 * Update mesh dimensions in Fortran commons
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (majgeo, MAJGEO)
(
 const cs_int_t   *ncel,    /* <-- number of cells */
 const cs_int_t   *ncelet,  /* <-- number of extended (real + ghost) cells */
 const cs_int_t   *nfac,    /* <-- number of interior faces */
 const cs_int_t   *nfabor,  /* <-- number of boundary faces */
 const cs_int_t   *nsom,    /* <-- number of vertices */
 const cs_int_t   *lndfac,  /* <-- interior face -> vertices array size */
 const cs_int_t   *lndfbr,  /* <-- boundary face -> vertices array size */
 const cs_int_t   *ncelbr,  /* <-- number of boundary cells */
 const cs_int_t   *ncelgb,  /* <-- global number of cells */
 const cs_int_t   *nfacgb,  /* <-- global number of interior faces */
 const cs_int_t   *nfbrgb,  /* <-- global number of boundary faces */
 const cs_int_t   *nsomgb,  /* <-- global number of vertices */
 const cs_int_t   *nthrdi,  /* <-- max threads per interior faces group */
 const cs_int_t   *nthrdb,  /* <-- max threads per boundary faces group */
 const cs_int_t   *ngrpi,   /* <-- number of interior face groups */
 const cs_int_t   *ngrpb,   /* <-- number of boundary face groups */
 const cs_int_t   *idxfi,   /* <-- interior face group/thread start/end ids */
 const cs_int_t   *idxfb,   /* <-- boundary face group/thread start/end ids */
 const cs_int_t    ifacel[],  /* <-- interior faces -> cells connectivity */
 const cs_int_t    ifabor[],  /* <-- boundary faces -> cells connectivity */
 const cs_int_t    ifmfbr[],  /* <-- boundary face families */
 const cs_int_t    ifmcel[],  /* <-- cell families */
 const cs_int_t    iprfml[],  /* <-- list of family properties */
 const cs_int_t    ipnfac[],  /* <-- interior faces -> vertices index */
 const cs_int_t    nodfac[],  /* <-- interior faces -> vertices connectivity */
 const cs_int_t    ipnfbr[],  /* <-- boundary faces -> vertices index */
 const cs_int_t    nodfbr[],  /* <-- boundary faces -> vertices connectivity */
 const cs_int_t    icelbr[],  /* <-- list of boundary cells */
 const cs_real_t  *volmin,    /* <-- minimum control volume */
 const cs_real_t  *volmax,    /* <-- maximum control volume */
 const cs_real_t  *voltot,    /* <-- total   control volume */
 const cs_real_t   xyzcen[],  /* <-- cell centers */
 const cs_real_t   surfac[],  /* <-- interior face surface vectors */
 const cs_real_t   surfbo[],  /* <-- boundary face surface vectors */
 const cs_real_t   cdgfac[],  /* <-- interior face centers */
 const cs_real_t   cdgfbo[],  /* <-- boundary face centers */
 const cs_real_t   xyznod[],  /* <-- vertex coordinates */
 const cs_real_t   volume[],  /* <-- cell volumes */
 const cs_real_t   surfan[],  /* <-- interior face surfaces */
 const cs_real_t   surfbn[],  /* <-- boundary face surfaces */
 const cs_real_t   dist[],    /* <-- distance IJ.Nij */
 const cs_real_t   distb[],   /* <-- likewise for border faces */
 const cs_real_t   pond[],    /* <-- weighting (Aij=pond Ai+(1-pond)Aj */
 const cs_real_t   dijpf[],   /* <-- vector I'J' */
 const cs_real_t   diipb[],   /* <-- likewise for border faces */
 const cs_real_t   dofij[]    /* <-- vector OF at interior faces */
);

/*----------------------------------------------------------------------------
 * Free Fortran allocated memory
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (memfin, MEMFIN) (void);

/*----------------------------------------------------------------------------
 * Mesh renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (numvec, NUMVEC)
(
 const cs_int_t  *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t  *ncel,     /* <-- number of local cells */
 const cs_int_t  *nfac,     /* <-- number of interior faces */
 const cs_int_t  *nfabor,   /* <-- number of boundary faces */
 const cs_int_t  *lregis,   /* <-- vector registor length */
 cs_int_t        *irveci,   /* <-> interior face vectorization indic. */
 cs_int_t        *irvecb,   /* <-> boundary face vectorization indic. */
 cs_int_t         ifacel[], /* <-- interior face->cell connectivity */
 cs_int_t         ifabor[], /* <-- boundary face->cell connectivity */
 cs_int_t         inumfi[], /* <-> interior faces renumbering (size: nfac) */
 cs_int_t         inumfb[], /* <-> boundary faces renumbering (size: nfabor) */
 cs_int_t         iworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_int_t         ismbs[]   /* --- work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * Test renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (tstvec, TSTVEC)
(
 const cs_int_t  *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t  *ncel,     /* <-- number of local cells */
 const cs_int_t  *nfac,     /* <-- number of interior faces */
 const cs_int_t  *nfabor,   /* <-- number of boundary faces */
 const cs_int_t   ifacel[], /* <-- interior face->cell connectivity */
 const cs_int_t   ifabor[], /* <-- boundary face->cell connectivity */
 cs_int_t         iworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_int_t         ismbs[],  /* --- work array, size: ncelet */
 cs_int_t         ismbv[],  /* --- work array, size: ncelet */
 cs_real_t        rworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_real_t        rsmbs[],  /* --- work array, size: ncelet */
 cs_real_t        rsmbv[]   /* --- work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * User function for modification of a post-processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (usmpst, USMPST)
(
 const cs_int_t  *nummai,    /* <-- number or post-processing mesh */
 const cs_int_t  *nvar,      /* <-- number of variables */
 const cs_int_t  *nscal,     /* <-- number of scalars */
 const cs_int_t  *nvlsta,    /* <-- number of statistical variables (lagr) */
 cs_int_t        *ncelps,    /* <-> number of post-processed cells */
 cs_int_t        *nfacps,    /* <-> number of post processed interior faces */
 cs_int_t        *nfbrps,    /* <-> number of post processed boundary faces */
 cs_int_t        *imodif,    /* <-> 1 if mesh is modified, 0 otherwise */
 const cs_int_t   itypps[3], /* <-- flag (0 or 1) for presence of cells, */
                             /*     interior faces, and boundary faces */
 cs_int_t         lstcel[],  /* <-> list of post-processed cells */
 cs_int_t         lstfac[],  /* <-> list of post-processed interior faces */
 cs_int_t         lstfbr[],  /* <-> list of post-processed boundary faces */
 const cs_real_t  dt[],      /* <-- local time step */
 const cs_real_t  rtpa[],    /* <-- cell variables at previous time step */
 const cs_real_t  rtp[],     /* <-- cell variables */
 const cs_real_t  propce[],  /* <-- cell physical properties */
 const cs_real_t  propfa[],  /* <-- interior face physical properties */
 const cs_real_t  propfb[],  /* <-- boundary face physical properties */
 const cs_real_t  coefa[],   /* <-- boundary conditions array */
 const cs_real_t  coefb[],   /* <-- boundary conditions array */
 const cs_real_t  statce[],  /* <-- cell statistics (Lagrangian) */
 cs_real_t        tracel[],  /* --- work array for output cells */
 cs_real_t        trafac[],  /* --- work array for output interior faces */
 cs_real_t        trafbr[]   /* --- work array for output boundary faces */
);

/*----------------------------------------------------------------------------
 * User function for output of variables on a post-processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (usvpst, USVPST)
(
 const cs_int_t  *nummai,    /* <-- number or post-processing mesh */
 const cs_int_t  *nvar,      /* <-- number of variables */
 const cs_int_t  *nscal,     /* <-- number of scalars */
 const cs_int_t  *nvlsta,    /* <-- number of statistical variables (lagr) */
 const cs_int_t  *ncelps,    /* <-- number of post-processed cells */
 const cs_int_t  *nfacps,    /* <-- number of post processed interior faces */
 const cs_int_t  *nfbrps,    /* <-- number of post processed boundary faces */
 const cs_int_t   itypps[3], /* <-- flag (0 or 1) for presence of cells, */
                             /*     interior faces, and boundary faces */
 const cs_int_t   lstcel[],  /* <-- list of post-processed cells */
 const cs_int_t   lstfac[],  /* <-- list of post-processed interior faces */
 const cs_int_t   lstfbr[],  /* <-- list of post-processed boundary faces */
 const cs_real_t  dt[],      /* <-- local time step */
 const cs_real_t  rtpa[],    /* <-- cell variables at previous time step */
 const cs_real_t  rtp[],     /* <-- cell variables */
 const cs_real_t  propce[],  /* <-- cell physical properties */
 const cs_real_t  propfa[],  /* <-- interior face physical properties */
 const cs_real_t  propfb[],  /* <-- boundary face physical properties */
 const cs_real_t  coefa[],   /* <-- boundary conditions array */
 const cs_real_t  coefb[],   /* <-- boundary conditions array */
 const cs_real_t  statce[],  /* <-- cell statistics (Lagrangian) */
 cs_real_t        tracel[],  /* --- work array for output cells */
 cs_real_t        trafac[],  /* --- work array for output interior faces */
 cs_real_t        trafbr[]   /* --- work array for output boundary faces */
);

/*----------------------------------------------------------------------------
 * User function to compute coarsening array for algebraic multigrid
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (ustmgr, USTMGR)
(
 const cs_int_t   *iappel,      /* <-- 1: initialization call
                                       2: computional call */
 const cs_int_t   *igr,         /* <-- new grid level (0 = base) */
 const cs_int_t   *isym,        /* <-- 1: symmetric; 2 nonsymmteric */
 const cs_int_t   *ncelf,       /* <-- number of cells in fine grid */
 const cs_int_t   *ncelfe,      /* <-- number of cells with halo in fine grid */
 const cs_int_t   *nfacf,       /* <-- number of faces in fine grid */
 const cs_int_t   *iwarnp,      /* <-- verbosity level */
 cs_int_t         *iusmgr,      /* --> 0: automatic method
                                       1: use this sub-routine */
 cs_int_t         *niw,         /* --> size of iw for call 2 */
 cs_int_t         *nrw,         /* --> size of rw for call 2 */
 const cs_int_t    ifacef[],    /* <-- fine grid face->cell connectivity */
 const cs_real_t   daf[],       /* <-- diagonal terms of fine grid */
 const cs_real_t   xaf[],       /* <-- extradiagonal terms of fine grid */
 const cs_real_t   surfaf[],    /* <-- fine grid face surface vectors */
 const cs_real_t   volumf[],    /* <-- fine grid cell volumes */
 const cs_real_t   xyzfin[],    /* <-- fine grid cell centers */
 cs_int_t          irscel[],    /* --> Fine -> coarse cell connectivity */
 cs_int_t          iw[],        /* --> work array of size niw (call 2) */
 cs_real_t         rw[]         /* --> work array of size nrw (call 2) */
);

/*============================================================================
 *  User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define mesh joinings.
 *----------------------------------------------------------------------------*/

void
cs_user_join(void);

/*----------------------------------------------------------------------------
 * Define mesh files to read and optional associated transformations.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_input(void);

/*----------------------------------------------------------------------------
 * Modifiy geometry and mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Define periodic faces.
 *----------------------------------------------------------------------------*/

void
cs_user_periodicity(void);

/*----------------------------------------------------------------------------
 * Define couplings with other instances of Code_Saturne.
 *----------------------------------------------------------------------------*/

void
cs_user_saturne_coupling(void);

/*----------------------------------------------------------------------------
 * Set user solver.
 *----------------------------------------------------------------------------*/

int
cs_user_solver_set(void);

/*----------------------------------------------------------------------------
 * Main call to user solver.
 *----------------------------------------------------------------------------*/

void
cs_user_solver(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Define couplings with SYRTHES code.
 *----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
