#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells_detection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Fortran function/subroutine prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Main Fortran subroutine
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (caltri, CALTRI)
(
 const cs_int_t   *iverif   /* <-- activate elementary tests */
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
 cs_real_t        trafbr[]   /* --- work array for output boundary faces */
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
 * Main Fortran options initialization
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (initi1, INITI1)
(
 const cs_int_t  *iverif          /* <-- Activate elementary tests */
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
 * User override of default frequency or calculation end based output.
 *
 * Fortran interface:
 *
 * subroutine pstusn (ntmabs, ntcabs, ttcabs)
 * *****************
 *
 * integer          ntmabs      : <-- : maximum time step number
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstusn, PSTUSN)
(
 const cs_int_t  *ntmabs,
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
);

/*----------------------------------------------------------------------------
 * Indicate if the variable considered is a component of a vector or tensor
 * in the presence of periodicity of rotation
 *
 * Fortran interface:
 *
 * subroutine pergra (ivar, ipvar)
 * *****************
 *
 * integer    ivar      : <-- : variable number
 * integer    idimtr    : --> : 0 if ivar does not match a vector or tensor
 *                      :     :   or there is no periodicity of rotation
 *                      :     : 1 for velocity, 2 for Reynolds stress
 *                      :     :   in case of periodicity of rotation
 * integer    irpvar    : --> : -1 if ivar does not match a vector or tensor
 *                      :     :    or there is no periodicity of rotation
 *                      :     : In presence of periodicity of rotation:
 *                      :     :  0 for iu, 1 for iv, 2 for iw
 *                      :     :  0 for ir11, 1 for ir22, 2 for ir33,
 *                      :     :  3 for ir12, 4 for ir13, 5 for ir23
 *----------------------------------------------------------------------------*/

void CS_PROCF (pergra, PERGRA)
(
 const cs_int_t  *ivar,
 const cs_int_t  *idimtr,
 const cs_int_t  *irpvar
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
 * Insert thin wall into a mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_thinwall(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Mesh smoothing.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_smoothe(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Tag bad cells within the mesh based on geometric criteria.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_bad_cells_tag(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Define periodic faces.
 *----------------------------------------------------------------------------*/

void
cs_user_periodicity(void);

/*----------------------------------------------------------------------------
 * Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void);

/*----------------------------------------------------------------------------
 * Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void);

/*----------------------------------------------------------------------------
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * parameters:
 *   nt_max_abs <-- maximum time step number
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs);

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
