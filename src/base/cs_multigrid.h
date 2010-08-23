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

#ifndef __CS_MULTIGRID_H__
#define __CS_MULTIGRID_H__

/*============================================================================
 * Multigrid solver.
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build a hierarchy of meshes starting from a fine mesh, for an
 * ACM (Additive Corrective Multigrid) method, grouping cells at
 * most 2 by 2.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmlga, CLMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of internal faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 cs_int_t         *iagmax,    /* <-> Maximum agglomeration count */
 const cs_int_t   *nagmax,    /* <-- Agglomeration count limit */
 const cs_int_t   *ncpost,    /* <-- If > 0, postprocess coarsening, using
                                     coarse cell numbers modulo ncpost */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 const cs_int_t   *ngrmax,    /* <-- Maximum number of grid levels */
 const cs_int_t   *ncegrm,    /* <-- Maximum local number of cells on
                                     coarsest grid */
 const cs_real_t  *rlxp1,     /* <-- P0/P1 relaxation parameter */
 const cs_real_t  *dam,       /* <-- Matrix diagonal */
 const cs_real_t  *xam        /* <-- Matrix extra-diagonal terms */
);

/*----------------------------------------------------------------------------
 * Destroy a hierarchy of meshes starting from a fine mesh, keeping
 * the corresponding system information for future calls.
 *----------------------------------------------------------------------------*/

void CS_PROCF(dsmlga, DSMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname      /* <-- variable name length */
);

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(resmgr, RESMGR)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *iresds,    /* <-- Descent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *iresas,    /* <-- Ascent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ireslp,    /* <-- Coarse Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ipol,      /* <-- Preconditioning polynomial degree
                                     (0: diagonal) */
 const cs_int_t   *ncymxp,    /* <-- Max number of cycles */
 const cs_int_t   *nitmds,    /* <-- Max number of iterations for descent */
 const cs_int_t   *nitmas,    /* <-- Max number of iterations for ascent */
 const cs_int_t   *nitmap,    /* <-- Max number of iterations for
                                     coarsest solution */
 const cs_int_t   *iinvpe,    /* <-- Indicator to cancel increments
                                     in rotational periodicity (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 cs_int_t         *ncyclf,    /* --> Number of cycles done */
 cs_int_t         *niterf,    /* --> Number of iterations done */
 const cs_real_t  *epsilp,    /* <-- Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* <-- Residue normalization */
 cs_real_t        *residu,    /* --> Final non normalized residue */
 const cs_int_t   *ifacel,    /* <-- Face -> cell connectivity  */
 const cs_real_t  *rhs,       /* <-- System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize multigrid solver API.
 *----------------------------------------------------------------------------*/

void
cs_multigrid_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize multigrid solver API.
 *----------------------------------------------------------------------------*/

void
cs_multigrid_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MULTIGRID_H__ */
