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

#ifndef __CS_SLES_H__
#define __CS_SLES_H__

/*============================================================================
 * Sparse Linear Equation Solvers
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solver types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SLES_PCG,       /* Preconditionned conjugate gradient */
  CS_SLES_JACOBI,    /* Jacobi */
  CS_SLES_BICGSTAB,  /* Bi-conjugate gradient stabilized */
  CS_SLES_N_TYPES    /* Number of implemented resolution algorithms */

} cs_sles_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(reslin, RESLIN)
(
 const char       *cname,     /* --> variable name */
 const cs_int_t   *lname,     /* --> variable name length */
 const cs_int_t   *ncelet,    /* --> Number of cells, halo included */
 const cs_int_t   *ncel,      /* --> Number of local cells */
 const cs_int_t   *nfac,      /* --> Number of faces */
 const cs_int_t   *isym,      /* --> Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ireslp,    /* --> Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ipol,      /* --> Preconditioning polynomial degree
                                     (0: diagonal) */
 const cs_int_t   *nitmap,    /* --> Number of max iterations */
 const cs_int_t   *iinvpe,    /* --> Indicator to cancel increments
                                     in rotational periodicty (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* --> Verbosity level */
 cs_int_t         *niterf,    /* <-- Number of iterations done */
 const cs_real_t  *epsilp,    /* --> Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* --> Residue normalization */
 cs_real_t        *residu,    /* <-- Final non normalized residue */
 const cs_int_t   *ifacel,    /* --> Face -> cell connectivity  */
 const cs_real_t  *dam,       /* --> Matrix diagonal */
 const cs_real_t  *xam,       /* --> Matrix extra-diagonal terms */
 const cs_real_t  *smbrp,     /* --> System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_finalize(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SLES_H__ */
