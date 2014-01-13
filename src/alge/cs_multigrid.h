#ifndef __CS_MULTIGRID_H__
#define __CS_MULTIGRID_H__

/*============================================================================
 * Multigrid solver.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include "cs_sles.h"

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
 * ACM (Additive Corrective Multigrid) method.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmlga, CLMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ibsize,    /* <-- Matrix block size */
 const cs_int_t   *iesize,    /* <-- Matrix extra diag block size */
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

/*----------------------------------------------------------------------------
 * Build a hierarchy of meshes starting from a fine mesh, for an
 * ACM (Additive Corrective Multigrid) method.
 *
 * parameters:
 *   var_name               <-- variable name
 *   verbosity              <-- verbosity level
 *   postprocess_block_size <-- if > 0, postprocess coarsening, using
 *                              coarse cell numbers modulo ncpost
 *   aggregation_limit      <-- maximum allowed fine cells per coarse cell
 *   n_max_levels           <-- maximum number of grid levels
 *   n_g_cells_min          <-- global number of cells on coarsest grid
 *                              under which no merging occurs
 *   p0p1_relax             <-- p0/p1 relaxation_parameter
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_build(const char       *var_name,
                   int               verbosity,
                   int               postprocess_block_size,
                   int               aggregation_limit,
                   int               n_max_levels,
                   cs_gnum_t         n_g_cells_min,
                   double            p0p1_relax,
                   bool              symmetric,
                   const int        *diag_block_size,
                   const int        *extra_diag_block_size,
                   const cs_real_t  *da,
                   const cs_real_t  *xa);

/*----------------------------------------------------------------------------
 * Destroy a hierarchy of meshes starting from a fine mesh, keeping
 * the corresponding system and postprocessing information for future calls.
 *
 * parameters:
 *   var_name <-- variable name
 *----------------------------------------------------------------------------*/

void
cs_multigrid_destroy(const char  *var_name);

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   var_name              <-- Variable name
 *   descent_smoother_type <-- Type of smoother for descent (PCG, Jacobi, ...)
 *   ascent_smoother_type  <-- Type of smoother for ascent (PCG, Jacobi, ...)
 *   coarse_solver_type    <-- Type of solver (PCG, Jacobi, ...)
 *   abort_on_divergence   <-- Call errorhandler if divergence is detected
 *   poly_degree           <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode         <-- Halo update option for rotational periodicity
 *   verbosity             <-- Verbosity level
 *   n_max_cycles          <-- Maximum number of cycles
 *   n_max_iter_descent    <-- Maximum nb. of iterations for descent phases
 *   n_max_iter_ascent     <-- Maximum nb. of iterations for ascent phases
 *   n_max_iter_coarse     <-- Maximum nb. of iterations for coarsest solution
 *   precision             <-- Precision limit
 *   r_norm                <-- Residue normalization
 *   n_cycles              --> Number of cycles
 *   n_iter                --> Number of iterations
 *   residue               <-> Residue
 *   rhs                   <-- Right hand side
 *   vx                    --> System solution
 *   aux_size              <-- Number of elements in aux_vectors
 *   aux_vectors           --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

int
cs_multigrid_solve(const char          *var_name,
                   cs_sles_type_t       descent_smoother_type,
                   cs_sles_type_t       ascent_smoother_type,
                   cs_sles_type_t       coarse_solver_type,
                   bool                 abort_on_divergence,
                   int                  poly_degree,
                   cs_halo_rotation_t   rotation_mode,
                   int                  verbosity,
                   int                  n_max_cycles,
                   int                  n_max_iter_descent,
                   int                  n_max_iter_ascent,
                   int                  n_max_iter_coarse,
                   double               precision,
                   double               r_norm,
                   int                 *n_cycles,
                   int                 *n_iter,
                   double              *residue,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MULTIGRID_H__ */
