#ifndef __CS_SLES_H__
#define __CS_SLES_H__

/*============================================================================
 * Sparse Linear Equation Solvers
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
#include "cs_halo_perio.h"
#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
  CS_SLES_PCG_SR,    /* Preconditionned conjugate gradient, single reduction*/
  CS_SLES_JACOBI,    /* Jacobi */
  CS_SLES_BICGSTAB,  /* Bi-conjugate gradient stabilized */
  CS_SLES_GMRES,     /* Generalized minimal residual */
  CS_SLES_N_TYPES    /* Number of resolution algorithms */

} cs_sles_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

extern const char *cs_sles_type_name[];

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(reslin, RESLIN)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ilved,     /* <-- Interleaved indicator  */
                              /*     1: interleaved; 2: not interleaved */
 const cs_int_t   *ibsize,    /* <-- Block size of element ii,ii */
 const cs_int_t   *iesize,    /* <-- Block size of element ij */
 const cs_int_t   *ireslp,    /* <-- Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab, 3: gmres,
                                     200: pcg_sr */
 const cs_int_t   *ipol,      /* <-- Preconditioning polynomial degree
                                     (0: diagonal; -1: non-preconditionned) */
 const cs_int_t   *nitmap,    /* <-- Number of max iterations */
 const cs_int_t   *iinvpe,    /* <-- Indicator to cancel increments
                                     in rotational periodicty (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 cs_int_t         *niterf,    /* --> Number of iterations done */
 const cs_real_t  *epsilp,    /* <-- Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* <-- Residue normalization */
 cs_real_t        *residu,    /* --> Final non normalized residue */
 const cs_int_t   *ifacel,    /* <-- Face -> cell connectivity  */
 const cs_real_t  *dam,       /* <-- Matrix diagonal */
 const cs_real_t  *xam,       /* <-- Matrix extra-diagonal terms */
 const cs_real_t  *smbrp,     /* <-- System right-hand side */
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

/*----------------------------------------------------------------------------
 * Update sparse linear equation solver API in case of mesh modification.
 *----------------------------------------------------------------------------*/

void
cs_sles_update_mesh(void);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Set MPI communicator for dot products.
 *----------------------------------------------------------------------------*/

void
cs_sles_set_mpi_reduce_comm(MPI_Comm comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Test if a general sparse linear system needs solving or if the right-hand
 * side is already zero within convergence criteria.
 *
 * The computed residue is also updated;
 *
 * parameters:
 *   var_name      <-- Variable name
 *   solver_name   <-- Name of solver
 *   n_rows        <-- Number of (non ghost) rows in rhs
 *   verbosity     <-- Verbosity level
 *   r_norm        <-- Residue normalization
 *   residue       <-> Residue
 *   rhs           <-- Right hand side
 *
 * returns:
 *   1 if solving is required, 0 if the rhs is already zero within tolerance
 *   criteria (precision of residue normalization)
 *----------------------------------------------------------------------------*/

int
cs_sles_needs_solving(const char        *var_name,
                      const char        *solver_name,
                      cs_int_t           n_rows,
                      int                verbosity,
                      double             r_norm,
                      double            *residue,
                      const cs_real_t   *rhs);

/*----------------------------------------------------------------------------
 * General sparse linear system resolution.
 *
 * parameters:
 *   var_name          <-- Variable name
 *   solver_type       <-- Type of solver (PCG, Jacobi, ...)
 *   update_stats      <-- Automatic solver statistics indicator
 *   symmetric         <-- Symmetric coefficients indicator
 *   a                 <-- Matrix
 *   poly_degree       <-- Preconditioning polynomial degree
 *                         (0: diagonal; -1: non-preconditioned)
 *   rotation_mode     <-- Halo update option for rotational periodicity
 *   verbosity         <-- Verbosity level
 *   n_max_iter        <-- Maximum number of iterations
 *   precision         <-- Precision limit
 *   r_norm            <-- Residue normalization
 *   n_iter            --> Number of iterations
 *   residue           <-> Residue
 *   rhs               <-- Right hand side
 *   vx                --> System solution
 *   aux_size          <-- Number of elements in aux_vectors
 *   aux_vectors       --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

int
cs_sles_solve(const char          *var_name,
              cs_sles_type_t       solver_type,
              bool                 update_stats,
              const cs_matrix_t   *a,
              int                  poly_degree,
              cs_halo_rotation_t   rotation_mode,
              int                  verbosity,
              int                  n_max_iter,
              double               precision,
              double               r_norm,
              int                 *n_iter,
              double              *residue,
              const cs_real_t     *rhs,
              cs_real_t           *vx,
              size_t               aux_size,
              void                *aux_vectors);

/*----------------------------------------------------------------------------
 * Output default post-processing data for failed system convergence.
 *
 * parameters:
 *   var_name         <-- Variable name
 *   mesh_id          <-- id of error output mesh, or 0 if none
 *   rotation_mode    <-- Halo update option for rotational periodicity
 *   a                <-- Linear equation matrix
 *   rhs              <-- Right hand side
 *   vx               <-> Current system solution
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_def(const char          *var_name,
                              int                  mesh_id,
                              cs_halo_rotation_t   rotation_mode,
                              const cs_matrix_t   *a,
                              const cs_real_t     *rhs,
                              cs_real_t           *vx);

/*----------------------------------------------------------------------------
 * Output post-processing variable for failed system convergence.
 *
 * parameters:
 *   var_name        <-- Variable name
 *   diag_block_size <-- Block size for diagonal
 *   mesh_id         <-- id of error output mesh, or 0 if none
 *   var             <-- Variable values
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_var(const char  *var_name,
                              int          mesh_id,
                              int          diag_block_size,
                              cs_real_t   *var);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_H__ */
