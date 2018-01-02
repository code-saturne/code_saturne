#ifndef __CS_GRID_H__
#define __CS_GRID_H__

/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

#include "cs_halo.h"
#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure associated with opaque grid object */

typedef struct _cs_grid_t cs_grid_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names for coarsening options */

extern const char *cs_grid_coarsening_type_name[];

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the multigrid layer
 * (cs_multigrid.c), not directly by the user, so they are no more
 * documented than private static functions)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create base grid by mapping from shared mesh values.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   n_cells               <-- Local number of cells
 *   n_cells_ext           <-- Local number of cells + ghost cells
 *   n_faces               <-- Local number of faces
 *   symmetric             <-- True if xam is symmetric, false otherwise
 *   diag_block_size       <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size <-- Block sizes for extra diagonal, or NULL
 *   face_cell             <-- Face -> cells connectivity
 *   halo                  <-- Halo structure associated with this level,
 *                             or NULL.
 *   cell_cen              <-- Cell center (size: 3.n_cells_ext)
 *   cell_vol              <-- Cell volume (size: n_cells_ext)
 *   face_normal           <-- Internal face normals (size: 3.n_faces)
 *   a                     <-- Associated matrix
 *   a_conv                <-- Associated matrix (convection)
 *   a_diff                <-- Associated matrix (diffusion)
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_shared(cs_lnum_t              n_cells,
                           cs_lnum_t              n_cells_ext,
                           cs_lnum_t              n_faces,
                           bool                   symmetric,
                           const int             *diag_block_size,
                           const int             *extra_diag_block_size,
                           const cs_lnum_2_t     *face_cell,
                           const cs_halo_t       *halo,
                           const cs_real_t       *cell_cen,
                           const cs_real_t       *cell_vol,
                           const cs_real_t       *face_normal,
                           const cs_matrix_t     *a,
                           const cs_matrix_t     *a_conv,
                           const cs_matrix_t     *a_diff);

/*----------------------------------------------------------------------------
 * Destroy a grid structure.
 *
 * parameters:
 *   grid <-> Pointer to grid structure pointer
 *----------------------------------------------------------------------------*/

void
cs_grid_destroy(cs_grid_t **grid);

/*----------------------------------------------------------------------------
 * Get grid information.
 *
 * parameters:
 *   g           <-- Grid structure
 *   level       --> Level in multigrid hierarchy (or NULL)
 *   symmetric   --> Symmetric matrix coefficients indicator (or NULL)
 *   db_size     --> Size of the diagonal block (or NULL)
 *   eb_size     --> Size of the extra diagonal block (or NULL)
 *   n_ranks     --> number of ranks with data (or NULL)
 *   n_cells     --> Number of local cells (or NULL)
 *   n_cells_ext --> Number of cells including ghosts (or NULL)
 *   n_faces     --> Number of faces (or NULL)
 *   n_g_cells   --> Number of global cells (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 bool             *symmetric,
                 int              *db_size,
                 int              *eb_size,
                 int              *n_ranks,
                 cs_lnum_t        *n_cells,
                 cs_lnum_t        *n_cells_ext,
                 cs_lnum_t        *n_faces,
                 cs_gnum_t        *n_g_cells);

/*----------------------------------------------------------------------------
 * Get number of cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of cells of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cells(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get number of extended (local + ghost) cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of extended cells of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cells_ext(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get maximum number of extended (local + ghost) cells corresponding to
 * a grid, both with and without merging between ranks
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   maximum number of extended cells of grid structure, with or without
 *   merging
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cells_max(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get global number of cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   global number of cells of grid structure
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_grid_get_n_g_cells(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get grid's associated matrix information.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   pointer to matrix structure
 *----------------------------------------------------------------------------*/

const cs_matrix_t *
cs_grid_get_matrix(const cs_grid_t  *g);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the MPI subcommunicator for a given grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
cs_grid_get_comm(const cs_grid_t  *g);

#endif

/*----------------------------------------------------------------------------
 * Create coarse grid from fine grid.
 *
 * parameters:
 *   f                    <-- Fine grid structure
 *   verbosity            <-- Verbosity level
 *   coarsening_type      <-- Coarsening traversal type:
 *                              0: algebraic with natural face traversal;
 *                              1: algebraic with face traveral by criteria;
 *                              2: algebraic with Hilbert face traversal
 *   aggregation_limit    <-- Maximum allowed fine cells per coarse cell
 *   relaxation_parameter <-- P0/P1 relaxation factor
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen(const cs_grid_t  *f,
                int               verbosity,
                int               coarsening_type,
                int               aggregation_limit,
                double            relaxation_parameter);

/*----------------------------------------------------------------------------
 * Compute coarse cell variable values from fine cell values
 *
 * parameters:
 *   f       <-- Fine grid structure
 *   c       <-- Fine grid structure
 *   f_var   <-- Variable defined on fine grid cells
 *   c_var   --> Variable defined on coarse grid cells
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

void
cs_grid_restrict_cell_var(const cs_grid_t  *f,
                          const cs_grid_t  *c,
                          const cs_real_t  *f_var,
                          cs_real_t        *c_var);

/*----------------------------------------------------------------------------
 * Compute fine cell integer values from coarse cell values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_num   --> Variable defined on coarse grid cells
 *   f_num   <-- Variable defined on fine grid cells
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_cell_num(const cs_grid_t  *c,
                         const cs_grid_t  *f,
                         int              *c_num,
                         int              *f_num);

/*----------------------------------------------------------------------------
 * Compute fine cell variable values from coarse cell values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_var   --> Variable defined on coarse grid cells
 *   f_var   <-- Variable defined on fine grid cells
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_cell_var(const cs_grid_t  *c,
                         const cs_grid_t  *f,
                         cs_real_t        *c_var,
                         cs_real_t        *f_var);

/*----------------------------------------------------------------------------
 * Project coarse grid cell numbers to base grid.
 *
 * If a global coarse grid cell number is larger than max_num, its
 * value modulo max_num is used.
 *
 * parameters:
 *   g             <-- Grid structure
 *   n_base_cells  <-- Number of cells in base grid
 *   max_num       <-- Values of c_cell_num = global_num % max_num
 *   c_cell_num    --> Global coarse cell number (modulo max_num)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_cell_num(const cs_grid_t  *g,
                         cs_lnum_t         n_base_cells,
                         int               max_num,
                         int               c_cell_num[]);

/*----------------------------------------------------------------------------
 * Project coarse grid cell rank to base grid.
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   f_cell_rank  --> Global coarse cell rank projected to fine cells
 *----------------------------------------------------------------------------*/

void
cs_grid_project_cell_rank(const cs_grid_t  *g,
                          cs_lnum_t         n_base_cells,
                          int               f_cell_rank[]);

/*----------------------------------------------------------------------------
 * Project variable from coarse grid to base grid
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   c_var        <-- Cell variable on coarse grid
 *   f_var        --> Cell variable projected to fine grid
 *----------------------------------------------------------------------------*/

void
cs_grid_project_var(const cs_grid_t  *g,
                    cs_lnum_t         n_base_cells,
                    const cs_real_t   c_var[],
                    cs_real_t         f_var[]);

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric and project it to base grid
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   diag_dom     --> Diagonal dominance metric (on fine grid)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_diag_dom(const cs_grid_t  *g,
                         cs_lnum_t         n_base_cells,
                         cs_real_t         diag_dom[]);

/*----------------------------------------------------------------------------
 * Return the merge_stride if merging is active.
 *
 * returns:
 *   grid merge stride if merging is active, 1 otherwise
 *----------------------------------------------------------------------------*/

int
cs_grid_get_merge_stride(void);

/*----------------------------------------------------------------------------
 * Finalize global info related to multigrid solvers
 *----------------------------------------------------------------------------*/

void
cs_grid_finalize(void);

/*----------------------------------------------------------------------------
 * Dump grid structure
 *
 * parameters:
 *   g <-- grid structure that should be dumped
 *----------------------------------------------------------------------------*/

void
cs_grid_dump(const cs_grid_t  *g);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query the global multigrid parameters for parallel grid merging.
 *
 * parameters:
 *   rank_stride          --> number of ranks over which merging
 *                            takes place, or NULL
 *   cells_mean_threshold --> mean number of cells under which merging
 *                            should be applied, or NULL
 *   cells_glob_threshold --> global number of cells under which merging
 *                            should be applied, or NULL
 *   min_ranks            --> number of active ranks under which
 *                            no merging takes place, or NULL
 *----------------------------------------------------------------------------*/

void
cs_grid_get_merge_options(int         *rank_stride,
                          int         *cells_mean_threshold,
                          cs_gnum_t   *cells_glob_threshold,
                          int         *min_ranks);

/*----------------------------------------------------------------------------
 * Set global multigrid parameters for parallel grid merging.
 *
 * parameters:
 *   rank_stride          <-- number of ranks over which merging
 *                            takes place
 *   cells_mean_threshold <-- mean number of cells under which merging
 *                            should be applied
 *   cells_glob_threshold <-- global number of cells under which merging
 *                            should be applied
 *   min_ranks            <-- number of active ranks under which
 *                            no merging takes place
 *----------------------------------------------------------------------------*/

void
cs_grid_set_merge_options(int         rank_stride,
                          int         cells_mean_threshold,
                          cs_gnum_t   cells_glob_threshold,
                          int         min_ranks);

/*----------------------------------------------------------------------------
 * Set matrix tuning behavior for multigrid coarse meshes.
 *
 * The finest mesh (level 0) is handled by the default tuning options,
 * so only coarser meshes are considered here.
 *
 * parameters:
 *   fill_type <-- associated matrix fill type
 *   max_level <-- maximum level for which tuning is active
 *----------------------------------------------------------------------------*/

void
cs_grid_set_matrix_tuning(cs_matrix_fill_type_t  fill_type,
                          int                    max_level);

/*----------------------------------------------------------------------------
 * Force matrix variant selection for multigrid coarse meshes.
 *
 * The finest mesh (level 0) is handled by the default options,
 * so only coarser meshes are considered here.
 *
 * parameters:
 *   fill_type <-- associated matrix fill type
 *   level     <-- level for which variant is assiged
 *   mv        <-- matrix variant to assign (NULL to unassign)
 *----------------------------------------------------------------------------*/

void
cs_grid_set_matrix_variant(cs_matrix_fill_type_t       fill_type,
                           int                         level,
                           const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Log the current settings for multigrid parallel merging.
 *----------------------------------------------------------------------------*/

void
cs_grid_log_merge_options(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GRID_H__ */
