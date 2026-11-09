#ifndef CS_GRID_PRIV_H
#define CS_GRID_PRIV_H

/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_base.h"

#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure associated with opaque grid object */

typedef struct _cs_grid_t cs_grid_t;

namespace cs {

/* Class associated with grid rank merging */
/*-----------------------------------------*/

#if defined(HAVE_MPI)

class grid_rank_merge {

public:
  int         rank_id{-1};       /* destination rank for this subset */
  int         rank_stride{0};    /* rank stride within this subset
                                    (equal to parent rank stride) */
  int         sub_root{-1};      /* source sub-root when merging */
  int         sub_rank{-1};      /* sub-rank when merging
                                    (0 <= sub-rank < merge_sub_size) */
  int         sub_size[2]{1,1};  /* number of merged ranks for this subset
                                    (actual and target)*/

  cs_lnum_t  *row_idx{nullptr};  /* start row_id for each sub-rank when
                                    merge_sub_rank = 0 (size: merge_size + 1) */

  void       *buf{nullptr};      /* additional buffer when merge
                                    send and receive overlap */
  size_t      buf_size{0};       /* size of merge buffer */

private:

  /* Prevent use of some non-relevant operators
     ------------------------------------------ */

  void * operator new[] (size_t);

public:

  /*--------------------------------------------------------------------------
   * Custom new and delete operators (implemented in cs_grid_rank_merge.cpp)
   *--------------------------------------------------------------------------*/

  void *
  operator new(std::size_t size);

  void
  operator delete(void *ptr);

  /*--------------------------------------------------------------------------
   * Constructors and destructors
   *--------------------------------------------------------------------------*/

  ~grid_rank_merge();

  /*--------------------------------------------------------------------------
   * Public memeber functions
   *--------------------------------------------------------------------------*/

  void
  dump();

};

#endif

} /* namespace cs */

/*----------------------------------------------------------------------------*/

/* Structure associated with grid */
/*--------------------------------*/

struct _cs_grid_t {

  int                 level;        /* Level in multigrid hierarchy */

  cs_alloc_mode_t     alloc_mode;   /* Preferred allocation mode */

  bool                conv_diff;    /* true if convection/diffusion case,
                                       false otherwise */
  bool                symmetric;    /* Symmetric matrix coefficients
                                       indicator */
  bool                use_faces;    /* True if face information is present */
  bool                use_xa;       /* True if face-based extra-diagonal
                                       are available */

  cs_lnum_t           db_size;      /* Block sizes for diagonal */
  cs_lnum_t           eb_size;      /* Block sizes for extra diagonal */

  cs_gnum_t           n_g_rows;     /* Global number of rows */

  cs_lnum_t           n_rows;       /* Local number of rows */
  cs_lnum_t           n_cols_ext;   /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */

  cs_lnum_t           n_elts_r[2];  /* Size of array used for restriction
                                       operations ({n_rows, n_cols_ext} when
                                       no grid merging has taken place) */

  /* Grid hierarchy information */

  const struct _cs_grid_t  *parent; /* Pointer to parent (finer) grid */

  /* Connectivity information */

  cs_lnum_t           n_faces;      /* Local number of faces */
  const cs_lnum_2_t  *face_cell;    /* Face -> cells connectivity (1 to n) */
  cs_lnum_2_t        *_face_cell;   /* Face -> cells connectivity
                                       (private array) */

  /* Restriction from parent to current level */

  cs_lnum_t          *coarse_row;   /* Fine -> coarse row connectivity;
                                       size: parent n_cols_ext
                                       -2: undetermined (initial value)
                                       -1: penalized
                                        {0, n-1}: standard mapping */
  cs_lnum_t          *coarse_face;  /* Fine -> coarse face connectivity
                                       (1 to n, signed:
                                       = 0 fine face inside coarse cell
                                       > 0 orientation same as parent
                                       < 0 orientation opposite as parent);
                                       size: parent n_faces */

  /* Face to cell data (if owner); uses same index as MSR matrix */

  cs_lnum_t    *cell_face;          /* Cell to faces adjacency */
  short int    *cell_face_sgn;      /* Cell to faces orientation */

  /* Geometric data */

  cs_real_t         relaxation;     /* P0/P1 relaxation parameter */
  const cs_real_3_t  *cell_cen;     /* Cell center (shared) */
  cs_real_3_t        *_cell_cen;    /* Cell center (private) */

  const cs_real_t    *cell_vol;     /* Cell volume (shared) */
  cs_real_t          *_cell_vol;    /* Cell volume (private) */

  cs_real_3_t        *face_normal;  /* Surface normal of internal faces.
                                       (L2 norm = face area).
                                       Not used for finest grid. */

  /* Parallel / periodic halo */

  const cs_halo_t  *halo;           /* Halo for this connectivity (shared) */
  cs_halo_t        *_halo;          /* Halo for this connectivity (private) */

  /* Matrix-related data */

  const cs_real_t  *da;             /* Diagonal (shared) */
  cs_real_t        *_da;            /* Diagonal (private) */

  const cs_real_t  *xa;             /* Extra-diagonal (shared) */
  cs_real_t        *_xa;            /* Extra-diagonal (private) */
  cs_real_t        *xa_conv;        /* Extra-diagonal (except level 0) */
  cs_real_t        *xa_diff;        /* Extra-diagonal (except level 0) */

  const cs_real_t  *xa0;            /* Symmetrized extra-diagonal (shared) */
  cs_real_t        *_xa0;           /* Symmetrized extra-diagonal (private) */
  cs_real_t        *xa0_diff;       /* Symmetrized extra-diagonal */

  cs_real_t        *xa0ij;

  cs_matrix_structure_t   *matrix_struct;  /* Associated matrix structure */
  const cs_matrix_t       *matrix;         /* Associated matrix (shared) */
  cs_matrix_t             *_matrix;        /* Associated matrix (private) */

#if defined(HAVE_MPI)

  /* Additional fields to allow merging grids */

  cs::grid_rank_merge  *rank_merge;

  int               n_ranks;           /* Number of active ranks */
  int               rank_stride;       /* Stride between active ranks
                                          <= 0 for inactive ranks) */
  MPI_Comm          comm;              /* Associated communicator */

#endif
};

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set factor to ensure diagonal dominance.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_set_diag_dom_clip_factor(double  factor);

/*----------------------------------------------------------------------------
 * Create base grid by mapping from shared mesh values.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   n_faces        <-- Local number of faces
 *   db_size        <-- Block sizes for diagonal
 *   eb_size        <-- Block sizes for extra-diagonal
 *   face_cell      <-- Face -> cells connectivity
 *   a              <-- Associated matrix
 *   conv_diff      <-- Convection-diffusion mode
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_shared(cs_lnum_t              n_faces,
                           cs_lnum_t              db_size,
                           cs_lnum_t              eb_size,
                           const cs_lnum_2_t     *face_cell,
                           const cs_matrix_t     *a,
                           bool                   conv_diff);

/*----------------------------------------------------------------------------
 * Create base grid by mapping from parent (possibly shared) matrix.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   a       <-- associated matrix
 *   n_ranks <-- number of active ranks (<= 1 to restrict to local values)
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_parent(const cs_matrix_t  *a,
                           int                 n_ranks);

/*----------------------------------------------------------------------------
 * Destroy a grid structure.
 *
 * parameters:
 *   grid <-> Pointer to grid structure pointer
 *----------------------------------------------------------------------------*/

void
cs_grid_destroy(cs_grid_t  **grid);

/*----------------------------------------------------------------------------
 * Free a grid structure's associated quantities.
 *
 * The quantities required to compute a coarser grid with relaxation from a
 * given grid are not needed after that stage, so may be freed.
 *
 * parameters:
 *   g <-> Pointer to grid structure
 *----------------------------------------------------------------------------*/

void
cs_grid_free_quantities(cs_grid_t *g);

/*----------------------------------------------------------------------------
 * Get grid information.
 *
 * parameters:
 *   g          <-- Grid structure
 *   level      --> Level in multigrid hierarchy (or NULL)
 *   symmetric  --> Symmetric matrix coefficients indicator (or NULL)
 *   db_size    --> Size of the diagonal block (or NULL)
 *   eb_size    --> Size of the extra diagonal block (or NULL)
 *   n_ranks    --> number of ranks with data (or NULL)
 *   n_rows     --> Number of local rows (or NULL)
 *   n_cols_ext --> Number of columns including ghosts (or NULL)
 *   n_entries  --> Number of entries (or NULL)
 *   n_g_rows   --> Number of global rows (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 bool             *symmetric,
                 cs_lnum_t        *db_size,
                 cs_lnum_t        *eb_size,
                 int              *n_ranks,
                 cs_lnum_t        *n_rows,
                 cs_lnum_t        *n_cols_ext,
                 cs_lnum_t        *n_entries,
                 cs_gnum_t        *n_g_rows);

/*----------------------------------------------------------------------------
 * Check if a grid has merged rank data at the current level.
 *
 * parameters:
 *   g          <-- Grid structure
 *
 * return:
 *   true if grid includes rank merge structure, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_grid_has_rank_merge(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get memory allocation mode corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   memory allocation typt
 *----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_grid_get_alloc_mode(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get number of rows corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of rows of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_rows(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get number of extended (local + ghost) columns corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of extended columns of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cols_ext(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get maximum number of extended (local + ghost) columns corresponding to
 * a grid, both with and without merging between ranks
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   maximum number of extended columns of grid structure, with or without
 *   merging
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cols_max(const cs_grid_t  *g);

/*----------------------------------------------------------------------------
 * Get global number of rows corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   global number of rows of grid structure
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_grid_get_n_g_rows(const cs_grid_t  *g);

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

/*----------------------------------------------------------------------------
 * Get the MPI subcommunicator for a given merge stride.
 *
 * parameters:
 *   parent       <-- parent MPI communicator
 *   merge_stride <-- associated merge stride
 *
 * returns:
 *   MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
cs_grid_get_comm_merge(MPI_Comm  parent,
                       int       merge_stride);

#endif // have_MPI

/*----------------------------------------------------------------------------
 * Create coarse grid from fine grid.
 *
 * parameters:
 *   f                          <-- Fine grid structure
 *   alloc_mode                 <-- Memory allocation mode
 *   coarsening_type            <-- Coarsening criteria type
 *   aggregation_limit          <-- Maximum allowed fine rows per coarse rows
 *   verbosity                  <-- Verbosity level
 *   merge_stride               <-- Associated merge stride
 *   merge_rows_mean_threshold  <-- mean number of rows under which
 *                                  merging should be applied
 *   merge_rows_glob_threshold  <-- global number of rows under which
 *                                  merging should be applied
 *   relaxation_parameter       <-- P0/P1 relaxation factor
 *   aggr_type                  --> aggregation type string
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen(cs_grid_t            *f,
                cs_alloc_mode_t       alloc_mode,
                cs_grid_coarsening_t  coarsening_type,
                int                   aggregation_limit,
                int                   verbosity,
                int                   merge_stride,
                int                   merge_rows_mean_threshold,
                cs_gnum_t             merge_rows_glob_threshold,
                double                relaxation_parameter,
                char                  aggr_type[64]);

/*----------------------------------------------------------------------------
 * Create coarse grid with only one row per rank from fine grid.
 *
 * parameters:
 *   f            <-- Fine grid structure
 *   alloc_mode   <-- Memory allocation mode
 *   merge_stride <-- Associated merge stride
 *   verbosity    <-- Verbosity level
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen_to_single(const cs_grid_t  *f,
                          cs_alloc_mode_t   alloc_mode,
                          int               merge_stride,
                          int               verbosity);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Merge grids on several ranks to a grid on a single rank
 *
 * The coarse grid must previously have been initialized with _coarse_init()
 * and its coarsening indicator determined (at least for the local cells;
 * it is extended to halo cells here if necessary).
 *
 * parameters:
 *   g           <-- Pointer to grid structure
 *   sub_size    <-- Associated merge subset size
 *   verbosity   <-- verbosity level
 *   rank_stride <-- destination rank stride, or -1 for default
 *----------------------------------------------------------------------------*/

void
cs_grid_merge_ranks(cs_grid_t  *g,
                    int         sub_size,
                    int         verbosity,
                    int         rank_stride = -1);

#endif // HAVE_MPI

/*----------------------------------------------------------------------------
 * Merge bottom grid over adjacent ranks, to reduce communicator size.
 *
 * If this would lead to a mean number of rows exceeding a given threshold,
 * no merging is done.
 *
 * parameters:
 *   g                <-- Grid structure
 *   verbosity        <-- Verbosity level
 *   n_max_ranks      <-- Maximum number of MPI ranks for bottom grid
 *   max_row_factor   <-- maximum acceptable mean ratio of merged rows
 *                        (per MPI rank) to finest rows.
 *----------------------------------------------------------------------------*/

void
cs_grid_merge_bottom(cs_grid_t      *g,
                     int             verbosity,
                     int             n_max_ranks,
                     float           max_row_factor);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*
 * \brief Gather coarse cell data in case of merged grids
 *
 * Though the grid is <const> here, its `merge_buf` member may be allocated
 * and used here.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_rank_merge_gather_row_data
(
  const cs_grid_t  *g,            /*!< grid structure */
  cs_datatype_t     data_type,    /*!< associated data type */
  size_t            data_stride,  /*!< associated data stride */
  void             *data          /*!< associated data buffer */
);

/*----------------------------------------------------------------------------*/
/*
 * \brief Scatter coarse cell data in case of merged grids
 *
 * Though the grid is <const> here, its `merge_buf` member may be allocated
 * and used here.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_rank_merge_scatter_row_data
(
  const cs_grid_t  *g,            /*!< grid structure */
  cs_datatype_t     data_type,    /*!< associated data type */
  size_t            data_stride,  /*!< associated data stride */
  void             *data          /*!< associated data buffer */
);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Project coarse grid row numbers to base grid.
 *
 * If a global coarse grid row number is larger than max_num, its
 * value modulo max_num is used.
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   max_num     <-- Values of c_row_num = global_num % max_num
 *   c_row_num   --> Global coarse row number (modulo max_num)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_row_num(const cs_grid_t  *g,
                        cs_lnum_t         n_base_rows,
                        int               max_num,
                        int               c_row_num[]);

/*----------------------------------------------------------------------------
 * Project coarse grid row rank to base grid.
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   f_row_rank  --> Global coarse row rank projected to fine rows
 *----------------------------------------------------------------------------*/

void
cs_grid_project_row_rank(const cs_grid_t  *g,
                         cs_lnum_t         n_base_rows,
                         int               f_row_rank[]);

/*----------------------------------------------------------------------------
 * Project variable from coarse grid to base grid
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   c_var       <-- Row variable on coarse grid
 *   f_var       --> Row variable projected to fine grid
 *----------------------------------------------------------------------------*/

void
cs_grid_project_var(const cs_grid_t  *g,
                    cs_lnum_t         n_base_rows,
                    const cs_real_t   c_var[],
                    cs_real_t         f_var[]);

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric and project it to base grid
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   diag_dom    --> Diagonal dominance metric (on fine grid)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_diag_dom(const cs_grid_t  *g,
                         cs_lnum_t         n_base_rows,
                         cs_real_t         diag_dom[]);

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

/*----------------------------------------------------------------------------
 * Compute coarse row variable values from fine row values
 *
 * parameters:
 *   ctx     <-> Reference to dispatch context
 *   f       <-- Fine grid structure
 *   c       <-- Fine grid structure
 *   f_var   <-- Variable defined on fine grid rows
 *   c_var   --> Variable defined on coarse grid rows
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

void
cs_grid_restrict_row_var(cs_dispatch_context  &ctx,
                         const cs_grid_t  *f,
                         const cs_grid_t  *c,
                         const cs_real_t  *f_var,
                         cs_real_t        *c_var);

/*----------------------------------------------------------------------------
 * Compute fine row variable values from coarse row values
 *
 * parameters:
 *   ctx       <-> Reference to dispatch context
 *   c         <-- Fine grid structure
 *   f         <-- Fine grid structure
 *   increment <-- if true, add value to f_var; otherwise, overwrite it
 *   c_var     <-- Variable defined on coarse grid rows
 *   f_var     <-> Variable defined on fine grid rows
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_row_var(cs_dispatch_context  &ctx,
                        const cs_grid_t      *c,
                        const cs_grid_t      *f,
                        bool                  increment,
                        cs_real_t            *c_var,
                        cs_real_t            *f_var);

/*----------------------------------------------------------------------------*/

#endif /* CS_GRID_PRIV_H */
